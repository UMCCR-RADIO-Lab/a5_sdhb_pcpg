#!/usr/bin/python3
import os
import sys
import argparse
import re
import math

parser = argparse.ArgumentParser()

parser.add_argument('--pileup', dest='pileup', help='mPileup file', required=True)
parser.add_argument('--outputValue', dest='output', help='The output value: VAF, AltCount, AltRefTotal, or Total', required=False, default="VAF")
parser.add_argument('--pileupSamples', dest='pileupSamples', help='list/order of samples in mPileup ', required=True, nargs='*')
parser.add_argument('--variantCalls', dest='variantCalls', help='Tab-sep variant calls (chr,pos, ref, alt)', required=True)
parser.add_argument('--minBaseQuality', dest='minBaseQuality', help='Minimum base quality to count read', required=False, type=int, default=13)
parser.add_argument('--OneForAll', dest='OneForAll', help='If one read is above minBaseQuality then count all reads', required=False, action="store_true")
parser.add_argument('--OneForAllMode', dest='OneForAllMode', help='Sample=Current sample must have HQ reads, Cohort=One of the samples must have HQ read', required=False, default="Cohort")
parser.add_argument('--debugLevel', dest='debugLevel', help='The desired debug output level', required=False, default=0)

args = parser.parse_args()

class BaseCall:
    
    def __init__(self, Base, quality_threshold=0):
        self.base = Base
        self.raw_fwd_support = 0
        self.raw_rvrs_support = 0
        self.qf_fwd_support = 0
        self.qf_rvrs_support = 0
        self.fwd_qualities = []
        self.rvrs_qualities = []
        self.quality_threshold = quality_threshold
        self.has_deletion = False
        self.has_insertion = False
 
    def ascii2numeric(self, ascii_bq):
        return ord(ascii_bq) - 33
 
    def add_support(self, orientation, base_quality):
        if orientation == "fwd":
            self.raw_fwd_support = self.raw_fwd_support + 1
            if self.ascii2numeric(base_quality) >= self.quality_threshold:
                self.qf_fwd_support = self.qf_fwd_support + 1
            self.fwd_qualities.append(base_quality)            
        
        if orientation == "rvrs":
            self.raw_rvrs_support = self.raw_rvrs_support + 1
            if self.ascii2numeric(base_quality) >= self.quality_threshold:
                self.qf_rvrs_support = self.qf_rvrs_support + 1
            self.rvrs_qualities.append(base_quality)

    def support(self, OneForAll=False, raw=False):
        #print(" -- Base:%s,rawf:%s,rawr:%s, qff:%s, qfr:%s -- " % (self.base, self.raw_fwd_support, self.raw_rvrs_support, self.qf_fwd_support, self.qf_rvrs_support))
        if raw:
            return self.raw_fwd_support + self.raw_rvrs_support
        else:
            if not OneForAll:
                return self.qf_fwd_support + self.qf_rvrs_support
            else:
                
                if self.qf_fwd_support > 0 or self.qf_rvrs_support > 0:
                    return self.raw_fwd_support + self.raw_rvrs_support
                else:
                    return 0
            
    def has_bidirectional(self):
        if self.qf_fwd_support > 0 and self.qf_rvrs_support > 0: return True
        return False

    def mean_quality(self):
        tq = 0
        for q in self.fwd_qualities + self. rvrs_qualities:
            nq = ord(q) - 33
            tq = tq + nq
        return round(tq/len(self.fwd_qualities + self. rvrs_qualities),2)
  
    def min_quality(self):
        qv = list(map(lambda x: ord(x) - 33, self.fwd_qualities + self.rvrs_qualities))
        qv.sort()
        return qv[0]

    def median_quality(self):
        qv = list(map(lambda x: ord(x) - 33, self.fwd_qualities + self.rvrs_qualities))
        if len(qv) > 2:
            qv.sort()
            median_index1 = math.floor((len(qv)+1)/2)
            median_index2 = math.ceil((len(qv)+1)/2)
            return ((qv[median_index1] + qv[median_index2])/2)
        elif len(qv) == 2: 
            return ((qv[0] + qv[1])/2)
        elif len(qv) == 1:
            return qv[0]

    def max_quality(self):
        qv = list(map(lambda x: ord(x) - 33, self.fwd_qualities + self.rvrs_qualities))
        qv.sort()
        return qv[len(qv)-1]

    def quality_summary(self):
        return "min:%s,mean:%s,med:%s,max:%s,rawF:%s,rawR:%s,qfF:%s,qfR:%s" % (self.min_quality(), self.mean_quality(),self.median_quality(), self.max_quality(), self.raw_fwd_support, self.raw_rvrs_support, self.qf_fwd_support, self.qf_rvrs_support)


def process_indels(bases, qualities):
    
    insertions = []
    deletions = []
    while indel := re.search("[.,ACGTNacgtn*#][+-][0-9]+", string=bases):
        match_start = indel.span()[0]
        match_type = bases[match_start+1:match_start+2]
        len_start = match_start + 2
        len_end = indel.span()[1]
        ins_len = int(bases[len_start:len_end])
        ins_seq = bases[(len_end):(len_end + ins_len)]
        
        if match_type == "+":
            insertions.append(ins_seq)
        else:
            deletions.append(ins_seq)

        bases=bases[0:len_start-2] + ("_" * ((len_end - match_start) + ins_len)) + bases[len_end + ins_len:]
        qualities=qualities[0:match_start] + ("_" * ((len_end - match_start) + ins_len)) + qualities[match_start + 1:]


    return ((bases, qualities, insertions, deletions))


colDefs = ["chr","pos", "ref"]
for s in args.pileupSamples:
    colDefs = colDefs + ["%s$total" % s, "%s$calls" % s, "%s$quality" % s]

VAFs = dict()
Totals = dict()
with open(args.pileup, 'r') as pileupfh:
        for row in pileupfh:
            fields = row.split("\t")
            Chr = fields[0]
            Pos = fields[1]
            Ref = fields[2]

            if (Chr not in VAFs.keys()):
                VAFs[Chr] = dict()
                Totals[Chr] = dict()

            if (Pos not in VAFs[Chr].keys()):
                VAFs[Chr][Pos] = dict()
                Totals[Chr][Pos] = dict()

            for i in range(1,len(fields)):
                if "total" in colDefs[i]:
                    s = colDefs[i].split("$")[0]
                    Totals[Chr][Pos][s] = fields[i]
                if "calls" in colDefs[i]:
                    s = colDefs[i].split("$")[0]
                    VAFs[Chr][Pos][s] = dict()
                    bases = fields[i]
                    bases = re.sub("[$><#]|\^.", "", bases)
                    qualities = fields[i+1].strip()

                    bases, qualities, insertions, deletions = process_indels(bases, qualities)
                    
                    
                    if(len(bases) != len(qualities)):
                        print("Base and quality string equality sanity check failed for %s %s Exiting. " % (Chr, Pos))
                        print(bases)
                        print(qualities)
                        exit(11)

                    for insertion in insertions:
                        ins_key = Ref + "/" + Ref + insertion
                        ins_key = ins_key.upper()
                        if (ins_key not in VAFs[Chr][Pos][s].keys()):
                            VAFs[Chr][Pos][s][ins_key] = BaseCall(insertion, args.minBaseQuality)
                            
                        if (insertion == insertion.upper()):
                            VAFs[Chr][Pos][s][ins_key].add_support("fwd", "F") #F=phread-33
                        else:
                            VAFs[Chr][Pos][s][ins_key].add_support("rvrs", "F")

                    for deletion in deletions:
                        del_key = Ref + deletion + "/" + Ref
                        del_key = del_key.upper()
                        if (del_key not in VAFs[Chr][Pos][s].keys()):
                            VAFs[Chr][Pos][s][del_key] = BaseCall(deletion, args.minBaseQuality)
                            
                        if (deletion == deletion.upper()):
                            VAFs[Chr][Pos][s][del_key].add_support("fwd", "F") #F=phread-33
                        else:

                            VAFs[Chr][Pos][s][del_key].add_support("rvrs", "F")

                    for q in range(0,len(bases)-1):
                        b=bases[q].upper()
                        
                        if b == "_": continue #spacer inserted for indel

                        if b in ".,":
                            sbs_key = Ref + "/" + Ref
                        else:
                            sbs_key = Ref + "/" + b

                        if (sbs_key not in VAFs[Chr][Pos][s].keys()):
                            VAFs[Chr][Pos][s][sbs_key] = BaseCall(b, args.minBaseQuality)
                            
                        if (b==bases[q]):
                            VAFs[Chr][Pos][s][sbs_key].add_support("fwd", qualities[q])
                        else:
                            VAFs[Chr][Pos][s][sbs_key].add_support("rvrs", qualities[q])


header="CHROM\tPOS\tREF\tALT\t"
if args.output.lower() == "vaf":
    header = header + '\t'.join([s+"_vaf" for s in args.pileupSamples])    
elif args.output.lower() == "altcount":
    header = header + '\t'.join([s+"_altCount" for s in args.pileupSamples])    
elif args.output.lower() == "altreftotal":
    header = header + '\t'.join(["%s_altCount\t%s_refCount\t%s_TotalCount" % (s,s,s) for s in args.pileupSamples])    
elif args.output.lower() == "total":
     header = header + '\t'.join([s+"_TotalCount" for s in args.pileupSamples])  
elif args.output.lower() == "all":
    header = header + '\t'.join(["%s_vaf\t%s_altCount\t%s_refCount\t%s_TotalCount" % (s,s,s,s) for s in args.pileupSamples])  
elif args.output.lower() == "allqc":
    header = header + '\t'.join(["%s_rawFreq\t%s_vaf\t%s_altCount\t%s_refCount\t%s_TotalCount\t%s_DepthQCFilt\t%s_QC" % (s,s,s,s,s,s,s) for s in args.pileupSamples])  
else:
    print ("Invalid output type")
    exit(808)

print (header)              

with open(args.variantCalls, 'r') as mutectfh:
    for row in mutectfh:
        if "ALT" in row: continue
        fields = row.split("\t")
        Chr = fields[0]
        Pos = fields[1]
        Ref = fields[2]
        Alt = fields[3].strip() 

        mut_key = Ref + "/" + Alt
        if len(Ref) == 1 and len(Alt) == 1:
            #mut_type="SBS"
            ref_key = Ref + "/" + Ref               
        elif len(Ref) > 1 and len(Alt) == 1:
            #mut_type="Deletion"
            ref_key = Alt + "/" + Alt                            
        elif len(Ref) == 1 and len(Alt) > 1:
            #mut_type="Insertion"
            ref_key = Ref + "/" + Ref             

        print("%s\t%s\t%s\t%s" % (Chr, Pos, Ref, Alt), end="")
        for s in args.pileupSamples:
            altCount_raw = 0
            altCount = 0
            try:
                altCount_raw = float(VAFs[Chr][Pos][s][mut_key].support(raw=True))
                if args.OneForAllMode == "Cohort":
                    for c in args.pileupSamples:
                        try:
                            if float(VAFs[Chr][Pos][c][mut_key].support(OneForAll=False)) > 0:
                                altCount = float(VAFs[Chr][Pos][s][mut_key].support(OneForAll=True, raw=True))
                                break
                        except KeyError:
                            continue
                else:
                    altCount = float(VAFs[Chr][Pos][s][mut_key].support(args.OneForAll))
            except: 
                if int(args.debugLevel) > 3:
                    print("Missing Key -defaulting AltCount to 0 \n")
                altCount= 0          
            try:
               key = Ref + "/" + Ref
               refCount = float(VAFs[Chr][Pos][s][ref_key].support(args.OneForAll))
            except: 
                refCount= 0       
            
            
            
            
            TotalCount_QC = 0
            try:                
                TotalCount =  float(Totals[Chr][Pos][s])
                for t in VAFs[Chr][Pos][s]:
                    TotalCount_QC = TotalCount_QC + VAFs[Chr][Pos][s][t].support(args.OneForAll)                
            except KeyError:
                if (int(args.debugLevel) > 3):
                    print("\nNo reads in %s at %s %s %s %s" % (s, Chr, Pos, Ref, Alt))
                TotalCount = 0

            try:
                vaf = round(altCount/TotalCount,3)
                raw_vaf = round(altCount_raw/TotalCount,3)
            except ZeroDivisionError:
                vaf="NA"
                raw_vaf="NA"

            if args.output.lower() == "vaf":
                print("\t%s" % vaf, end="")
            elif args.output.lower() == "altcount":
                print("\t%s" % int(altCount), end="")
            elif args.output.lower() == "altreftotal":
                print("\t%s\t%s\t%s" % (int(altCount), int(refCount), int(TotalCount)), end="")
            elif args.output == "Total":
                print("\t%s" % (int(TotalCount)), end="")
            elif args.output.lower() == "all":
                print("\t%s\t%s\t%s\t%s" % (vaf, int(altCount), int(refCount), int(TotalCount)), end="")
            elif args.output.lower() == "allqc":
                try:
                    qc = VAFs[Chr][Pos][s][Alt].quality_summary()
                except KeyError:
                    qc = ""
                print("\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (raw_vaf, vaf, int(altCount), int(refCount), int(TotalCount), int(TotalCount_QC), qc), end="")                
            else:
                print ("Invalid output type")
                exit(808)
            
        print("")


