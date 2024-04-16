#!/bin/python3

#This script merges all the first pass splice junction output from STAR across all A5 samples for use in generating a custom genome 

from pathlib import Path
from math import floor

out_file="/g/data/pq08/projects/ppgl/a5/wts/reference/star_second_pass_genome_neb/A5_pooled_SJ_neb.out.tab"

files_processed=0

splice_dict = dict()

for path in Path('/g/data/pq08/projects/ppgl/a5/wts/results/star/neb').rglob('*_SJ.out.tab'):
    
    files_processed += 1
    with open(path) as fh:
        for entry in fh:
            chromosome, intron_start, intron_end, strand, intron_motif, is_annotated, n_unq_mapped, n_multi_mapped, max_overhang = entry.strip().split("\t") 

            splice_key=chromosome + "&" + intron_start + "&" + intron_end + "&" + strand + "&" + intron_motif + "&" + is_annotated
            

            if splice_key not in splice_dict.keys():
                splice_dict[splice_key] = [0,0,0]

            splice_dict[splice_key][0] += int(n_unq_mapped)
            splice_dict[splice_key][1] += int(n_multi_mapped)
            splice_dict[splice_key][2] += int(max_overhang)
    
    print("Processed: " + str(files_processed))

with open(out_file,"w") as fh:
    for splice_key in splice_dict.keys():
        chromosome, intron_start, intron_end, strand, intron_motif, is_annotated = splice_key.split("&") 
        n_unq_mapped =  splice_dict[splice_key][0]
        n_multi_mapped =  splice_dict[splice_key][1]
        max_overhang  =  splice_dict[splice_key][2]
        if n_unq_mapped > 10:
            print("\t".join([chromosome, intron_start, intron_end, strand, intron_motif, is_annotated, str(n_unq_mapped), str(n_multi_mapped), str(floor(int(max_overhang)/files_processed))]),file=fh)
    

