#########################################
# This script creates position files    #
# from previously called variants       #
# for running mpileup against normal    #
# samples to create a variant blacklist #
#                                       #
# Author: Aidan Flynn                   #
# Date: 22/05/2023                      #
#########################################


setwd("/g/data/pq08/projects/ppgl/")

################
# Dependencies #
################

library(BSgenome.Hsapiens.UCSC.hg38)
library(furrr)

################
# Data loaders #
################

source("./a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_somatic_variants(somatic_vcf_dir = NULL, quickload = T)

#############
# Functions #
#############

#This function determines the DNA sequence after deletion 
# and expands the deletion to a series of single base REF/ALT calls
expand_deletion <- function (seqnames, start, end, ALT, REF)
{
  del_len <- nchar(REF)-1
  
  padding=max(del_len*10, 100)
  padding=min(padding, 300)
  
  ref_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, seqnames, 
                    start=start, 
                    end=start + padding, 
                    width=NA,
                    strand="+", as.character=TRUE)
  
  mutant_read_seq <- ref_seq
  stringr::str_sub(string = mutant_read_seq, start = 2, end = del_len+1) <- ""
  
  ref_seq <- stringr::str_sub(string = ref_seq, start = 1, end=nchar(mutant_read_seq))
  
  covered_pos <- start:(start  + padding - del_len)
  
  return(data.frame(seqnames=seqnames, 
                    start=covered_pos, 
                    end=covered_pos, 
                    ALT=list_c(stringr::str_split(string = mutant_read_seq, pattern = "")), 
                    REF=list_c(stringr::str_split(string = ref_seq, pattern = "")),
                    id=paste(seqnames,start,end,ALT, sep="_")) )
}

#This function determines the DNA sequence after insertion 
# and expands the insertion to a series of single base REF/ALT calls
expand_insertion <- function (seqnames, start, end, ALT, REF)
{
  ins_len <- nchar(ALT)-1
  
  padding=max(ins_len*10, 100)
  padding=min(padding, 300)
  
  ref_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, seqnames, 
                    start=start, 
                    end=start  + padding, 
                    width=NA,
                    strand="+", as.character=TRUE)
  
  mutant_read_seq <- paste0(ALT,stringr::str_sub(string = ref_seq, start = 2, end=-(nchar(ALT))))
  
  covered_pos <- start:(start  + padding)
  
  return(data.frame(seqnames=seqnames, 
                    start=covered_pos, 
                    end=covered_pos, 
                    ALT=list_c(stringr::str_split(string = mutant_read_seq, pattern = "")), 
                    REF=list_c(stringr::str_split(string = ref_seq, pattern = "")),
                    id=paste(seqnames,start,end,ALT, sep="_")) 
  )
}


########################
# Subset Variant Types #
########################

snv <- a5_somatic_variants %>% ungroup() %>% 
  dplyr::select(seqnames, start, end, REF, ALT) %>% 
  filter(nchar(ALT)==1, nchar(REF)==1, ) %>% 
  distinct()

del <- a5_somatic_variants %>% ungroup() %>% 
  filter(nchar(REF) > 1, nchar(ALT) == 1) %>% 
  dplyr::select(seqnames, start, end, REF, ALT) %>% 
  distinct()

ins <- a5_somatic_variants %>% ungroup() %>% 
  filter(nchar(REF) == 1, nchar(ALT) > 1) %>% 
  dplyr::select(seqnames, start, end, REF, ALT) %>% 
  distinct()


#################
# Expand InDels #
#################

# Here I expand the indels to single base records
# eg 
#seqnames start end REF ALT
#chr1     1     4   A   ATCG
# becomes
#chr1     1     1   A   A
#chr1     2     2   C   T
#chr1     3     3   G   C
#chr1     4     4   C   G

plan(strategy = "multisession", workers = 7)

del_expanded <- furrr::future_pmap(.l = del, expand_deletion)
del_expanded <- bind_rows(del_expanded)
del_expanded_use <- del_expanded %>% filter(ALT != REF) %>% group_by(id) %>% slice_min(n=5, order_by = start)

ins_expanded <- furrr::future_pmap(.l = ins, expand_insertion)
ins_expanded <- bind_rows(ins_expanded)
ins_expanded_use <- del_expanded %>% filter(ALT != REF) %>% group_by(id) %>% slice_min(n=5, order_by = start)

######################
# Recombine Variants #
######################

var_pos <- bind_rows(snv %>% mutate(id=paste(seqnames,start,end,ALT, sep="_")), 
                     del_expanded_use, 
                     ins_expanded_use)

# postemp <- read.delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/pileups/pos.temp") 
# colnames(postemp) <- c("seqnames", "start")
# var_pos <- var_pos %>% anti_join(postemp) 

var_by_chrom <- var_pos %>% 
  group_by(seqnames) %>% 
  group_split() %>% 
  purrr::set_names(purrr::map_chr(., ~.x$seqnames[1]))

out_dir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/positions"

items_per_thread=2000
for (chr in names(var_by_chrom))
{
  if(nrow(var_by_chrom[[chr]]) > items_per_thread)
  {
    nbreaks <- ceiling(nrow(var_by_chrom[[chr]]) / items_per_thread)
    breakfactor <- rep(1:nbreaks, each=items_per_thread)
    breakfactor <- breakfactor[1:nrow(var_by_chrom[[chr]])]
    
    subsets <- split(var_by_chrom[[chr]], breakfactor)
    
    for (s in 1:nbreaks)
    {
      write.table(x = subsets[[s]], 
                  sep="\t", 
                  file = file.path(out_dir,paste0(chr,"_part",s,".txt")),
                  row.names=F,
                  col.names = F,
                  quote=F)
    }
    
  } else {
    write.table(x = var_by_chrom[[chr]], 
                sep="\t", 
                file = file.path(out_dir,paste0(chr,"_part",1,".txt")),
                row.names=F,
                col.names = F,
                quote=F)
  } 
}
