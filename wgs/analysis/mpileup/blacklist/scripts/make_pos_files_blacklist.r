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
# Data loaders #
################

source("./a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_somatic_variants(somatic_vcf_dir = NULL, quickload = T)


########################
# Print Position Files #
########################

var_pos <- a5_somatic_variants %>% ungroup() %>% 
  dplyr::select(seqnames, start, end, REF, ALT) %>% 
  mutate(type=case_when(
  nchar(ALT)==1 & nchar(REF)==1 ~ "SBS",
  nchar(REF) > 1 & nchar(ALT) == 1 ~ "DEL",
  nchar(REF) == 1 & nchar(ALT) > 1 ~ "INS")) %>% 
  distinct()

postemp <- read.delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/blacklist/positions/pos.done") 
colnames(postemp) <- c("seqnames", "start")
var_pos <- var_pos %>% anti_join(postemp) 

var_by_chrom <- var_pos %>% 
  group_by(seqnames) %>% 
  group_split() %>% 
  purrr::set_names(purrr::map_chr(., ~.x$seqnames[1]))

out_dir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/blacklist/positions"

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
