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

##################
# Paired Samples #
##################

Samples <- list(
  E122=c("E122-1","E122-2"),
  E128=c("E128-1","E128-2"),
  E132=c("E132-1","E132-2"),
  E136=c("E136-1","E136-2"),
  E143=c("E143-1","E143-2","E143-3"),
  E146=c("E146-1","E146-2"),
  E158=c("E158-1","E158-2"),
  E159=c("E159-1","E159-2","E159-3","E159-4"),
  E166=c("E166-1","E166-2"),
  E167=c("E167-1","E167-2"),
  E169=c("E169-1","E169-2"),
  E225=c("E225-1","E225-2"),
  E229=c("E229-1","E229-2")
)


########################
# Print Position Files #
########################

out_dir="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/paired_samples/positions/"

for (patient in names(Samples))
{
  var_pos <- a5_somatic_variants %>% 
    filter(A5_ID %in% Samples[[patient]]) %>% ungroup() %>% 
    dplyr::select(seqnames, start, end, REF, ALT)  %>% 
    mutate(type=case_when(
      nchar(ALT)==1 & nchar(REF)==1 ~ "SBS",
      nchar(REF) > 1 & nchar(ALT) == 1 ~ "DEL",
      nchar(REF) == 1 & nchar(ALT) > 1 ~ "INS")) %>% 
    distinct()
  
  write.table(x = var_pos, 
              sep="\t", 
              file = file.path(out_dir,paste0(patient,".txt")),
              row.names=F,
              col.names = F,
              quote=F)
}





  