# isofox splicing analysis 
rm(list=ls())

library(tidyverse)
library(Homo.sapiens)

# isofox fisher exact results 
# NOTES:
# FetProb = Fisher exact test probabilty = (P value i think)
# ExpVal = Expected Value 

atrx_vs_tert <- read_csv("Data/bulk RNA-seq/isofox/fisher_exact/atrx_vs_tert/isofox.alt_sj_cohort_compare.csv")

# ENSG are hard to interpret so I will append gene symbols on the end
genes <- AnnotationDbi::select(Homo.sapiens,
                               keys=atrx_vs_tert$GeneId,
                               columns=c("SYMBOL", "ENSEMBL"),
                               keytype="ENSEMBL") %>%
  mutate(ENSEMBL_SYMBOL = paste0(ENSEMBL,"_",SYMBOL)) %>%
  mutate(ENSEMBL_SYMBOL = gsub("_NA$", "", x = ENSEMBL_SYMBOL)) %>% 
  distinct()

atrx_vs_tert <- atrx_vs_tert %>%
  dplyr::rename("ENSEMBL" = GeneId) %>%
  left_join(genes, by = "ENSEMBL")

colnames(atrx_vs_tert) <- colnames(atrx_vs_tert) %>%
  str_replace(pattern = "CohortA", replacement = "ATRX") %>% 
  str_replace(pattern = "CohortB", replacement = "TERT")

# ----
# Calculate the FDR-adjusted P-values
# ----

# benjamini hochberg
# TODO: I could set number of comparisons to the number of observed alt-sjs in the cohort
# I think this was what charles was suggesting but will double check this...
# the number from the isofox output is:
# 1279364 
# I could use this with bonferroni method to adjust the p-values?


# use BH method to calculate adjusted P-values
# TODO: not sure if this is legit - need to check with Charles...
atrx_vs_tert <- atrx_vs_tert %>%
  mutate(p.adjust = p.adjust(p = atrx_vs_tert$FetProb, method = "BH"))

# ----
# look at suspected ATRX splice-site mutations 
# ----

atrx_vs_tert <- atrx_vs_tert %>% 
  dplyr::select(SYMBOL, ENSEMBL, p.adjust, FetProb, everything()) %>% 
  arrange(FetProb) %>%
  write_csv("Data/bulk RNA-seq/isofox/fisher_exact/atrx_vs_tert/isofox_tert_vs_atrx_fisher_exact.csv")
