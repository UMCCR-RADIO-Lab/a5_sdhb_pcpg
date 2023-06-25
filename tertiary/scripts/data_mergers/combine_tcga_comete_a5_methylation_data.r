######################################################################
# Script for loading and merging TCGA/A5/Flynn methylation data sets #
# Author: Aidan Flynn                                                #
# Date: 30/11/2022                                                   #
# Languages: R                                                       #
######################################################################

old_wd <- getwd()
setwd("/g/data/pq08/projects/ppgl")

library(limma)
library(edgeR)
library(tidyverse)

#######################
# Import data loaders #
#######################

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
source("./public_data/data_loaders/comete_methylation_dataloader.r")
source("./public_data/data_loaders/tcga_methylation_dataloader.r")
source("./a5/methylation/scripts/data_loaders/a5_methylation_dataloader.r")

######################
# Colours and Themes #
######################

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

##################
# Load data sets #
##################

if(!exists("a5_anno")){
  data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.org", use_cache = T)
}

message("Loading COMETE methylation array data ...")
data_loader_comete_methylation(remove_27k_availaible_on_450k = T, only_zethoven_annotated = T)

message("Loading TCGA methylation array data ...")
data_loader_tcga_methylation(quickload = T, normalisation = "functional")

tcga_methylation_450k_beta <- 
  data.frame(getBeta(tcga_methylation_450k_filtered), check.names = F) %>% 
  tibble::rownames_to_column("Probe")

tcga_methylation_450k_mval <-  
  data.frame(getM(tcga_methylation_450k_filtered), check.names = F) %>% 
  tibble::rownames_to_column("Probe")

message("Loading A5 methylation array data ...")
data_loader_a5_methylation_array(quickload = T, remove_excluded_samples=T, normalisation = "functional")

a5_methylation_beta <- 
  data.frame(getBeta(a5_methylation_filtered), check.names = F) %>% 
  tibble::rownames_to_column("Probe")

a5_methylation_mval <-  
  data.frame(getM(a5_methylation_filtered), check.names = F) %>% 
  tibble::rownames_to_column("Probe")

#####################
# Combine data sets #
#####################
message("Merging 27K datasets ...")
meth_tcga_comete_a5_27k_mval <- comete_methylation_27k_mval %>% 
  inner_join(comete_methylation_450k_mval, by="Probe") %>% 
  inner_join(tcga_methylation_450k_mval, by="Probe") %>% 
  inner_join(a5_methylation_mval, by="Probe")

meth_tcga_comete_a5_27k_beta <- comete_methylation_27k_beta %>% 
  inner_join(comete_methylation_450k_beta, by="Probe") %>% 
  inner_join(tcga_methylation_450k_beta, by="Probe") %>% 
  inner_join(a5_methylation_beta, by="Probe")

message("Merging 450K datasets ...")
meth_tcga_comete_a5_450k_mval <- comete_methylation_450k_mval %>% 
  inner_join(tcga_methylation_450k_mval, by="Probe") %>% 
  inner_join(a5_methylation_mval, by="Probe")

meth_tcga_comete_a5_450k_beta <- comete_methylation_450k_beta %>% 
  inner_join(tcga_methylation_450k_beta, by="Probe") %>% 
  inner_join(a5_methylation_beta, by="Probe")

#Reformat A5 annotation to match zethoven annotation
message("Preparing annotation ...")
a5_anno_zethoven_format <- 
  a5_to_zethoven_anno_format(a5_anno %>% 
                               filter(A5_ID %in% a5_methylation_targets$Sample_Name,
                                      Exclude!="Y")) %>% 
  inner_join(a5_anno %>% dplyr::select(A5_ID, methylation_batch), by=c("Sample"="A5_ID")) %>% 
  mutate(methylation_batch=paste("A5", methylation_batch, sep="-")) 
  

#Combine sample annotation
meth_tcga_comete_a5_anno <- bind_rows(comete_methylation_27k_anno %>% mutate(across(everything(),as.character)),
                           comete_methylation_450k_anno %>% mutate(across(everything(),as.character)),
                           tcga_methylation_450k_anno %>% mutate(across(everything(),as.character))) %>% 
  dplyr::select(Sample, Cluster, Dataset, Sex, Castro_Vega_miRNA_Classification, Malignancy, Genotype, methylation_batch, new_naming) %>%
  bind_rows(a5_anno_zethoven_format) %>% 
  mutate(Malignancy=gsub("Malignant","Metastatic",Malignancy), 
         Malignancy=gsub("Benign","Non-Metastatic",Malignancy)) 

#27k sample annotation
meth_tcga_comete_a5_27k_anno <- meth_tcga_comete_a5_anno %>% filter(Sample %in% colnames(meth_tcga_comete_a5_27k_beta))

#450k sample annotation
meth_tcga_comete_a5_450k_anno <- meth_tcga_comete_a5_anno %>% filter(Sample %in% colnames(meth_tcga_comete_a5_450k_beta))  

rm(meth_tcga_comete_a5_anno)

#Matching columns/samples sanity check
if (length(setdiff(colnames(meth_tcga_comete_a5_27k_beta)[-1], meth_tcga_comete_a5_27k_anno$Sample)) > 0 | 
    length(setdiff(meth_tcga_comete_a5_27k_anno$Sample, colnames(meth_tcga_comete_a5_27k_beta)[-1])) > 0 |
    length(setdiff(colnames(meth_tcga_comete_a5_450k_beta)[-1], meth_tcga_comete_a5_450k_anno$Sample)) > 0 | 
    length(setdiff(meth_tcga_comete_a5_450k_anno$Sample, colnames(meth_tcga_comete_a5_450k_beta)[-1])) > 0)
    { stop("combine_tcga_comete_a5_methylation_data: Annotation sample list and data columns do not match") }

#Reorder annotation to match data
meth_tcga_comete_a5_27k_anno <- meth_tcga_comete_a5_27k_anno[match(colnames(meth_tcga_comete_a5_27k_beta)[-1], meth_tcga_comete_a5_27k_anno$Sample),]
meth_tcga_comete_a5_450k_anno <- meth_tcga_comete_a5_450k_anno[match(colnames(meth_tcga_comete_a5_450k_beta)[-1], meth_tcga_comete_a5_450k_anno$Sample),]

# Sanity check for ordering
if(!(all(colnames(meth_tcga_comete_a5_27k_beta)[-1]==meth_tcga_comete_a5_27k_anno$Sample) & 
     all(colnames(meth_tcga_comete_a5_450k_beta)[-1]==meth_tcga_comete_a5_450k_anno$Sample))) {
  stop("Annotation ordering sanity check failed")
}

#harmonise m-value and beta value data-frames
meth_tcga_comete_a5_27k_mval <- meth_tcga_comete_a5_27k_mval[,colnames(meth_tcga_comete_a5_27k_beta)]
meth_tcga_comete_a5_450k_mval <- meth_tcga_comete_a5_450k_mval[,colnames(meth_tcga_comete_a5_450k_beta)]

#Convert to matrix objects
mval_all_27k.matrix <- as.matrix(meth_tcga_comete_a5_27k_mval[,-1])
mval_all_450k.matrix <- as.matrix(meth_tcga_comete_a5_450k_mval[,-1])
rownames(mval_all_27k.matrix) <- meth_tcga_comete_a5_27k_mval$Probe
rownames(mval_all_450k.matrix) <- meth_tcga_comete_a5_450k_mval$Probe

message("Performing dataset/batch correction ...")
#Create design matrices
design_27k <- model.matrix(~0 + Genotype, data = meth_tcga_comete_a5_27k_anno)
colnames(design_27k) <- gsub("Genotype| ","", colnames(design_27k))
rownames(design_27k) <- colnames(mval_all_27k.matrix)
design_450k <- model.matrix(~0 + Genotype, data = meth_tcga_comete_a5_450k_anno)
colnames(design_450k) <- gsub("Genotype| ","", colnames(design_450k))
rownames(design_450k) <- colnames(mval_all_450k.matrix)

# Remove the library specific batch effects
meth_tcga_comete_a5_27k_mval.batch_removed <-
  removeBatchEffect(
    mval_all_27k.matrix,
    batch = meth_tcga_comete_a5_27k_anno$methylation_batch,
    batch2 = meth_tcga_comete_a5_27k_anno$Sex,
    design = design_27k
  )
meth_tcga_comete_a5_27k_mval.batch_removed <- meth_tcga_comete_a5_27k_mval.batch_removed %>% data.frame(check.names = F) %>% filter(!if_any(.fns=is.infinite))

meth_tcga_comete_a5_450k_mval.batch_removed <-
  removeBatchEffect(
    mval_all_450k.matrix,
    batch = meth_tcga_comete_a5_450k_anno$methylation_batch,
    batch2 = meth_tcga_comete_a5_450k_anno$Sex,
    design = design_450k
  )
meth_tcga_comete_a5_450k_mval.batch_removed <- meth_tcga_comete_a5_450k_mval.batch_removed %>% data.frame(check.names = F) %>% filter(!if_any(.fns=is.infinite))

message("Loaded objects into global environment:
          meth_tcga_comete_a5_27k_mval (raw M-values for 27k probes - ",ncol(meth_tcga_comete_a5_27k_mval)-1," arrays),
          meth_tcga_comete_a5_450k_mval (raw M-values for 450K probes - ",ncol(meth_tcga_comete_a5_450k_mval)-1," arrays),
          meth_tcga_comete_a5_27k_mval.batch_removed (batch corrected M-values for 27k probes),
          meth_tcga_comete_a5_450k_mval.batch_removed (batch corrected M-values for 450K probes)
          meth_tcga_comete_a5_27k_anno (sample annotation 27k probeset),
          meth_tcga_comete_a5_450k_anno (sample annotation 450K probeset)")

############
# Clean up #
############
rm(mval_all_27k.matrix)
rm(mval_all_450k.matrix)
rm(design_27k)
rm(design_450k)
rm(a5_anno_zethoven_format)
setwd(old_wd)
rm(old_wd)