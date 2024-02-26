#########################################################################################
# Script for loading and merging TCGA/A5/Flynn Whole Transcriptome Sequencing data sets #
# Author: Aidan Flynn                                                                   #
# Date: 30/11/2022                                                                      #
# Languages: R                                                                          #
#########################################################################################

old_wd <- getwd()
setwd("/g/data/pq08/projects/ppgl")

library(tidyverse)
library(SummarizedExperiment)
library(limma)
library(edgeR)

#######################
# Import data loaders #
#######################

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
source("./a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
source("./public_data/data_loaders/tcga_wts_dataloader.r")
source("./public_data/data_loaders/flynn_wts_dataloader.r")

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#############
# Load data #
#############

#Load A5 Samples
htseq_outs <- "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/truseq/gene/"
data_loader_a5_wts_counts(htseq_outs)

#load TCGA count data
data_loader_tcga_wts_counts(remove_normals = F)

# Get the counts that Magnus sent from 'The genomic landscape of phaeochromocytoma' (Flynn et al)
data_loader_flynn_wts_counts()



##############
# Merge data #
##############

wts_tcga_flynn_a5_counts <- flynn_wts_counts %>% 
  inner_join(tcga_wts_counts) %>% 
  inner_join(a5_wts_counts %>% 
               tibble::rownames_to_column("ENSGID_Symbol") %>% 
               separate(col=ENSGID_Symbol, into=c("ensembl_gene_id","Symbol"), sep="_", extra="merge") %>% 
               mutate(ensembl_gene_id=gsub("(ENSG[0-9]+)[.][0-9]+","\\1",ensembl_gene_id)) %>%
               filter(!grepl("PAR_Y_", Symbol)) %>% dplyr::select(-Symbol)) 

####################
# Merge Annotation #
####################

if(!exists("a5_anno")) {
  data_loader_a5_clinical_anno(use_cache = T) }

a5_anno_zethoven_format <- a5_to_zethoven_anno_format(a5_anno)

keep_a5_samples <- colnames(a5_wts_dge_list$qc_ok)

a5_anno_zethoven_format <- a5_anno_zethoven_format %>% filter(Sample %in% keep_a5_samples)

wts_tcga_flynn_a5_anno <- bind_rows(tcga_wts_anno, flynn_wts_anno, a5_anno_zethoven_format) %>% 
  mutate(Malignancy=gsub("Malignant","Metastatic",Malignancy), 
         Malignancy=gsub("Benign","Non-Metastatic",Malignancy))

# Keep only samples with both annotation and RNA-seq data
available_samples <- intersect(wts_tcga_flynn_a5_anno$Sample,colnames(wts_tcga_flynn_a5_counts)[-1])
wts_tcga_flynn_a5_counts <- wts_tcga_flynn_a5_counts %>% dplyr::select(ensembl_gene_id, !!available_samples)

wts_tcga_flynn_a5_anno <- wts_tcga_flynn_a5_anno %>% 
  filter(Sample %in% intersect(colnames(wts_tcga_flynn_a5_counts), 
                                wts_tcga_flynn_a5_anno$Sample))

wts_tcga_flynn_a5_anno <- wts_tcga_flynn_a5_anno %>%  
  mutate(ATRX_mutation = ifelse(is.na(ATRX_mutation) | ATRX_mutation == "0", 
                                FALSE, 
                                TRUE)) %>% 
  mutate(ATRX_mutation = ifelse(Sample %in% c("V-PH-22T"," V-PH-23T", a5_anno %>%  filter(TERT_ATRX_Mutation == "ATRX") %>%  pull(A5_ID)), 
                                TRUE, 
                                ATRX_mutation)) 

wts_tcga_flynn_a5_anno <- wts_tcga_flynn_a5_anno %>%  
  mutate(TERT_mutation = ifelse(Sample %in% c("V-PH-03T", "V-PH-04T", a5_anno %>%  filter(TERT_ATRX_Mutation == "TERT") %>%  pull(A5_ID)), 
                                TRUE, 
                                FALSE)) 

####################
# Batch Correction #
####################

#Convert to matrix
wts_tcga_flynn_a5_counts_mat <- 
  wts_tcga_flynn_a5_counts %>%
  dplyr::select(-ensembl_gene_id) %>%
  as.matrix()

# Reset rownames
rownames(wts_tcga_flynn_a5_counts_mat) <- wts_tcga_flynn_a5_counts$ensembl_gene_id

#Convert to DGE object
wts_tcga_flynn_a5_dge <- DGEList(wts_tcga_flynn_a5_counts_mat)

# Filter lowly expressed genes
keep.exprs <- filterByExpr(wts_tcga_flynn_a5_dge, group = wts_tcga_flynn_a5_anno$Genotype)
wts_tcga_flynn_a5_dge <- wts_tcga_flynn_a5_dge[keep.exprs,, keep.lib.sizes=FALSE]

# Apply TMM normalisation to the DGElist
wts_tcga_flynn_a5_dge <- calcNormFactors(wts_tcga_flynn_a5_dge)
# Convert the TMM counts to log2 CPMs
wts_tcga_flynn_a5_lcpm <- edgeR::cpm(wts_tcga_flynn_a5_dge, log = T)

design_genotype <- model.matrix(~0 + Genotype, data = wts_tcga_flynn_a5_anno)
colnames(design_genotype) <- gsub("Genotype| ","", colnames(design_genotype))
rownames(design_genotype) <- rownames(wts_tcga_flynn_a5_dge$samples)

# Remove the library specific batch effects
wts_tcga_flynn_a5_lcpm.batch_removed <- removeBatchEffect(wts_tcga_flynn_a5_lcpm, batch=wts_tcga_flynn_a5_anno$Dataset,design = design_genotype)

message("Loaded objects into global environment:
          wts_tcga_flynn_a5_counts (raw counts dataframe for merged dataset),
          wts_tcga_flynn_a5_lcpm.batch_removed (batch corrected log2-cpm values)
          wts_tcga_flynn_a5_anno (sample annotation)")

############
# Clean up #
############

rm(wts_tcga_flynn_a5_counts_mat)
rm(a5_anno_zethoven_format)
rm(keep_a5_samples)
setwd(old_wd)
rm(old_wd)