#####################################################################
# Script for loading and merging TCGA/A5/COMETE small-RNA data sets #
# Author: Aidan Flynn                                               #
# Date: 30/11/2022                                                  #
# Languages: R                                                      #
#####################################################################

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
source("./public_data/data_loaders/comete_smallrna_dataloader.r")
source("./public_data/data_loaders/tcga_smallrna_dataloader.r")
source("./a5/small_rna/scripts/data_loaders/a5_smallrna_seq_dataloader.r")

######################
# Colours and Themes #
######################

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#############
# Load data #
#############

data_loader_tcga_smallrna(genome_version="hg38")
data_loader_comete_smallrna(genome_version="hg38", only_zethoven_annotated=T, remove_outlier_samples=T)
data_loader_a5_smallrna(genome_version="hg38")

if(!exists("a5_anno")) {
  data_loader_a5_clinical_anno(use_cache = T) }

##############
# Merge data #
##############

# Sanity check labels are the same
if(
sum(tcga_smallrna_read_counts$Geneid == a5_smallrna_read_counts$Geneid) != nrow(a5_smallrna_read_counts) |
sum(comete_smallrna_read_counts$Geneid == a5_smallrna_read_counts$Geneid) != nrow(a5_smallrna_read_counts)) {
  stop("Sanity check for row counts and ordering across datasets failed")
}

a5_keep_samples <- colnames(a5_smallrna_dge_list$qc_ok)

# Combine the counts based on geneid
smallrna_tcga_comete_a5_counts <- 
  tcga_smallrna_read_counts %>% 
  left_join(a5_smallrna_read_counts %>% dplyr::select(Geneid,all_of(a5_keep_samples))) %>% 
  left_join(comete_smallrna_read_counts)


smallrna_tcga_comete_a5_counts_mat <- as.matrix(smallrna_tcga_comete_a5_counts[,-1])
rownames(smallrna_tcga_comete_a5_counts_mat) <- as.matrix(smallrna_tcga_comete_a5_counts$Geneid)


a5_anno_zethoven_format <- a5_to_zethoven_anno_format(a5_anno)
  
smallrna_tcga_comete_a5_anno <- bind_rows(comete_smallrna_anno %>% mutate(across(everything(),as.character)),
                        tcga_smallrna_anno %>% mutate(across(everything(),as.character))) %>% 
  dplyr::select(Sample, Dataset, Sex, Castro_Vega_miRNA_Classification, Malignancy, Genotype, new_naming) %>%
  bind_rows(a5_anno_zethoven_format) %>% 
  mutate(Malignancy=gsub("Malignant","Metastatic",Malignancy), 
         Malignancy=gsub("Benign","Non-Metastatic",Malignancy)) %>% 
  filter(Sample %in% colnames(smallrna_tcga_comete_a5_counts_mat))

smallrna_tcga_comete_a5_dge <- DGEList(smallrna_tcga_comete_a5_counts_mat[,smallrna_tcga_comete_a5_anno$Sample])

# Filter lowly expressed genes
keep.exprs <- filterByExpr(smallrna_tcga_comete_a5_dge, group = smallrna_tcga_comete_a5_anno$new_naming)
smallrna_tcga_comete_a5_dge <- smallrna_tcga_comete_a5_dge[keep.exprs,, keep.lib.sizes=FALSE]

# Apply TMM normalisation to the DGElist
smallrna_tcga_comete_a5_dge <- calcNormFactors(smallrna_tcga_comete_a5_dge)
# Convert the TMM counts to log2 CPMs
smallrna_tcga_comete_a5_lcpm <- edgeR::cpm(smallrna_tcga_comete_a5_dge, log = T)

# Sanity check for ordering
if(sum(smallrna_tcga_comete_a5_anno$Sample == rownames(smallrna_tcga_comete_a5_dge$samples)) != 
   length(rownames(smallrna_tcga_comete_a5_dge$samples)))
{
  stop("Sanity check for annotation ordering failed - combine_tcga_comete_a5_smallrna_seq")  
}

design <- model.matrix(~0 + Genotype, data = smallrna_tcga_comete_a5_anno)
colnames(design) <- gsub("Genotype| ","", colnames(design))
rownames(design) <- rownames(smallrna_tcga_comete_a5_dge$samples)

# Remove the library specific batch effects
smallrna_tcga_comete_a5_lcpm.batch_removed <- removeBatchEffect(smallrna_tcga_comete_a5_lcpm, 
                                                                batch=smallrna_tcga_comete_a5_anno$Dataset, 
                                                                batch2 = smallrna_tcga_comete_a5_anno$Sex, 
                                                                design = design)

message("Loaded objects into global environment: 
          smallrna_tcga_comete_a5_dge (pre-batch correction feature counts as a DGEList object),
          smallrna_tcga_comete_a5_lcpm.batch_removed (batch corrected log2 cpm values)
          smallrna_tcga_comete_a5_anno (sample annotation),")

############
# Clean up #
############

rm(smallrna_tcga_comete_a5_counts)
rm(smallrna_tcga_comete_a5_counts_mat)
rm(smallrna_tcga_comete_a5_lcpm)
rm(a5_anno_zethoven_format)
setwd(old_wd)
rm(old_wd)