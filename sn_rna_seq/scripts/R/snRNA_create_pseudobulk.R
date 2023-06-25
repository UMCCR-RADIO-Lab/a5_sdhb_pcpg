# functions for pseudobulking a Seurat object 
# REFERENCES: https://github.com/hbctraining/scRNA-seq_online/blob/939b70563e57dca000681c0b8e3e75edd2053c5a/lessons/pseudobulk_DESeq2_scrnaseq.md

rm(list=ls())

library(Seurat)
library(Matrix.utils)
library(tidyverse)

make_pseudobulk_matrix <- function(obj, grp){
  #' pseudobulking is basically calculating the sum of counts within each group
  #' @param obj a seurat object
  #' @param grp character vector containing the names of metadata columns to group the cells by for aggregation
  #' @returns matrix with genes as rows and groups as columns 
  
  # get the specified groups from the metadata 
  grp_df <- rna@meta.data[,grp]
  # Aggregate across cluster-sample groups
  pb <- aggregate.Matrix(t(rna@assays$RNA@counts), 
                         groupings = grp_df, fun = "sum") 
  pb <- data.matrix(t(pb))
}

rna <- readRDS("Data/snRNA-seq/processed/snRNA_cs_regressed_annotated.RDS")

# the groups that I want to aggregate by are sample and cell type
pb <- make_pseudobulk_matrix(rna, grp = c('cell_type', 'A5.ID'))

# write.csv(pb, "Data/snRNA-seq/processed/snRNA_pseudobulk.csv")
