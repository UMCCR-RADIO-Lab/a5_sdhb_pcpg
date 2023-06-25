#----
# Chromaffin Cell Differentiation marker expression
#----

rm(list=ls())
library(Seurat)
library(ComplexHeatmap)
library(tidyverse)
library(circlize)
library(ggsci)
library(scales)
library(pals)
set.seed(1234)

seurat_obj <- readRDS("Data/snRNA-seq/processed/snRNA_cs_regressed_annotated_doublets_removed.RDS")
seurat_obj <- subset(seurat_obj, subset = cell_type_integration == "Tumour" |
                       cell_type_integration == "Chromaffin_cell", downsample = 800)
metadata <- seurat_obj@meta.data
metadata$genotype <- recode(metadata$orig.ident, 
                            "E156_1" = "TERT",
                              "E123_1" = "TERT",
                              "E174"= "TERT_NF1",
                              "E146_1" = "TERT_subclone",
                              "E197_1" = "ATRX",
                              "E166_1" = "ATRX",
                              "E019" = "ATRX",
                              "E018" = "ATRX_RET",
                              "PGL1" = "Unknown_secondary_driver",
                              "PGL3" = "Unknown_secondary_driver",
                              "NAM025" = "Normal",
                              "NAM021" = "Normal")
metadata_ordered <- metadata[order(metadata$genotype, metadata$cell_type),]

cc_gene_table <- read_csv("Data/chromaffin_differentiation_genes.csv")
cc_genes <- cc_gene_table$Gene # chromaffin cell differentiation genelist

# get scaled counts from the seurat object
counts <- as.matrix(seurat_obj[["RNA"]]@data)
# subset genes of interest and order counts by the metadata
counts <- counts[cc_genes, rownames(metadata_ordered)]
counts_scaled <-  t(scale(t(counts)))
dim(counts_scaled)

# make annotation colourmap 
cell_types <- as.character(unique(metadata$cell_type))
cell_type_colours <- pal_d3("category20")(length(cell_types))
names(cell_type_colours) <- cell_types

genotypes <- as.character(unique(metadata$genotype))
genotype_colours <- pal_igv()(length(genotypes))
names(genotype_colours) <- genotypes

annotation_colors <- list("cell_type" = cell_type_colours,
                          "genotype" = genotype_colours)
anno_df <- metadata_ordered %>%
  dplyr::select(cell_type, genotype)

# heat map annotation 
ha <-  HeatmapAnnotation(df = anno_df,
                         show_annotation_name = TRUE,
                         col = annotation_colors)


col_fun <-  colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

Heatmap(counts_scaled,
        col=col_fun,
        cluster_rows = T,
        cluster_columns = F,
        column_order = NULL,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_row_names = TRUE,
        show_column_names = FALSE,
        use_raster = TRUE,
        bottom_annotation = NULL,
        top_annotation = ha)

