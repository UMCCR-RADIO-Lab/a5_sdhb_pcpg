# Projection of neoplastic cells on to the Jansky foetal diffusion map
rm(list=ls())

library(Seurat)
library(tidyverse)
library(future)
library(patchwork)
library(ggsci)
library(scales)
source("Scripts/colour_palettes/A5_colour_palettes.R")

set.seed(1234)

# ----
# make neo-only seurat objects for projection
# ----

rna <- readRDS("Data/snRNA-seq/processed/snRNA_cs_regressed_annotated.RDS")

#TODO: remove unneccesary metadata to make the files smaller?

rna@meta.data <- rna@meta.data %>%
  select(barcode, cell_type, A5.ID, nCount_RNA, nFeature_RNA, percent.mt, metastasis_driver)

# Filter out the non-neoplastic cells
rna <- subset(rna, subset = cell_type %in% c("Tumor", "Chromaffin cells"))

# split into individual samples
obj.list <- SplitObject(rna, split.by = "A5.ID")

sample_names <- names(obj.list)

for (i in 1:length(obj.list)){
  sample_name <- sample_names[i]
  sample_obj <- obj.list[[i]]
  # get the counts
  print(sample_name)
  # save each counts matrix as an RDS
  neo_counts <- data.matrix(sample_obj@assays$RNA@counts)
  saveRDS(neo_counts, file.path("Data/snRNA-seq/processed/neo_only_counts", paste0(sample_name, "-neo-only-counts.RDS")))
  # save individual seurat objects as RDS
  saveRDS(sample_obj, file.path("Data/snRNA-seq/processed/neo_only_seuratobj", paste0(sample_name, "-neo-only-seuratobj.RDS")))
  }

