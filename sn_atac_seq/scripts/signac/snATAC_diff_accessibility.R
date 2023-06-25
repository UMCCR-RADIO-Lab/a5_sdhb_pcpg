
library(Seurat)
library(SeuratWrappers)
library(Signac)
library(future)
library(rtracklayer)
library(tidyverse)

set.seed(1234)

# setup parrallelisation
plan("multisession", workers = 8) # use 8 cpu cores
options(future.globals.maxSize = 200 * 1024 ^ 3) # 200GB RAM
plan()

setwd("/data/gpfs/projects/punim0648/Projects/Blake/A5_R_project")

atac <- readRDS("Data/snATAC-seq/processed/snATAC_AGG_02_processed.RDS")

# compare TERT vs ATRX mutants 
DefaultAssay(atac) <- "MACS2"
da_peaks <- FindMarkers(
  object = atac,
  ident.1 = "TERT",
  ident.2 = "ATRX",
  group.by = "metastasis_driver",
  min.pct = 0,
  logfc.threshold = 0.25,
  test.use = "LR",
  latent.vars = "peak_region_fragments")
da_peaks <- da_peaks %>%
  rownames_to_column(var = "peak")
# find gene closest to each peak
closest_features <- ClosestFeature(atac,
                                   regions = da_peaks$peak,
                                   annotation = Annotation(atac))
# merge the da peaks and closest features into a single dataframe
da_peaks <- closest_features %>%
  rename(peak = query_region, closest_gene = gene_name) %>%
  dplyr::select(peak, closest_gene, distance, gene_biotype) %>% 
  left_join(y = da_peaks, by = "peak")

write.csv(da_peaks, "results/seurat_da_peaks/tert_vs_atrx.csv")