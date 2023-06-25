# Paradifference figures 
rm(list=ls())
library(Seurat)
library(SeuratWrappers)
library(Signac)
library(future)
library(rtracklayer)
library(tidyverse)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)
library(viridis)
library(ggsci)
library(scales)
# library(clustree)
# library(JASPAR2020)
# library(TFBSTools)
# library(motifmatchr)
set.seed(1234)

atac <- readRDS("Data/snATAC-seq/processed/snATAC_hg38_filtered_macs2.RDS")
atac <- subset(atac, subset = cluster != "Tumor")
# coverage plots 
DefaultAssay(atac) <- "MACS2"
# set the order order of coverage tracks

cell_types <- c(
  "E123-1","E146-1","E156-1","E188-1","E200-1", "E201-1", # TERT
  "E140-1", "E166-1", "E197-1", "E225-1", # ATRX
  "E171-1", # unknown met driver
  "Chromaffin cells",
  "Adrenocortical cells",
  "SCLCs",
  "Fibroblasts",
  "Endothelial cells",
  "Lymphocytes",
  "Myeloid cells")

covplot_colours <-c(pal_d3("category20")(20)[c(2,2,2,2,2,2,1,1,1,1,3)], 
                   rep("lightgrey",7))

levels(atac) <- cell_types


# TODO: get the hg38 positions for these
AFF4_breakpoint <- StringToGRanges("chr5-132292963-132292984")
AFF4_breakpoint$color <- "red"

TERT_breakpoint <- StringToGRanges("chr5-1295481-1295503")
TERT_breakpoint$color <- "red"

TERT_covplot <- CoveragePlot(atac, region = TERT_breakpoint,
                             extend.upstream = 45000, extend.downstream = 5000,
                             region.highlight = TERT_breakpoint,
                             peaks = FALSE) &
  scale_fill_manual(values = covplot_colours)

AFF4_covplot <- CoveragePlot(atac, region = AFF4_breakpoint,
                             extend.upstream = 80000, extend.downstream = 10000,
                             region.highlight = AFF4_breakpoint,
                             peaks = FALSE) &
  scale_fill_manual(values = covplot_colours)

tert_fusion_covplot <- TERT_covplot | AFF4_covplot
tert_fusion_covplot

# ggsave(tert_fusion_covplot, filename = "Figures/tert_fusion_covplot.pdf",
#        height = unit(6, "cm"),
#        width = unit(12, "cm"))

