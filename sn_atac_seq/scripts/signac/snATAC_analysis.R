#TODO: rerun the analysis for the hg38 aligned data

# TODO: 1. fix merging of the updated clinical data in the seurat object metadata, 
#    including TERT and ATRX mutation Status these are all <NA> at the moment 

# TODO: 4. get the TERT vs Normal differential peaks and intersect these with TERT v ATRX gene list

setwd("/media/gadi/projects/ppgl/a5/")

rm(list=ls())

library(Seurat)
library(Signac)
library(rtracklayer)
library(tidyverse)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)

source("./sample_annotation/scripts/data_loaders/blake_a5_colour_palettes.r")
set.seed(1234)

# library(clustree)
# library(JASPAR2020)
# library(TFBSTools)
# library(motifmatchr)

#----
# read in the aggregated data and create a seurat object 
#----

counts <- Read10X_h5(filename = "./sn_atac_seq/analysis/cellranger/aggregated_hg38/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "./sn_atac_seq/analysis/cellranger/aggregated_hg38/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)
fragments <- "/sn_atac_seq/analysis/cellranger/aggregated_hg38/outs/fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = "hg38",
  fragments = fragments,
  min.cells = 10,
  min.features = 200
)
atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "cellranger",
  meta.data = metadata
)

#----
# Organise gene  
#----

# extract gene annotations from EnsDb
# get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
Annotation(atac) <- annotation

#----
# label barcodes according to sample of origin 
#----

atac_samples <- c("1" = "E123-1",
                  "2" = "E140-1",
                  "3" = "E146-1",
                  "4" = "E156-1",
                  "5" = "E166-1",
                  "6" = "E171-1",
                  "7" = "E188-1",
                  "8" = "E197-1",
                  "9" = "E200-1",
                  "10" = "E201-1",
                  "11" = "E225-1",
                  "12" = "NAM018")

#Aggregation tables contain the original sample identities (number appended to each cell barcode corresponde to the row number of the aggregation table)
#located in Data/cellranger_aggregation_tables
#aggregation_table <- read.csv("Data/cellranger_aggregation_tables/all_samples_recalled_peaks.csv")
#transfer original sample ID info into the metadata
#extract the original Identity from the Barcode name and store in metadata
Barcodes <- colnames(atac)
atac[["sequencing_library"]] <- substr(Barcodes, start = 18, stop=19) #number corresponds to library in Cell_Numbers
#transfer original sample ID info into the metadata
atac[["sampleID"]] <- atac$sequencing_library %>%
  recode(!!!atac_samples)
Idents(atac) <- atac$sampleID

#----
# Add the clinical data to the seurat object metadata 
#----

# read in the clinical data 
clindata <- read_csv("Data/A5_full_clinical_and_genomic_version_2.csv") %>% 
  dplyr::filter(!is.na(`A5 ID`))
# filter to retain only the snATAC-seq samples
clindata <- clindata %>%
  dplyr::filter(`A5 ID` %in% atac_samples) %>%
  mutate(sampleID = `A5 ID`)

# add the clinical data to the metadata
atac_md <- atac[[]]
atac_md <- atac_md %>%
  rownames_to_column(var = "barcode") %>%
  left_join(clindata, by = "sampleID")
rownames(atac_md) <- atac_md$barcode
atac@meta.data <- atac_md

#----
# QC metric computation 
#----

# compute nucleosome signal score per cell
atac <- NucleosomeSignal(object = atac)
# compute TSS enrichment score per cell
atac <- TSSEnrichment(atac, fast = F)
# add blacklist ratio and fraction of reads in peaks
atac$pct_reads_in_peaks <- atac$peak_region_fragments / atac$passed_filters * 100
# TODO: remove blacklist ratio as this is no longer calculated by cellranger-atac v2
atac$blacklist_ratio <- atac$blacklist_region_fragments / atac$passed_filters 
# blacklist_ratio calculation was changed from peak_region_fragments to passed_filters as I want to avoid bias towards the more common cell types
# TSS Plot
atac$high.tss <- ifelse(atac$TSS.enrichment > 2, "High", "Low")
TSSPlot(atac, group.by = "high.tss") + NoLegend()
# Nucleosome plot 
atac$nucleosome_group <- ifelse(atac$nucleosome_signal > 4, "NS > 4", "NS < 4")
FragmentHistogram(object = atac, group.by = "nucleosome_group")
# QC plot
VlnPlot(
  object = atac,
  features = c("pct_reads_in_peaks", "peak_region_fragments","passed_filters",
               "TSS.enrichment", "blacklist_ratio", "nucleosome_signal"),
  pt.size = 0,
  ncol = 3,
  group.by = "sampleID")

# vln plot the qc cutoffs 
vln_1 <- VlnPlot(
  object = atac,
  group.by = "sampleID", 
  features = "passed_filters", 
  pt.size = 0) +
  geom_hline(yintercept=5000, linetype = "dashed", colour = "red") +
  NoLegend() +
  labs(title = "fragments", x = NULL)+ 
  scale_fill_manual(values = sample_colours)

# NOTE: Blacklist ratio is not supported by cellranger v2.0
# there is a bug in cellranger atac aggr that is making blacklist_region_fragments the same as TSS_fragments 

vln_2 <- VlnPlot(
  object = atac,
  group.by = "sampleID", 
  features = "blacklist_ratio", 
  pt.size = 0) +
  geom_hline(yintercept=0.005, linetype = "dashed", colour = "red") +
  NoLegend() +
  labs(title = "blacklist ratio", x = NULL)+ 
  scale_fill_manual(values = sample_colours)

vln_3 <- VlnPlot(
  object = atac,
  group.by = "sampleID", 
  features = "TSS.enrichment", 
  pt.size = 0) +
  geom_hline(yintercept=2, linetype = "dashed", colour = "red") +
  NoLegend() +
  labs(title = "TSS enrichment score", x = NULL)+ 
  scale_fill_manual(values = sample_colours)
# individual fragment histograms
Frag_hist <- FragmentHistogram(atac) + 
  scale_fill_manual(values = sample_colours)

# Make a QC plot for fragment length of all samples on the same plot 
# using geom_smooth

reads <- Signac:::MultiGetReadsInRegion(atac, region = "chr1-1-2000000")

# add sample names to the fragments dataframe
fragments_df <- atac@meta.data %>%
  dplyr::select(sampleID, sequencing_library) %>% 
  rownames_to_column(var = "cell") %>% 
  left_join(reads, by = "cell")

fragment_dist <- ggplot(fragments_df) +
  geom_density(aes(x = length, colour = sampleID)) +
  theme_classic() + 
  scale_color_manual(values = sample_colours)

# Nice QC plot:
vln_1 +
  # vln_2 +
  vln_3 + # remove this because blacklist ratio is incorrect from cellranger bug
  fragment_dist + plot_layout(ncol = 2, guides="collect") &
  theme(text = element_text(size=10),
        axis.text = element_text(size=10))

#----
# Filter out poor quality cells 
#----

# quantile(atac$TSS.enrichment, 0.05) 
# quantile(atac$blacklist_ratio, 0.99) 
# quantile(atac$passed_filters, 0.05) 
# quantile(atac$passed_filters, 0.95)\

table(atac$passed_filters > 5000 & atac$TSS.enrichment > 2) # passed filters is the fragment number
# 3274 low quality cells are removed 

atac$remove <-  if_else(atac$passed_filters > 5000 &
                          atac$TSS.enrichment > 2,
                        "keep","remove")
# compare the cells i'm keeping vs removing
FragmentHistogram(object = atac, group.by = "remove")+
  TSSPlot(atac, group.by = "remove") + NoLegend()

table(atac$remove)
# removing 3274 cells and keeping 21995 cells

# saveRDS(atac, "Data/snATAC-seq/temp/snATAC_hg38_raw.RDS")

# I want to filter based on total reads (passed_filters) and TSS enrichment
# avoiding filtering based on the cellranger peaks because I'm calling peaks with MACS2 later on
# reads within peaks could bias the counts of under-represented cell types

atac <- subset(
  atac,
  remove == "keep")

# saveRDS(atac, "Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")

#----
# Signac pipeline
#----

# atac <- readRDS("Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")
# not absolutely necessary to do this pipe for both cellranger and MACS2 assays, can remove this step later
# will keep it here for now so I can check the umaps to make sure everything ok 

atac <- FindTopFeatures(atac, min.cutoff = 10)
# Latent semantic indexing
atac <- RunTFIDF(atac)
atac <- RunSVD(atac, reduction.key = "LSIcellranger_", reduction.name = "lsi.cellranger")
DepthCor(atac, reduction = 'lsi.cellranger')

# the first component is correlated to sequencing depth, therefore it is removed from downstream analysis 
atac <- RunUMAP(atac, reduction = "lsi.cellranger", dims = 2:30, reduction.key = "UMAPcellranger_",reduction.name = "umap.cellranger")
atac <- FindNeighbors(atac, reduction = "lsi.cellranger", dims = 2:50)
atac <- FindClusters(atac, algorithm = 3, resolution = 1.5)
#saveRDS(atac, "Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")

# which resolution to use for clustering?
# use clustree?

#----
# generate gene-activity matrix 
#----

#compute gene activity
gene.activities <- GeneActivity(atac)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
atac[["RNA"]] <- CreateAssayObject(counts = gene.activities)
atac <- NormalizeData(
  object = atac,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = median(atac$nCount_RNA)
)
DefaultAssay(atac) <- "RNA"
# saveRDS(atac, "Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")

#----
# annotate ATAC nuclei with snRNA-seq reference 
#----

#atac <- readRDS("Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")
DefaultAssay(atac) <- "RNA" # make sure correct assay is set

rna <- readRDS("Data/snRNA-seq/processed/snRNA_cs_regressed_annotated.RDS")
Idents(rna) <- rna$cell_type # make sure correct cell identities are set
# perform label transfer 
transfer.anchors <- FindTransferAnchors(
  reference = rna,
  query = atac,
  reduction = "cca"
)
predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = rna$cell_type,
  weight.reduction = atac[["lsi.cellranger"]],
  dims = 2:30
)
atac <- AddMetaData(object = atac, metadata = predicted.labels)
# create a metadata column for describing individual tumour samples and predicted cell types 
atac$neoplastic_cell_sample <- case_when(
  atac$sampleID == "NAM018" ~ "Normal",
  atac$predicted.id == "Tumor"~ atac$sampleID,
  TRUE ~ "Normal")

# make a column for the normal cell type or tumor sample ID
# this annotation is used for peak calling, so that it is performed
# individually for each cell type and tumor sample

atac$cluster <- case_when(
  atac$sampleID == "NAM018" & atac$predicted.id == "Tumor" ~ "Chromaffin cells", # this corrects the 1 normal cell was classified as tumor
  atac$sampleID == "NAM018" ~ atac$predicted.id,
  atac$predicted.id == "Tumor"~ atac$sampleID,
  TRUE ~ atac$predicted.id)

table(atac$neoplastic_cell_sample)
Idents(atac) <- atac$predicted.id

DimPlot(atac)
DimPlot(atac, group.by="cluster")
DimPlot(atac,group.by = "predicted.id")
# saveRDS(atac, "Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")

rm(rna)

#----
# UMAP plots 
#----

# atac <- readRDS("Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")

# plot samples
atac_umap_sample <- DimPlot(atac, group.by = "sampleID", pt.size = 0.5, label = T) +
  ggtitle("ATAC sample") + NoLegend()
# plot clusters
atac_umap_cluster <- DimPlot(atac, group.by = "seurat_clusters", pt.size = 0.5, label = T) +
  ggtitle("ATAC cluster") + NoLegend()
# plot cell types 
atac_umap_predictedid <- DimPlot(atac, group.by = "predicted.id", pt.size = 0.5) +
  ggtitle("ATAC predicted cell type")
atac_umap_celltype <- DimPlot(atac, group.by = "neoplastic_cell_sample", cols = c(sample_colours, Normal="grey"), pt.size = 0.5, label = F) +
  ggtitle("Tumour samples (neoplastic cells)")
atac_umap_sample + atac_umap_cluster + atac_umap_celltype

FeaturePlot(atac, features = c("prediction.score.max")) + scale_color_viridis() + atac_umap_celltype
VlnPlot(atac, features = c("prediction.score.max"), group.by = "seurat_clusters")
hist(atac$prediction.score.max)
abline(v = 0.5, col = "red")
table(atac$prediction.score.max > 0.5)
# filter out the low-scoring cells 
# filtered out 68 cells 
atac <- subset(atac,
               subset = prediction.score.max > 0.5) 
# saveRDS(atac, "Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")

#----
# Gene Activity plots 
#----

# atac <- readRDS("Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")
DefaultAssay(atac) <- "RNA"
# Marker gene expression DotPlot
fibroblast_markers <- c("FAP", "PDGFRB", "PDGFRA", "ACTA2", "COL1A1")
mono_macro_markers <-c("MSR1", "CD163", "CCL3")
chromaffin_markers <- c("PNMT", "TH", "DBH","CHGA", "CHGB")
adrenocortical_markers <- c("STAR", "CYP11B1", "CYP11A1")
endothelial_markers <- c("EPAS1", "FLT1")
sustentacular_markers <- c("CDH19", "SOX10", "S100B", "VIM")
lymphocyte_markers <- c("CD2", "CD3E" , "MS4A1")
# reorder to be more aesthetic
marker_genes <- c(
  adrenocortical_markers,
  chromaffin_markers,
  endothelial_markers,
  fibroblast_markers,
  lymphocyte_markers,
  mono_macro_markers,
  sustentacular_markers)
# make a dotplot for marker gene activity in the atac cell types  
atac_dp_celltype <- DotPlot(
  atac,
  features = rev(marker_genes),
  group.by = "predicted.id",
  cols = c("grey", "red")) +
  scale_colour_distiller(palette = "RdYlBu")+
  ggtitle("snATAC-seq gene activity") + theme(axis.text.x = element_text(angle = 90, hjust=0.9, vjust=0.5))
atac_dp_celltype
# make a dotplot for marker gene activity in the atac clusters
atac_dp_cluster <- DotPlot(
  atac,
  features = rev(marker_genes),
  group.by = "seurat_clusters",
  cols = c("grey", "red"))+
  scale_colour_distiller(palette = "RdYlBu")+
  ggtitle("snATAC-seq gene activity") + theme(axis.text.x = element_text(angle = 90, hjust=0.9, vjust=0.5))

atac_dp_cluster
FeaturePlot(atac, assay = "RNA", features = c("TERT"))

#----
# Peak calling
#----

# TODO: check and make sure that the clusters that are being used for peak calling are good

rm(list = ls())
gc()

atac <- readRDS("Data/snATAC-seq/temp/snATAC_hg38_filtered.RDS")
DefaultAssay(atac) <- "cellranger"




################## repeated from above ################## 

atac$cluster <- case_when(
  atac$sampleID == "NAM018" & atac$predicted.id == "Tumor" ~ "Chromaffin cells", # this corrects the 1 normal cell was classified as tumor
  atac$sampleID == "NAM018" ~ atac$predicted.id,
  atac$predicted.id == "Tumor"~ atac$sampleID,
  TRUE ~ atac$predicted.id)

FeaturePlot(atac, features = "TERT", order = TRUE)
DimPlot(atac, group.by = "cluster")

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"
genome(annotation) <- "hg38"
fragments <- "Data/snATAC-seq/aggregated/hg38/A5_snATAC_AGG_hg38/outs/fragments.tsv.gz"

#########################################################

peaks <- CallPeaks(
  object = atac,
  group.by = "cluster",
  cleanup = FALSE,
  macs2.path = "/data/gpfs/projects/punim0648/Projects/Blake/conda/miniconda3/envs/macs2env/bin/macs2"
)

export.bed(peaks, "Data/snATAC-seq/temp/macs2_peaks.bed")
#remove unmapped scaffold peaks
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
# remove peaks in blacklist regions
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38, invert = TRUE)
write_csv(data.frame(peaks), "Data/snATAC-seq/temp/macs2_peaks.bed")
################################################################
# check that TERT promoter peak was called as a sanity check
closest_genes <- ClosestFeature(object = atac, regions = peaks)
closest_genes %>%
  dplyr::filter(gene_name == "TERT")

TERT_peak <- StringToGRanges(regions = "chr5-1294927-1295454",
                             sep = c("-", "-"))
TERT_peaks <- subsetByOverlaps(x = peaks, ranges = TERT_peak) 
TERT_peak_string <- "chr5-1294927-1295454"
# TERT peak was called in E156-1, E200-1, E188-1 
################################################################

# I want to remove peaks exclusive to the lower-quality samples
peaks <- peaks[peaks$peak_called_in != "E146_1"]
peaks <- peaks[peaks$peak_called_in != "E225-1"]
peaks <- peaks[peaks$peak_called_in != "E171-1"]
peaks <- peaks[peaks$peak_called_in != "E140-1"]

# count reads in macs2 peaks 
peakcounts <- FeatureMatrix(
  fragments = Fragments(atac),
  features = peaks,
  cells = colnames(atac))

# add a the feature matrix as a new assay  
atac[["MACS2"]] <- CreateChromatinAssay(
  counts = peakcounts,
  min.cells = 10,
  genome = "hg38",
  fragments = fragments,
  annotation = annotation
)
# saveRDS(atac, "Data/snATAC-seq/processed/snATAC_hg38_filtered_macs2.RDS")

atac <- readRDS("Data/snATAC-seq/processed/snATAC_hg38_filtered_macs2.RDS")

DefaultAssay(atac) <- "MACS2"
atac <- FindTopFeatures(atac, min.cutoff = 10)
atac <- RunTFIDF(atac)
atac <- RunSVD(atac, reduction.key = "LSIMACS2_", reduction.name = "lsi.MACS2")
atac <- RunUMAP(atac, reduction = "lsi.MACS2", dims = 2:20, reduction.key = "UMAPMACS2_", reduction.name = "umap.MACS2")
atac <- FindNeighbors(atac, reduction = "lsi.MACS2", dims = 2:20)
atac <- FindClusters(atac, algorithm = 3)

# saveRDS(atac, "Data/snATAC-seq/processed/snATAC_hg38_filtered_macs2.RDS")

# ----
# UMAP plots
# ----

atac <- readRDS("Data/snATAC-seq/processed/snATAC_hg38_filtered_macs2.RDS")

DimPlot(atac, group.by='predicted.id')
DimPlot(atac, group.by = "predicted.id", reduction = "umap.cellranger") +
  DimPlot(atac, group.by = "predicted.id", reduction = "umap.MACS2") +
  DimPlot(atac, group.by = "sampleID", reduction = "umap.MACS2")

atac$cell_type_abbreviated <- recode(atac$predicted.id,
                                     "Chromaffin cells"= "CCs",
                                     "Myeloid cells"= "MCs",
                                     "Lymphocytes"= "LCs",
                                     "Adrenocortical cells"= "ACs",
                                     "Fibroblasts"= "FCs",
                                     "Endothelial cells"= "ECs")

DimPlot(atac, group.by = "cell_type_abbreviated", reduction = "umap.MACS2")
DimPlot(atac, group.by = "TERT_or_ATRX", reduction = "umap.MACS2") + 
  DimPlot(atac, group.by = "TERT_or_ATRX", reduction = "umap.MACS2")


  
#----
# Differential Accessibility Analysis
#----

atac$TERT_or_ATRX <- case_when(
  atac$sampleID == "NAM018" | atac$predicted.id != "Tumor" ~ "Normal",
  atac$sampleID %in% c("E123-1","E146-1","E156-1","E188-1","E200-1","E201-1") ~ "TERT",
  atac$sampleID %in% c("E140-1", "E166-1", "E197-1", "E225-1") ~ "ATRX",
  atac$sampleID == "E171-1" ~ "Unknown met driver")

atac$TERT_or_ATRX_or_lowquality <- case_when(
  atac$sampleID == "NAM018" & atac$cell_type_abbreviated == "CCs" ~ "Normal CCs",
  atac$sampleID == "NAM018" | atac$predicted.id != "Tumor" ~ "Normal other",
  atac$sampleID %in% c("E156-1","E188-1","E200-1","E201-1") ~ "TERT",
  atac$sampleID %in% c("E166-1", "E197-1") ~ "ATRX",
  atac$sampleID %in% c("E123-1", "E146-1", "E140-1", "E225-1","E171-1") ~ "Low quality"
)

Idents(atac) <- atac$TERT_or_ATRX

# TERT vs ATRX mutants 
DefaultAssay(atac) <- "MACS2"

# compare TERT mutant tumors to ATRX mutant tumors (neoplastic cells only)
da_peaks <- FindMarkers(
  object = atac,
  ident.1 = "TERT",
  ident.2 = "ATRX",
  min.pct = 0,
  logfc.threshold = 0.1,
  test.use = "LR",
  latent.vars = "peak_region_fragments"
)
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
  left_join(y = da_peaks) %>%
  dplyr::filter(p_val_adj < 0.05)

# write.csv(da_peaks, "results/seurat_da_peaks/tert_vs_atrx.csv")

da_peaks %>% dplyr::filter(closest_gene == "TWIST1" & p_val_adj < 0.05)
da_peaks %>% dplyr::filter(closest_gene == "SNAI2" & p_val_adj < 0.05)
da_peaks %>% dplyr::filter(closest_gene == "TERT" & p_val_adj < 0.05)

tert_peaks <- da_peaks %>%
  dplyr::filter(
    p_val_adj < 0.05 & 
      avg_log2FC > 0
  )

tert_promoter_peaks <- tert_peaks %>%
  dplyr::filter(distance < 2000 &
                  p_val_adj < 0.01)

dim(tert_promoter_peaks)  

atrx_peaks <-  da_peaks %>%
  dplyr::filter(
    p_val_adj < 0.05 & 
      avg_log2FC < 0
  )

# TODO: repeat DAA for TERT vs Normal Chromaffin Cells

#----
# Coverage plots 
#----

DefaultAssay(atac) <- "MACS2"

Idents(atac) <- atac$cluster

# set the order order of coverage tracks
levels(atac) <- c(
  "E123-1","E146-1","E156-1","E188-1","E200-1", "E201-1", # TERT
  
  "E140-1", "E166-1", "E197-1", "E225-1", # ATRX
  
  "E171-1",
  
  "Chromaffin cells",
  "Adrenocortical cells",
  "SCLCs",
  "Fibroblasts",
  "Endothelial cells",
  "Lymphocytes",
  "Myeloid cells")

da_peaks <- read.csv("results/seurat_da_peaks/tert_vs_atrx.csv")
tert_da_peaks_ranges <- da_peaks %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  pull(peak) %>% 
  StringToGRanges()

# TERT Coverage Plot
CoveragePlot(
  object = atac,
  region = "chr5-1251335-1301433",
  ranges = tert_da_peaks_ranges
)

# SNAI2 Coverage Plot
CoveragePlot(
  object = atac,
  region = "chr8-48917048-48922886",
  extend.upstream = 30000,
  extend.downstream = 10000,
  ranges = tert_da_peaks_ranges,
  ranges.title = "TERT DARs"
)

#----
# Signac Motif analysis 
#----
# TODO: rerun below code on the new hg38 data 

# Overrepresented motifs in the DARs  
atac <- readRDS("Data/snATAC-seq/processed/snATAC_AGG_01_processed.RDS")
DefaultAssay(atac) <- "MACS2"
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE)
)
# add motif information
atac <- AddMotifs(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)
da_peaks <- read.csv("results/seurat_da_peaks/tert_vs_atrx.csv")
top_da_peaks <- da_peaks %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(avg_log2FC < -1 | avg_log2FC > 1) %>% 
  pull(X)
top_da_peaks

# find peaks open in TERT or ATRX cells 
open.peaks <- AccessiblePeaks(atac, idents = c(tert_mutant, atrx_mutant))
# match the overall GC content in the peak set
meta.feature <- GetAssayData(atac, assay = "MACS2", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top_da_peaks, ],
  n = 50000
)
# Find enriched motifs within the TERT deferentially accessible peaks 
enriched.motifs <- FindMotifs(
  object = atac,
  features = top_da_peaks,
  background = peaks.matched
)
head(enriched.motifs, 50)
MotifPlot(
  object = atac,
  motifs = head(rownames(enriched.motifs),75)
)

enriched.motifs %>%
  dplyr::filter(pvalue < 0.05) %>% # filter significant motifs
  # arrange(desc(fold.enrichment)) %>% # order by fold enrichment 
  head(75)

#----
# ChromVAR motif analysis 
#----

# use ChromVAR to compute motif activities per cell
atac <- RunChromVAR(
  object = atac,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
#saveRDS(atac, "Data/snATAC-seq/processed/snATAC_AGG_01_processed.RDS")
atac <- readRDS("Data/snATAC-seq/processed/snATAC_AGG_01_processed.RDS")
DefaultAssay(atac) <- "chromvar"
differential.motifs <- FindMarkers(
  object = atac,
  ident.1 = tert_mutant,
  ident.2 = atrx_mutant,
  only.pos = TRUE,
  test.use = "LR",
  latent.vars = "passed_filters"
)

#write.csv(differential.motifs, "results/chromVAR_motifs/tert_vs_atrx_motifs.csv")
differential.motifs <- read_csv("results/chromVAR_motifs/tert_vs_atrx_motifs.csv")
top_motifs <- differential.motifs %>%
  arrange(desc(avg_log2FC)) %>% # arrange motifs by fold change 
  pull(X1)
# plot top motifs by fold change 
MotifPlot(
  object = atac,
  motifs = head((top_motifs), 75),
  assay = "MACS2"
)

# look at the activity of SNAI2
cell_types_dp <- DimPlot(
  atac,
  reduction = "umap.MACS2",
  group.by = "cluster",
  label = T, label.size = 3, repel = T,
  pt.size = 0.1)

snai2_motif_fp <- FeaturePlot(
  object = atac,
  features = "MA0745.2",
  min.cutoff = "q10",
  max.cutoff = "q90",
  pt.size = 0.1, label = F, order = T,
  reduction = "umap.MACS2"
) + ggtitle("SNAI2 (MA0745.2) motif activity")

# snai2 motif does not appear to have higher accessibility in the TERT-mutants as expected

snai2_activity_fp <- FeaturePlot(
  object = atac,
  features = "SNAI2",
  min.cutoff = "q10",
  max.cutoff = "q90",
  pt.size = 0.1, label = F, order = T,
  reduction = "umap.MACS2"
) + ggtitle("SNAI2 gene activity ")

cell_types_dp + snai2_motif_fp + snai2_activity_fp + plot_layout(ncol = 3)

klf4_motif_fp <- FeaturePlot(
  object = atac,
  features = "MA0039.4",
  min.cutoff = "q10",
  max.cutoff = "q90",
  pt.size = 0.1, label = F, order = T,
  reduction = "umap.MACS2"
) + ggtitle("KLF4 (MA0039.4) motif activity")

klf4_activity_fp <- FeaturePlot(
  object = atac,
  features = "KLF4",
  min.cutoff = "q10",
  max.cutoff = "q90",
  pt.size = 0.1, label = F, order = T,
  reduction = "umap.MACS2"
) + ggtitle("KLF4 gene activity")

cell_types_dp + klf4_motif_fp + klf4_activity_fp + plot_layout(ncol = 3)

VlnPlot(atac, features = head(differential.motifs$X1)) # these results seem biased towards E188-1?
VlnPlot(atac, features = head(enriched.motifs$motif))
VlnPlot(atac, features = "nCount_MACS2")

#----
# Transcription factor footprinting 
#----

atac <- readRDS("Data/snATAC-seq/processed/snATAC_AGG_01_processed.RDS")
DefaultAssay(atac) <- "MACS2"

# gather the footprinting information for sets of motifs
atac <- Footprint(
  object = atac,
  motif.name = c("SNAI2", "EBF1", "EBF3", "MAZ", "SP1", "KLF4", "TWIST1", "ZEB1", "JUNB"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# plot the footprint data for each group of cells
footprint_plot <- PlotFootprint(atac, features = c("SNAI2", "EBF1", "EBF3", "MAZ", "SP1", "KLF4", "TWIST1", "ZEB1", "JUNB"), label.idents = F)
footprint_plot + plot_layout(guides = 'collect')

junb_footprint <- PlotFootprint(atac, features = "JUNB",label.top = 20)
junb_footprint

#----
# investigate motifs
#----

# look at the <- to see which TFs are expressed 
rna <- readRDS("Data/snRNA-seq/processed/snRNA_cs_regressed_annotated_doublets_removed.RDS")
levels(rna) <- c(
  "E123_1",
  "E156_1",
  "E146_1",
  "E174", 
  "E019",
  "E166_1",
  "E197_1", 
  "E018",
  "E025",
  "PGL1", 
  "PGL3", 
  "Chromaffin_cell", 
  "Sustentacular_cell",
  "Adrenocortical_cell",
  "Endothelial_cell",
  "Fibroblast",
  "Monocyte_Macrophage",
  "B_cell",
  "T_cell")

#----
# cell types histogram 
#----

atac <- readRDS("Data/snATAC-seq/processed/snATAC_AGG_01_processed.RDS")

# make a stacked bar plot of the sample-cell types
plot.data <- atac@meta.data %>%
  rownames_to_column("barcode") %>%
  dplyr::select(sampleID, predicted.id) %>% 
  as_tibble() %>% 
  group_by(sampleID) %>% 
  mutate(count_per_sample = n()) %>% 
  ungroup() %>% 
  group_by(sampleID, predicted.id) %>% 
  mutate(count_per_type_sample = n()) %>% 
  distinct()
plot.data

cell_prop_plot <- ggplot(
  plot.data,
  aes(x = sampleID,
      y = count_per_type_sample,
      fill = predicted.id)) +
  geom_bar(position="fill", stat="identity") +
  ggtitle("Cell type proportion by sample (snATAC)") +
  ylab("Proportion of sample") +
  theme_classic()

cell_prop_plot 

#----
# Cicero co-accessible networks 
#----

levels(atac) <- rev(levels(atac))
levels(rna) <- rev(levels(rna))
DefaultAssay(atac) <- "RNA"
a <- VlnPlot(atac, features = "LMX1B", pt.size = 0.2, stack = F) + ggtitle("LMX1B Gene Activity snATAC") 
b <- VlnPlot(rna, features = "LMX1B", pt.size = 0.2, stack = F) + ggtitle("LMX1B Gene Expression snRNA")
a+b

emt_genes <- c("LOXL2", "TWIST1", "TCF3", "MMP2", "MMP1", "KRT19", "CDH2", "ZEB1", "SNAI1", "SNAI2")

FeaturePlot(rna, features = "GATA4", order = T, label = T)
VlnPlot(rna, features = "GATA4") + NoLegend()

