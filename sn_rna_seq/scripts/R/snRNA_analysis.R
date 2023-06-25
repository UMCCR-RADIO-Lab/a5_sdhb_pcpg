# script for pre-processing A5 scRNA-seq
rm(list=ls())

library(Seurat)
library(tidyverse)
library(future)
library(patchwork)
library(ggsci)
library(scales)
library(ggthemes)
source("Scripts/colour_palettes/A5_colour_palettes.R")

set.seed(1234)

#####
# PREPROCESSING AND QC
#####

#----
# load in and merge datasets 
#----

# input data:
# scrublet tables - Data/snRNA-seq/scrublet/
# snRNA-seq - Data/snRNA-seq/

make_seurat <- function(sample_name, matrix_dir, scrublet_csv){
  #' read seurat object from 10x output and add scrublet predictions to metadata 
  sample.data <- Read10X(data.dir = matrix_dir)
  sample.seurat <- CreateSeuratObject(counts = sample.data, project = sample_name, min.cells = 3, min.features = 200)
  # add scrublet results to metadata
  sample.scrublet <- read.csv(scrublet_csv, sep = ",")
  rownames(sample.scrublet) <- colnames(sample.data)
  sample.seurat <- AddMetaData(sample.seurat, metadata = sample.scrublet)
  return(sample.seurat)
}

# make a list with all the sample names 
samples_all <- list.dirs(path = "Data/snRNA-seq/hg38-counts/", full.names = F, recursive = F)
# loop through the samples, load them into seurat objects 
for (i in 1:length(samples_all)){
  sample_name <- (samples_all[i])
  matrix_dir <- file.path("Data/snRNA-seq/hg38-counts/", sample_name, "/outs/filtered_feature_bc_matrix/")
  scrublet_csv <- file.path("Data/snRNA-seq/scrublet/Scrublet_results_hg38/", paste0(sample_name, "_scrublet_output_table.csv"))
  seurat_obj <- make_seurat(sample_name = sample_name,
                            matrix_dir = matrix_dir,
                            scrublet_csv = scrublet_csv)
  assign(make.names(sample_name), seurat_obj)
}

# merge all of the samples into one object
rna <- merge(
  E019,
  y = c(E123.1,
        E140.1,
        E143.1,
        E146.1,
        E156.1,
        E166.1,
        E171.1,
        E197.1,
        E225.1,
        NAM021,
        NAM025,
        NPG103,
        P018.PGL1,
        P018.PGL3),
  add.cell.ids = samples_all,
  project = "A5_single_nuclei")

#save the raw counts matrix as a csv so for scMatch cell type analysis
#counts <- GetAssayData(object = rna, slot = "counts")
#write.csv(counts, "Data/scMatch/input/rna_all_rawcounts_mat.csv", quote = T)

# clear up memory:
rm(E123.1, E140.1, E143.1, E146.1, E156.1, E166.1, E171.1, E197.1, E225.1, NAM021, NAM025, NPG103, P018.PGL1, P018.PGL3)

#----
# Calculate QC metrics
#----

# calculate % mitochondrial gene expression
rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^MT-")
# calculate log10 genes detected per umi 
rna[["log10.genes.per.umi"]] <- log10(rna$nFeature_RNA) / log10(rna$nCount_RNA)

# ----
# organise the sample-level metadata, put it in the Seurat Obj
# ----

# read in the clinical data
a5_clindata <- read_csv("Data/A5_full_clinical_and_genomic_version_2.csv")
# filter to retain only the snRNA-seq samples
a5_clindata <- a5_clindata %>%
  dplyr::filter(`A5 ID` %in% c(samples_all, "E019-1")) 

# get available metadata for the pheo-atlas paper samples
singlecell_metadata <- read_tsv("Data/metadata_sc_samples.tsv") 

# harmonize the important single cell metadata fields with A5 metadata  
singlecell_metadata <- singlecell_metadata %>%
  dplyr::rename("A5 ID" = Sample) %>% 
  filter(`A5 ID` %in% c("P018-PGL1", "P018-PGL3", "E240", "E243")) %>% 
  dplyr::select(`A5 ID`, Gender, Location, Malignancy) %>% 
  mutate(Gender = recode(Gender,"M" = "male")) %>% 
  mutate(`Metastatic state of tumour sample` = if_else(
    `A5 ID` %in% c("E240", "E243"),
    "Normal", "Primary")) %>%  # these have malignancy = "benign" so they must all be primaries 
  mutate(`Tumour metastasised` = if_else(
    `A5 ID` %in% c("E240", "E243"),
    "Normal", "No")) %>% 
  dplyr::select(-Malignancy, -Location)
  
# join the A5 and singlecell metadata
# add the clinical data to the metadata, retaining the rownames 
metadata <- rna@meta.data
metadata <- metadata %>%
  # add cell barcode column 
  rownames_to_column(var = "barcode") %>%
  # add a sample column 
  mutate(`A5 ID` = recode(orig.ident,
                         "NAM021" = "E240-1",
                         "NAM025" = "E243-1",
                         "E019" = "E019-1")) %>% 
  # join the clinical data to the metadata 
  left_join(a5_clindata, by = "A5 ID") %>% 
  # TODO: fix the merging issue duplicating columns with same name
  left_join(singlecell_metadata, by = "A5 ID")%>% 
  # add a column describing the library preparation chemistry that was used 
  mutate(chemistry = if_else(`A5 ID` %in% c("E140-1",
                                            "E143-1",
                                            "E171-1",
                                            "E225-1"),
                             "3prime", "5prime"))

# make the metadata colnames compatible with seurat commands
compatiblenames <- setNames(colnames(metadata), make.names(colnames(metadata)))
metadata <- metadata %>%
  rename(!!!compatiblenames)

rownames(metadata) <- metadata$barcode
rna@meta.data <- metadata

#----
# read in scMatch annotations 
#----

# TODO: re-run this on the new data?
# #read in the top scMatch annotations for each cell, add annotations to metadata
# scMatch_annos <- read.csv("Data/snRNA-seq/scMatch/results/annotation_result_keep_all_genes/human_Spearman_top_ann.csv")
# #make barcodes the same
# scMatch_annos$cell <- gsub('.', '-',scMatch_annos$cell, 
#                            fixed = TRUE) 
# #check that cell order in the scMatch is still correct
# table(colnames(rna)==(scMatch_annos$cell))
# names(scMatch_annos$cell.type) <- (scMatch_annos$cell)
# rna <- AddMetaData(rna, metadata = scMatch_annos$cell.type, col.name = 'scMatch_predicted_cell_type')

#----
# Make QC plots 
#----

# Visualize the number of cell counts per sample
ncells_bar <- metadata %>% 
  ggplot(aes(x=A5.ID, fill=A5.ID)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells") + 
  scale_fill_manual(values = sample_colours)
# ggsave(ncells_bar, "Figures/QC_plots/ncells_per_sample_barplot.pdf")

# Visualize the number UMIs/transcripts per cell
nCount_density <- metadata %>% 
  ggplot(aes(color=A5.ID, x=nCount_RNA, fill= A5.ID)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 750, linetype = "dashed", colour = "red") + 
  scale_colour_manual(values = sample_colours) + 
  scale_fill_manual(values = sample_colours) +
  theme(legend.position = "none")

# Visualize the distribution of genes detected per cell via histogram
nfeatures_density <- metadata %>% 
  ggplot(aes(color=A5.ID, x=nFeature_RNA, fill= A5.ID)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() +
  ylab("Cell density") +
  geom_vline(xintercept = 500, linetype = "dashed", colour = "red") + 
  scale_colour_manual(values = sample_colours) + 
  scale_fill_manual(values = sample_colours) +
  theme(legend.position = "none")

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
complexity_density <- metadata %>%
  ggplot(aes(x=log10.genes.per.umi, color = A5.ID, fill=A5.ID)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  scale_colour_manual(values = sample_colours) + 
  scale_fill_manual(values = sample_colours)

nfeatures_density + nCount_density + complexity_density 

# visualise percentage mitochondrial, number of counts and number of genes as violin plots 
nfeature_vln <- VlnPlot(rna, features = "nFeature_RNA", group.by = "A5.ID", pt.size = 0) +
  geom_hline(yintercept = c(500, 6000), linetype = "dashed", colour = "red") +
  scale_fill_manual(values = sample_colours) + 
  labs(title = "nFeature_RNA", x = NULL) 
ncount_vln <- VlnPlot(rna, features = "nCount_RNA", group.by = "A5.ID", pt.size = 0) +
  geom_hline(yintercept = 750, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = sample_colours) + 
  lims(y = c(0 , 25000)) +
  labs(title = "nCount_RNA", x = NULL) 
mt_vln <- VlnPlot(rna, features = "percent.mt", group.by = "A5.ID", pt.size = 0) +
  geom_hline(yintercept = 10, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = sample_colours) + 
  lims(y = c(0 , 40)) +
  labs(title = "percent.mt", x = NULL)

nfeature_vln + ncount_vln + mt_vln +
  plot_layout(ncol = 3) &
  theme(legend.position = "none")

featurescatter1 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'predicted_doublet')
featurescatter2 <- FeatureScatter(rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'predicted_doublet')

featurescatter1 + featurescatter2

# investigate data quality
quantile(rna$nCount_RNA, c(0.005, 0.05, 0.5, 0.95, 0.995))
table(rna$nCount_RNA < 16000)
table(rna$nCount_RNA > 1000)

quantile(rna$nFeature_RNA, c(0.005, 0.05, 0.5, 0.95, 0.995))
table(rna$nFeature_RNA < 6000)
table(rna$nFeature_RNA > 500)

quantile(rna$percent.mt, c(0.005, 0.05, 0.5, 0.95, 0.975, 0.995))
table(rna$percent.mt < 10)

#saveRDS(rna, "Data/snRNA-seq/temp/snRNA_raw.RDS")

#----
# Remove poor-quality cells 
#----

# removing all predicted doublets,
# discarding NPG103 because it has low library quality and no chromaffin cells,
# discarding E019-1 because of low library quality
rna <- readRDS("Data/snRNA-seq/temp/snRNA_raw.RDS") # 92676 cells before filtering (only counting cells with 300+ counts)
rna <- subset(
  rna,
  subset = nFeature_RNA > 500 &
    nCount_RNA > 750 &
    percent.mt < 10 & 
    A5.ID != "NPG103" &
    A5.ID != "E019-1" &
    predicted_doublet == "False")

# ____ cells remain after filtering

# saveRDS(rna, "Data/snRNA-seq/temp/snRNA_filtered.RDS")

#####
# ANALYSIS
#####

#----
# Regress unwanted sources of variation incl. cell cycle
#----

rna <- readRDS("Data/snRNA-seq/temp/snRNA_filtered.RDS")
# Normalisation and Scaling- required for cell cycle scoring
rna <- NormalizeData(rna)
rna <- FindVariableFeatures(rna, selection.method = "vst", nfeatures = 2000)
rna <- ScaleData(rna, features = rownames(rna))

# make list of genes for s and g2 phases
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#assign phase based on phase-specific gene expression
rna <- CellCycleScoring(rna, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# saveRDS(rna, "Data/snRNA-seq/temp/snRNA_filtered.RDS")
rna <- readRDS("Data/snRNA-seq/temp/snRNA_filtered.RDS")

#To visualise the cell-cycle-related variance, run PCA using the s and g2 gene lists
#rna <- RunPCA(rna, features = c(s.genes, g2m.genes)) 
#looks like some of the cells are have phase-specific clustering in the PCA

#rna <- RunPCA(rna)
# plot2 <- DimPlot(rna, group.by = "Phase") + ggtitle("Phase - before regression") #Clearly the G2M cells are on the right side of PC1
# plot3 <- FeaturePlot(rna, features = "nCount_RNA") # high UMI cells appear to coalesce
# plot4 <- FeaturePlot(rna, features = "percent.mt")

# regress out cell cyle phase, percentage mitochrondrial rna, sequencing depth 
#Log normalise and scale: this will create a normalised and regressed dataset that I can use for downstream analysis
rna <- ScaleData(rna, vars.to.regress = c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt"), features = rownames(rna))

# Principal component analysis and visualisation
rna <- RunPCA(rna)
plot5 <- DimPlot(rna, group.by = "Phase") + ggtitle("Phase - after regression") 
plot6 <- FeaturePlot(rna, features = "nCount_RNA", order = TRUE)
plot7 <- FeaturePlot(rna, features = "percent.mt", order = TRUE)

# plot PCA before and after 
# (plot2+plot3+plot4) /
#   (plot5+plot6+plot7) 

(plot5+plot6+plot7) 

# saveRDS(rna, file = "Data/snRNA-seq/temp/snRNA_cs_regressed.RDS")

#----
# Dimensionality reduction 
#----

rna <- readRDS(file = "Data/snRNA-seq/temp/snRNA_cs_regressed.RDS")        
#PCA with variable genes, clustering and dimensionality reduction
rna$chemistry <- if_else(A5.ID %in% c("E140-1",
                                        "E143-1",
                                        "E171-1",
                                        "E225-1"),
                         "3prime", "5prime")

ElbowPlot(rna, ndims = 50)

rna <- RunUMAP(rna, dims = 1:25)
rna <- FindNeighbors(rna, dims = 1:25)
rna <- FindClusters(rna, resolution = 1) 

# plot umap by sample and cluster to see what she looks like
clusters_dimplot <- DimPlot(rna, group.by = 'seurat_clusters', label = TRUE) + NoLegend() +ggtitle('Seurat Cluster')
samples_dimplot <- DimPlot(rna, group.by = 'A5.ID', cols = sample_colours, label = T) +ggtitle('Sample')
genotypes_dimplot <- DimPlot(rna, group.by = 'TERT_or_ATRX', cols = ) +ggtitle('Genotype')
chemistry_dimplot <- DimPlot(rna, group.by = 'chemistry') + ggtitle('Chemistry')

samples_dimplot +
  clusters_dimplot +
  genotypes_dimplot +
  chemistry_dimplot +
  plot_layout(ncol = 2)

# saveRDS(rna, file = "Data/snRNA-seq/temp/snRNA_cs_regressed.RDS")

#----
# Plots 
#----

#rna <- readRDS(file = "Data/snRNA-seq/temp/snRNA_cs_regressed.RDS")

# check quality control metrics across the clusters 
VlnPlot(rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",  'doublet_score'), ncol = 4, pt.size=0, group.by = 'seurat_clusters')
FeaturePlot(rna, features = c('percent.mt', 'nFeature_RNA', 'nCount_RNA', 'doublet_score'))
DimPlot(rna, group.by = 'Phase')

fibroblast_markers <- c("FAP", "PDGFRB", "PDGFRA", "ACTA2", "COL1A1")
mono_macro_markers <-c("MSR1", "CD163", "CCL3")
chromaffin_markers <- c("PNMT", "TH", "DBH","CHGA", "CHGB")
adrenocortical_markers <- c("STAR", "CYP11B1", "CYP11A1")
endothelial_markers <- c("FLT1", "EPAS1")
sustentacular_markers <- c("CDH19", "SOX10", "S100B", "VIM")
lymphocyte_markers <- c("CD2", "CD3E" , "MS4A1")

marker_genes_dp <- DotPlot(
  rna,
  group.by = "seurat_clusters",
  features = rev(c(
    chromaffin_markers,
    adrenocortical_markers,
    sustentacular_markers,
    mono_macro_markers,
    endothelial_markers,
    lymphocyte_markers,
    fibroblast_markers))
) +
  RotatedAxis()

VlnPlot(rna, features= "doublet_score")
clusters_dimplot / samples_dimplot
marker_genes_dp

#----
# DGE for all clusters
#----

# to help annotate, identify top DEGs for each cluster 
# cluster_markers <- FindAllMarkers(rna,
#                                   only.pos = TRUE,
#                                   min.pct = 0.25,
#                                   logfc.threshold = 0.25)

#write_tsv(cluster_markers, file = "results/differential_gene_expression/cluster_markers.tsv")
#cluster_markers <- read_tsv("results/differential_gene_expression/cluster_markers.tsv")
# get top 10 significant genes for each cluster 
# top_5 <- cluster_markers %>% 
#   filter(p_val_adj < 0.05) %>% 
#   group_by(cluster) %>% 
#   slice_max(n = 5, order_by = avg_log2FC)
# 
# top_5 %>% filter(cluster == 27)

#----
# Annotate clusters
#----

# TODO: change this to just annotate tumour clusters as "Tumour"
# can ten make another metadata column where they are annotated per sample of origin

#annotate the clusters according to cell type:
#based these annotations on the expression of marker genes (dotPlot) as well as scMatch predictions

# all cell types -  c("Tumour", "Chromaffin cells", "Adrenocortical cells", "Endothelial cells", "Fibroblasts", "SCLCs", "Myeloid cells", "Lymphocytes")

new.cluster.ids <- c(
  "Tumor", #0
  "Tumor", #1
  "Tumor", #2
  "Tumor", #3
  "Tumor", #4
  "Tumor", #5 
  "Tumor", #6
  "Tumor", #7
  "Tumor", #8
  "Tumor", #9
  "Tumor", #10
  "Chromaffin cells", #11 
  "Adrenocortical cells", #12
  "Myeloid cells", #13
  "SCLCs", #14
  "Fibroblasts",#15
  "Adrenocortical cells", #16
  "Endothelial cells", #17
  "Tumor", #18
  "Lymphocytes", #19
  "Endothelial cells", #20
  "SCLCs", #21
  "Tumor", #22
  "Tumor", #23
  "Fibroblasts", #24 
  "Tumor", #25
  "Tumor", #26
  "Myeloid cells" #27
  )

Idents(rna) <- rna$seurat_clusters
names(x = new.cluster.ids) <- levels(x = rna)
rna <- RenameIdents(object = rna, new.cluster.ids)
rna$cell_type <- Idents(rna)

DimPlot(rna, group.by = "cell_type")

# clusters labelled according to scMatch results and these marker genes

marker_genes_dp <- DotPlot(
  rna,
  features = rev(c(
    chromaffin_markers,
    adrenocortical_markers,
    sustentacular_markers,
    mono_macro_markers,
    endothelial_markers,
    lymphocyte_markers,
    fibroblast_markers))
) +
  RotatedAxis()

annotated_dimplot <- DimPlot(rna, cols = pal_genotypes_cell_types, label=T) + NoLegend()

(annotated_dimplot + samples_dimplot) /
  marker_genes_dp

# add a metadata columnn describing the metastasis driver genotype
metadata <- rna@meta.data
metadata <- metadata %>% 
  mutate(metastasis_driver = case_when(
  A5.ID %in% c("E240-1", "E243-1") ~ "Normal",
  cell_type == "Tumor" ~ Assumed.driver.of.metastais, 
  TRUE ~ "Normal"))
rownames(metadata) <- metadata$barcode
rna@meta.data <- metadata

DimPlot(rna, group.by = "metastasis_driver") + DimPlot(rna, group.by = "cell_type", cols = pal_genotypes_cell_types)

# TODO: add P018-PGL1/3 met driver information. I believe this should be "unknown" ...

# saveRDS(rna, file = "Data/snRNA-seq/processed/snRNA_cs_regressed_annotated.RDS")

# ----
# cell type proportions stacked bar plot
# ----

rna <- readRDS("Data/snRNA-seq/processed/snRNA_cs_regressed_annotated.RDS")

# add a metadata columnn describing the metastasis driver genotype
metadata <- rna@meta.data
metadata <- metadata %>% 
  mutate(sample_metastasis_driver = case_when(
    A5.ID %in% c("E240-1", "E243-1") ~ "Normal",
    A5.ID %in% c("P018-PGL1", "P018-PGL3") ~ "Unknown",
    TRUE ~ Assumed.driver.of.metastais))
rownames(metadata) <- metadata$barcode
rna@meta.data <- metadata

cell_type_frequencies <- rna@meta.data %>%
  dplyr::select(A5.ID, cell_type, sample_metastasis_driver) %>%
  group_by(A5.ID, cell_type, sample_metastasis_driver) %>%
  summarise(cell_type_count = n()) %>% # cell count in groups
  group_by(A5.ID) %>% 
  mutate(total = sum(cell_type_count)) %>%
  mutate(cell_type_prop = cell_type_count / total)

cell_types_stacked_bar <- ggplot(cell_type_frequencies,
                                  aes(x = A5.ID,
                                      y = cell_type_prop,
                                      fill = cell_type)) +
  geom_col(position = 'fill') +
  theme_classic() +
  facet_grid(~sample_metastasis_driver, scales = "free_x", space = "free_x")+
  labs(y = "Fraction") +
  scale_fill_manual(values = cell_type_colours) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank())

cell_types_stacked_bar +
  labs(colour = "") +
  plot_layout(guides = 'collect') 

# ----
# Schwann-cell-like cells in metastatic vs nonmetastatic tumours
# ----

# prior observations have suggested that metastatic PCPG have lower numbers of associated SCLCs



#----
# Look at differences in cell cycle phases between the genotypes 
#----

cell_cycle_phase <- rna@meta.data %>%
  # Only plot the chromaffin / neoplastic cells
  dplyr::filter(cell_type %in% c("Tumor","Chromaffin cells")) %>% 
  dplyr::select(A5.ID, Phase, sample_metastasis_driver) %>%
  group_by(A5.ID, Phase, sample_metastasis_driver) %>%
  summarise(phase_count = n()) %>% # cell count in groups
  group_by(A5.ID) %>% 
  mutate(total = sum(phase_count)) %>%
  mutate(phase_prop = phase_count / total)

cell_cycle_phase_stacked_bar <- ggplot(cell_cycle_phase,
                                 aes(x = A5.ID,
                                     y = phase_prop,
                                     fill = Phase)) +
  geom_col(position = 'fill') +
  theme_classic() +
  facet_grid(~sample_metastasis_driver, scales = "free_x", space = "free_x")+
  labs(y = "Fraction") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  labs(title = "neoplastic and chromaffin cells phases")

cell_cycle_phase_stacked_bar +
  labs(colour = "") +
  plot_layout(guides = 'collect') 

#----
# make a cells file for each tumor (neoplastic cells only)
#----

# these will be used to make neo-only bams for telomerehunter
# get barcodes corresponding to the neoplastic cells for each sample 
cell_tables <- rna[[c("barcode", "cell_type", "A5.ID")]] %>% 
  tibble() %>% 
  filter(cell_type %in% c("Tumor", "Chromaffin cells")) %>% 
  dplyr::select(-cell_type) %>% 
  mutate(barcode = str_extract(barcode, pattern="[^_]*$")) %>% 
  group_by(A5.ID) %>% 
  group_split()

all_ids <- unique(rna$A5.ID)

for (i in 1:length(all_ids)){
  sample_name <- cell_tables[[i]] %>% pull(A5.ID) %>% unique()
  write_tsv(cell_tables[[i]],
            file.path("results/snRNA-seq/neo_only_bams/",
                      paste0(sample_name, "_neo_only_cells.tsv")),
            col_names = FALSE)
}
