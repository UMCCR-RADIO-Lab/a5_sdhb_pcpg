# make heatmap for the bulk DEGs, show which are expressed in normal cell populations
rm(list=ls())

library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggsci)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
source("Scripts/Functions/dotplot_functions_A5.R")
# source("Analysis/load_signature_genesets_jansky.R")

# ----
# read in and organise the single cell PCPG data
# ----

# TODO: change to the 'doublets removed' version for this until I have a properly processed dataset
sn_pcpg <- readRDS("Data/snRNA-seq/processed/snRNA_cs_regressed_annotated.RDS")

# # remove the samples not belonging to the A5 cohort 
# sn_pcpg <- subset(sn_pcpg,
#                  subset = cell_type %in% c("E018", "E025", "E174", "PGL1","PGL3"),
#                  invert = TRUE)

sn_pcpg$cell_type_integration <- factor(sn_pcpg$cell_type_integration,
                                       levels = all_cell_types)

# make a metadata column describing tumour subtype and cell type
sn_md <- sn_pcpg@meta.data
sn_md <- sn_md %>%
  mutate(cell_type_tumour_genotype = recode(cell_type,
                                                    "E156-1" = "TERT",
                                                    "E123-1" = "TERT",
                                                    "E146-1" = "TERT",
                                                    "E143-1" = "TERT",
                                                    "E140-1" = "ATRX",
                                                    "E225-1" = "ATRX",
                                                    "E166-1" = "ATRX",
                                                    "E197-1" = "ATRX",
                                                    "E171-1" = "Unknown met driver")) 
sn_pcpg@meta.data <- sn_md
sn_pcpg$cell_type_tumour_genotype <- factor(sn_pcpg$cell_type_tumour_genotype,
                                           levels = c(
                                             "ATRX",
                                             "TERT",
                                             # "TERT (s)",
                                             "Unknown met driver",
                                             all_cell_types[1:length(all_cell_types)]))
# ----
# read in and organise the fetal single cell
# ----

sn_fetal <- readRDS("Data/published_datasets/Jansky_2021/adrenal_medulla_seurat.RDS")

sn_fetal$cell_type <- Idents(object = sn_fetal)

# make colour palette
colours_fetal <- setNames(object = discrete_palette, nm = as.character(unique(Idents(sn_fetal))))

# # remove the SCP cell types
# sn_fetal <- subset(sn_fetal,
#                                subset = cell_type %in% c("SCPs",
#                                                          "cycling SCPs",
#                                                          "late SCPs" ), 
#                                invert = TRUE)

fetal_md <- sn_fetal@meta.data

fetal_md <- fetal_md %>%
  mutate(cell_type_plot = recode(cell_type,
                                 "connecting Chromaffin cells" = "Connecting progenitor cells",
                                 "Chromaffin cells" = "Early chromaffin cells",
                                 "late Chromaffin cells" = "Late chromaffin cells",
                                 "Neuroblasts" = "Early neuroblasts",
                                 "cycling Neuroblasts" = "Cycling neuroblasts",
                                 "late Neuroblasts" = "Late neuroblasts",
                                 "SCPs" = "Early SCPs",
                                 "cycling SCPs"= "Cycling SCPs",
                                 "late SCPs" = "Late SCPs"))
sn_fetal@meta.data <- fetal_md

# ----
# organise the genes to plot 
# ----

tert_degs_bulk <- read_tsv("results/differential_gene_expression/bulk_RNA_seq/tertvswt_dge.tsv") %>%
  mutate(genotype = "TERT") %>% 
  filter(adj.P.Val < 0.05)

atrx_degs_bulk <- read_tsv("results/differential_gene_expression/bulk_RNA_seq/atrxvswt_dge.tsv") %>%
  mutate(genotype = "ATRX") %>% 
  filter(adj.P.Val < 0.05)

degs_bulk_all <- tert_degs_bulk %>% 
  bind_rows(atrx_degs_bulk)

# TODO: update singlecell with new ENSEMBL annotations
# FOR NOW: Remove genes that are not also present in the singlecell assay 
degs_bulk_all <- degs_bulk_all %>%
  dplyr::filter(gene_name %in% rownames(sn_pcpg)) %>%
  # note - no noncoding rna genes are in the jansky data, this step removes some LINCRNAs from the top list
  dplyr::filter(gene_name %in% rownames(sn_fetal)) 

# get the top 20 genes that are up in each genotype
top_degs_up <- degs_bulk_all %>% 
  group_by(genotype) %>% 
  slice_max(order_by = logFC, n = 20) %>% 
  ungroup()
top_degs_down <- degs_bulk_all %>% 
  group_by(genotype) %>% 
  slice_min(order_by = logFC, n = 20) %>% 
  ungroup()

genes_plot_df <- bind_rows(top_degs_up, top_degs_down) %>%
  mutate(direction = ifelse(logFC > 0, "up", "down")) %>% 
  arrange(genotype, direction, adj.P.Val) %>% 
  mutate(genotype_direction = paste0(genotype, " ",  direction)) %>% 
  mutate(genotype_direction = factor(genotype_direction, levels = c("ATRX up", "ATRX down", "TERT up", "TERT down")))

# remove all duplicated genes
genes_plot_df <- genes_plot_df[!duplicated(genes_plot_df$gene_name),]

# ----
# read in and organise the bulk data
# ----

# specify plotting order
# genotypes.lvl <- c("Normal", "SDHB", "SDHA", "SDHD", "VHL", "EPAS1", "FH", "HRAS", "NF1", "RET", "TMEM127", "MAX", "H3F3A", "Unknown", "MAML3")
# subgroups.lvl <- c("Normal", "C1Ai","C1Aii", "C1Bi", "C1Bii", "C2A", "C2Bi", "C2Bii")

# read in the bulk rna-seq metadata
bulk_metadata <- read_csv("Data/A5_full_clinical_and_genomic_version_2.csv") 
# bottom row has lot of NA's - this removes that row
bulk_metadata <- bulk_metadata[!is.na(bulk_metadata$`A5 ID`),]

bulk_qc <- read_csv("Data/bulk RNA-seq/RNA_seq_key_QC_stats.csv") %>%
  dplyr::rename(`A5 ID` = A5_ID) %>%
  dplyr::filter(!duplicated(`A5 ID`))

# Remove a sample with adrenocortical contamination,  another that is low quality, one that only has WGS, and E166-2
bulk_metadata <- bulk_metadata %>%
  dplyr::filter(!`A5 ID` %in% c("E154-1", "E144-1", "E233-1", "E166-2")) %>%
  mutate(`Metastatic state of tumour sample` = replace(`Metastatic state of tumour sample`, `Metastatic state of tumour sample` == "Metastatic", "Metastasis")) %>%
  # Join on the RNA-Seq specific QC
  left_join(bulk_qc) %>%
  # Add in patient info for MDS later on
  mutate(Patient = gsub("-.*", "", `A5 ID`)) %>%
  mutate(Patient_plot = replace(Patient,!(duplicated(Patient)| duplicated(Patient,fromLast = T)),NA))
# make column describing genotype
bulk_metadata <- bulk_metadata %>%
  mutate(metastasis_driver = recode(TERT_or_ATRX, 
                                    "Non_met_primary" = "Non-met primary",
                                    "Unknown_met_driver_primary" = "other/unknown",
                                    "Short_follow_up_primary" = "other/unknown",
                                    "Unknown_met_driver_met" = "other/unknown"))
bulk_metadata$metastasis_driver <- factor(bulk_metadata$metastasis_driver, levels = all_genotypes)
# make column for ALT status 
bulk_metadata <- bulk_metadata %>%
  mutate(ALT = if_else(`Difference in mean tumour and blood telomere lengths` > 1, "ALT", "No ALT"))

# order the bulk metadata according to mutation status 
bulk_metadata <- bulk_metadata %>%
  arrange(metastasis_driver, ALT, `A5 ID`, Gender)

# specify order I want the subgroups to be plotted
# bulk_metadata$Subgroup4 <- factor(bulk_metadata$Subgroup4, levels = c(subgroups.lvl))

# Read in bulk data (batch-normalised expression matrix)
bulk_rna <- read_csv("Data/bulk RNA-seq/RNA-Seq genewise log2 CPMs.csv")
colnames(bulk_rna) <- sub(colnames(bulk_rna), pattern = "\\.", replacement = "-")

# remove metadata for samples not in the bulk dataset
bulk_metadata <- bulk_metadata[bulk_metadata$`A5 ID` %in% colnames(bulk_rna),]

# remove the samples from the bulk that arent in WGS cohort 
bulk_rna <- bulk_rna %>%
  dplyr::select(all_of(c("SYMBOL", "ENSEMBL",
                         bulk_metadata$`A5 ID`))) %>% 
  # replace NA symbols with ENSG IDs
  dplyr::mutate("Gene" = if_else(is.na(SYMBOL), ENSEMBL, SYMBOL))
#TODO change the duplicated symbols for ENSEMBL IDS
# remove duplicated symbols
bulk_rna <- bulk_rna[!duplicated(bulk_rna$Gene),]

# Get the batch normed expression from the bulk
bulk_plot_mat <- as.matrix(dplyr::select(bulk_rna, -ENSEMBL, -Gene, -SYMBOL))
# Add genes as rownames
rownames(bulk_plot_mat) <- bulk_rna$Gene
# Z score transform
bulk_plot_mat <- t(scale(t(bulk_plot_mat)))

# filter out the genes that arent in the bulk (1 LncRNA gene)
genes_plot_df <- genes_plot_df[genes_plot_df$gene_name %in% rownames(bulk_plot_mat), ]
# subset the gene expression matrix to just include the genes i want to plot
bulk_heatmap_deg <- bulk_plot_mat[genes_plot_df$gene_name, ]

# check that the order of the heatmap columns hasn't changed 
table(bulk_metadata$`A5 ID` == colnames(bulk_plot_mat))

# ----
# Draw the bulk heatmap
# ----

# make an annotation describing TERT or ATRX mutation status
top_anno_bulk <- HeatmapAnnotation(
  "genotype" = bulk_metadata$metastasis_driver,
  "ALT" = bulk_metadata$ALT, 
  "Sex" = bulk_metadata$Gender, 
  show_annotation_name = c("genotype" = FALSE,
                           "ALT" = FALSE,
                           "Sex" = FALSE),
  annotation_legend_param = list(
    "genotype" = list(legend_height = unit(3, "cm")),
    "ALT" = list(legend_height = unit(3, "cm"))),
  col = list("genotype" = genotype_cols,
             "ALT" = ALT_cols,
             "Sex" = gender_cols),
  show_legend = c("genotype" = TRUE))

# plot the heatmap
bulk_hm <- Heatmap(bulk_heatmap_deg,
                   row_split = genes_plot_df$genotype_direction,
                   column_split = bulk_metadata$metastasis_driver,
                   width = ncol(bulk_heatmap_deg)*unit(1, "mm"),
                   height = ncol(bulk_heatmap_deg)*unit(5, "mm"),
                   row_gap = unit(2, "mm"),
                   column_gap = unit(2, "mm"),
                   row_title_rot = 0,
                   top_annotation = top_anno_bulk,
                   column_title_rot = 90,
                   cluster_columns = FALSE,
                   #show_column_dend = FALSE,
                   cluster_column_slices = FALSE,
                   cluster_rows = FALSE,
                   #show_row_dend = FALSE,
                   cluster_row_slices = FALSE,
                   show_row_names = TRUE,
                   row_names_side = "left",
                   show_column_names  = FALSE,
                   row_names_gp = gpar(fontface = "italic"),
                   heatmap_legend_param  = hm_legend_params)

# ----
# Draw the single cell PCPG Dotplot
# ----

genes_plot_df$gene_name %in% rownames(sn_pcpg)

pcpg_dot = HeatmapDotPlot.Seurat(sn_pcpg,
                                      features = genes_plot_df$gene_name,
                                      gene_grouping = genes_plot_df$genotype,
                                      aggr.by = "cell_type_tumour_genotype",
                                      cluster_row_slices = FALSE,
                                      cluster_columns = FALSE,
                                      cluster_rows = FALSE,
                                      assay = "RNA",
                                      row_names_side = "left",
                                      column_names_side = "top",
                                      slot = "data",
                                      col = c("blue", "grey", "red"),
                                      column_title_rot = 90,
                                      # left_annotation = side_anno_deg,
                                      row_names_gp = gpar(fontface = "italic"),
                                      heatmap_legend_param  = hm_legend_params)


# ----
# Draw the single cell FETAL Dotplot
# ----

genes_plot_fetal <- genes_plot_df[genes_plot_df$gene_name %in% rownames(sn_fetal),]

fetal_dot = HeatmapDotPlot.Seurat(sn_fetal,
                                       features = genes_plot_fetal$gene_name,
                                       gene_grouping = genes_plot_fetal$genotype,
                                       aggr.by = "cell_type",
                                       cluster_row_slices = FALSE,
                                       cluster_columns = FALSE,
                                       cluster_rows = FALSE,
                                       assay = "RNA",
                                       row_names_side = "left",
                                       column_names_side = "top",
                                       slot = "data",
                                       col = c("blue", "grey", "red"),
                                       column_title_rot = 90,
                                       # left_annotation = side_anno_deg,
                                       row_names_gp = gpar(fontface = "italic"),
                                       heatmap_legend_param  = hm_legend_params)

# dev.off()
# pdf(file = "Figures/heatmap_bulk_degs.pdf",
#     height = 24, width = 14)
# draw(bulk_hm + pcpg_dot + fetal_dot, annotation_legend_list=hm.legend,
#      merge_legend = FALSE,
#      heatmap_legend_side = "bottom",
#      annotation_legend_side = "right",
#      ht_gap = unit(c(10), "mm"))
# dev.off()

