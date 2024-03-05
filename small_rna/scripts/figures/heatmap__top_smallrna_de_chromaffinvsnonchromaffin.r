library(patchwork)


setwd("/g/data/pq08/projects/ppgl")

######################
# Plotting Functions #
######################

# source("/g/data/pq08/projects/ppgl/a5/sn_rna_seq/scripts/r_scripts/sn_rna_dotplot_functions.r")

source("/g/data/pq08/projects/ppgl/a5/methylation/scripts/helpers/plot_methylation.r")

################
# Data Loaders #
################

#Small RNA
source("./a5/small_rna/scripts/differential_expression/a5_smallrna_differential_expression.r")

#############################
# Heatmap Genes - Small RNA #
#############################

genes_per_group <- 100
gene_select <- "topbottom_logfc" #"top_abslogfc" "topbottom_logfc" "pval_logfc_direction"
adj_pval_cutoff <- 0.05
lfc_cutoff <- 1

get_sig_gene_ids <- . %>% filter(adj.P.Val < adj_pval_cutoff, abs(logFC) > lfc_cutoff) %>% pull(Gene)

GOI_lfc_pipe <- . %>% filter(adj.P.Val < adj_pval_cutoff, abs(logFC) > lfc_cutoff) %>% 
  mutate(direction = ifelse(logFC > 0, "Up", "Down")) %>% 
  mutate(Symbol = gsub("ENSG[0-9]+([.][0-9]+)?_(.+)","\\2", Gene)) %>% 
  #filter(Symbol %in% rownames(a5_snrna)) %>% 
  mutate(rank=row_number()) %>% 
  dplyr::select(Gene, Symbol, adj.P.Val, logFC, rank, direction) %>% 
  {
    if(gene_select == "top_abslogfc")
    { arrange(., desc(abs(logFC))) }
    else if(gene_select == "topbottom_logfc")
    { arrange(., logFC) } 
    else if(gene_select == "pval_logfc_direction") {
      arrange(., adj.P.Val)
    }
  }

GOI_slice_pipe <-  . %>% {
  if(gene_select == "top_abslogfc" & nrow(.) > genes_per_group) { 
    slice_min(., n = genes_per_group, order_by = rank) 
  }
  else if(gene_select == "topbottom_logfc" & nrow(.) > genes_per_group) { 
    bind_rows(slice_max(., n = genes_per_group/2, order_by = logFC), 
              slice_min(., n = genes_per_group/2, order_by = logFC)) 
  }
  else if(gene_select == "pval_logfc_direction" & nrow(.) > genes_per_group) { 
    group_by(., direction) %>% 
      slice_min(n = genes_per_group/2, order_by = adj.P.Val) %>% 
      ungroup()
  } else ( . )
}


GOI_plot <- smallrna_top_tables[["Non_chromaffin_vs_Chromaffin"]][["Non_chromaffin_vs_Chromaffin"]] %>% 
  GOI_lfc_pipe %>%  
  GOI_slice_pipe

#######################
# Heatmap - Small RNA #
#######################

a5_anno_use <- a5_anno %>% 
  mutate(primary_location_plotting = case_when(
    A5_ID == "E185-1" ~ "Head and neck", 
    A5_ID == "E128-1" ~ "Extraadrenal (abdominal/thoracic)",
    primary_location_plotting == "Extraadrenal (bladder)" ~ "Extraadrenal (abdominal/thoracic)",
    differential_group_anatomy == "Thoracic_non_chromaffin" ~ "Thoracic (non-chromaffin)", 
    .default = primary_location_plotting)) %>%
  filter(A5_ID %in% colnames(a5_smallrna_lcpm_list[["SDHB"]])) %>% 
  arrange(differential_group_anatomy, cell_of_origin) %>% 
  mutate(A5_ID = factor(A5_ID, levels=.$A5_ID))

sample.order <- levels(a5_anno_use$A5_ID)

plot.data <- a5_smallrna_lcpm_list[["SDHB"]][GOI_plot$Gene, sample.order] %>% 
  data.frame(check.names = F) %>% 
  tibble::rownames_to_column("Gene") 
rownames(plot.data) <-plot.data$Gene

hm.annot.data.col <- a5_anno_use %>% 
  dplyr::select(A5_ID, differential_group_anatomy, 
                TERT_ATRX_Mutation, differential_group_purity, 
                differential_group_anatomy, Catecholamine_profile, 
                differential_group_sampletype)

hm.annot.data.col <- hm.annot.data.col[match(sample.order, hm.annot.data.col$A5_ID),]

hm.annot.col = HeatmapAnnotation(
  df = hm.annot.data.col[, 2, drop=F],
  which = "column",
  show_legend = TRUE,
  col = list(
    "differential_group_anatomy" = location_cols#,
    #"differential_group_sampletype" = specimen_type_cols,
    #"TERT_ATRX_Mutation" = driver_cols,
    # "primary_location_plotting" = location_cols,
    # "Catecholamine_profile" = setNames(
    #   RColorBrewer::brewer.pal(
    #     n = length(unique(hm.annot.data.col$Catecholamine_profile)),
    #     name = "Set1"),
    #   unique(hm.annot.data.col$Catecholamine_profile))
  ))


mir_chr_anno <- smallrna_feature_annotation %>% 
  dplyr::select(Geneid, Chr) %>% separate_rows(Chr, sep=";") %>% 
  distinct() %>% 
  group_by(Geneid) %>% 
  summarise(Chr=ifelse(n() >1, "Multiple", Chr)) %>% 
  dplyr::rename(Gene=Geneid)

plot.data <- plot.data %>% left_join(mir_chr_anno)

chrom_cols<-c(RColorBrewer::brewer.pal(n=8, name="Set1"),RColorBrewer::brewer.pal(n=8, name="Set2"),RColorBrewer::brewer.pal(n=length(unique(plot.data$Chr))-16, name="Set3"))

hm.annot.row = HeatmapAnnotation(Chromosome = plot.data$Chr, 
                                 which="row", 
                                 show_legend = TRUE,
                                 col=list("Chromosome"=setNames(
                                   chrom_cols,
                                   unique(plot.data$Chr))[1:length(unique(plot.data$Chr))])
)

rownames(plot.data) <-  plot.data$Gene

hm <- Heatmap(t(scale(t(plot.data %>% dplyr::select(-Gene, -Chr)))),
              #col = col_fun,
              #row_split = plot.data$source,
              column_split = hm.annot.data.col$differential_group_anatomy,
              #width = ncol(bulk_heatmap_deg)*unit(2, "mm"),
              #height = ncol(bulk_heatmap_deg)*unit(4, "mm"),
              row_gap = unit(2, "mm"),
              column_gap = unit(2, "mm"),
              row_title_rot = 0,
              top_annotation = hm.annot.col,
              left_annotation = hm.annot.row,
              column_title_rot = 90,
              cluster_columns = FALSE,
              #show_column_dend = FALSE,
              cluster_column_slices = FALSE,
              cluster_rows = TRUE,
              #show_row_dend = FALSE,
              cluster_row_slices = FALSE,
              show_row_names = TRUE,
              row_names_side = "left",
              show_column_names  = TRUE,
              row_names_gp = gpar(fontsize = 5, fontface = "italic"),
              #heatmap_legend_param  = hm_legend_params
)

#############
# Dot Plots #
#############

hm_snrna_dotplot <- HeatmapDotPlot.Seurat(object = a5_snrna,
                                          features = plot.data$Symbol,
                                          aggr.by = "cell_type",
                                          #gene_grouping = GOI.plot$source,
                                          show_row_names = FALSE,
                                          cell.size = 0.15
)

hm_snrna_dotplot_genotype <- HeatmapDotPlot.Seurat(object = subset(x = a5_snrna, subset = cell_type == "Tumor"),
                                                   features = plot.data$Symbol,
                                                   aggr.by = "TERT_ATRX_Mutation",
                                                   #gene_grouping = GOI.plot$source,
                                                   show_row_names = FALSE,
                                                   cell.size = 0.15
)

plot_list[[contrast]][[mode]]  <- hm_bulk_gex + hm_snrna_dotplot_genotype + hm_snrna_dotplot


# saveRDS(plot_list[[contrast]][[mode]], "/g/data/pq08/projects/ppgl/a5/tertiary/results/snrna_wts_crossreference/heatmaplist_20240115.rds")
# 
# plot_config <-
#   as_tibble(
#     matrix(
#       c(
#         "Metastatic_All_vs_NonMetPri_WT", 6,9.5,
#         "TERT_PriMet_vs_NonMetPri_WT", 6,9.5,
#         "ATRX_PriMet_vs_NonMetPri_WT", 6,9.5,
#         "ATRX_All_vs_TERT_All", 6,9.5,
#         "ATRX_TERT_Met_vs_NMP", 6,9.5,
#         "ATRXvsTERT_ATRXvsNMP", 6,9.5,
#         "ATRXvsTERT_TERTvsNMP", 4,9.5
#       ),
#       byrow = T,
#       ncol = 3,
#       dimnames = list(c(),c("contrast","height","width"))))
# 
# purrr::pwalk(.l = plot_config,
#              .f = \(contrast, height, width) {
#                pdf(file = file.path("/g/data/pq08/projects/ppgl/a5/tertiary/results/snrna_wts_crossreference/", paste0("snrna_wts_heatmap_",contrast, "_20240115.pdf")),
#                    height = as.numeric(height),
#                    width =  as.numeric(width))
#                print(plot_list[[contrast]]$contrast_only_extended)
#                dev.off()
#              }
# )
