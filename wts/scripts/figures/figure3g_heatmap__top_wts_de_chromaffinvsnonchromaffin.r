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

# source("./a5/sn_rna_seq/scripts/data_loaders/a5_snrna_dataloader.r")
# data_loader_a5_snrna()
# a5_snrna <- snrna_annotate_cell_types(a5_snrna)

#load methylation array data
quickload_diff_meth = T; quickload_gsea <- T; quickload_dmr <- T;
source("./a5/methylation/scripts/a5_methylation_analysis.r")

##Differential Expression 
#WTS
output_qc <- F; output_tables <- F; output_plots <- F;
source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")

#######################
# Heatmap Genes - WTS #
#######################

genes_per_group <- 60
gene_select <- "topbottom_logfc" #"top_abslogfc" "pval_logfc_direction"
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




GOI <- wts_top_tables[["Non_chromaffin_vs_Chromaffin"]][["Non_chromaffin_vs_Chromaffin"]] %>% 
  filter(gene_biotype == "protein_coding") %>% 
  GOI_lfc_pipe %>%  
  GOI_slice_pipe


GOI.plot <- GOI %>% 
  mutate(Symbol = gsub("ENSG[0-9]+([.][0-9]+)?_(.+)","\\2", Gene),
         Symbol_Label=ifelse(grepl("^A[CL][0-9]", Symbol), 
                             stringr::str_extract(string = Gene, pattern = "ENSG[0-9]+"), 
                             Symbol)) %>% 
  dplyr::select(Gene, Symbol, Symbol_Label)

##################
# DMR Membership #
##################

dmr_genes <- dmr_lists[["hn"]][["Parasympathetic_vs_Sympathetic"]] %>% 
  GenomicRanges::as.data.frame() %>% 
  filter(Fisher < 0.05) %>% 
  separate_rows(overlapping.genes, sep=", ") %>% 
  filter(!is.na(overlapping.genes)) %>% 
  pull(overlapping.genes)

GOI.plot$In_DMR <- "No"
for (i in 1:nrow(GOI.plot))
{
  if(GOI.plot$Symbol[i] %in% dmr_genes) { GOI.plot$In_DMR[i] <- "Yes"} 
}


#################
# Heatmap - WTS #
#################

a5_anno_use <- a5_anno %>% 
  mutate(primary_location_plotting = case_when(
    A5_ID == "E185-1" ~ "Head and neck", 
    A5_ID == "E128-1" ~ "Extraadrenal (abdominal/thoracic)",
    primary_location_plotting == "Extraadrenal (bladder)" ~ "Extraadrenal (abdominal/thoracic)",
    differential_group_anatomy == "Thoracic_non_chromaffin" ~ "Thoracic (non-chromaffin)", 
    .default = primary_location_plotting)) %>%
  filter(A5_ID %in% colnames(a5_wts_lcpm_list[["SDHB"]])) %>% 
  arrange(differential_group_anatomy, cell_of_origin) %>% 
  mutate(A5_ID = factor(A5_ID, levels=.$A5_ID))

sample.order <- levels(a5_anno_use$A5_ID)

contrast_membership <- contrastdesign_to_memberlist("Non_chromaffin_vs_Chromaffin", contrast_matrix_hn, design_matrix_hn, contrast_name_sep="_vs_")

plot.data <- GOI.plot %>% 
  inner_join(a5_wts_lcpm_list[["SDHB"]][GOI.plot$Gene,sample.order] %>% 
               data.frame(check.names = F) %>%  
               tibble::rownames_to_column("Gene")) %>% inner_join(ensid_to_biotype) %>% 
  data.frame(check.names = F) %>% 
  distinct()


plot.data.EnsIds <- gsub("(ENSG[0-9]+)([.][0-9]+)?_.+","\\1", plot.data$Gene)
data_loader_ensgid_to_chr(plot.data.EnsIds, use_cache=TRUE)
ensgid_to_chr <- data.frame(Gene=plot.data$Gene, ensembl_gene_id=plot.data.EnsIds) %>% left_join(ensgid_to_chr) %>% dplyr::select(-ensembl_gene_id, -start_position, -end_position)
plot.data <- plot.data %>% left_join(ensgid_to_chr)
plot.data$chromosome_name <- paste0("chr", plot.data$chromosome_name)

rownames(plot.data) <- plot.data$Symbol_Label

hm.annot.data.col <- a5_anno_use %>% 
  dplyr::select(A5_ID, primary_location_plotting)


hm.annot.col = HeatmapAnnotation(
  df = hm.annot.data.col[, -1, drop=F],
  which = "column",
  show_legend = TRUE,
  col = list("primary_location_plotting" = location_cols))


chrom_cols<-c(RColorBrewer::brewer.pal(n=8, name="Set1"),
              RColorBrewer::brewer.pal(n=8, name="Set2"),
              RColorBrewer::brewer.pal(n=length(unique(plot.data$chromosome_name))-16, name="Set3"))[1:length(unique(plot.data$chromosome_name))]

hm.annot.row = HeatmapAnnotation(df = plot.data[,c( "In_DMR", "chromosome_name"), drop=F], #"source",
                                 which="row", 
                                 show_legend = TRUE,
                                 col=list(
                                   #"source"=setNames(RColorBrewer::brewer.pal(n=length(unique(plot.data$source)), name="Set3")[1:length(unique(plot.data$source))], unique(plot.data$source)),
                                   In_DMR=c("Yes"="Black", "No"="lightgrey"),
                                   gene_biotype=bio_type_cols,
                                   "chromosome_name"=setNames(
                                     chrom_cols,
                                     unique(plot.data$chromosome_name)) 
                                   ))



hm_bulk_gex <- Heatmap(t(scale(t(plot.data %>% dplyr::select(-Gene, -Symbol, -gene_biotype, -In_DMR, -Symbol_Label, -chromosome_name)))),
                       #col = col_fun,
                       #row_split = plot.data$source,
                       column_split = hm.annot.data.col$primary_location_plotting,
                       #width = ncol(bulk_heatmap_deg)*unit(2, "mm"),
                       #height = ncol(bulk_heatmap_deg)*unit(4, "mm"),
                       row_gap = unit(2, "mm"),
                       column_gap = unit(2, "mm"),
                       row_title_rot = 90,
                       top_annotation = hm.annot.col,
                       left_annotation = hm.annot.row,
                       column_title_rot = 90,
                       cluster_columns = FALSE,
                       #show_column_dend = FALSE,
                       cluster_column_slices = FALSE,
                       cluster_rows = TRUE,
                       #show_row_dend = FALSE,
                       cluster_row_slices = FALSE,
                       #row_labels = plot.data$Symbol_Label,
                       show_row_names = TRUE,
                       row_names_side = "left",
                       show_column_names  = TRUE,
                       row_names_gp = gpar(fontsize = 7, fontface = "italic"),
                       #heatmap_legend_param  = hm_legend_params
)


# #############
# # Dot Plots #
# #############
# 
# hm_snrna_dotplot <- HeatmapDotPlot.Seurat(object = a5_snrna,
#                                           features = plot.data$Symbol,
#                                           aggr.by = "cell_type",
#                                           #gene_grouping = GOI.plot$source,
#                                           show_row_names = FALSE,
#                                           cell.size = 0.15
# )
# 
# hm_snrna_dotplot_genotype <- HeatmapDotPlot.Seurat(object = subset(x = a5_snrna, subset = cell_type == "Tumor"),
#                                                    features = plot.data$Symbol,
#                                                    aggr.by = "TERT_ATRX_Mutation",
#                                                    #gene_grouping = GOI.plot$source,
#                                                    show_row_names = FALSE,
#                                                    cell.size = 0.15
# )
# 
# plot_list[[contrast]][[mode]]  <- hm_bulk_gex + hm_snrna_dotplot_genotype + hm_snrna_dotplot
# 

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
