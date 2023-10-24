library(patchwork)


setwd("/g/data/pq08/projects/ppgl")

######################
# Plotting Functions #
######################

source("/g/data/pq08/projects/ppgl/a5/sn_rna_seq/scripts/r_scripts/sn_rna_dotplot_functions.r")

source("/g/data/pq08/projects/ppgl/a5/methylation/scripts/helpers/plot_methylation.r")

################
# Data Loaders #
################

source("./a5/sn_rna_seq/scripts/data_loaders/a5_snrna_dataloader.r")
data_loader_a5_snrna()
a5_snrna <- snrna_annotate_cell_types(a5_snrna)

#load methylation array data
quickload_diff_meth = T; quickload_gsea <- T; quickload_dmr <- T;
source("./a5/methylation/scripts/a5_methylation_analysis_v2.r")

#Differential Expression - WTS
source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")

#Differential Expression - top de gene correlations
quickload_correlation_matrix = T
source("./a5/wts/scripts/differential_expression/wts_de_toptable_correlation_matrix.r")
mki67_cor <- c(cor_mat_rho$TERT_PriMet_vs_NonMetPri_WT["ENSG00000148773.14_MKI67",],
               cor_mat_rho$ATRX_PriMet_vs_NonMetPri_WT["ENSG00000148773.14_MKI67",],
               cor_mat_rho$Metastasis_All_vs_NonMetPri_WT["ENSG00000148773.14_MKI67",])
mki67_cor <- mki67_cor[!duplicated(names(mki67_cor))]

##############
# Annotation #
##############

#Add diff group info and fill in P018-PGL1, P018-PGL3, E240-1, and E243-1
md <- a5_snrna@meta.data %>% 
  left_join(a5_anno %>% dplyr::select(A5_ID, starts_with("differential"))) %>% 
  mutate(
    differential_group_sampletype_strict= as.character(differential_group_sampletype_strict),
    differential_group_anatomy=case_when(
      orig.ident %in% c("P018-PGL1", "P018-PGL3") ~ "Abdominal_Thoracic",
      orig.ident %in% c("E240-1", "E243-1" ) ~ "Normal",
      TRUE ~ differential_group_anatomy),
    differential_group_sampletype = case_when(
      orig.ident %in% c("P018-PGL1", "P018-PGL3") ~ "Non-metastatic primary",
      orig.ident %in% c("E240-1", "E243-1" ) ~ "Normal",
      TRUE ~ differential_group_sampletype),
    differential_group_sampletype_strict = case_when(
      orig.ident %in% c("P018-PGL1", "P018-PGL3") ~ "Non-metastatic primary",
      orig.ident %in% c("E240-1", "E243-1" ) ~ "Normal",
      TRUE ~ differential_group_sampletype_strict),
    differential_group_genotype = case_when(
      orig.ident %in% c("P018-PGL1", "P018-PGL3") ~ "WT_VAFOK",
      orig.ident %in% c("E240-1", "E243-1" ) ~ "WT_VAFOK",
      TRUE ~ differential_group_genotype),
    differential_group_purity = case_when(
      orig.ident %in% c("P018-PGL1", "P018-PGL3") ~ "PurityOK",
      orig.ident %in% c("E240-1", "E243-1" ) ~ "PurityOK",
      TRUE ~ differential_group_purity),
    TERT_ATRX_Mutation = case_when(
      orig.ident %in% c("P018-PGL1", "P018-PGL3") ~ "WT",
      orig.ident %in% c("E240-1", "E243-1" ) ~ "WT",
      TRUE ~ TERT_ATRX_Mutation))

cell_count_breaks <- c(0,5,10,20,50,100,500,50000)
cell_count_labels <- c("0-5", "5-10", "10-20", "20-50", "50-100", "100-500", "500+")

md <- md %>% 
  left_join(a5_snrna@meta.data %>% 
              group_by(orig.ident, cell_type) %>% 
              dplyr::count() %>% 
              mutate(n_cells_sample_celltype = cut(n,  
                                                   breaks=cell_count_breaks, 
                                                   labels=cell_count_labels)))
#left join removes rownames
rownames(md)  <- rownames(a5_snrna@meta.data)
a5_snrna@meta.data <- md


#################
# Heatmap Genes #
#################

genes_per_group <- 30
gene_select <- "topbottom_logfc" #"top_abslogfc"
adj_pval_cutoff <- 0.01
mki67_correlation_cutoff <- 0.3

GOI_lfc_pipe <- . %>% filter(adj.P.Val < adj_pval_cutoff) %>% 
  arrange(desc(abs(logFC))) %>% 
  #slice_head(n = 200) %>% 
  mutate(Symbol = gsub("ENSG[0-9]+([.][0-9]+)?_(.+)","\\2", Gene)) %>% 
  filter(Symbol %in% rownames(a5_snrna)) %>% 
  mutate(rank=row_number()) %>% 
  dplyr::select(Gene, Symbol, adj.P.Val, logFC, rank) %>% 
  {
    if(gene_select == "top_abslogfc")
    { . }
    else if(gene_select == "topbottom_logfc")
    { arrange(., logFC) }
  } %>% 
  mutate(mki67_correlation=mki67_cor[Gene])

GOI_slice_pipe <-  . %>% {
  if(gene_select == "top_abslogfc") { slice_min(., n = genes_per_group, order_by = rank) 
  }
  else if(gene_select == "topbottom_logfc") { 
    bind_rows(slice_max(., n = genes_per_group/2, order_by = logFC), 
              slice_min(., n = genes_per_group/2, order_by = logFC)) 
  }
}

GOI.met <- 
  wts_top_tables[["genosampletype"]][["Metastasis_All_vs_NonMetPri_WT"]] %>% 
  GOI_lfc_pipe %>%  
  GOI_slice_pipe %>% 
  mutate(source ="Metastasis_All_vs_NonMetPri_WT")
          

GOI.tert <- wts_top_tables[["genosampletype"]][["TERT_All_vs_NonTERT"]] %>% 
  GOI_lfc_pipe %>% 
  GOI_slice_pipe %>% 
  mutate(source="TERT_All_vs_NonTERT")
          

GOI.atrx <-
  wts_top_tables[["genosampletype"]][["ATRX_All_vs_NonATRX"]] %>% 
  GOI_lfc_pipe %>%  
  GOI_slice_pipe %>% 
  mutate(source="ATRX_All_vs_NonATRX")
  

GOI <- bind_rows(GOI.tert %>% dplyr::select(Gene, source),
                 GOI.atrx %>% dplyr::select(Gene, source),
                 GOI.met %>% dplyr::select(Gene, source)) #%>% filter(abs(mki67_correlation) < mki67_correlation_cutoff)

GOI$Gene[duplicated(GOI$Gene) | duplicated(GOI$Gene, fromLast=T)] <- paste0(GOI$Gene[duplicated(GOI$Gene) | duplicated(GOI$Gene, fromLast=T)], "*")
GOI <- GOI[!duplicated(GOI$Gene),]

# GOI.common <- GOI %>% 
#   group_by(Gene) %>% 
#   summarise(n=n(), rank=mean(rank), abs_logFC=mean(abs(logFC)), logFC=mean(logFC)) %>% 
#   filter(n==3) %>% 
#   GOI_slice_pipe %>% 
#   mutate(source="Common") %>% 
#   dplyr::select(Gene, rank, source)

# GOI.unique <- GOI %>% 
#   group_by(Gene) %>% 
#   mutate(n=n()) %>% 
#   filter(n==1) %>% 
#   group_by(source) %>% 
#   GOI_slice_pipe %>% 
#   dplyr::select(Gene, rank, source)


# GOI.plot <- bind_rows(GOI.common, GOI.unique) %>% distinct()

GOI.plot <- GOI %>% 
  mutate(Symbol = gsub("ENSG[0-9]+([.][0-9]+)?_(.+)","\\2", Gene),
         Symbol_Label=ifelse(grepl("^A[CL][0-9]", Symbol), 
                             stringr::str_extract(string = Gene, pattern = "ENSG[0-9]+"), 
                             Symbol),
         Symbol=gsub("[*]$","", Symbol),
         Gene=gsub("[*]$","", Gene))

GOI.plot$source <- factor(GOI.plot$source, 
                     levels = c(#"Common",
                                "Metastasis_All_vs_NonMetPri_WT",
                                "TERT_All_vs_NonTERT",
                                "ATRX_All_vs_NonATRX"))

##################
# DMR Membership #
##################

#Gene must be in a DMR in all comparisons to be flagged as DMR in common category (T=Require all, F=Require any)
strict_common_dmr = T

dmr_genes <- purrr::map(.x = dmr_lists[["genosampletype"]], 
           .f =  function (dmrs) {
             dmrs %>% GenomicRanges::as.data.frame() %>% 
               filter(Fisher < 0.05) %>% 
               separate_rows(overlapping.genes, sep=", ") %>% 
               filter(!is.na(overlapping.genes)) %>% 
               pull(overlapping.genes) })

GOI.plot$In_DMR <- "No"
for (i in 1:nrow(GOI.plot))
{
  if(GOI.plot$source[[i]] == "Common") {
    if (strict_common_dmr) {
    if(GOI.plot$Symbol[i] %in% purrr::reduce(dmr_genes,intersect)) { GOI.plot$In_DMR[i] <- "Yes"} }
    else {
    if(GOI.plot$Symbol[i] %in% unlist(dmr_genes)) { GOI.plot$In_DMR[i] <- "Yes"}  }
  } else {
    if(GOI.plot$Symbol[i] %in% dmr_genes[[GOI.plot$source[[i]]]]) { GOI.plot$In_DMR[i] <- "Yes"} 
  }
}

############
# WTS Expr #
############

sample.order <- a5_anno %>% 
  filter(A5_ID %in% colnames(a5_wts_lcpm_list[["SDHB_abdothoracic"]])) %>% 
  arrange(TERT_ATRX_Mutation, desc(differential_group_sampletype_strict)) %>% 
  pull(A5_ID)

plot.data <- GOI.plot %>% 
  inner_join(a5_wts_lcpm_list[["SDHB_abdothoracic"]][GOI.plot$Gene,sample.order] %>% 
               data.frame(check.names = F) %>%  
               tibble::rownames_to_column("Gene")) %>% inner_join(ensid_to_biotype) %>% 
  data.frame(check.names = F)
rownames(plot.data) <- plot.data$Symbol_Label

plot.data <- plot.data %>% mutate(mki67_correlation=mki67_cor[Gene])

hm.annot.data.col <- a5_anno %>% 
  filter(A5_ID %in% sample.order) %>% 
  mutate(A5_ID=factor(A5_ID, levels=sample.order)) %>% 
  arrange(A5_ID) %>% 
  dplyr::select(A5_ID, TERT_ATRX_Mutation, 
                differential_group_sampletype_strict)


hm.annot.col = HeatmapAnnotation(
  df = hm.annot.data.col[, -1],
  which = "column",
  show_legend = TRUE,
  col = list(
    "differential_group_sampletype_strict" = sampletype_strict_cols,
    "TERT_ATRX_Mutation" = driver_cols
  )
)

hm.annot.row = HeatmapAnnotation(df = plot.data[,c("source", "In_DMR","gene_biotype", "mki67_correlation"), drop=F], 
                                 which="row", 
                                 show_legend = TRUE,
                                 col=list("source"=setNames(
                                   RColorBrewer::brewer.pal(n=length(unique(plot.data$source)), name="Set3"),
                                   unique(plot.data$source)),
                                   In_DMR=c("Yes"="Black", "No"="lightgrey"),
                                 gene_biotype=setNames(
                                   RColorBrewer::brewer.pal(n=length(unique(plot.data$gene_biotype)), name="Set1"),
                                   unique(plot.data$gene_biotype)),
                                 mki67_correlation=colorRamp2(c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1), c("blue", "lightblue","white", "white", "white", "orange", "red"))))



hm_bulk_gex <- Heatmap(t(scale(t(plot.data %>% dplyr::select(-Gene, -Symbol,-source,-gene_biotype, -In_DMR, -Symbol_Label, -mki67_correlation)))),
        #col = col_fun,
        row_split = plot.data$source,
        #column_split = a5_anno.wts.nohn_noex$genotype_groups,
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
        #row_labels = plot.data$Symbol_Label,
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
                      features = GOI.plot$Symbol, 
                      aggr.by = "cell_type",
                      gene_grouping = GOI.plot$source,
                      show_row_names = FALSE
                      )

hm_snrna_dotplot_genotype <- HeatmapDotPlot.Seurat(object = subset(x = a5_snrna, subset = cell_type == "Tumor"),
                                          features = GOI.plot$Symbol, 
                                          aggr.by = "TERT_ATRX_Mutation",
                                          gene_grouping = GOI.plot$source,
                                          show_row_names = FALSE
)

hm_bulk_gex + hm_snrna_dotplot_genotype + hm_snrna_dotplot
