library(patchwork)


setwd("/g/data/pq08/projects/ppgl")

################
# Data Loaders #
################

#load clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#WTS
source("./a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
htseq_outs <- "./a5/wts/analysis/htseq/truseq/gene"
data_loader_a5_wts_counts(count_file_dir=htseq_outs)

#Differential Expression - smallRNA
source("./a5/small_rna/scripts/differential_expression/a5_smallrna_differential_expression.r")


###########
# Heatmap #
###########

sample.order <- a5_anno %>% filter(Exclude=="N") %>% 
  arrange(TERT_ATRX_Mutation, desc(differential_group_sampletype_strict)) %>% pull(A5_ID)
sample.order <- intersect(sample.order, colnames(a5_wts_dge_list[["SDHB_abdothoracic"]]))
missing <- setdiff(sample.order, colnames(a5_smallrna_lcpm_list[["SDHB_abdothoracic"]]))
sample.order.complete <- sample.order
sample.order <- intersect(sample.order, colnames(a5_smallrna_lcpm_list[["SDHB_abdothoracic"]]))

genes_per_group <- 20

GOI_pipe <- . %>% filter(adj.P.Val < 0.05) %>% 
  arrange(desc(abs(logFC))) %>% 
  #slice_head(n = 200) %>% 
  mutate(rank=row_number()) %>% 
  dplyr::select(Gene, adj.P.Val, logFC, rank)

GOI.met <- 
  smallrna_top_tables[["genosampletype"]][["Metastasis_All_vs_NonMetPri_WT"]] %>% 
  GOI_pipe %>% 
  mutate(source ="Metastasis_All_vs_NonMetPri_WT")


GOI.tert <- 
  smallrna_top_tables[["genosampletype"]][["TERT_All_vs_NonTERT"]] %>% 
  GOI_pipe %>%  
  mutate(source="TERT_All_vs_NonTERT")


GOI.atrx <-
  smallrna_top_tables[["genosampletype"]][["ATRX_All_vs_NonATRX"]] %>% 
  GOI_pipe %>% 
  mutate(source="ATRX_All_vs_NonATRX")


GOI <- bind_rows(GOI.met, GOI.tert, GOI.atrx) 

# GOI.unique <- GOI %>% 
#   group_by(Gene) %>% 
#   mutate(n=n()) %>% 
#   filter(n==1) %>% 
#   group_by(source) %>% 
#   slice_min(n = genes_per_group, order_by = rank) %>% 
#   dplyr::select(Gene, rank, source)


GOI.plot <- GOI

GOI.plot$source <- factor(GOI.plot$source, 
                          levels = c("Metastasis_All_vs_NonMetPri_WT",
                                     "TERT_All_vs_NonTERT",
                                     "ATRX_All_vs_NonATRX"))

########
# Expr #
########

plot.data <- GOI.plot %>% 
  inner_join(a5_smallrna_lcpm_list[["SDHB_abdothoracic"]][GOI.plot$Gene,sample.order] %>% 
               as_tibble(rownames = "Gene")) %>% mutate(gene_biotype="miRNA") %>% 
  data.frame(check.names = F) %>% 
  distinct()
rownames(plot.data) <- gsub("[.]","-",make.names(plot.data$Gene, unique = T))
plot.data[,missing] <- NA
plot.data <- plot.data[,c("source", "gene_biotype", sample.order.complete)]

hm.annot.data.col <- a5_anno %>% 
  filter(A5_ID %in% sample.order) %>% 
  dplyr::select(A5_ID, TERT_ATRX_Mutation, 
                differential_group_sampletype_strict)

hm.annot.data.col <- hm.annot.data.col[match(sample.order.complete, hm.annot.data.col$A5_ID),]

hm.annot.col = HeatmapAnnotation(
  df = hm.annot.data.col[, -1],
  which = "column",
  show_legend = TRUE,
  col = list(
    "differential_group_sampletype_strict" = sampletype_strict_cols,
    "TERT_ATRX_Mutation" = driver_cols
  )
)

hm.annot.row = HeatmapAnnotation(df = plot.data[,c("source", "gene_biotype"), drop=F], 
                                 which="row", 
                                 show_legend = TRUE,
                                 col=list("source"=setNames(
                                   RColorBrewer::brewer.pal(n=10, name="Set3")[1:length(unique(plot.data$source))],
                                   unique(plot.data$source)),
                                   gene_biotype=c("miRNA"="#637b62ff")))



hm_bulk_smallrna <- Heatmap(t(scale(t(plot.data %>% dplyr::select(-source,-gene_biotype)))),
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
                       show_row_names = TRUE,
                       row_names_side = "left",
                       show_column_names  = TRUE,
                       row_names_gp = gpar(fontsize = 5, fontface = "italic"),
                       #heatmap_legend_param  = hm_legend_params
)


hm_bulk_smallrna
