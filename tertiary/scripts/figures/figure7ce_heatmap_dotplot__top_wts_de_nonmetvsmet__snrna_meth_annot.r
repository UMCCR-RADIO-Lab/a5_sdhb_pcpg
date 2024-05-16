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
source("./a5/methylation/scripts/a5_methylation_analysis.r")

#Differential Expression - WTS
output_qc <- F; output_tables <- F; output_plots <- F;
source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")

#Differential Expression - top de gene correlations
quickload_correlation_matrix = T
source("./a5/wts/scripts/differential_expression/wts_de_toptable_correlation_matrix.r")
mki67_cor <- c(cor_mat_rho$TERT_PriMet_vs_NonMetPri_WT["ENSG00000148773.14_MKI67",],
               cor_mat_rho$ATRX_PriMet_vs_NonMetPri_WT["ENSG00000148773.14_MKI67",],
               cor_mat_rho$Metastatic_All_vs_NonMetPri_WT["ENSG00000148773.14_MKI67",])
mki67_cor <- mki67_cor[!duplicated(names(mki67_cor))]

#Assays
gs4_auth("aidan.flynn@umccr-radio-lab.page")
a5_assays <- read_sheet("1hnXdXI29KvvuLxsaTBID1-EbE7mTSk6bG05HcSFfhgo", 
                        sheet="AssaysPerformed", 
                        col_types = "c")

snrna_samples <- a5_assays %>%  filter(`snRNA-Seq` == 1) %>%  pull(A5_ID)

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



#########################################
# Remove non-A5 samples for publication #
#########################################

pub_keep <- !(a5_snrna@meta.data$orig.ident %in% c("E018", "P018-PGL1", "P018-PGL3"))
a5_snrna@meta.data$publication_include <- pub_keep
a5_snrna <- subset(a5_snrna, publication_include)

#################
# Heatmap Genes #
#################

genes_per_group <- 30
gene_select <- "pval_logfc_direction" #"top_abslogfc" "topbottom_logfc"
adj_pval_cutoff <- 0.05
lfc_cutoff <- 1
mki67_correlation_cutoff <- 0.3

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
  } %>% 
  mutate(mki67_correlation=mki67_cor[Gene])

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


GOI <- list()
for(contrast in c("Metastatic_All_vs_NonMetPri_WT", "TERT_PriMet_vs_NonMetPri_WT", "ATRX_PriMet_vs_NonMetPri_WT", "ATRX_All_vs_TERT_All"))# "TERT_All_vs_NonTERT", "ATRX_All_vs_NonATRX",
{
  GOI[[contrast]]  <-
    wts_top_tables[["genosampletype"]][[contrast]] %>% 
    filter(gene_biotype == "protein_coding") %>% 
    GOI_lfc_pipe %>%  
    GOI_slice_pipe %>% 
    mutate(source = contrast)
}

GOI[["ATRX_TERT_Met_vs_NMP"]]  <-
  wts_top_tables[["genosampletype"]][["Metastatic_All_vs_NonMetPri_WT"]] %>% 
  filter(gene_biotype == "protein_coding") %>% 
  filter(Gene %in% (wts_top_tables[["genosampletype"]][["TERT_PriMet_vs_NonMetPri_WT"]] %>% get_sig_gene_ids), 
         Gene %in% (wts_top_tables[["genosampletype"]][["ATRX_PriMet_vs_NonMetPri_WT"]] %>% get_sig_gene_ids)) %>% 
  GOI_lfc_pipe %>%  
  GOI_slice_pipe %>% 
  mutate(source = "ATRX_TERT_Met_vs_NMP")

GOI[["ATRXvsTERT_ATRXvsNMP"]]  <-
  wts_top_tables[["genosampletype"]][["ATRX_All_vs_TERT_All"]] %>% 
  filter(gene_biotype == "protein_coding") %>% 
  filter(Gene %in% (wts_top_tables[["genosampletype"]][["ATRX_PriMet_vs_NonMetPri_WT"]] %>% get_sig_gene_ids)) %>% 
  GOI_lfc_pipe %>%  
  GOI_slice_pipe %>% 
  mutate(source = "ATRXvsTERT_ATRXvsNMP")

GOI[["ATRXvsTERT_TERTvsNMP"]]  <-
  wts_top_tables[["genosampletype"]][["ATRX_All_vs_TERT_All"]] %>% 
  filter(gene_biotype == "protein_coding") %>% 
  filter(Gene %in% (wts_top_tables[["genosampletype"]][["TERT_PriMet_vs_NonMetPri_WT"]] %>% get_sig_gene_ids)) %>% 
  GOI_lfc_pipe %>%  
  GOI_slice_pipe %>% 
  mutate(source = "ATRXvsTERT_TERTvsNMP")

GOI <- purrr::map(GOI, .f = function(x) { x %>% dplyr::select(Gene, source) }) %>%  bind_rows()

#GOI$Gene[duplicated(GOI$Gene) | duplicated(GOI$Gene, fromLast=T)] <- paste0(GOI$Gene[duplicated(GOI$Gene) | duplicated(GOI$Gene, fromLast=T)], "*")
#GOI <- GOI[!duplicated(GOI$Gene),]


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
                             Symbol))#,
#Symbol=gsub("[*]$","", Symbol),
#Gene=gsub("[*]$","", Gene))

# GOI.plot$source <- factor(GOI.plot$source, 
#                      levels = c(#"Common",
#                                 "Metastatic_All_vs_NonMetPri_WT",
#                                 #"TERT_All_vs_NonTERT",
#                                 #"ATRX_All_vs_NonATRX",
#                                 "ATRX_PriMet_vs_NonMetPri_WT",
#                                 "TERT_PriMet_vs_NonMetPri_WT"#,
#                                 #"ATRX_All_vs_TERT_All"
#                                 ))



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
  } else if(GOI.plot$source[[i]] == "ATRX_TERT_Met_vs_NMP") {
    if(GOI.plot$Symbol[i] %in% dmr_genes[["Metastatic_All_vs_NonMetPri_WT"]]) { GOI.plot$In_DMR[i] <- "Yes"} 
  } else if(GOI.plot$source[[i]] %in% c("ATRXvsTERT_ATRXvsNMP","ATRXvsTERT_TERTvsNMP")) {
    if(GOI.plot$Symbol[i] %in% dmr_genes[["ATRX_PriMet_vs_NonMetPri_WT"]]) { GOI.plot$In_DMR[i] <- "Yes"} 
  } else if(GOI.plot$source[[i]] %in% c("ATRXvsTERT_TERTvsNMP","ATRXvsTERT_TERTvsNMP")) {
    if(GOI.plot$Symbol[i] %in% dmr_genes[["TERT_PriMet_vs_NonMetPri_WT"]]) { GOI.plot$In_DMR[i] <- "Yes"} 
  } else {
    if(GOI.plot$Symbol[i] %in% dmr_genes[[GOI.plot$source[[i]]]]) { GOI.plot$In_DMR[i] <- "Yes"} 
  }
}


############
# Heatmaps #
############

plot_list <- list()

for (contrast in c("ATRXvsTERT_ATRXvsNMP", "ATRXvsTERT_TERTvsNMP"))#unique(GOI.plot$source)
{
  
  plot_list[[contrast]] <- list()
  
  for (mode in c("contrast_only_extended"))#"all_sample", "contrast_only", 
  {
    plot_list[[contrast]][[mode]] <- list()
    
    
    message("Processing - ", contrast, mode)
    
    GOI.plot.current <- GOI.plot %>%  filter(source == contrast)
    
    sample.order <- a5_anno %>% 
      filter(A5_ID %in% colnames(a5_wts_lcpm_list[["SDHB_abdothoracic"]])) %>% 
      arrange(TERT_ATRX_Mutation, desc(differential_group_sampletype_strict)) %>% 
      pull(A5_ID)
    
    contrast_group <- case_when(contrast == "ATRX_TERT_Met_vs_NMP" ~ "Metastatic_All_vs_NonMetPri_WT",
                                contrast %in% c("ATRXvsTERT_ATRXvsNMP", "ATRXvsTERT_TERTvsNMP") ~ "ATRX_All_vs_TERT_All",
                                TRUE ~ contrast)
    contrast_membership <- contrastdesign_to_memberlist(contrast_group, contrast_matrix_genosampletype, design_matrix_genosampletype, contrast_name_sep="_vs_")
    if (mode == "contrast_only") {
      keep_samples <- contrast_membership %>% filter(group != "non_participant") 
      
      sample.order <- keep_samples %>% arrange(group) %>% pull(A5_ID)
      
    } else if(mode == "contrast_only_extended") {
      keep_samples_participating <- contrast_membership %>% filter(group != "non_participant") %>%  pull(A5_ID)
      
      keep_samples_extended <- a5_anno %>% 
        filter((differential_group_sampletype_strict %in%  c("Metastatic primary", "Metastasis" ) & TERT_ATRX_Mutation == "WT") |
                 TERT_ATRX_Mutation != "WT")  %>%  pull(A5_ID)
      keep_samples_extended <- c(keep_samples_extended, snrna_samples)
      
      
      set.seed(100)
      keep_samples_extended_variable <- 
        contrast_membership %>% 
        filter(group == "non_participant") %>% 
        slice_sample(prop = 0.2)  %>%  pull(A5_ID)
      
      keep_samples <- c(keep_samples_participating, keep_samples_extended, keep_samples_extended_variable)
      
      sample.order <- intersect(sample.order, keep_samples)
    }
    
    plot.data <- GOI.plot.current %>% 
      inner_join(a5_wts_lcpm_list[["SDHB_abdothoracic"]][GOI.plot.current$Gene,sample.order] %>% 
                   data.frame(check.names = F) %>%  
                   tibble::rownames_to_column("Gene")) %>% inner_join(ensid_to_biotype) %>% 
      data.frame(check.names = F) %>% 
      distinct()
    rownames(plot.data) <- plot.data$Symbol_Label
    
    plot.data <- plot.data %>% mutate(mki67_correlation=mki67_cor[Gene])
    
    ATRXvsTERT_sig_genes <- wts_top_tables$genosampletype$ATRX_All_vs_TERT_All %>% filter(adj.P.Val < 0.05, abs(logFC) > 0.5) %>%  pull(Gene)
    plot.data <- plot.data %>% mutate(InTERTvsATRX = ifelse(Gene %in% ATRXvsTERT_sig_genes, "Yes", "No"))
    
    plot.data <- plot.data %>%  mutate(gene_biotype = ifelse(grepl("[Pp]seudogene", gene_biotype), "pseudogene", gene_biotype))
    
    
    
    hm.annot.data.col <- a5_anno %>% 
      filter(A5_ID %in% sample.order) %>% 
      mutate(A5_ID=factor(A5_ID, levels=sample.order)) %>% 
      arrange(A5_ID) %>% 
      dplyr::select(A5_ID, TERT_ATRX_Mutation, 
                    differential_group_sampletype_strict) %>% 
      left_join(contrast_membership) %>% 
      dplyr::rename(contrast_group=group) %>% 
      mutate(snrna = A5_ID %in% snrna_samples)
    
    
    hm.annot.col = HeatmapAnnotation(
      df = hm.annot.data.col[, -1],
      which = "column",
      show_legend = TRUE,
      col = list(
        "differential_group_sampletype_strict" = sampletype_strict_cols,
        "TERT_ATRX_Mutation" = driver_cols,
        contrast_group = setNames(c("red", "black", "grey"), c(setdiff(unique(hm.annot.data.col$contrast_group),"non_participant"),"non_participant"))
      )
    )
    
    hm.annot.row = HeatmapAnnotation(df = plot.data[,c( "In_DMR","gene_biotype", "mki67_correlation", "InTERTvsATRX"), drop=F], #"source",
                                     which="row", 
                                     show_legend = TRUE,
                                     col=list(
                                       #"source"=setNames(RColorBrewer::brewer.pal(n=length(unique(plot.data$source)), name="Set3")[1:length(unique(plot.data$source))], unique(plot.data$source)),
                                       In_DMR=c("Yes"="Black", "No"="lightgrey"),
                                       gene_biotype=bio_type_cols,
                                       mki67_correlation=colorRamp2(c(-1, -0.5, -0.25, 0, 0.25, 0.5, 1), c("blue", "lightblue","white", "white", "white", "orange", "red")),
                                       InTERTvsATRX=c("Yes"="Black", "No"="lightgrey")))
    
    
    
    hm_bulk_gex <- Heatmap(t(scale(t(plot.data %>% dplyr::select(-Gene, -Symbol,-source,-gene_biotype, -In_DMR, -Symbol_Label, -mki67_correlation,-InTERTvsATRX)))),
                           row_title=paste("contrast:", contrast, " __ samples:", mode),
                           #col = col_fun,
                           #row_split = plot.data$source,
                           column_split = hm.annot.data.col$TERT_ATRX_Mutation,
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
                           row_names_gp = gpar(fontsize = 5, fontface = "italic"),
                           #heatmap_legend_param  = hm_legend_params
    )
    
    
    #############
    # Dot Plots #
    #############
    
    hm_snrna_dotplot <- HeatmapDotPlot.Seurat(object = a5_snrna,
                                              features = plot.data$Symbol, 
                                              aggr.by = "cell_type",
                                              #gene_grouping = GOI.plot.current$source,
                                              show_row_names = FALSE,
                                              cell.size = 0.15
    )
    
    hm_snrna_dotplot_genotype <- HeatmapDotPlot.Seurat(object = subset(x = a5_snrna, subset = cell_type == "Tumor"),
                                                       features = plot.data$Symbol, 
                                                       aggr.by = "TERT_ATRX_Mutation",
                                                       #gene_grouping = GOI.plot.current$source,
                                                       show_row_names = FALSE,
                                                       cell.size = 0.15
    )
    
    plot_list[[contrast]][[mode]]  <- hm_bulk_gex + hm_snrna_dotplot_genotype + hm_snrna_dotplot
    
  }
  
  
}


saveRDS(plot_list[[contrast]][[mode]], "/g/data/pq08/projects/ppgl/a5/tertiary/results/snrna_wts_crossreference/heatmaplist_20240115.rds")

plot_config <-
  as_tibble(
    matrix(
      c(
        "Metastatic_All_vs_NonMetPri_WT", 6,9.5,
        "TERT_PriMet_vs_NonMetPri_WT", 6,9.5,
        "ATRX_PriMet_vs_NonMetPri_WT", 6,9.5,
        "ATRX_All_vs_TERT_All", 6,9.5,
        "ATRX_TERT_Met_vs_NMP", 6,9.5,
        "ATRXvsTERT_ATRXvsNMP", 6,9.5,
        "ATRXvsTERT_TERTvsNMP", 4,9.5
      ),
      byrow = T,
      ncol = 3,
      dimnames = list(c(),c("contrast","height","width"))))

purrr::pwalk(.l = plot_config,
             .f = \(contrast, height, width) {
               pdf(file = file.path("/g/data/pq08/projects/ppgl/a5/tertiary/results/snrna_wts_crossreference/", paste0("snrna_wts_heatmap_",contrast, "_20240115.pdf")),
                   height = as.numeric(height),
                   width =  as.numeric(width))
               print(plot_list[[contrast]]$contrast_only_extended)
               dev.off()
             }
)
