#This script has been formatted for both sourcing and converting to markdown with knitr::spin
#'---
#' title: "Differential expression analysis of A5 WTS data"
#' author: "Aidan Flynn"
#' date: "28/04/2023"
#' ---

#+ knitr_options, echo=FALSE
#################
# Knitr options #
#################

knitr::opts_knit$set(root.dir = "/g/data/pq08/projects/ppgl/") 
knitr::opts_chunk$set(echo = FALSE, message=FALSE, warning = F)

#+ config, echo=FALSE
#Config flags to control output of plots and tables
output_qc <<- FALSE
output_tables <<- FALSE
output_plots <<- FALSE
setwd('/g/data/pq08/projects/ppgl/')

#+ dependencies, echo=F
################
# Dependencies #
################

library(tidyverse)
library(limma)
library(edgeR)
library(ggrepel)
library(patchwork)
library(umap)
library(ComplexHeatmap)
library(ggplot2)
#library(biomaRt)
#library(org.Hs.eg.db)

theme_set(theme_bw())

#+ dataloaders, echo=FALSE

#######################
# Import data loaders #
#######################

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
source("./a5/wts/scripts/data_loaders/a5_wts_dataloader.r")

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#######################
# Import data helpers #
#######################

source("./a5/sample_annotation/scripts/data_loaders/a5_contrast_groups.r")

#############
# Load data #
#############

data_loader_a5_clinical_anno("aidan.flynn@umccr-radio-lab.page", use_cache = T)

# Read in the htseq counts ----
htseq_outs <- "./a5/wts/analysis/htseq/truseq/gene"
data_loader_a5_wts_counts(count_file_dir=htseq_outs)


#get gene biotypes
EnsIds <- gsub("(ENSG[0-9]+)([.][0-9]+)?_.+","\\1", rownames(a5_wts_counts))
data_loader_get_biotypes(EnsIds, use_cache = T)
ensid_to_biotype <- data.frame(Gene=rownames(a5_wts_counts), ensembl_gene_id=EnsIds) %>% left_join(ensid_to_biotype) %>% dplyr::select(-ensembl_gene_id)

#Gene symbol lookup table
all_genes <- unique(unlist(purrr::map(.x = a5_wts_dge_list, .f = rownames)))
symbol_lookup <- setNames(all_genes, 
                          gsub("ENSG[0-9]+([.][0-9]+)?_(.+)",
                               "\\2", 
                               all_genes))

#############
# Gene sets #
#############

#+ genesets, echo=TRUE 

marker_genes <- list(
  HIF = c("VEGFA", "EPAS1", "HIF3A", "ISM1", "HK2", "SLC2A1", 
          "EGLN3", "EGLN1", "BNIP3", "PGK1", "CA12", "HK1", 
          "LDHA", "ARNT", "ABL2"),
  ECS = c("FLT1", "DLL4", "HEY1", "HEY2", "ANGPT2", "JAG1"),
  Fibroblast = c("MCAM", "IL1RAPL1", "ACTA2", "PDGFRB", "FAP", 
                 "INHBA", "POSTN", "FMO2", "ACSM3", "MFAP5", "PECAM1"),
  MP_M1 = c("CD80", "CD86", "HLA-DRA", "HLA-DRB1", "TNF", 
            "CCL2", "CCL3", "CCL4", "CCL8", "NOS2"),
  MP_M2 = c("MRC1", "CD163", "ARG1"),
  MP_Pan = c("CD68", "CD14"),
  MP_TA = c("PLXDC1", "PLAU"),
  Differentiation = c("SLC18A1","ARG2","CHGB","CHGA","TRIB3","DBH","DDC","LEF1","CARTPT","SLC18A2","TH","PENK","ALK","ASCL1","ATF3","CDH9","ERBB4","FOS","GATA3","HAND2","ISL1","JUN","JUNB","NPY","NTRK3","PHOX2A","PHOX2B","PNMT","RET","SLC6A2","SOX11","TFAP2B"),
  Catecholamin_Biosynth=c("TH", "QDPR", "DDC", "DBH", "SLC18A1", "PNMT", "SLC6A2", "MAOA", "COMT"),
  Calsina_met_markers = c("C5orf49", "DNASE1L3", "HMGB3", "RRM2", "CDT1", "TTC9", "CCNB2", "IQGAP3", "FAM83H", "CDK1", "PBK", "NETO2", "CSMD2", "PMAIP1"))


marker_genes_df <- data.frame(Symbol=unlist(marker_genes), 
                              source=rep(x = names(marker_genes), times=unlist(lapply(marker_genes, length))))

#+ TMM, echo=FALSE

#######
# TMM #
#######

#Convert the TMM counts to log2 CPMs
a5_wts_lcpm_list <- lapply(a5_wts_dge_list, cpm, log = T)

#+ MDS_header, eval=output_qc, echo=FALSE, results='asis'

#############
# MDS Plots #
#############

knitr::knit_expand(text="# MDS Plots") %>% cat()

#+ MDS_code, eval=output_qc, echo=FALSE

if(output_qc)
{ 
  mds_list <- map(a5_wts_lcpm_list, plotMDS, plot=FALSE)
  
  annotate_mds <- function (mds, annotation)
  {
    tibble(leading_logFC_dim1 = mds$x,
           leading_logFC_dim2 = mds$y,
           A5_ID = colnames(mds$distance.matrix),
           Dim = "1+2") %>% 
      left_join(annotation) %>%
      mutate(Primary_Location_Simplified=gsub("_abdominal|_thoracic|_bladder|_cardiac|_[Ll]eft|_[Rr]ight",
                                              "",
                                              Primary_Location_Simplified)) %>% 
      mutate(`sample_purity`=as.numeric(`sample_purity`)) 
  }
  
  create_dimplot_nested <- function (colour_variable,
                                     plot_prototype) {
    plot_prototype +
      geom_point(aes(colour = !!sym(colour_variable))) +
      ggtitle(colour_variable)
  }
  
  create_dimplot <- function (mds_tbl, colour_by, dim1_col_name, dim2_col_name) {
    mds_ggplot <-
      ggplot(data = mds_tbl, 
             mapping = aes(x = !!sym(dim1_col_name), 
                           y = !!sym(dim2_col_name)))
    
    plot_list <- map(
      .x = colour_by,
      .f = \(x) create_dimplot_nested(colour_variable = x,
                                      plot_prototype = mds_ggplot))
    
    plot_list[[length(plot_list) + 1]] <-
      mds_ggplot + geom_text(aes(label = A5_ID), 
                             size = 3)
    
    return(plot_list)
  }
  
  mds_tbl_list <- lapply(mds_list,
                         annotate_mds, 
                         annotation=a5_anno)
  
  # Features to annotate
  colour_by <- c("Primary_Location_Simplified", "sample_purity", 
                 "TERT_ATRX_Mutation", "differential_group_anatomy", 
                 "Gender", "differential_group_sampletype", 
                 "differential_group")
  
  
  # Create Plots
  mds_plot_list <- map(
    .x = mds_tbl_list,
    .f = \(x) create_dimplot(x, 
                             colour_by=colour_by, 
                             dim1_col_name="leading_logFC_dim1", 
                             dim2_col_name="leading_logFC_dim2")
  )
}

#+ MDS_allsamp_print_header, eval=output_qc, echo=FALSE, results='asis'
knitr::knit_expand(text="\n\n### MDS Plots with all samples\n\n") %>% cat()

#+ MDS_all_samp_print, eval=output_qc, echo=FALSE, fig.width=12, fig.height=40
if(output_qc)
{
  #pdf(file = "./a5/wts/results/plots/mds/mds_wts_allsamples.pdf", height = 30, width=15)
  Reduce(f = "+", x=mds_plot_list[["all"]] ) + plot_layout(ncol=1) 
  #dev.off()
}  

#+ MDS_sdhb_header, eval=output_qc, echo=FALSE, results='asis'
knitr::knit_expand(text="\n\n### MDS Plots with Head and Neck samples excluded\n\n") %>% cat()

#+ MDS_sdhb_print, eval=output_qc, echo=FALSE, fig.width=12, fig.height=40
if(output_qc)
{ #pdf(file = "./a5/wts/results/plots/mds/mds_wts_sdhbonly.pdf", height = 30, width=15)
  Reduce(f = "+", x=mds_plot_list[["SDHB"]] ) + plot_layout(ncol=1) 
  #dev.off()
} 

#+ MDS_sdhbat_header, eval=output_qc, echo=FALSE, results='asis'
knitr::knit_expand(text="\n\n### MDS Plots with abdo-thoracic cases only\n\n") %>% cat()

#+ MDS_sdhbat_print, eval=output_qc, echo=FALSE, fig.width=12, fig.height=40
if(output_qc)
{  
  #pdf(file = "./a5/wts/results/plots/mds/mds_wts_sdhbabdothoracic_only.pdf", height = 30, width=15)
  Reduce(f = "+", x=mds_plot_list[["SDHB_abdothoracic"]] ) + plot_layout(ncol=1) 
  #dev.off()
}

#+ umap_header, eval=output_qc, echo=FALSE, results='asis'

##############
# UMAP Plots #
##############

knitr::knit_expand(text="\n\n### UMAPS\n\n") %>% cat()

#+ umap_code, eval=output_qc, echo=FALSE, fig.width=12, fig.height=40
if(output_qc)
{ 
  set.seed(10)
  umap_config <- umap.defaults
  umap_config$n_neighbors=10
  #umap_config$spread=3
  umap_list <- map(a5_wts_lcpm_list, \(lcpm) umap(t(lcpm), config = umap_config)) 
  
  annotate_umap <- function (umap_result, annotation)
  {
    umap_layout <- data.frame(umap_result$layout)
    return(
      tibble(UMAP1 = umap_layout$X1,
             UMAP2 = umap_layout$X2,
             A5_ID = rownames(umap_layout),
             Dim = "1+2") %>% 
        left_join(annotation) %>%
        mutate(Primary_Location_Simplified=gsub("_abdominal|_thoracic|_bladder|_cardiac|_[Ll]eft|_[Rr]ight",
                                                "",
                                                Primary_Location_Simplified)) %>% 
        mutate(`sample_purity`=as.numeric(`sample_purity`)) 
    )
  }
  
  #Apply annotation
  umap_annotated_list <- map(umap_list, \(x) annotate_umap(umap_result = x, annotation=a5_anno))
  
  # write.table(umap_annotated_list$All, "./a5/wts/results/differential_gene_expression/tables/umap_coord_allsample_seed42_nn10_hg38_tmm_log2cpm.tsv", sep="\t", row.names=F)
  # write.table(umap_annotated_list$NoHN, "./a5/wts/results/differential_gene_expression/tables/umap_coord_nohn_seed42_nn10_hg38_tmm_log2cpm.tsv", sep="\t", row.names=F)
  # write.table(umap_annotated_list$NoHN_Core, "./a5/wts/results/differential_gene_expression/tables/umap_coord_noexnohn_seed42_nn10_hg38_tmm_log2cpm.tsv", sep="\t", row.names=F)
  
  colour_by <- c("Primary_Location_Simplified", "sample_purity", "TERT_ATRX_Mutation", 
                 "differential_group_anatomy", "Gender", "differential_group_sampletype", 
                 "differential_group")
  
  
  # Create Plots
  umap_plot_list <- map(
    .x = umap_annotated_list,
    .f = \(x) create_dimplot(x, 
                             colour_by=colour_by, 
                             dim1_col_name="UMAP1", 
                             dim2_col_name="UMAP2")
  )
  
}

#+ umap_allsamp_header, eval=output_qc, echo=FALSE, results='asis'
knitr::knit_expand(text="\n\n### UMAP - All Samples\n\n") %>% cat()

#+ umap_allsamp_print, eval=output_qc, echo=FALSE, fig.width=12, fig.height=40
if(output_qc)
{  
  #pdf(file = "./a5/wts/results/plots/umap/umap_wts_allsamples.pdf", height = 30, width=15)
  Reduce(f = "+", x=umap_plot_list[["all"]] ) + plot_layout(ncol=1) 
  #dev.off()
}

#+ umap_sdhb_header, eval=output_qc, echo=FALSE, results='asis'
knitr::knit_expand(text="\n\n### UMAP - Head and neck and abdothoracic\n\n") %>% cat()

#+ umap_sdhb_print, eval=output_qc, echo=FALSE, fig.width=12, fig.height=40
if(output_qc)
{  
  #pdf(file = "./a5/wts/results/plots/umap/umap_wts_sdhbonly_abdothoracic_hn.pdf", height = 30, width=15)
  Reduce(f = "+", x=umap_plot_list[["SDHB"]] ) + plot_layout(ncol=1) 
  #dev.off()
}

#+ umap_at_header, eval=output_qc, echo=FALSE, results='asis'
knitr::knit_expand(text="\n\n### UMAP - Abdothoracic Only\n\n") %>% cat()

#+ umap_at_print, eval=output_qc, echo=FALSE, results='asis', fig.width=12, fig.height=40
if(output_qc)
{  
  #pdf(file = "./a5/wts/results/plots/umap/umap_wts_sdhbonly_abdothoracic.pdf", height = 30, width=15)
  Reduce(f = "+", x=umap_plot_list[["SDHB_abdothoracic"]] ) + plot_layout(ncol=1) 
  #dev.off()
}

############################
############################
#                          #
# Differential Expression  #
#                          #
############################
############################

#' # Differential Expression

#lists to store DE limma fit objects and top tables
wts_de_fits <- list()
wts_top_tables <- list()


#####################
# Contrast Matrices #
#####################

#Generates contrast_matrix_hn and design_matrix_hn
make_hn_vs_abdominothoracic_contrasts(sample_anno = a5_anno %>% 
                                        filter(A5_ID %in% colnames(a5_wts_lcpm_list$SDHB)))

#Generates contrast_matrix_genosampletype and design_matrix_genosampletype
make_genotype_sampletype_contrasts(sample_anno = a5_anno %>% 
                                     filter(A5_ID %in% colnames(a5_wts_lcpm_list$SDHB_abdothoracic)))


#########################
# Head/Neck Vs AdboThor #
#########################

#' ## Differential expression between head and neck and abdominothoracic

#' Here we perform differential expression between those samples marked abdominothoracic and those marked head and neck in there clinical annotation.  Samples with UMAP clustering incongruent with their clinical annotation or where tumour location was not specified were excluded from contrasts (but retained within dataset for LIMMAs variance/dispersion calculations).  
#' Excluded samples: 
{{toString(a5_anno %>% filter(differential_group_anatomy=="Ambiguous") %>% pull(A5_ID))}}
#' 

count_contrast_members(contrast_matrix_hn, design_matrix_hn) %>% 
  knitr::kable(caption = "Contrast group member counts")

#Create working verson of counts and sample annotation
counts_dge <- a5_wts_dge_list[["SDHB"]]
counts_anno <- a5_anno %>% filter(A5_ID %in% colnames(a5_wts_dge_list[["SDHB"]]))
counts_anno <- counts_anno[match(colnames(counts_dge), counts_anno$A5_ID),]

#Run voom
v <- voom(counts_dge, design_matrix_hn, plot=TRUE)

# estimate correlation between samples from the same patient
file_hashes <- purrr::map(list(a5_wts_dge_list[["SDHB"]]$counts, 
                               design_matrix_hn, 
                               counts_anno), 
                          digest::digest, algo = "md5")
suffix <- paste(purrr::map(file_hashes, stringr::str_sub, start=25, end=32), collapse = "_")
checkpoint_file <- paste0("./a5/wts/quickload_checkpoints/dupcor_",suffix,".rds")
if(file.exists(checkpoint_file))
{
  dupcor <- readRDS(checkpoint_file)
} else {
  
  message("Checkpoint file not found. Running duplicateCorrelation ... ")
  dupcor <- duplicateCorrelation(v, 
                                 design_matrix_hn, 
                                 block=counts_anno$`Patient ID`)
  saveRDS(object = dupcor, 
          file = checkpoint_file)
  
}

# differential expression analysis, keeping track of samples that are correlated as they come from the same patient
vfit <- lmFit(v, design_matrix_hn, block=counts_anno$`Patient ID`, correlation=dupcor$consensus)
vfit <- contrasts.fit(vfit, contrasts=contrast_matrix_hn)
efit <- eBayes(vfit)

plotSA(efit, main="Final model: Mean-variance trend")

#Store fit
wts_de_fits[["Parasympathetic_vs_Sympathetic"]] <- efit

# get top tables 
wts_top_tables[["Parasympathetic_vs_Sympathetic"]] <- list()
for (contrast in dimnames(contrast_matrix_hn)$Contrasts)
{
  wts_top_tables[["Parasympathetic_vs_Sympathetic"]][[contrast]] <- 
    topTable(efit, 
             coef = contrast, 
             number = Inf, 
             sort.by = "P") %>% 
    rownames_to_column(var = "Gene") %>% 
    mutate(gene_symbol = if_else(grepl(x = Gene, pattern = "_"), 
                                 str_extract(Gene, pattern = '[^_]*$'),
                                 Gene)) %>% 
    left_join(ensid_to_biotype %>% dplyr::select(Gene, gene_biotype), by=c("Gene")) %>% 
    dplyr::relocate(gene_biotype, .after=Gene)
}



#+ de_hnabdo_resulttable_header, eval=output_tables, echo=FALSE, results='asis'  
knitr::knit_expand(text="\n\n# DE Result Tables\n\n") %>% cat()

#########
# Result Tables 
#########

#+ de_hnabdo_resulttable_summary_code, eval=output_tables, echo=FALSE
if(output_tables)
{
  
  #DE summary
  decideTests(efit, lfc = 0.5, 
              adjust.method = "BH", 
              p.value = 0.05) %>% 
    summary() %>% 
    knitr::kable(caption = "DE gene counts at adj.p<0.05 and LFC>0.5")
}

#+ de_hnabdo_summary_header, eval=output_tables, echo=FALSE, results='asis'  
knitr::knit_expand(text="\n\n# Top 1000 DE genes - HN vs abdo/thoracic\n\n") %>% cat()

#+ de_hnabdo_summary_print, eval=output_tables, echo=FALSE
if(output_tables)
{
  DT::datatable(data = wts_top_tables[["Parasympathetic_vs_Sympathetic"]][["Parasympathetic_vs_Sympathetic"]] %>% 
                  filter(abs(logFC)>1, 
                         adj.P.Val<0.05) %>% 
                  slice_min(n = 1000, 
                            order_by = adj.P.Val) %>% 
                  mutate(across(.cols = c(logFC, AveExpr, t, 
                                          P.Value, adj.P.Val, B), 
                                round, 4)), 
                options = list(pageLength = 20, scrollX = T))
  
}

#+ de_hnabdo_volcano_header, eval=output_plots, echo=FALSE, results='asis' 

#########
# Result Plots 
#########

  #########
  # Volcano Plot
  #########

knitr::knit_expand(text="\n\n# Volcano Plot - H/N vs abdo/thoracic\n\n") %>% cat()

#+ de_hnabdo_volcano_code, eval=output_plots, fig.width=7, fig.height=5
if(output_plots){
  plot_volcano <- function(tt, n_label){
    
    genes_to_label <- tt %>%
      slice_min(adj.P.Val, n = n_label) %>% 
      pull(Gene)
    
    tt_plot <- tt %>% #flag the most significant genes to label
      mutate(significant = if_else(adj.P.Val < 0.05, "Adj. P Value < 0.05", "not significant")) %>% 
      mutate(label = if_else(Gene %in% genes_to_label, TRUE, FALSE))
    
    volcano <- ggplot(tt_plot, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(col=significant)) +
      scale_color_manual(values=c("red", "black")) +
      geom_text_repel(
        data = . %>% filter(label == TRUE),
        box.padding = unit(0.5, "lines"),
        point.padding = unit(0.1, "lines"),
        segment.color = "grey",
        max.overlaps = Inf,
        size = 3,
        aes(label=gene_symbol))
    
    return(volcano)
    
  }
  
  plot_volcano(wts_top_tables[["Parasympathetic_vs_Sympathetic"]][["Parasympathetic_vs_Sympathetic"]], 60) + ggtitle("Head and Neck vs Abdo-thoracic")
  
}


#+ de_hnabdo_boxplots_header, eval=output_plots, echo=FALSE, results='asis' 


  #########
  # Box Plots
  #########
  
  knitr::knit_expand(text="\n\n# Boxplots of Top 20 DE genes - H/N vs abdo/thoracic\n\n") %>% cat()
  
#+ de_hnabdo_boxplots_code, eval=output_plots, echo=FALSE
  if(output_plots){
  gene_list <- wts_top_tables[["Parasympathetic_vs_Sympathetic"]][["Parasympathetic_vs_Sympathetic"]] %>% 
    filter(Gene %in% (ensid_to_biotype %>% 
                        filter(gene_biotype=="protein_coding") %>% 
                        pull(Gene))) %>% 
    filter(abs(logFC)>1, adj.P.Val<0.05) %>% 
    slice_min(n = 20, order_by = adj.P.Val) %>% 
    pull(Gene)
  
  extra_genes <- c(marker_genes$Catecholamin_Biosynth[c(1,4,6,7)], "LEF1", "CHGA")
  
  gene_list <- c(gene_list, symbol_lookup[extra_genes])
  
  plot_data <- a5_wts_lcpm_list[["SDHB"]][gene_list,,drop=F] %>% 
    data.frame(check.names = F) %>% 
    tibble::rownames_to_column("Gene") %>% 
    pivot_longer(cols=-Gene, 
                 names_to="A5_ID", 
                 values_to="log2_cpm") %>% 
    inner_join(a5_anno %>% 
                 dplyr::select(A5_ID, Major_Cluster, 
                               Primary_Location_Simplified)) %>% 
    mutate(Primary_Location_Simplified=gsub("_abdominal|_thoracic|_bladder|_cardiac|_[Ll]eft|_[Rr]ight",
                                            "",
                                            Primary_Location_Simplified),
           Major_Cluster_simple=gsub(" [(].+|Parasympathetic/",
                                     "",
                                     Major_Cluster),
           Gene=gsub("ENSG[0-9.]+_","",Gene))
  #%>% filter(Primary_Location_Simplified != "Unspecified")
  }
#+ de_hnabdo_boxplots_print, eval=output_plots, fig.height=21, fig.width=18
if(output_plots){
  ggplot(plot_data, aes(x=Primary_Location_Simplified, y=log2_cpm)) + 
    geom_boxplot(outlier.alpha = 0) + 
    #geom_jitter(mapping=aes(color=Concordant), width = 0.1) + 
    geom_jitter(width = 0.1, size=0.5) + 
    scale_color_manual(values=c("Black", "Red", "Orange")) +
    facet_wrap("Gene", scales="free_y") + 
    theme_bw() + 
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
}

#+ de_hnabdo_heatmap_header, eval=output_plots, echo=FALSE, results='asis' 

  #########
  # Heatmap
  #########
  

  knitr::knit_expand(text="\n\n# Annotated heatmap of top 100 genes - H/N vs abdo/thoracic\n\n") %>% cat()    
  
#+ de_hnabdo_heatmap_print, eval=output_plots, echo=FALSE, fig.height=20, fig.width=15
  
if(output_plots){
  genes_per_group <- 100
  
  GOI <- wts_top_tables[["Parasympathetic_vs_Sympathetic"]][["Parasympathetic_vs_Sympathetic"]] %>% filter(adj.P.Val < 0.05) %>% arrange(desc(abs(logFC))) %>% slice_head(n = genes_per_group) %>% dplyr::pull(Gene)
  
  extra_genes <- c(marker_genes$Catecholamin_Biosynth[c(1,4,6,7)], "LEF1", "CHGA")
  
  GOI <- c(GOI, symbol_lookup[extra_genes])
  
  plot.data <- a5_wts_lcpm_list[["SDHB"]][GOI,] %>% 
    data.frame(check.names = F) %>% 
    tibble::rownames_to_column("Gene") %>% 
    inner_join(ensid_to_biotype)
  rownames(plot.data) <-plot.data$Gene
  
  sample.order <- colnames(plot.data)[grepl("^E[0-9]+", colnames(plot.data))]
  
  hm.annot.data.col <- a5_anno %>% 
    filter(A5_ID %in% sample.order) %>% 
    dplyr::select(A5_ID, Primary_Location_Simplified, 
                  TERT_ATRX_Mutation, differential_group_purity, 
                  differential_group_anatomy, Catecholamine_profile, 
                  differential_group_sampletype) %>%
    mutate(Primary_Location_Simplified=gsub("_abdominal|_thoracic|_bladder|_cardiac|_[Ll]eft|_[Rr]ight",
                                            "",
                                            Primary_Location_Simplified))
  
  hm.annot.data.col <- hm.annot.data.col[match(sample.order, hm.annot.data.col$A5_ID),]
  
  hm.annot.col = HeatmapAnnotation(
    df = hm.annot.data.col[, -1],
    which = "column",
    show_legend = TRUE,
    col = list(
      "differential_group_sampletype" = specimen_type_cols,
      "differential_group_anatomy" = setNames(
        RColorBrewer::brewer.pal(
          n = length(unique(hm.annot.data.col$differential_group_anatomy)),
          name = "Set1"),unique(hm.annot.data.col$differential_group_anatomy)),
      "TERT_ATRX_Mutation" = driver_cols,
      "Primary_Location_Simplified" = setNames(
        RColorBrewer::brewer.pal(
          n = length(unique(hm.annot.data.col$Primary_Location_Simplified)),
          name = "Set1"),
        unique(hm.annot.data.col$Primary_Location_Simplified)),
      "Catecholamine_profile" = setNames(
        RColorBrewer::brewer.pal(
          n = length(unique(hm.annot.data.col$Catecholamine_profile)),
          name = "Set1"),
        unique(hm.annot.data.col$Catecholamine_profile))
    ))
  
  plot.data.EnsIds <- gsub("(ENSG[0-9]+)([.][0-9]+)?_.+","\\1", rownames(plot.data))
  data_loader_ensgid_to_chr(plot.data.EnsIds, use_cache=TRUE)
  ensgid_to_chr <- data.frame(Gene=rownames(plot.data), ensembl_gene_id=plot.data.EnsIds) %>% left_join(ensgid_to_chr) %>% dplyr::select(-ensembl_gene_id, -start_position, -end_position)
  plot.data <- plot.data %>% left_join(ensgid_to_chr)
  plot.data$chromosome_name <- paste0("chr", plot.data$chromosome_name)
  
  chrom_cols<-c(RColorBrewer::brewer.pal(n=8, name="Set1"),RColorBrewer::brewer.pal(n=8, name="Set2"),RColorBrewer::brewer.pal(n=length(unique(plot.data$chromosome_name))-16, name="Set3"))
  
  if(length(unique(plot.data$gene_biotype)) < 13)
  {
    bio_type_cols <- setNames(
      RColorBrewer::brewer.pal(n = length(unique(
        plot.data$gene_biotype
      )), name = "Set3"),
      unique(plot.data$gene_biotype)
    )
  } else {
    bio_type_cols <- setNames(
      c(
        RColorBrewer::brewer.pal(n = 12, name = "Set3"),
        RColorBrewer::brewer.pal(n = (length(
          unique(plot.data$gene_biotype)
        ) - 12), name = "Set1")
      ),
      unique(plot.data$gene_biotype)
    )
  } 
  
  hm.annot.row = HeatmapAnnotation(df = plot.data[,c("gene_biotype","chromosome_name")], 
                                   which="row", 
                                   show_legend = TRUE,
                                   col=list("gene_biotype"=bio_type_cols,
                                            "chromosome_name"=setNames(
                                              chrom_cols,
                                              unique(plot.data$chromosome_name)))
  )
  
  rownames(plot.data) <-  gsub("ENSG[0-9]+([.][0-9]+)?_(.+)","\\2", plot.data$Gene)
  
  hm <- Heatmap(t(scale(t(plot.data %>% dplyr::select(-Gene,-gene_biotype, -chromosome_name)))),
                #col = col_fun,
                #row_split = plot.data$source,
                #column_split = a5_anno.wts.nohn_noex$genotype_groups,
                #width = ncol(bulk_heatmap_deg)*unit(2, "mm"),
                #height = ncol(bulk_heatmap_deg)*unit(4, "mm"),
                row_gap = unit(2, "mm"),
                column_gap = unit(2, "mm"),
                row_title_rot = 0,
                top_annotation = hm.annot.col,
                left_annotation = hm.annot.row,
                column_title_rot = 90,
                cluster_columns = TRUE,
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
  
  draw(hm, heatmap_legend_side="right", annotation_legend_side="bottom",
       legend_grouping = "original", padding = unit(c(2, 2, 2, 40), "mm"))
}

###########################
# Genotype/Sample-type DE #
###########################

#' ## Differential Expression - ATRX/TERT/Met

#+ AT_DE_header, echo=FALSE, results='asis'

knitr::knit_expand(text="Make a design matrix with biological groups:\n\n* {{paste(differential_group_levels, collapse='  \n* ')}}") %>% cat()
#'

#' Andrew Note: To account for sample-specific effects, the patient ID was used as a blocking variable (as suggested here: <https://stat.ethz.ch/pipermail/bioconductor/2013-July/054087.html>

#+ AT_DE_contrasts, echo=FALSE

count_contrast_members(contrast_matrix_genosampletype, design_matrix_genosampletype) %>% 
  knitr::kable(caption = "Contrast group member counts")

useful_contrasts <- c("TERT_PriMet_vs_NonMetPri_WT",
                      "ATRX_PriMet_vs_NonMetPri_WT",
                      "Metastasis_All_vs_NonMetPri_WT")

#+ AT_DE_contrasts_print, echo=FALSE

#' Using contrast groups:
#' 
#' * 
{{paste(useful_contrasts, collapse='  \n* ')}}
#'

#+ AT_DE_code, echo=FALSE

counts_dge <- a5_wts_dge_list$SDHB_abdothoracic

counts_anno <- a5_anno %>% filter(A5_ID %in% colnames(a5_wts_dge_list[["SDHB_abdothoracic"]]))
counts_anno <- counts_anno[match(colnames(counts_dge), counts_anno$A5_ID),]

# voom transform

par(mfrow=c(1,2))
v <- voom(counts_dge, design_matrix_genosampletype, plot=TRUE)

# estimate correlation between samples from the same patient  
file_hashes <- purrr::map(list(a5_wts_dge_list[["SDHB_abdothoracic"]]$counts, 
                               design_matrix_genosampletype, 
                               counts_anno), 
                          digest::digest, algo = "md5")
suffix <- paste(purrr::map(file_hashes, stringr::str_sub, start=25, end=32), collapse = "_")
checkpoint_file <- paste0("./a5/wts/quickload_checkpoints/dupcor_",suffix,".rds")
if(file.exists(checkpoint_file))
{
  dupcor <- readRDS(file = checkpoint_file)
} else {
  dupcor <- duplicateCorrelation(v, design_matrix_genosampletype, block=counts_anno$`Patient ID`)
  saveRDS(object = dupcor, file = checkpoint_file)
}


# differential expression analysis, keeping track of samples that are correlated as they come from the same patient
vfit <- lmFit(v, 
              design_matrix_genosampletype, 
              block=counts_anno$`Patient ID`, 
              correlation=dupcor$consensus)
vfit <- contrasts.fit(vfit, contrasts=contrast_matrix_genosampletype)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

wts_de_fits[["genosampletype"]] <- efit

sum.fit <- decideTests(efit, lfc = 0.5, adjust.method = "BH", p.value = 0.05)
summary(sum.fit) %>% knitr::kable(caption = "DE gene counts at adj.p<0.05 and LFC>0.5")

wts_top_tables[["genosampletype"]] <- list()
for (contrast in dimnames(contrast_matrix_genosampletype)$Contrasts)
{
  # get top ATRX signature genes 
  wts_top_tables[["genosampletype"]][[contrast]] <- 
    topTable(efit, 
             coef = contrast, 
             number = Inf, 
             sort.by = "P") %>% 
    rownames_to_column(var = "Gene") %>% 
    mutate(gene_symbol = if_else(grepl(x = Gene, pattern = "_"), 
                                 str_extract(Gene, pattern = '[^_]*$'),
                                 Gene)) %>% 
    left_join(ensid_to_biotype %>% dplyr::select(Gene, gene_biotype), by=c("Gene")) %>% 
    dplyr::relocate(gene_biotype, .after=Gene)
}
# Get the top TERT signature genes 


# #write_wts_top_tables_genotype
# 
# for (contrast in dimnames(contrast_matrix_genosampletype)$Contrasts)
# {
#   wts_top_tables[["genosampletype"]][[contrast]] %>%
#     write_tsv(paste0("./a5/wts/results/differential_gene_expression/",
#                      contrast, 
#                      "_tt_dge_nohn_noex_hg38.tsv"))
# }



#+ AT_DE_toptable_print_header, eval=output_tables, echo=FALSE, results='asis'

#########
# Result Tables 
#########

knitr::knit_expand(text="\n\n### Top 1000 DE - All Metastases vs. All Non-Metastatic Primaries\n\n") %>% cat()

knitr::knit_expand(text="\nFiltered: Absolute-logFC >1, adj.P.Val < 0.05") %>% cat()

#+ AT_DE_toptable_print_code, eval=output_tables, echo=FALSE
if(output_tables)
{
  contrast = "Metastasis_All_vs_NonMetPri_WT"
  DT::datatable(data = wts_top_tables[["genosampletype"]][[contrast]] %>% filter(abs(logFC)>1, adj.P.Val<0.05) %>% slice_min(n = 1000, order_by = adj.P.Val) %>% 
                  mutate(across(.cols = c(logFC, AveExpr, t, 
                                          P.Value, adj.P.Val, B), round, 4)), 
                options = list(pageLength = 20, scrollX = T))
}

#+ AT_DE_ATRXTERT_toptables_print_header, eval=output_tables, echo=FALSE, results='asis'

knitr::knit_expand(text="\n\n### Top 1000 DE - All TERT mutants vs all ATRX mutants\n\n") %>% cat()

knitr::knit_expand(text="\nFiltered: Absolute-logFC >1, adj.P.Val < 0.05") %>% cat()

#+ AT_DE_ATRXTERT_toptables_print, eval=output_tables, echo=FALSE
if(output_tables){
  contrast = "ATRX_All_vs_TERT_All"
  DT::datatable(data = wts_top_tables[["genosampletype"]][[contrast]] %>% filter(abs(logFC)>1, adj.P.Val<0.05) %>% slice_min(n = 1000, order_by = adj.P.Val) %>% 
                  mutate(across(.cols = c(logFC, AveExpr, t, 
                                          P.Value, adj.P.Val, B), round, 4)), 
                options = list(pageLength = 20, scrollX = T))
}
#+ AT_DE_ATRXNonmet_toptables_header, eval=output_tables, echo=FALSE, results='asis'

knitr::knit_expand(text="\n\n### Top 1000 DE - All ATRX vs. All Non-Metastatic Primaries\n\n") %>% cat()

knitr::knit_expand(text="\nFiltered: Absolute-logFC >1, adj.P.Val < 0.05") %>% cat()

#+ AT_DE_ATRXNonmet_toptables_print, eval=output_tables, echo=FALSE
if(output_tables){
  contrast = "ATRX_PriMet_vs_NonMetPri_WT"
  DT::datatable(data = wts_top_tables[["genosampletype"]][[contrast]] %>% filter(abs(logFC)>1, adj.P.Val<0.05) %>% slice_min(n = 1000, order_by = adj.P.Val) %>% 
                  mutate(across(.cols = c(logFC, AveExpr, t, 
                                          P.Value, adj.P.Val, B), round, 4)), 
                options = list(pageLength = 20, scrollX = T))
}
#+ AT_DE_TERTNonmet_toptables_header, eval=output_tables, echo=FALSE, results='asis'

knitr::knit_expand(text="\n\n### Top 1000 DE - All TERT mutatants vs. All Non-Metastatic Primaries\n\n") %>% cat()

knitr::knit_expand(text="\nFiltered: Absolute-logFC >1, adj.P.Val < 0.05") %>% cat()

#+ AT_DE_TERTNonmet_toptables_print, eval=output_tables, echo=FALSE
if(output_tables){
  contrast = "TERT_PriMet_vs_NonMetPri_WT"
  DT::datatable(data = wts_top_tables[["genosampletype"]][[contrast]] %>% filter(abs(logFC)>1, adj.P.Val<0.05) %>% slice_min(n = 1000, order_by = adj.P.Val) %>% 
                  mutate(across(.cols = c(logFC, AveExpr, t, 
                                          P.Value, adj.P.Val, B), round, 4)), 
                options = list(pageLength = 20, scrollX = T))
}

#+ AT_DE_TERTATRX_venn_header, eval=output_plots, echo=FALSE, results='asis'

#########
# Result Plots 
#########

knitr::knit_expand(text="\n\n### TERT and ATRX vs Non-Metastatic Primaries significant gene overlap\n\n") %>% cat()

#+ AT_DE_TERTATRX_venn_print_header, eval=output_plots, echo=FALSE
if(output_plots) {
  vennDiagram(sum.fit[,c("ATRX_PriMet_vs_NonMetPri_WT","TERT_PriMet_vs_NonMetPri_WT")], circle.col=c("turquoise", "salmon"))
}
#+ AT_DE_volcano_header, eval=output_plots, echo=FALSE, results='asis'

knitr::knit_expand(text="\n\n### Volcano plot\n\n") %>% cat() 

knitr::knit_expand(text="\nTop 20 genes labelled (sorted by adj.P.Val)") %>% cat()

#+ AT_DE_volcano_code, eval=output_plots, fig.width=14, fig.height=25
if(output_plots){
  plot_volcano <- function(tt, n_label){
    
    genes_to_label <- tt %>%
      slice_min(adj.P.Val, n = n_label) %>% 
      pull(Gene)
    
    tt_plot <- tt %>% #flag the most significant genes to label
      mutate(significant = if_else(adj.P.Val < 0.05, "Adj. P Value < 0.05", "not significant")) %>% 
      mutate(label = if_else(Gene %in% genes_to_label, TRUE, FALSE))
    
    volcano <- ggplot(tt_plot, aes(logFC, -log10(adj.P.Val))) +
      geom_point(aes(col=significant)) +
      scale_color_manual(values=c("red", "black")) +
      geom_text_repel(
        data = . %>% filter(label == TRUE),
        box.padding = unit(0.5, "lines"),
        point.padding = unit(0.1, "lines"),
        segment.color = "grey",
        max.overlaps = Inf,
        size = 3,
        aes(label=gene_symbol))
    
    return(volcano)
    
  }
  
  tert <- plot_volcano(wts_top_tables[["genosampletype"]][["TERT_PriMet_vs_NonMetPri_WT"]], 20) + ggtitle("TERT Pri/Met Vs WT Non-metastatic primaries")
  atrx <- plot_volcano(wts_top_tables[["genosampletype"]][["ATRX_PriMet_vs_NonMetPri_WT"]], 20) + ggtitle("ATRX Pri/Met Vs WT Non-metastatic primaries")
  tert_vs_atrx <- plot_volcano(wts_top_tables[["genosampletype"]][["ATRX_All_vs_TERT_All"]], 20) + ggtitle("TERT Vs ATRX")
  met_vs_nonmet <- plot_volcano(wts_top_tables[["genosampletype"]][["Metastasis_All_vs_NonMetPri_WT"]], 20) + ggtitle("All metastases vs non-metastatic primaries")
  
  tert + atrx + tert_vs_atrx + met_vs_nonmet + plot_layout(ncol=1)
}

#+ AT_DE_BoxPlots_header, eval=output_plots, echo=FALSE, results='asis'

knitr::knit_expand(text="\n\n### Box plots of top genes\n\n") %>% cat()

knitr::knit_expand(text="\nGenes filtered for 'adj.P.Val < 0.05' and sorted by Absolute-logFC. Top 10 from each contrast shown.") %>% cat() 

#+ AT_DE_BoxPlots_print, eval=output_plots, echo=FALSE, fig.height=21, fig.width=18
if(output_plots){
  tert_top <- wts_top_tables[["genosampletype"]][["TERT_PriMet_vs_NonMetPri_WT"]] %>% 
    filter(adj.P.Val < 0.05) %>% 
    arrange(desc(abs(logFC))) %>%  
    slice_head(n = 50) %>% pull(Gene)
  
  atrx_top <- wts_top_tables[["genosampletype"]][["ATRX_PriMet_vs_NonMetPri_WT"]] %>% 
    filter(adj.P.Val < 0.05) %>% 
    arrange(desc(abs(logFC))) %>%  
    slice_head(n = 50) %>% 
    pull(Gene)
  
  met_top <- wts_top_tables[["genosampletype"]][["Metastasis_All_vs_NonMetPri_WT"]] %>% 
    filter(adj.P.Val < 0.05) %>% 
    arrange(desc(abs(logFC))) %>%  
    slice_head(n = 50) %>% 
    pull(Gene)
  
  GOI <- c(met_top[1:10], tert_top[1:10], atrx_top[1:10])
  
  top_logcpm <- a5_wts_lcpm_list[["SDHB_abdothoracic"]][GOI,, drop=F] %>% 
    data.frame(check.names = F) %>% 
    tibble::rownames_to_column(var = "Gene") %>% 
    pivot_longer(cols=-Gene, names_to = "A5_ID", values_to = "log2cpm") %>% 
    left_join(ensid_to_biotype)
  
  top_logcpm <- top_logcpm %>% 
    left_join(a5_anno %>% 
                dplyr::select(A5_ID, is_primary_or_met,
                              TERT_ATRX_Mutation, differential_group_purity, 
                              differential_group_sampletype, differential_group))
  
  top_logcpm$differential_group <- factor(top_logcpm$differential_group, levels=differential_group_levels)
  
  plot.data <- top_logcpm %>% filter(gene_biotype=="protein_coding") 
  
  ggplot(mapping=aes(x=paste(TERT_ATRX_Mutation, differential_group_sampletype, sep="_"), 
                     y=log2cpm)) + 
    geom_boxplot(data=plot.data %>% filter(differential_group_purity=="PurityOK")) + 
    geom_jitter(data=plot.data, 
                mapping = aes(color=differential_group,
                              alpha=differential_group_purity), 
                width=0.2)  + 
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) + 
    scale_color_manual(values = differential_group_colors) +
    scale_alpha_manual(values = c(0.2,1)) + 
    facet_wrap("Gene", scale="free_y")
}

#+ AT_DE_heatmap_header, eval=output_plots, echo=FALSE, results='asis'

knitr::knit_expand(text="\n\n### Heatmap of top genes\n\n") %>% cat()

knitr::knit_expand(text="\nGenes filtered for 'adj.P.Val < 0.05' and sorted by Absolute-logFC. Top 50 from each contrast shown. Some groups have less because duplicate genes are removed.") %>% cat() 

#+ AT_DE_heatmap_print, eval=output_plots, echo=FALSE, fig.height=18, fig.width=18
if(output_plots){
  genes_per_group <- 50
  
  GOI <- bind_rows(
    wts_top_tables[["genosampletype"]][["TERT_PriMet_vs_NonMetPri_WT"]] %>% 
      filter(adj.P.Val < 0.05) %>% 
      arrange(desc(abs(logFC))) %>%  
      slice_head(n = genes_per_group) %>%  
      mutate(source="TERTvsNonMet") %>% 
      dplyr::select(source, Gene),
    wts_top_tables[["genosampletype"]][["ATRX_PriMet_vs_NonMetPri_WT"]] %>% 
      filter(adj.P.Val < 0.05) %>% 
      arrange(desc(abs(logFC))) %>% 
      slice_head(n = genes_per_group) %>% 
      mutate(source="ATRXvsNonMet") %>% 
      dplyr::select(source, Gene),
    wts_top_tables[["genosampletype"]][["Metastasis_All_vs_NonMetPri_WT"]] %>% 
      filter(adj.P.Val < 0.05) %>% 
      arrange(desc(abs(logFC))) %>% 
      slice_head(n = genes_per_group) %>% 
      mutate(source="MetvsNonMet") %>% dplyr::select(source, Gene)
  )
  GOI <- GOI[!duplicated(GOI$Gene),]
  
  plot.data <- GOI %>% 
    inner_join(a5_wts_lcpm_list[["SDHB_abdothoracic"]][GOI$Gene,] %>% 
                 data.frame(check.names = F) %>%  
                 tibble::rownames_to_column("Gene")) %>% inner_join(ensid_to_biotype)
  rownames(plot.data) <-plot.data$Gene
  
  sample.order <- colnames(plot.data)[grepl("^E[0-9]+", colnames(plot.data))]
  hm.annot.data.col <- a5_anno %>% 
    filter(A5_ID %in% sample.order) %>% 
    dplyr::select(A5_ID, TERT_ATRX_Mutation, 
                  differential_group_purity, 
                  differential_group_sampletype)
  
  hm.annot.data.col <- hm.annot.data.col[match(sample.order, hm.annot.data.col$A5_ID),]
  
  hm.annot.col = HeatmapAnnotation(
    df = hm.annot.data.col[, -1],
    which = "column",
    show_legend = TRUE,
    col = list(
      "differential_group_sampletype" = specimen_type_cols,
      "TERT_ATRX_Mutation" = driver_cols,
      "differential_group_purity" = setNames(c("black", "grey"),
                                             unique(hm.annot.data.col$differential_group_purity))
    )
  )
  
  plot.data.EnsIds <- gsub("(ENSG[0-9]+)([.][0-9]+)?_.+","\\1",
                           rownames(plot.data))
  data_loader_ensgid_to_chr(plot.data.EnsIds, use_cache=TRUE)
  ensgid_to_chr <- data.frame(Gene=rownames(plot.data), ensembl_gene_id=plot.data.EnsIds) %>% left_join(ensgid_to_chr) %>% dplyr::select(-ensembl_gene_id, -start_position, -end_position)
  plot.data <- plot.data %>% left_join(ensgid_to_chr)
  plot.data$chromosome_name <- paste0("chr", plot.data$chromosome_name)
  
  chrom_cols<-c(RColorBrewer::brewer.pal(n=8, name="Set1"),RColorBrewer::brewer.pal(n=8, name="Set2"),RColorBrewer::brewer.pal(n=length(unique(plot.data$chromosome_name))-16, name="Set3"))
  
  hm.annot.row = HeatmapAnnotation(df = plot.data[,c("source","gene_biotype","chromosome_name")], 
                                   which="row", 
                                   show_legend = TRUE,
                                   col=list("source"=setNames(
                                     RColorBrewer::brewer.pal(n=length(unique(plot.data$source)), name="Set3"),
                                     unique(plot.data$source)), "gene_biotype"=setNames(
                                       RColorBrewer::brewer.pal(n=length(unique(plot.data$gene_biotype)), name="Set1"),
                                       unique(plot.data$gene_biotype)),
                                     "chromosome_name"=setNames(
                                       chrom_cols,
                                       unique(plot.data$chromosome_name)))
  )
  
  rownames(plot.data) <-  gsub("ENSG[0-9]+([.][0-9]+)?_(.+)","\\2", plot.data$Gene)
  
  Heatmap(t(scale(t(plot.data %>% dplyr::select(-Gene,-source,-gene_biotype, -chromosome_name)))),
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
          cluster_columns = TRUE,
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
  
}
