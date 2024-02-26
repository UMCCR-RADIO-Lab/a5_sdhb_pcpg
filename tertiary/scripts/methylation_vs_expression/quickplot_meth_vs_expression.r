setwd("/g/data/pq08/projects/ppgl/a5/")

library(tidyverse)
library(ggplot2)
library(patchwork)

###########
# Imports #
###########

#Modified version of gometh from missmethyl to use an offline cache for KEGG pathways
source("./methylation/scripts/go_meth_offline.r")

#Plotting helper functions
source("./methylation/scripts/helpers/plot_methylation.r")

###############
# DataLoaders #
###############

#load clinical annotation
source("./sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#load methylation array data
source("./methylation/scripts/data_loaders/a5_methylation_dataloader.r")
data_loader_a5_methylation_array(quickload = T, output_qc = F, normalisation="functional")

#load WTS data
source("/g/data/pq08/projects/ppgl/a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
htseq_outs <- "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/truseq/gene"
data_loader_a5_wts_counts(count_file_dir=htseq_outs)

######################
# Compute M/B-values #
######################

# calculate M-values for statistical analysis
m_vals <- getM(a5_methylation_filtered)
b_vals <- getBeta(a5_methylation_filtered)

#########################
# Compute logcpm values #
#########################


a5_wts_lcpm_list <- lapply(a5_wts_dge_list, cpm, log = T)

plot_methylation_vs_expr(#gene_symbol = "MSH6",
                         probes = c("cg00539935", "cg12755110", "cg16080759", "cg16975576", "cg02269272", "cg05628496"),
                         #samples_to_label = c("E169-1", "E169-2"),
                         label_outliers = F,
                         sample_annotation = (a5_anno %>% dplyr::select(Sample_ID=A5_ID,differential_group_sampletype_strict, Gender, TERT_ATRX_Mutation, differential_group_anatomy)), 
                         b_vals = b_vals, 
                         summarise_probes_by_region = F,
                         m_vals = m_vals,
                         log2_cpm =a5_wts_lcpm_list$SDHB, 
                         plot_mode = "beta", 
                         array_type = "EPIC", 
                         colour_column = "differential_group_sampletype_strict", 
                         grouping_column = "differential_group_sampletype_strict",
                         shape_column = "TERT_ATRX_Mutation",
                         )

plot_methylation(gene_symbol = "ATRX", 
                 show_mean_summarised_only = F,
                 sample_annotation = (a5_anno %>% dplyr::select(Sample_ID=A5_ID,TERT_ATRX_Mutation, Gender)), 
                 b_vals = b_vals, 
                 m_vals = m_vals, 
                 plot_mode = "beta", 
                 array_type = "EPIC", 
                 colour_column = "TERT_ATRX_Mutation", 
                 grouping_column = "TERT_ATRX_Mutation",
                 shape_column = "Gender")
