setwd("/g/data/pq08/projects/ppgl/a5")

library(tidyverse)
library(ggplot2)
library(patchwork)

###########
# Imports #
###########

#Plotting helper functions
source("./methylation/scripts/helpers/plot_methylation.r")

###############
# DataLoaders #
###############

#load clinical annotation
source("./sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#load methylation array data
#Methylation
source("./tertiary/scripts/data_mergers/combine_tcga_comete_a5_methylation_data.r")

#load WTS data
source("./tertiary/scripts/data_mergers/combine_tcga_flynn_a5_wts_data.r")

##############
# Sample Set #
##############

sample_anno <- wts_tcga_flynn_a5_anno %>% 
  filter(Sample %in% intersect(wts_tcga_flynn_a5_anno$Sample, meth_tcga_comete_a5_450k_anno$Sample)) %>%
  mutate(Malignancy=replace_na(Malignancy, "Unknown")) %>%
  dplyr::rename("subtype"="new_naming")

######################
# Compute M/B-values #
######################

# calculate M-values for statistical analysis
m_vals <- meth_tcga_comete_a5_450k_mval.batch_removed[, intersect(colnames(meth_tcga_comete_a5_450k_mval.batch_removed), sample_anno$Sample)]
b_vals <- lumi::m2beta(meth_tcga_comete_a5_450k_mval.batch_removed[, intersect(colnames(meth_tcga_comete_a5_450k_mval.batch_removed), sample_anno$Sample)])

#################
# Relabel genes #
#################

ensgid_symbol_lookup <- data.frame(ensgid_symbol=rownames(a5_wts_counts)) %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), remove=F, extra = "merge", sep="_") %>% 
  mutate(ensgid=gsub("[.].+$", "", ensgid))

rownames(wts_tcga_flynn_a5_lcpm.batch_removed) <- ensgid_symbol_lookup$symbol[match(rownames(wts_tcga_flynn_a5_lcpm.batch_removed), ensgid_symbol_lookup$ensgid)]


############
# Plotting #
############

plot_methylation_vs_expr(#gene_symbol = "MGMT",
                         probes= c("cg00539935", "cg12755110", "cg16080759", "cg16975576", "cg02269272", "cg05628496"),
                         samples_to_label = c("E169-1", "E169-2", "E167-2"),
                         label_outliers = F,
                         sample_annotation = (sample_anno %>% dplyr::select(Sample_ID=Sample, subtype, Sex, Malignancy, Genotype, PCPG, ATRX_mutation)), 
                         b_vals = b_vals, 
                         summarise_probes_by_region = F,
                         m_vals = m_vals,
                         log2_cpm =wts_tcga_flynn_a5_lcpm.batch_removed, 
                         plot_mode = "beta", 
                         array_type = "450K", 
                         colour_column = "subtype", 
                         grouping_column = "subtype",
                         shape_column = "Sex", 
                         colour_scale = subtype_cols 
)

plot_methylation(gene_symbol = "SLC6A2", 
                 #region = "chr10:131265043-131265625",
                 show_mean_summarised_only = F,
                 label_outliers = F,
                 sample_annotation = (sample_anno %>% dplyr::select(Sample_ID=Sample, Dataset, subtype, Malignancy, Genotype, PCPG, ATRX_mutation)), 
                 b_vals = b_vals, 
                 m_vals = m_vals, 
                 plot_mode = "m", 
                 array_type = "450K", 
                 colour_column = "subtype", 
                 grouping_column = "Malignancy",
                 colour_scale = subtype_cols,
                 shape_column = "Dataset")  + coord_cartesian(ylim=c(-6,4))
