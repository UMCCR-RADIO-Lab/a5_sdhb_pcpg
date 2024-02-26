setwd("/g/data/pq08/projects/ppgl/a5/")

library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggpubr)

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
source("./methylation/scripts/data_loaders/a5_methylation_dataloader.r")
data_loader_a5_methylation_array(quickload = T, output_qc = F, normalisation="functional")

#load WTS data
source("/g/data/pq08/projects/ppgl/a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
htseq_outs <- "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/truseq/gene"
data_loader_a5_wts_counts(count_file_dir=htseq_outs)

source("./sample_annotation/scripts/data_loaders/a5_color_scheme.r")

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

sample_anno <- a5_anno %>% dplyr::select(Sample_ID=A5_ID,differential_group_sampletype_strict, Gender, TERT_ATRX_Mutation, Primary_Location_Simplified) %>%  
  mutate(Primary_Location_Simplified=dplyr::recode(Primary_Location_Simplified,
                                                   "Extraadrenal_abdominal"="Extraadrenal (abdominal/thoracic)",
                                                   "Extraadrenal_thoracic"="Extraadrenal (abdominal/thoracic)",
                                                   "Extraadrenal_thoracic_cardiac"="Extraadrenal (abdominal/thoracic)",
                                                   "Extraadrenal_bladder"="Extraadrenal (bladder)",
                                                   "Extraadrenal_thoracic_mediastinum"="Extraadrenal (mediastinum)",
                                                   "Adrenal_left"="Adrenal",
                                                   "Adrenal_right"="Adrenal",
                                                   "Head_neck"="Head and neck"),
         Primary_Location_Simplified=factor(Primary_Location_Simplified,
                                            levels=c("Extraadrenal (abdominal/thoracic)",
                                                     "Extraadrenal (bladder)",
                                                     "Extraadrenal (mediastinum)",
                                                     "Adrenal",
                                                     "Head and neck",
                                                     "Unspecified")))



######
# TH #
######

gg_th <- plot_methylation_vs_expr(gene_symbol = "TH",
                                            probes = c("cg08573687", "cg12158055", "cg26253759"),
                                            label_outliers = F,
                                            sample_annotation = sample_anno,
                                            b_vals = b_vals, 
                                            summarise_probes_by_region = F,
                                            m_vals = m_vals,
                                            log2_cpm =a5_wts_lcpm_list$SDHB, 
                                            plot_mode = "m", 
                                            array_type = "EPIC", 
                                            colour_column = "Primary_Location_Simplified", 
                                            colour_scale = location_cols,
                                            grouping_column = "Primary_Location_Simplified",
                                            shape_column = "Gender"
)


gg_th <- gg_th + 
  geom_smooth(method='lm') + 
  stat_cor(label.x = -1, label.y = 15, method = "spearman") +
  coord_cartesian() + facet_wrap("probe_id", scales="free_x")



#######
# DBH #
#######

gg_dbh <- plot_methylation_vs_expr(gene_symbol = "DBH",
                                  probes = c("cg08902605", "cg02928015", "cg07824742"),
                                  label_outliers = F,
                                  sample_annotation = sample_anno,
                                  b_vals = b_vals, 
                                  summarise_probes_by_region = F,
                                  m_vals = m_vals,
                                  log2_cpm =a5_wts_lcpm_list$SDHB, 
                                  plot_mode = "m", 
                                  array_type = "EPIC", 
                                  colour_column = "Primary_Location_Simplified", 
                                  colour_scale = location_cols,
                                  grouping_column = "Primary_Location_Simplified",
                                  shape_column = "Gender"
)

gg_dbh <- gg_dbh + 
  geom_smooth(method='lm') + 
  stat_cor(label.x = -1, label.y = 15, method = "spearman") +
  coord_cartesian() + facet_wrap("probe_id", scales="free_x")


##########
# SLC6A2 #
##########

gg_slc6a2 <- plot_methylation_vs_expr(gene_symbol = "SLC6A2",
                         probes = c("cg16629702","cg09746736","cg03226000"), #,"cg01277542","cg14112935","cg07158132","cg02693870","cg02652273","cg11282495","cg26668679","cg27047406","cg24261673","cg10362591","cg11413593"),
                         #samples_to_label = c("E169-1", "E169-2"),
                         label_outliers = F,
                         sample_annotation = sample_anno,
                         b_vals = b_vals, 
                         summarise_probes_by_region = F,
                         m_vals = m_vals,
                         log2_cpm =a5_wts_lcpm_list$SDHB, 
                         plot_mode = "m", 
                         array_type = "EPIC", 
                         colour_column = "Primary_Location_Simplified", 
                         colour_scale = location_cols,
                         grouping_column = "Primary_Location_Simplified",
                         shape_column = "Gender",
                         export_plot_data = "slc6a_meth_expr"
)


gg_slc6a2 <- gg_slc6a2 + 
  geom_smooth(method='lm') + 
  stat_cor(label.x = -2, method = "spearman") +
  coord_cartesian() + facet_wrap("probe_id", scales="free_x")

##########
# MGMT #
##########

gg_mgmt <- plot_methylation_vs_expr(gene_symbol = "MGMT",
                                      probes = c("cg12434587","cg12981137"),
                                      samples_to_label = c("E169-1", "E169-2", "E167-2"),
                                      label_outliers = F,
                                      sample_annotation = sample_anno,
                                      b_vals = b_vals, 
                                      summarise_probes_by_region = F,
                                      m_vals = m_vals,
                                      log2_cpm =a5_wts_lcpm_list$SDHB, 
                                      plot_mode = "m", 
                                      array_type = "EPIC", 
                                      colour_column = "Primary_Location_Simplified", 
                                      colour_scale = location_cols,
                                      grouping_column = "Primary_Location_Simplified",
                                      shape_column = "Gender",
                                      export_plot_data = "mgmt_meth_expr"
)


gg_mgmt <- gg_mgmt + 
  geom_smooth(method='lm') + 
  stat_cor(label.x = -2, method = "spearman") +
  coord_cartesian() + facet_wrap("probe_id", scales="free_x")
