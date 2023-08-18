###########################################
# This script generates a jitter-box plot #
# of Tumour Mutation Burden versus        #
# clinical outcome/category               #
#                                         #
# Author: Aidan Flynn                     #
# Date: 17/08/2023                        #
###########################################

library(ggplot2)

################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

############
# Plotting #
############

plot_data <- a5_anno %>% 
         mutate(wgs_tmb=as.numeric(wgs_tmb),
                differential_group_sampletype_strict=factor(differential_group_sampletype_strict,
                                                            levels=c("Metastasis",
                                                                      "Metastatic primary",
                                                                      "Metastatic local recurrence",
                                                                      "Local recurrence (metastasis present)",
                                                                      "Primary (metastasis present)",
                                                                      "Non-metastatic local recurrence",
                                                                      "Primary (short follow up)",
                                                                      "Non-metastatic primary"))) %>% 
         filter(A5_ID != "E167-1")
         
  
ggplot(plot_data, 
       aes(x=differential_group_sampletype_strict, y=wgs_tmb)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(mapping = aes(color=TERT_ATRX_Mutation),width = 0.2) + 
  scale_color_manual(values = driver_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1,hjust=1)) +
  ylab("TMB") +
  xlab("")
