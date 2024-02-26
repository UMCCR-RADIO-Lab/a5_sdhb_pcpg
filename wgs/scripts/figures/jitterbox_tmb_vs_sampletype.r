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
                                                                      "Local recurrence (metastasis reported)",
                                                                      "Primary (metastasis reported)",
                                                                      "Non-metastatic local recurrence",
                                                                      "Primary (short follow up)",
                                                                      "Non-metastatic primary"))) %>% 
         filter(A5_ID != "E167-1") %>% 
  mutate(met_nonmet=recode(as.character(differential_group_sampletype_strict),
                                        "Metastasis"="Metastatic",
                                        "Metastatic primary"="Metastatic",
                                        "Metastatic local recurrence"="Metastatic",
                                        "Local recurrence (metastasis reported)"="Other",
                                        "Primary (metastasis reported)"="Other",
                                        "Non-metastatic local recurrence"="Other",
                                        "Primary (short follow up)"="Other",
                                        "Non-metastatic primary"="Non-metastatic"),
         met_nonmet=factor(met_nonmet, levels=c("Metastatic", "Non-metastatic", "Other")))
         
##Detailed categories  
ggplot(plot_data, 
       aes(x=differential_group_sampletype_strict, y=wgs_tmb)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(mapping = aes(color=TERT_ATRX_Mutation),width = 0.2) + 
  scale_color_manual(values = driver_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1,hjust=1)) +
  ylab("TMB") +
  xlab("")



##Broad categories
t_test_data <- plot_data %>%  
  group_by(`Patient ID`, met_nonmet) %>% 
  summarise(tmb=mean(wgs_tmb)) %>% 
  filter(met_nonmet != "Other")

t_result <- t.test(tmb~met_nonmet,  t_test_data, alternative="greater")


pj=position_jitter(width = 0.2,seed = 10)
ggplot(plot_data %>% arrange(met_nonmet, `Patient ID`), 
       aes(x=met_nonmet, y=wgs_tmb)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(mapping = aes(color=differential_group_sampletype_strict, 
                            shape=TERT_ATRX_Mutation),
              position=pj) + 
  geom_path(mapping=aes(group=`Patient ID`), position=pj, linetype=2, alpha=0.3) +
  scale_color_manual(values = sampletype_strict_cols) +
  scale_shape_manual(values = c(TERT=8,ATRX=4,WT=16)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1,hjust=1)) +
  ylab("TMB") +
  xlab("")


