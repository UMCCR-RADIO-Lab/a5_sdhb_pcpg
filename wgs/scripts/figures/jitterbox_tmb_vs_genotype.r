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
  mutate(wgs_tmb = as.numeric(wgs_tmb),
         differential_group_sampletype_strict = factor(differential_group_sampletype_strict,
                                                     levels=c("Metastasis",
                                                              "Metastatic primary",
                                                              "Primary (metastasis reported)",
                                                              "Non-metastatic local recurrence",
                                                              "Primary (short follow up)",
                                                              "Non-metastatic primary"))) %>%   
  arrange(TERT_ATRX_Mutation, `Patient ID`, differential_group_sampletype_strict)
  
         
##t-test

t_test_data <- plot_data %>% filter(cell_of_origin == "Chromaffin") %>%   
  group_by(`Patient ID`, TERT_ATRX_Mutation) %>% 
  filter(A5_ID != "E167-1") %>% 
  summarise(wgs_tmb=mean(wgs_tmb)) 

t_result_atrx_vs_tert <- t.test(wgs_tmb~TERT_ATRX_Mutation,  
                                t_test_data %>% filter(TERT_ATRX_Mutation != "WT"), alternative = c("greater"))


t_result_atrx_vs_wt <- t.test(wgs_tmb~TERT_ATRX_Mutation,  
                              t_test_data %>% filter(TERT_ATRX_Mutation != "TERT"), alternative = c("greater"))

t_result_tert_vs_wt <- t.test(wgs_tmb~TERT_ATRX_Mutation,  
                              t_test_data %>% filter(TERT_ATRX_Mutation != "ATRX"), alternative = c("greater"))


jp = position_jitter(width = 0.2, seed=25)
gg_bottom <- ggplot(plot_data, 
                    aes(x=TERT_ATRX_Mutation, y=wgs_tmb)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(mapping = aes(color=differential_group_sampletype_strict),position = jp) + 
  geom_path(mapping=aes(group=`Patient ID`), position = jp, linetype=2) +
  scale_color_manual(values = sampletype_strict_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1,hjust=1)) +
  ylab("TMB") +
  xlab("") + 
  facet_grid(~cell_of_origin, space="free_x", scales = "free_x") +
  coord_cartesian(ylim=c(0,5))


gg_top <- ggplot(plot_data, 
                 aes(x=TERT_ATRX_Mutation, y=wgs_tmb)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(mapping = aes(color=differential_group_sampletype_strict),position = jp) + 
  geom_path(mapping=aes(group=`Patient ID`), position = jp, linetype=2) +
  geom_segment(data=tibble(y=c(745,753,748),
                           yend=c(745,753,748),
                           x=c(0.6,1.6,0.6), 
                           xend=c(2.4,3.4,3.4),
                           cell_of_origin = factor("Chromaffin", levels=levels(plot_data$cell_of_origin))),
               mapping=aes(x=x,y=y,xend=xend,yend=yend, color=NULL,fill=NULL,group=NULL)) +
  geom_text(data=tibble(y=c(746,754,749),
                        x=c(1.4,2.6,2), 
                        label=c(paste0("p=",(round(t_result_atrx_vs_tert$p.value,4))), 
                                paste0("p=",round(t_result_tert_vs_wt$p.value,6)),
                                paste0("p=",as.character(round(t_result_atrx_vs_wt$p.value,6)))),
                        cell_of_origin = factor("Chromaffin", levels=levels(plot_data$cell_of_origin))),
            mapping=aes(x=x,y=y,label=label, color=NULL,fill=NULL,group=NULL)) +
  scale_color_manual(values = sampletype_strict_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1,hjust=1)) +
  ylab("TMB") +
  xlab("") + 
  facet_grid(~cell_of_origin, space="free_x", scales = "free_x") +
  coord_cartesian(ylim=c(735,755)) + 
  scale_y_continuous(breaks = c(735,740,745,750,755))

gg_top + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(6,6,0,6)) + 
  gg_bottom + theme(strip.background = element_blank(), strip.text = element_blank(), plot.margin = margin(0,6,6,6)) + 
  plot_layout(nrow = 2, guides = "collect", heights = c(1,3))
