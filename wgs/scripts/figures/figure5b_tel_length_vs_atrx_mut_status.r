################################
# Plot of telomere length vs   #
# atrx expression highlighting #
# c-circle and TERT/ATRX       #
# mutation status              #
#                              #  
# Author: Aidan Flynn          #
# Date: 31/07/2023             #
################################

library(ggplot2)
library(dplyr)
library(tidyr)

#################
# Clinical Data #
#################

if(!exists("a5_anno"))
{
  source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
  data_loader_a5_clinical_anno("aidan.flynn@umccr-radio-lab.page", use_cache = T)
}

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")


############
# Plotting #
############

plot_data <- a5_anno %>%  
              dplyr::select(A5_ID, `Patient ID`, TERT_ATRX_Mutation,
                            telhunter_log2_telcontentratio, cell_of_origin , 
                            differential_group_sampletype_strict, c_circle_result) %>% 
              mutate(c_circle_result = recode(c_circle_result, "Positive - Low"="Positive"),
                     c_circle_result=factor(c_circle_result, levels=c("No data","Negative", "Positive"))) %>% 
  arrange(TERT_ATRX_Mutation, c_circle_result, differential_group_sampletype_strict) %>% 
  mutate(telhunter_log2_telcontentratio = as.numeric(telhunter_log2_telcontentratio))

t_test_data <- plot_data %>% filter(cell_of_origin == "Chromaffin") %>%   
  group_by(`Patient ID`, TERT_ATRX_Mutation) %>% 
  summarise(telhunter_log2_telcontentratio=mean(telhunter_log2_telcontentratio)) 

t_result_atrx_vs_tert <- t.test(telhunter_log2_telcontentratio~TERT_ATRX_Mutation,  
                                t_test_data %>% filter(TERT_ATRX_Mutation != "WT"))


t_result_atrx_vs_wt <- t.test(telhunter_log2_telcontentratio~TERT_ATRX_Mutation,  
                              t_test_data %>% filter(TERT_ATRX_Mutation != "TERT"))

t_result_tert_vs_wt <- t.test(telhunter_log2_telcontentratio~TERT_ATRX_Mutation,  
                                t_test_data %>% filter(TERT_ATRX_Mutation != "ATRX"))


  
pj=position_jitter(width=0.2, seed=10)  
ggplot(plot_data ,aes(y=as.numeric(telhunter_log2_telcontentratio), 
                      x=TERT_ATRX_Mutation, 
                      color=c_circle_result,
                      #fill=differential_group_sampletype_strict, 
                      group=`Patient ID`)) + 
  geom_boxplot(aes(color=NULL,fill=NULL, group=NULL), color="lightgrey",outlier.alpha = 0) +
  geom_point(position=pj, size=3, stroke=1, shape=20) + 
  geom_path(position=pj, alpha=0.3, size=0.5, color="black", linetype=2) +
  geom_segment(data=tibble(y=c(1.7,1.85,2),
                           yend=c(1.7,1.85,2),
                           x=c(0.6,1.6,0.6), 
                           xend=c(2.4,3.4,3.4), 
                           cell_of_origin = "Chromaffin"),
               mapping=aes(x=x,y=y,xend=xend,yend=yend, color=NULL,fill=NULL,group=NULL)) +
  geom_text(data=tibble(y=c(1.78,1.91,2.08),
                        x=c(1.4,2.6,2), 
                        label=c(paste0("p=",(round(t_result_atrx_vs_tert$p.value,4))), 
                                paste0("p=",round(t_result_tert_vs_wt$p.value,3)),
                                paste0("p=",as.character(round(t_result_atrx_vs_wt$p.value,6)))),
                        cell_of_origin = "Chromaffin"),
            mapping=aes(x=x,y=y,label=label, color=NULL,fill=NULL,group=NULL)) +
  #geom_text(aes(label=A5_ID)) + 
  #scale_fill_manual(values = sampletype_strict_cols) +
  scale_color_manual(values = c("Positive"="red", 
                                "Weakly Positive"="orange",
                                "Negative"="black", 
                                "No data"="Grey")) +
  # scale_shape_manual(values = c("WT"=21,
  #                               "TERT"=24,
  #                               "ATRX"=25)) +
  #scale_shape_manual(values = c(8,18,15,17,19,4)) +
  # scale_size_manual(values=c(1.3,2.2,2.2)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(shape=21)),
         color = guide_legend(override.aes = list(shape=1, linetype=c(0,0,0)))) +
  ylab("Tumour/Normal Telomere Content Ratio (log2)") +
  facet_grid(~cell_of_origin, scales="free_x", space="free_x") +
  xlab("")


t_test_data <- plot_data %>% filter(cell_of_origin == "Chromaffin") %>%   
  group_by(`Patient ID`, c_circle_result) %>% 
  summarise(telhunter_log2_telcontentratio=mean(telhunter_log2_telcontentratio)) 

t_result_pos_vs_neg <- t.test(telhunter_log2_telcontentratio~c_circle_result,  
                                t_test_data %>% filter(c_circle_result != "No data"))

pj=position_jitter(width=0.2, seed=10)  
ggplot(plot_data %>% 
         filter(c_circle_result != "No data"),
       aes(y=as.numeric(telhunter_log2_telcontentratio))) + 
  geom_boxplot(aes(x=c_circle_result), color="lightgrey",outlier.alpha = 0) +
  geom_point(aes(x=c_circle_result, color=TERT_ATRX_Mutation), position=pj, size=3, stroke=1.5, shape=20) + 
  geom_path(aes(x=c_circle_result, group=`Patient ID`), position=pj, alpha=0.3, size=0.5, color="black", linetype=2) +
  geom_segment(data=tibble(y=c(1.7),yend=c(1.7),x=c(0.6,1.6), xend=c(2.4), cell_of_origin = "Chromaffin"),
               mapping=aes(x=x,y=y,xend=xend,yend=yend, color=NULL,fill=NULL,group=NULL)) +
  geom_text(data=tibble(y=c(1.78),
                        x=c(1.4),
                        label=c(paste0("p=",as.character(round(t_result_pos_vs_neg$p.value,6)))),
                        cell_of_origin = "Chromaffin"),
            mapping=aes(x=x,y=y,label=label, color=NULL,fill=NULL,group=NULL)) +
  #geom_text(aes(label=A5_ID)) + 
  #scale_fill_manual(values = sampletype_strict_cols) +
  scale_color_manual(values = c(ATRX = driver_cols[['ATRX']], TERT = driver_cols[['TERT']], WT = "grey")) +
  # scale_shape_manual(values = c("WT"=21,
  #                               "TERT"=24,
  #                               "ATRX"=25)) +
  #scale_shape_manual(values = c(8,18,15,17,19,4)) +
  # scale_size_manual(values=c(1.3,2.2,2.2)) +
  theme_bw() +
  #theme(legend.position = "bottom") +
  # guides(fill = guide_legend(override.aes = list(shape=21, linetype=c(0,0,0), color="white")),
  #        color = guide_legend(override.aes = list(shape=1, linetype=c(0,0,0)))) +
  ylab("Tumour/Normal Telomere Content Ratio (log2)") +
  facet_grid(~cell_of_origin, scales="free_x", space="free_x") +
  xlab("C-circle status")
