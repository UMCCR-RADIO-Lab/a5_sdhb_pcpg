################################################
# This script produces a plot of tumour/normal #
# ratios for each TVR (that is significant by  #
# ANOVA) versus C-Circle status.               #
#                                              #
# Author: Aidan Flynn                          #
# Date: 04/08/2023                             #
################################################

library(dplyr)
library(tidyr)
library(ggplot2)

#################
# Clinical Data #
#################

if(!exists("a5_anno"))
{
  source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
  data_loader_a5_clinical_anno("aidan.flynn@umccr-radio-lab.page", use_cache = T)
}

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#Simplify c-circle Positive annotation
a5_anno <- a5_anno %>%  mutate(c_circle_result = recode(c_circle_result, "Positive - Low"="Positive"))

#########################################
# Telomere Hunter Normalized TVR counts #
#########################################

a5th.norm_tvr <- read.delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/telomerehunter/combined/telomerehunter_a5_normalized_tvr_counts.tsv")


########
# TVRs normalised by all reads
########

#Find sig diff TVRs with ANOVA

sigtest_data <- a5th.norm_tvr %>% mutate(PID=gsub("T0","",PID)) %>%  
  mutate(log2_ratio_count_norm_by_all_reads=log2(Count_norm_by_all_reads_T/Count_norm_by_all_reads_C)) %>% 
  dplyr::select(PID, Pattern, log2_ratio_count_norm_by_all_reads, log2_ratio_count_norm_by_intratel_reads) %>% 
  inner_join(a5_anno %>% 
               dplyr::select(A5_ID, TERT_ATRX_Mutation, c_circle_result) %>% 
               dplyr::rename(PID=A5_ID)) %>% dplyr::select(-PID) %>% 
  filter(c_circle_result %in% c("Negative", "Positive"))

sigtest_data$Pattern <- as.factor(sigtest_data$Pattern)
sigtest_data$TERT_ATRX_Mutation <- as.factor(sigtest_data$TERT_ATRX_Mutation)
sigtest_data$c_circle_result <- as.factor(sigtest_data$c_circle_result)

sig_patterns <- data.frame(Pattern=vector(mode = "character"), t_pval_total=vector(mode="double"), t_pval_tel=vector(mode="double"), t_pval_tel_lt=vector(mode="double"))
for (p in levels(sigtest_data$Pattern))
{
  
  t_result_totalreads <- t.test(formula = log2_ratio_count_norm_by_all_reads ~ c_circle_result, 
                                data = sigtest_data %>% filter(Pattern==p, !is.infinite(log2_ratio_count_norm_by_all_reads)))
  
  t_result_telreads <- t.test(formula = log2_ratio_count_norm_by_intratel_reads ~ c_circle_result, 
                              data = sigtest_data %>% filter(Pattern==p, !is.infinite(log2_ratio_count_norm_by_intratel_reads)))
  
  t_result_telreads_lt <- t.test(formula = log2_ratio_count_norm_by_intratel_reads ~ c_circle_result, 
                                 data = sigtest_data %>% filter(Pattern==p, !is.infinite(log2_ratio_count_norm_by_intratel_reads)),alternative="greater")
  
  
  sig_patterns <- sig_patterns %>% 
    tibble::add_row(data.frame(Pattern=p, 
                               t_pval_total=t_result_totalreads$p.value, 
                               t_pval_tel=t_result_telreads$p.value,
                               t_pval_tel_lt=t_result_telreads_lt$p.value))
  
}

sig_patterns$t_pval_total_adj <- p.adjust(sig_patterns$t_pval_total, method = "BH")
sig_patterns$t_pval_tel_adj <- p.adjust(sig_patterns$t_pval_tel, method = "BH")
sig_patterns$t_pval_tel_lt_adj <- p.adjust(sig_patterns$t_pval_tel_lt, method = "BH")
sig_patterns <- sig_patterns %>%  filter(t_pval_total < 0.05)

plot.data <-  a5th.norm_tvr %>%  
  #filter(Pattern %in% sig_patterns) %>%  
  mutate(PID=gsub("T0","",PID)) %>%  
  mutate(log2_ratio_count_norm_by_all_reads=log2(Count_norm_by_all_reads_T/Count_norm_by_all_reads_C)) %>%
  inner_join(a5_anno %>% 
               dplyr::select(A5_ID, `Patient ID`, TERT_ATRX_Mutation, 
                             Primary_Location_Simplified, is_primary_or_met, 
                             tumour_metastasised, c_circle_result, telhunter_log2_telcontentratio) %>% 
               dplyr::rename(PID=A5_ID, 
                             Primary_Met=is_primary_or_met, 
                             Met_Case=tumour_metastasised)) %>% 
  filter(c_circle_result %in% c("Negative", "Positive"))
# %>% filter(Pattern %in% sig_patterns)

plot.data$telhunter_log2_telcontentratio <- as.numeric(plot.data$telhunter_log2_telcontentratio)
plot.data$c_circle_result <- factor(as.character(plot.data$c_circle_result), levels = c("Negative", "Positive")) #"No data",

plot.data.long <- plot.data %>%  
  mutate(PID=factor(as.character(PID), 
                    levels= (plot.data %>% arrange(c_circle_result, telhunter_log2_telcontentratio, TERT_ATRX_Mutation) %>% pull(PID) %>% unique()))) %>%
  dplyr::select(PID, Pattern, log2_ratio_count_norm_by_intratel_reads, log2_ratio_count_norm_by_all_reads, 
                c_circle_result, telhunter_log2_telcontentratio, TERT_ATRX_Mutation) %>% 
  pivot_longer(cols = c(log2_ratio_count_norm_by_intratel_reads, log2_ratio_count_norm_by_all_reads), 
               names_to="NormMethod", values_to="log2_ratio_count") %>% mutate(NormMethod=gsub("log2_ratio_count_norm_by_","",NormMethod))
plot.data.long$NormMethod <- factor(as.character(plot.data.long$NormMethod), levels=c("all_reads","intratel_reads"))

POI <- sig_patterns$Pattern

## TVR ratio vs C-Circle result
ggplot(
  plot.data.long %>% filter(Pattern %in% POI) %>% 
    mutate(NormMethod=recode(NormMethod, "all_reads"="Total Reads","intratel_reads"="Telomeric Reads")) %>% 
    arrange(c_circle_result, telhunter_log2_telcontentratio)
) +
  geom_boxplot(mapping = aes(
    x = interaction(c_circle_result, NormMethod),
    y = log2_ratio_count,
    fill=NormMethod
  ), outlier.alpha = 0) +
  geom_point(
    aes(x = interaction(c_circle_result, NormMethod),
        y = log2_ratio_count,
        color = telhunter_log2_telcontentratio > 0.5,
        fill=NormMethod), 
    position = position_jitter(width=0.15, 
                               height=0, 
                               seed=42),
    alpha=0.6,
    size=0.8) +
  #geom_text(aes(x = interaction(c_circle_result, NormMethod), y = log2_ratio_count,label=PID)) +
  scale_color_manual(values = c(ColorPalette[["DarkGrey2"]],
                                ColorPalette[["DarkRed2"]]),
                     name="T/N Tel. Ratio > 0.5") +
  scale_fill_manual(values = c(ColorPalette[["LightBlue3"]],
                               ColorPalette[["LightGreen1"]]),
                    name="Normalisation") +
  #scale_shape_manual(values = c(16,8,17)) +
  scale_x_discrete(labels=rep(levels(plot.data.long$c_circle_result),2)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)) +
  ylab("log2_ratio_count") +
  xlab("C-Circle Result") + 
  facet_wrap("Pattern", scale="free_y", ncol = 5) +
  geom_vline(xintercept = 2.5, linetype=2)
