library(ggplot2)
library(patchwork)
library(tidyr)

setwd("/g/data/pq08/projects/ppgl/")

################
# Data Loaders #
################

if(!exists("a5_anno"))
{
  source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
  data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)
}

if(!exists("a5_seg_keep"))
{
  source("./a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
  data_loader_cna_segs()
}

if(!exists("ColorPalette"))
{
  source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")
}

source("./a5/wgs/scripts/pcnt_genome_altered/percent_genome_altered.r")
source("./a5/wgs/scripts/figures/cna_frequency_plot.r")


############
# Plotting #
############

nox <-  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
noy <-  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
zero_margin <- theme(plot.margin = margin(0,0,0,0))
no_strip <- theme(strip.background = element_blank(), strip.text.y = element_blank())

plot_data <- a5_anno %>%  
  inner_join(genome_altered_stats) %>% 
  dplyr::select(A5_ID, `Patient ID`, 
                differential_group_sampletype_strict, 
                cell_of_origin,
                `sample_ploidy`, pcnt_altered, 
                structural_variant_count,
                TERT_ATRX_Mutation) 

plot_data <- plot_data %>% 
  mutate(WGD=ifelse(as.numeric(`sample_ploidy`) > 3.2,"Yes","No"),
         WGD=factor(as.character(WGD), levels=c("Yes","No")),
         clinical_class = case_match(as.character(differential_group_sampletype_strict),
                                     "Non-metastatic primary" ~ "Non-metastatic primary/Primary (short follow up)" , 
                                     "Primary (short follow up)" ~ "Non-metastatic primary/Primary (short follow up)",  
                                     "Primary (metastasis reported)" ~ "Primary (metastasis reported)",
                                     "Local recurrence (metastasis reported)" ~ "Primary (metastasis reported)",
                                     "Metastatic primary" ~ "Metastatic primary/Metastasis",
                                     "Metastatic local recurrence" ~ "Metastatic primary/Metastasis",
                                     "Metastasis" ~ "Metastatic primary/Metastasis",
                                     .default = differential_group_sampletype_strict),
         clinical_class = factor(clinical_class, 
                                 levels=c("Non-metastatic primary/Primary (short follow up)",
                                          "Primary (metastasis reported)",
                                          "Metastatic primary/Metastasis")),
         TERT_ATRX_Mutation = factor(TERT_ATRX_Mutation, levels=c("WT","ATRX","TERT"))) %>% 
  inner_join(genome_altered_stats) %>%  
  arrange(cell_of_origin, TERT_ATRX_Mutation, clinical_class,  pcnt_altered) 

SampleOrder.A5_ID <- plot_data$A5_ID 
  
plot_data$A5_ID <- factor(as.character(plot_data$A5_ID), levels=SampleOrder.A5_ID)



gg_cna_seg <- ggplot(a5_seg_keep_wgd_norm %>%  
                       inner_join(plot_data %>%  dplyr::select(A5_ID, clinical_class, cell_of_origin)) %>% 
                       filter(Class %in% c("Gain","Loss", "CNLOH", "Chromothripsis", "Subclonal Loss")) %>%  
                       mutate(A5_ID = factor(A5_ID, levels=SampleOrder.A5_ID)), 
                   aes(y=A5_ID, yend=A5_ID, x=start_offset, xend=end_offset, color=Class)) + 
  geom_segment(linewidth=1) +
  geom_vline(data = chr_offsets, mapping=aes(xintercept=offset), linewidth=0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1,hjust=1), panel.grid = element_blank()) +
  scale_x_continuous(breaks = chr_offsets$offset+(chr_offsets$chr_size/2), labels = chr_offsets$chromosome, expand=c(0,0)) +
  scale_color_manual(values=cna_palette) + 
  ylab("Copy Number") + 
  facet_grid(rows="cell_of_origin", scale = "free_y", space="free") 

gg_sampletype_rug <- ggplot(plot_data, aes(x="clinical_class", y=A5_ID, fill=clinical_class)) +
  geom_tile(height=0.8) +
  scale_fill_manual(values=c("Non-metastatic primary/Primary (short follow up)"="white",
                             "Primary (metastasis reported)"="grey",
                             "Metastatic primary/Metastasis"="black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1),
        panel.grid = element_blank())  + 
  facet_grid(rows="cell_of_origin", scale = "free_y", space="free", switch = "y") 

gg_tert_atrx_rug <- ggplot(plot_data, aes(x="TERT_ATRX", y=A5_ID, fill=TERT_ATRX_Mutation)) +
  geom_tile(height=0.8) +
  scale_fill_manual(values = driver_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1),
        panel.grid = element_blank())  + 
  facet_grid(rows="cell_of_origin", scale = "free_y", space="free", switch = "y") 

gg_wgd_rug <- ggplot(plot_data, aes(x="WGD", y=A5_ID, fill=WGD)) +
  geom_tile(height=0.8) +
  scale_fill_manual(values=c("Yes"="black", "No"="white")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1),
        panel.grid = element_blank())  + 
  facet_grid(rows="cell_of_origin", scale = "free_y", space="free") 

gg_pcnt_altered <- ggplot(plot_data, aes(x=pcnt_altered, y=A5_ID)) +
  geom_col(width = 0.7) + theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1,hjust=1))  + 
  facet_grid(rows="cell_of_origin", scale = "free_y", space="free") 

gg_sv_count <- ggplot(plot_data, aes(x=structural_variant_count, y=A5_ID)) +
  geom_col(width = 0.7) + theme_bw() + 
  scale_x_log10(minor_breaks=c(seq(1,9,1), seq(10,90,10), seq(100,900,100), seq(1000,10000,1000))) +
  theme(axis.text.x = element_text(angle=45, vjust = 1,hjust=1))  + 
  facet_grid(rows="cell_of_origin", scale = "free_y", space="free") 

gg_cna_freq_anatomy_flip + nox + theme(plot.margin = margin(6,0,0,0)) +
gg_sampletype_rug + noy + zero_margin  +
  gg_tert_atrx_rug + noy + zero_margin + no_strip +
  gg_wgd_rug + noy + zero_margin + no_strip +
  gg_cna_seg +noy + theme(plot.margin = margin(0,0,6,0)) + no_strip +
  gg_pcnt_altered + noy + zero_margin + no_strip +
  gg_sv_count + noy  + zero_margin + no_strip +
  plot_layout(guides="collect", widths = c(1,1,1,30,2,2), heights=c(1,3), design = "###A##\nBCDEFG")

