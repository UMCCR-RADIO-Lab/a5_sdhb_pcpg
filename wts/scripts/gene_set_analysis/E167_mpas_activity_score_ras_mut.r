library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

################
# Data Loaders #
################

#clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
if(!exists("a5_anno")) {
  data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)
}

#WTS
source("./a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
if (!exists("a5_wts_dge_list")) {
  htseq_outs <- "./a5/wts/analysis/htseq/neb/gene"
  data_loader_a5_wts_counts(count_file_dir=htseq_outs)
}

logcpm <- cpm(a5_wts_dge_list$qc_ok,log = T)

#MAPK Pathway Activity Score (MPAS) 
MPAS_genes <- c("SPRY2", "SPRY4", "ETV4", "ETV5", "DUSP4", "DUSP6", "CCND1", "EPHA2", "EPHA4")

MPAS_genes_expr <- logcpm %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), sep="_", extra = "merge") %>% 
  filter(symbol %in% MPAS_genes) %>%
  pivot_longer(cols = c(-ensgid,-symbol), names_to = "A5_ID", values_to = "log2_cpm") %>% 
  group_by(symbol) %>% 
  mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm))/sd(log2_cpm)) %>% 
  mutate(color=ifelse(A5_ID %in% c("E124-1", "E167-1","E167-2"),A5_ID, "Other"),
         label=ifelse(A5_ID %in% c("E124-1", "E167-1","E167-2"),A5_ID, NA))

MPAS_score <- MPAS_genes_expr %>% 
  group_by(A5_ID) %>% 
  summarise(`MPAS Score`=sum(log2_cpm_z)) %>%  
  arrange(desc(`MPAS Score`))  %>% 
  mutate(color=ifelse(A5_ID %in% c("E124-1", "E167-1","E167-2"), A5_ID, "Other"),
         label=ifelse(A5_ID %in% c("E124-1", "E167-1","E167-2"), A5_ID, NA))


col_scale <- scale_color_manual(values = c("E124-1" = "green", 
                                             "E167-1" = "red",
                                             "E167-2" = "blue", 
                                             "Other" = "grey"))

gg_mpas_components <- 
  ggplot(MPAS_genes_expr, 
         aes(x=symbol,
             y=log2_cpm_z, 
             color=color)) + 
  geom_jitter(width = 0.2) + 
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Gene Expression (z- scaled log2 CPM)") 
  
gg_mpas_score <- 
  ggplot(MPAS_score, 
         aes(x="MPAS Score",
             y=`MPAS Score`, 
             color=color)) + 
  geom_jitter(width = 0.2) + 
  geom_text(mapping = aes(label=label), color="black") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) 

gg_mpas_components +
  gg_mpas_score +
  plot_layout(widths = c(5,1), 
              guides="collect") &
  scale_color_manual(values = c("E124-1" = "green", 
                                "E167-1" = "red",
                                "E167-2" = "blue", 
                                "Other" = "grey")) &
  theme(axis.title.x = element_blank())

