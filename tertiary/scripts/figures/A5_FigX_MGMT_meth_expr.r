###########################################
###########################################
##                                       ##
##   E167 Treatment signature / MGMT     ##
##                                       ##
##   Script for producing a scatter      ##
##   plots showing MGMT methylation      ##
##   and expression highlighting         ##
##   E167                                ##
##                                       ##
###########################################
###########################################


setwd("/g/data/pq08/projects/ppgl")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

###############################
# Import data loaders and run #
###############################

#Methylation
source("./a5/methylation/scripts/data_loaders/a5_methylation_dataloader.r")
data_loader_a5_methylation_array(quickload = T)

#WTS
source("./a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
htseq_outs <- "./a5/wts/analysis/htseq/truseq/gene"
data_loader_a5_wts_counts(count_file_dir=htseq_outs)

#clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

b_vals <- getBeta(a5_methylation_filtered)

##################
# MGMT meth/expr #
##################

MGMT_promoter_probes <- c("cg12434587","cg12981137")

probe_beta <- b_vals %>% 
  as_tibble(rownames = "Probe") %>% 
  filter(Probe %in% MGMT_promoter_probes) %>% 
  pivot_longer(cols=-Probe, 
               names_to = "A5_ID", 
               values_to = "Beta") %>% 
  mutate(point_color=ifelse(grepl("167", A5_ID), "Red", "Black"), 
         Label=ifelse(grepl("167", A5_ID), A5_ID, NA))

mgmt_expr <- cpm(a5_wts_dge_list$SDHB, log = T) %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  filter(ensgid_symbol == "ENSG00000170430.10_MGMT") %>% 
  pivot_longer(cols=-ensgid_symbol, 
               names_to = "A5_ID", 
               values_to = "log2cpm") %>% 
  mutate(point_color=ifelse(grepl("167", A5_ID), "Red", "Black"), 
         Label=ifelse(grepl("167", A5_ID), A5_ID, NA))

pj <- position_jitter(width=0.2, seed=100)
gg_meth <- ggplot(probe_beta, 
       aes(x=Probe, y=Beta, color=point_color, label=Label)) + 
  geom_point(position = pj) + 
  geom_text(mapping=aes(y=Beta+0.02),position = pj, color="black") +
  scale_color_identity() + 
  theme_bw() +
  xlab("Probe") +
  ylab("Beta")

gg_expr <- ggplot(mgmt_expr, 
       aes(x="MGMT", y=log2cpm, color=point_color, label=Label)) + 
  geom_point(position = pj) + 
  geom_text_repel(mapping=aes(y=log2cpm),position = pj, color="black") +
  scale_color_identity() + 
  theme_bw() +
  xlab("") +
  ylab("CPM (log2)")

gg_meth + gg_expr + plot_layout(nrow=1, widths = c(2,1))

