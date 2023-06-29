###########################################
###########################################
##                                       ##
##   E167/E169 MGMT                      ##
##   Methylation/Expression              ##
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
               values_to = "Beta") 

mgmt_expr <- cpm(a5_wts_dge_list$SDHB, log = T) %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  filter(ensgid_symbol == "ENSG00000170430.10_MGMT") %>% 
  pivot_longer(cols=-ensgid_symbol, 
               names_to = "A5_ID", 
               values_to = "log2cpm")

mgmt_plots <- list()
for (patient in c("E167", "E169"))
{
  
  probe_beta <- probe_beta %>% 
  mutate(point_color=ifelse(grepl(patient, A5_ID), "Red", "Black"), 
         Label=ifelse(test = grepl(patient, A5_ID), 
                      yes = dplyr::recode(A5_ID, !!!setNames(a5_anno$PublicationID, a5_anno$A5_ID)), 
                      no = NA))

  mgmt_expr <- mgmt_expr %>% 
  mutate(point_color=ifelse(grepl(patient, A5_ID), "Red", "Black"), 
         Label=ifelse(test = grepl(patient, A5_ID), 
                      yes = dplyr::recode(A5_ID, !!!setNames(a5_anno$PublicationID, a5_anno$A5_ID)), 
                      no = NA))

  pj <- position_jitter(width=0.2, seed=100)
  gg_meth <-
    ggplot(probe_beta, 
           aes(x=Probe, y=Beta, 
               color=point_color, 
               label=Label)) + 
    geom_point(position = pj) + 
    geom_text(mapping=aes(y=Beta+0.02),position = pj, color="black") +
    scale_color_identity() + 
    theme_bw() +
    xlab("Probe") +
    ylab("Beta")
  
  gg_expr <- 
    ggplot(mgmt_expr, 
         aes(x="MGMT", y=log2cpm, color=point_color, label=Label)) + 
    geom_point(position = pj) + 
    geom_text_repel(mapping=aes(y=log2cpm),position = pj, color="black") +
    scale_color_identity() + 
    theme_bw() +
    xlab("") +
    ylab("CPM (log2)")
  
  mgmt_plots[[patient]] <- 
    gg_meth + 
    gg_expr + 
    plot_layout(nrow=1, widths = c(2,1))

  ggsave(filename = paste0("./a5/tertiary/results/plots/methylation_vs_expression/",patient,"_MGMT_expr_meth.pdf"),
        plot = mgmt_plots[[patient]],device = pdf(), width = 5, height = 9,units = "in")
}

