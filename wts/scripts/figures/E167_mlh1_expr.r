###########################################
###########################################
##                                       ##
##   E167 MLH1 Expression                ##
##                                       ##
##   Script for producing a scatter      ##
##   plot showing MLH1 expression        ##
##   highlighting E167                   ##
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

#WTS
source("./a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
htseq_outs <- "./a5/wts/analysis/htseq/neb/gene"
if (!exists("a5_wts_dge_list"))
{
  data_loader_a5_wts_counts(count_file_dir=htseq_outs)
}

#clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
if (!exists("a5_anno"))
{
  data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)
}

#############
# MLH1 expr #
#############

mlh1_expr <- cpm(a5_wts_dge_list$SDHB, log = T) %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  filter(ensgid_symbol == "ENSG00000076242.16_MLH1") %>% 
  pivot_longer(cols=-ensgid_symbol, 
               names_to = "A5_ID", 
               values_to = "log2cpm")

patient = "E167"

mlh1_expr <- mlh1_expr %>% 
    mutate(point_color=ifelse(grepl(patient, A5_ID), "black", "grey"), 
           Label=ifelse(test = grepl(patient, A5_ID), 
                        yes = dplyr::recode(A5_ID, !!!setNames(a5_anno$PublicationID, a5_anno$A5_ID)), 
                        no = NA))
  
pj <- position_jitter(width=0.2, seed=100)
gg_mlh1_expr <- 
  ggplot(mlh1_expr, 
         aes(x="MLH1", y=log2cpm, color=point_color, label=Label)) + 
  geom_point(position = pj) + 
  geom_text_repel(mapping=aes(y=log2cpm),position = pj, color="black") +
  scale_color_identity() + 
  theme_bw() +
  xlab("") +
  ylab("CPM (log2)")


  

