###########################################
###########################################
##                                       ##
##   E167/E169 mlh Expression           ##
##                                       ##
##   Script for producing a scatter      ##
##   plots showing mlh expression       ##
##                                       ##
##  Author: Aidan Flynn                  ##
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
if (!exists("a5_wts_dge_list")) {
  htseq_outs <- "./a5/wts/analysis/htseq/neb/gene"
  data_loader_a5_wts_counts(count_file_dir=htseq_outs)
}

#clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
if(!exists("a5_anno")) {
  data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)
}

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

##############
# Expression #
##############

expr <- cpm(a5_wts_dge_list$SDHB, log = T) %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  filter(ensgid_symbol %in% c("ENSG00000170430.10_MGMT", "ENSG00000076242.16_MLH1")) %>% 
  pivot_longer(cols=-ensgid_symbol, 
               names_to = "A5_ID", 
               values_to = "log2cpm")


expr <- expr %>% 
  left_join(a5_anno %>%  dplyr::select(A5_ID, `Patient ID`, PublicationID)) %>% 
  mutate(point_color=case_match(PublicationID, 
                                "E169-M1" ~ ColorPalette[["LightOrange1"]],
                                "E169-M2" ~ ColorPalette[["LightRed3"]],
                                "E167-M1" ~ ColorPalette[["LightGreen1"]],
                                "E167-M2" ~ ColorPalette[["DarkGreen1"]],
                                .default = "black"), 
         Label=ifelse(test = PublicationID %in% c("E169-M1", "E169-M2", "E167-M1", "E167-M2"), 
                      yes = dplyr::recode(A5_ID, !!!setNames(a5_anno$PublicationID, a5_anno$A5_ID)), 
                      no = NA),
         line_alpha = ifelse(!is.na(Label), 0.5, 0),
         symbol=gsub("ENSG.+_","",ensgid_symbol))


pj <- position_jitter(width=0.2, seed=100)

gg_expr <- list()
for (gene in c("ENSG00000170430.10_MGMT", "ENSG00000076242.16_MLH1"))
{
gg_expr[[gene]] <- 
  ggplot(expr %>%  filter(ensgid_symbol==gene), 
         aes(x=symbol, y=log2cpm, color=point_color, label=Label)) + 
  geom_point(position = pj) + 
  geom_line(aes(group=`Patient ID`, alpha= line_alpha),
            position = pj,
            linetype=2, 
            color="black") +
  geom_text_repel(mapping=aes(y=log2cpm),position = pj, color="black") +
  scale_color_identity() + 
  scale_alpha_identity() +
  theme_bw() +
  xlab("") +
  ylab("CPM (log2)")
}

