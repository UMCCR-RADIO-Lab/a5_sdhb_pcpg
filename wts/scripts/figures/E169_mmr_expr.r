###########################################
###########################################
##                                       ##
##   E169 MMR Expression                 ##
##                                       ##
##   Script for producing a scatter      ##
##   plot showing MMR gene expression    ##
##   highlighting E169                   ##
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
htseq_outs <- "./a5/wts/analysis/htseq/truseq/gene"
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

dna_repair_genes <- read.delim("/g/data/pq08/reference/gene_lists/dna_damage_repair_pearl_etal_s3.txt")

#Gene symbol lookup table
all_genes <- unique(unlist(purrr::map(.x = a5_wts_dge_list, .f = rownames)))
symbol_lookup <- setNames(all_genes, 
                          gsub("ENSG[0-9]+([.][0-9]+)?_(.+)",
                               "\\2", 
                               all_genes))

mmr_genes <- dna_repair_genes %>%  
  filter(Pathway.2 %in% c("MutL homologs", "Mismatch and loop recognition factors")) %>% 
  pull(Gene.ID)

mmr_genes <- symbol_lookup[mmr_genes]

mmr_expr <- cpm(a5_wts_dge_list$SDHB, log = T) %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  filter(ensgid_symbol %in% mmr_genes) %>% 
  pivot_longer(cols=-ensgid_symbol, 
               names_to = "A5_ID", 
               values_to = "log2_cpm") %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), extra = "merge", se="_") %>% 
  group_by(symbol) %>% mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm))/sd(log2_cpm)) 

patients = c("E167","E169")

mmr_expr <- mmr_expr %>% 
    mutate(point_color=ifelse(gsub("-.","",A5_ID) %in% patients, A5_ID, "Other"),
           PublicationID = dplyr::recode(A5_ID, !!!setNames(a5_anno$PublicationID, a5_anno$A5_ID)),
           Label=ifelse(test = gsub("-.","",A5_ID) %in% patients, 
                        yes = dplyr::recode(A5_ID, !!!setNames(a5_anno$PublicationID, a5_anno$A5_ID)), 
                        no = NA)) %>%  arrange(desc(point_color))
  
pj <- position_jitter(width=0.1, seed=100)
gg_mmr_expr <- 
  ggplot(mmr_expr, 
         aes(x=symbol, y=log2_cpm_z, color=Label, label=Label)) + 
  geom_point(position = pj) + 
  #geom_text_repel(mapping=aes(y=log2_cpm_z),position = pj, color="black") +
  scale_color_manual(values = c("E167-M1"="#cbc262","E167-M2"="#af6884","E169-M1"="#637b62","E169-M2"="#eb7b54"), na.value ="#d6d6d4") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  xlab("") +
  ylab("CPM (log2, z-scaled)")


  

