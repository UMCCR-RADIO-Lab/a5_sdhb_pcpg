setwd("/g/data/pq08/projects/ppgl/a5/")

library(tidyverse)
library(ggplot2)
library(patchwork)

###########
# Imports #
###########

#Modified version of gometh from missmethyl to use an offline cache for KEGG pathways
source("./methylation/scripts/go_meth_offline.r")

#Plotting helper functions
source("./methylation/scripts/helpers/plot_methylation.r")

###############
# DataLoaders #
###############

#load clinical annotation
source("./sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#Color schemes
source("./sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#load methylation array data
source("./methylation/scripts/data_loaders/a5_methylation_dataloader.r")
data_loader_a5_methylation_array(quickload = T, output_qc = F, normalisation="functional")

#load WTS data
source("/g/data/pq08/projects/ppgl/a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
htseq_outs <- "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/neb/gene"
data_loader_a5_wts_counts(count_file_dir=htseq_outs)

####################
# Fetch probe data #
####################

b_vals <- getBeta(a5_methylation_filtered)

probes =c("cg10896616", "cg23250191", "cg11625005", "cg10767223")

probe_annotation <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

region_probes <- data.frame(probe_annotation[probe_annotation$Name %in% probes, c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group")])

region_probes <- region_probes %>% 
  tidyr::separate_rows(UCSC_RefGene_Name, UCSC_RefGene_Group, sep=";") %>% 
  distinct()

gene_symbol <- unique(region_probes$UCSC_RefGene_Name)

probe_data <- b_vals[region_probes$Name, ,drop=F] %>% data.frame(check.names = F) %>% 
  tibble::rownames_to_column("probe_id") %>% 
  pivot_longer(cols=-probe_id, names_to = "A5_ID", values_to = "b_val")

#########################
# Fetch expression data #
#########################

a5_wts_lcpm_list <- lapply(a5_wts_dge_list, cpm, log = T) 

expr_data <- a5_wts_lcpm_list[["SDHB"]] %>%
  data.frame(check.names = F) %>% 
  tibble::rownames_to_column("ensgid_symbol") %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), extra = "merge", sep="_") %>% 
  dplyr::filter(symbol==gene_symbol) %>% 
  dplyr::select(-ensgid) %>%
  pivot_longer(cols=-symbol, 
               names_to = "A5_ID", 
               values_to = "log2cpm")
  

#####################
# Prepare plot data #
#####################

plot.data <- probe_data %>%  
  inner_join(region_probes %>% dplyr::select(probe_id=Name, chr, pos, UCSC_RefGene_Name, UCSC_RefGene_Group)) %>%  
  inner_join(a5_anno)

plot.data <- plot.data %>% inner_join(expr_data, by=c("UCSC_RefGene_Name"="symbol","A5_ID"="A5_ID"))

samples_to_label = c()
plot.data <- plot.data %>% mutate(log2cpm_z=(log2cpm-mean(log2cpm))/sd(log2cpm),
                                  label=ifelse(A5_ID %in% samples_to_label, A5_ID, NA))


plot.data <- plot.data %>%  
  mutate(probe_id=paste0(probe_id,"(",UCSC_RefGene_Group,")")) %>% 
  arrange(pos) %>% mutate(probe_id = factor(probe_id, levels=unique(.$probe_id)))

plot.data <- plot.data %>% mutate(z_class = cut(log2cpm_z, breaks = c(-10,1,10), label=c("Z < 1","Z > 1")))

plot.data <- plot.data %>%  mutate(TERT_ATRX_Mutation=ifelse(A5_ID == "E171-1", "WT", TERT_ATRX_Mutation))

########
# Plot #
########

# gg_meth <- ggplot(plot.data, aes(y=b_val, x=z_class, color=TERT_ATRX_Mutation)) +
#   geom_boxplot(aes(color=NULL)) +
#   geom_jitter(width = 0.2) +
#   #geom_text(aes(label=A5_ID)) + 
#   scale_color_manual(values=driver_cols) +
#   guides(color=guide_legend(override.aes = list(shape=19, label="", linewidth=c(0,0,0)))) +
#   theme_bw() + 
#   #theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1)) + 
#   ylab("Methylation (beta)") +
#   xlab("TERT expression (log2 CPM, z-score)") +
#   facet_wrap("probe_id", ncol=5) 

gg_meth_vs_expr <- ggplot(plot.data %>% filter(probe_id %in% c("cg11625005(TSS1500)", "cg10767223(TSS1500)")), 
                   aes(x=b_val, y=log2cpm, color=TERT_ATRX_Mutation)) +
  geom_point() +
  geom_text(aes(label=A5_ID)) + 
  scale_color_manual(values=driver_cols) +
  guides(color=guide_legend(override.aes = list(shape=19, label="", linewidth=c(0,0,0)))) +
  theme_bw() + 
  #theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1)) + 
  xlab("Methylation (beta)") +
  ylab("TERT expression (log2 CPM, z-score)") +
  facet_wrap("probe_id", ncol=1) 

gg_meth_density <- ggplot(plot.data %>% filter(probe_id %in% c("cg11625005(TSS1500)", "cg10767223(TSS1500)")), 
       aes(x=b_val, color=TERT_ATRX_Mutation)) +
  geom_density(adjust=0.8) +
  #geom_text(aes(label=A5_ID)) + 
  scale_color_manual(values=driver_cols) +
  theme_bw() + 
  #theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1)) + 
  facet_wrap("probe_id", ncol=1) 

gg_meth_vs_expr + gg_meth_density + plot_layout(ncol = 1, heights=c(1.5,0.5))
