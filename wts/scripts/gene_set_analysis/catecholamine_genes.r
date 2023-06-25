#######################
# Import data loaders #
#######################

setwd("/g/data/pq08/projects/ppgl/")

source("./a5/wts/scripts/data_loaders/a5_wts_dataloader.r")

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#############
# Load data #
#############

htseq_outs <- "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/truseq/gene/"
data_loader_a5_wts_counts(htseq_outs)

#########
# Plots #
#########

plot_data <- cpm(a5_wts_dge_list$qc_ok,log = T) %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), se="_", extra = "merge") %>% 
  filter(symbol %in% c("TH","DBH","PNMT")) %>% #, "SLC6A2"
  pivot_longer(cols = c(-ensgid,-symbol), names_to = "A5_ID", values_to = "log2_cpm") %>% 
  group_by(symbol) %>% mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm))/sd(log2_cpm)) %>% 
  inner_join(a5_anno %>% 
               filter(!Exclude_RNA) %>% 
               dplyr::select(A5_ID, differential_group_anatomy, Catecholamine_profile)) %>% 
  mutate(symbol=factor(symbol, levels=rev(c("TH","DBH","PNMT", "SLC6A2")))) %>% 
  arrange(symbol, differential_group_anatomy, A5_ID) %>% 
  mutate(label=ifelse(A5_ID %in% c("E152-1","E186-1","E182-1","E188-1",
                                   "E225-1","E225-2", 
                                   "E125-1","E133-1", 
                                   "E146-1", "E146-2", "E160-1",
                                   "E141-1","E200-1"), A5_ID,NA))


ggplot(plot_data %>% filter(Catecholamine_profile %in% c("Norepinephrine", "Non secreting", "Dopamine")) %>% 
         mutate(Catecholamine_profile=factor(Catecholamine_profile, levels=c("Non secreting", "Dopamine", "Norepinephrine"))), 
       aes(y=symbol, 
           x=log2_cpm_z, 
           color=differential_group_anatomy,
           label=label,
           group=A5_ID)) + 
  geom_point() +
  geom_path(alpha=0.5,linetype="22") +
  #ggrepel::geom_text_repel(nudge_y = 0.2) +       
  #geom_text(nudge_y = 0.2) +
  facet_wrap("Catecholamine_profile", ncol=3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme_bw()


# 
# E188/E225 - low levels cats (2-4xULN), low level TH (large tumours ~6cm, so maybe just low output) 
# E146 - low TH, 14cm tumour plus mets, 2-6xULN, so again probably low output
# E125/133: marked non-secreting, elevated TH, limited anno for E125, E133 notes:
#   Plasma fractionated catecholamines – mild elevations of Epi 61 (5-58) and dopamine 39 (5-23). Methoxytyramine: <20 pg/ml (UL: <20 pg/mL) 
# E186 - marked non-secreting, elevated TH, limited anno
# E165 - incorrect annotation? Marked as epi but blood levels are high for norepi, expression PNMT low, I've changed to norepi 
# E160 - case we discussed yesterday, cats probably from para-aortic second primary
# E152 - marked as dopamine, low TH, high DBH, sketchy clinical notes:
# "September 24, 2013 (outside NIH): Urine catecholamines were within normal limits. Twenty-four hour urine dopamine was 205.6 mcg/24 hr (0-498), 24 hr urine norepinephrine was 10.6 mcg/24 hr (12.1-85.5), and 24 hr urine epinephrine was 4.3 mcg/24 hr (1.7-22.4). Labs drawn at NIH before operation on 11/6/13 found elevated Dopamine >25 pg/mL (range: 0-24). All other levels within normal limits. No record of methoxytyramine at NIH."


cat_path_genes <- tibble(ensgid=c("ENSG00000180176","ENSG00000123454","ENSG00000141744"), 
symbol=c("TH","DBH","PNMT"))

plot_data <- wts_tcga_flynn_a5_lcpm.batch_removed %>% 
  as_tibble(rownames = "ensgid") %>% 
  inner_join(cat_path_genes) %>%
  pivot_longer(cols = c(-ensgid,-symbol), names_to = "Sample", values_to = "log2_cpm") %>% 
  group_by(symbol) %>% mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm))/sd(log2_cpm)) %>% 
  inner_join(wts_tcga_flynn_a5_anno %>% 
               dplyr::select(Sample, new_naming, Dataset)) %>% 
  mutate(symbol=factor(symbol, levels=rev(c("TH","DBH","PNMT", "SLC6A2")))) %>% 
  mutate(subtype=dplyr::recode(new_naming,
                               "A5 - Extraadrenal"="C1A1 (SDHx)",
                               "A5 - NF1"="C2A (Kinase)",
                               "A5 - Adrenal"="C1A1 (SDHx)",
                               "A5 - Head_neck"="C1A2 (SDHx-HN)",
                               "A5 - VHL"="C1B1 (VHL)",
                               "A5 - Extraadrenal_aortic"="C1A2 (SDHx-HN)",
                               "A5 - Unspecified"="C1A1 (SDHx)",
                               "C2B2 (MAML)"="C2B2 (MAML3)")) %>% 
  arrange(symbol, Sample) %>% 
  filter(subtype %in% c("C1A1 (SDHx)", "C1A2 (SDHx-HN)","C2A (Kinase)"), Dataset %in% c("TCGA", "A5"))
                                   
ggplot(plot_data, 
       aes(y=symbol, 
           x=log2_cpm_z, 
           color=subtype,
           group=Sample)) + 
  geom_point() +
  geom_path(alpha=0.5,linetype="22") +
  scale_color_manual(values = subtype_cols) +
  #ggrepel::geom_text_repel(nudge_y = 0.2) +       
  #geom_text(nudge_y = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme_bw()
