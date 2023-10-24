#######################################
# Script to generate expression plots #
# of genes involved in catecholamine  # 
# biosyntheses annotated with         #
# clinical catecholamine profile      #
#                                     #
# Author: Aidan Flynn                 #
# Date: 14/08/2023                    #
#######################################

setwd("/g/data/pq08/projects/ppgl/")

#############
# Load data #
#############

source("./a5/tertiary/scripts/data_mergers/combine_tcga_flynn_a5_wts_data.r")

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#########
# Plots #
#########

ensgid_to_hgnc <- ensgid_to_hgnc_from_biomart(ens_gids = rownames(wts_tcga_flynn_a5_lcpm.batch_removed),
                                              use_cache = T) 

plot_data <- wts_tcga_flynn_a5_lcpm.batch_removed %>% 
  as_tibble(rownames = "ensembl_gene_id") %>% 
  inner_join(ensgid_to_hgnc) %>% 
  filter(hgnc_symbol %in% c("PAH", "QDPR", "TH", "DDC" , "DBH","SLC18A2","PNMT", "SLC6A2")) %>% #, 
  pivot_longer(cols = c(-ensembl_gene_id,-hgnc_symbol), names_to = "Sample", values_to = "log2_cpm") %>% 
  group_by(hgnc_symbol) %>% mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm))/sd(log2_cpm)) %>% 
  inner_join(a5_anno %>% 
               filter(!Exclude_RNA) %>% 
               dplyr::select(A5_ID, Primary_Location_Simplified, Catecholamine_profile, Catecholamine_class, differential_group_sampletype_strict, TERT_ATRX_Mutation),
             by=c("Sample"="A5_ID")) %>% 
  dplyr::rename("A5_ID"="Sample") %>% 
  mutate(hgnc_symbol=factor(hgnc_symbol, levels=rev(c("QDPR", "PAH", "TH", "DDC" , "DBH", "PNMT", "SLC18A2", "SLC6A2")))) %>% 
  arrange(hgnc_symbol, Primary_Location_Simplified, A5_ID) %>% 
  mutate(label=ifelse(A5_ID %in% c("E152-1","E186-1","E182-1","E188-1",
                                   "E225-1","E225-2", 
                                   "E125-1","E133-1", 
                                   "E146-1", "E146-2", "E160-1", "E165-1", "E152-1",
                                   "E141-1","E200-1"), A5_ID,NA)) %>% 
  mutate(Primary_Location_Simplified=gsub("_abdominal|_thoracic|_bladder|_cardiac|_[Ll]eft|_[Rr]ight",
                                         "",
                                         Primary_Location_Simplified)) %>% 
  filter(Catecholamine_profile != "No data") %>% 
  mutate(Catecholamine_profile=recode(Catecholamine_profile, 
         "Norepinephrine (minor dopamine)" = "Norepinephrine",
         "Norepinephrine and dopamine" = "Norepinephrine"))


ggplot(plot_data #%>%  
         # mutate(Catecholamine_profile=factor(Catecholamine_profile, 
         #                                     levels=c("Non secreting", "Dopamine", "Norepinephrine", "Mixed")))
       , 
       aes(y=hgnc_symbol, 
           x=log2_cpm_z, 
           color=Primary_Location_Simplified,
           shape=Catecholamine_class,
           label=label,
           group=A5_ID)) + 
  geom_point(position=position_jitter(width = 0, height = 0.1)) +
  geom_path(alpha=0.5,linetype="22") +
  scale_color_manual(values=location_cols) +
  scale_shape_manual(values=c("<3xULN"=6, ">3xULN"=2,"No data"=4)) +
  ggrepel::geom_text_repel(nudge_y = 0.2) +       
  facet_wrap("Catecholamine_profile", ncol=4) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme_bw() +
  coord_cartesian(xlim = c(-4,4))


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


##########################
# Super-set PNMT Example #
##########################

# cat_path_genes <- tibble(ensgid=c("ENSG00000180176", "ENSG00000132437","ENSG00000123454","ENSG00000141744"), 
# symbol=c("TH", "DDC","DBH","PNMT"))
# 
# plot_data <- wts_tcga_flynn_a5_lcpm.batch_removed %>% 
#   as_tibble(rownames = "ensgid") %>% 
#   inner_join(cat_path_genes) %>%
#   pivot_longer(cols = c(-ensgid,-symbol), names_to = "Sample", values_to = "log2_cpm") %>% 
#   group_by(symbol) %>% mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm))/sd(log2_cpm)) %>% 
#   inner_join(wts_tcga_flynn_a5_anno %>% 
#                dplyr::select(Sample, new_naming, Dataset)) %>% 
#   mutate(symbol=factor(symbol, levels=rev(c("TH","DDC","DBH","PNMT", "SLC6A2")))) %>% 
#   mutate(subtype=dplyr::recode(new_naming,
#                                "A5 - Extraadrenal"="C1A1 (SDHx)",
#                                "A5 - NF1"="C2A (Kinase)",
#                                "A5 - Adrenal"="C1A1 (SDHx)",
#                                "A5 - Head_neck"="C1A2 (SDHx-HN)",
#                                "A5 - VHL"="C1B1 (VHL)",
#                                "A5 - Extraadrenal_aortic"="C1A2 (SDHx-HN)",
#                                "A5 - Unspecified"="C1A1 (SDHx)",
#                                "C2B2 (MAML)"="C2B2 (MAML3)")) %>% 
#   arrange(symbol, Sample) %>% 
#   filter(subtype %in% c("C1A1 (SDHx)", "C1A2 (SDHx-HN)","C2A (Kinase)"), Dataset %in% c("TCGA", "A5"))
#                                    
# ggplot(plot_data, 
#        aes(y=symbol, 
#            x=log2_cpm_z, 
#            color=subtype,
#            group=Sample)) + 
#   geom_point() +
#   geom_path(alpha=0.5,linetype="22") +
#   scale_color_manual(values = subtype_cols) +
#   #ggrepel::geom_text_repel(nudge_y = 0.2) +       
#   #geom_text(nudge_y = 0.2) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   theme_bw()
