############################################
# Import data loaders and run data mergers #
############################################

setwd("/g/data/pq08/projects/ppgl")

#WTS
source("./a5/tertiary/scripts/data_mergers/combine_tcga_flynn_a5_wts_data.r")
source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

############################
# Define genes of interest #
############################

hgnc_symbols <- c("RPRM", "DRG2")
#"ENSG00000086205"

##########################
# Fetch Ensembl gene ids #
##########################

hgnc_to_ensgid <- hgnc_to_ensgid_from_biomart(hgnc_symbols = hgnc_symbols, use_cache = T, update_cache = F)

#######################
# Annotate expression #
#######################
  
exp <- wts_tcga_flynn_a5_lcpm.batch_removed[hgnc_to_ensgid$ensembl_gene_id,,drop=F] %>% 
as_tibble(rownames = "ensembl_gene_id") %>% 
pivot_longer(cols = -ensembl_gene_id, 
             names_to = "Sample", 
             values_to = "log2_cpm") %>%
  group_by(ensembl_gene_id) %>% 
  mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm, na.rm=T))/sd(log2_cpm, na.rm=T)) %>% 
inner_join(wts_tcga_flynn_a5_anno) %>% 
  left_join(hgnc_to_ensgid)

exp <- exp %>%
  mutate(Malignancy=replace_na(Malignancy, "Unknown")) %>% 
  dplyr::rename("subtype"="new_naming") 

########
# Plot #
########

ggplot(exp, aes(x=subtype, y=log2_cpm)) +
  geom_boxplot() +
  #geom_jitter(mapping=aes(color=subtype, shape=Dataset), width = 0.2) +
  geom_jitter(mapping=aes(color=Malignancy, shape=Dataset), width = 0.2) +
  #scale_color_manual(values=subtype_cols) +
  facet_wrap("hgnc_symbol") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=1))


######
# A5 #
######

a5_wts_lcpm_list <- lapply(a5_wts_dge_list, cpm, log = T)



#Gene symbol lookup table
all_genes <- unique(unlist(purrr::map(.x = a5_wts_dge_list, .f = rownames)))
symbol_lookup <- setNames(all_genes, 
                          gsub("ENSG[0-9]+([.][0-9]+)?_(.+)",
                               "\\2", 
                               all_genes))



exp <- a5_wts_lcpm_list[["SDHB"]][symbol_lookup[hgnc_symbols],,drop=F] %>% 
  as_tibble(rownames = "ensembl_gene_id_symbol") %>% 
  separate(ensembl_gene_id_symbol, into=c("ensembl_gene_id", "hgnc_symbol"), sep="_", extra = "merge") %>% 
  pivot_longer(cols = c(-ensembl_gene_id,-hgnc_symbol), 
               names_to = "A5_ID", 
               values_to = "log2_cpm") %>% 
  group_by(ensembl_gene_id) %>% 
  mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm, na.rm=T))/sd(log2_cpm, na.rm=T)) %>% 
  inner_join(a5_anno) %>% 
  mutate(Primary_Location_Simplified=gsub("_left|_right","", Primary_Location_Simplified),
         label=ifelse(abs(log2_cpm_z) > 2, A5_ID, NA),
         label=ifelse(A5_ID %in% c("E167-1","E167-2","E169-1","E225-2"), A5_ID, label))
  


jp <- position_jitter(width = 0.2, seed=10)
ggplot(exp %>% arrange(differential_group_sampletype_strict, `Patient ID`), 
       aes(x=TERT_ATRX_Mutation, y=log2_cpm, group=`Patient ID`)) +
  geom_boxplot(mapping=aes(group=NULL)) +
  geom_path(alpha=0.2, position=jp, linetype=2) +
  geom_point(mapping=aes(colour=differential_group_sampletype_strict), position=jp) +
  ggrepel::geom_text_repel(mapping=aes(label=A5_ID), position=jp, size=4) +
  facet_wrap("hgnc_symbol") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=1))


# jp <- position_jitter(width = 0.3, seed=10)
# ggplot(exp %>% arrange(differential_group_sampletype_strict, `Patient ID`) %>% 
#          mutate(resect_post_treatment = case_when(
#            grepl("CVD", Resection_post_dna_damaging_treatment) & grepl("MIBG", Resection_post_dna_damaging_treatment) ~ "CVD/MIBG",
#            grepl("CVD", Resection_post_dna_damaging_treatment) ~ "CVD",
#            grepl("MIBG", Resection_post_dna_damaging_treatment) ~ "MIBG",
#            grepl("Lutate", Resection_post_dna_damaging_treatment) ~ "Lutate",
#            grepl("TMZ", Resection_post_dna_damaging_treatment) ~ "TMZ",
#            grepl("Uncertain", Resection_pre_dna_damaging_treatment) == "Uncertain" ~ "No Data",
#            Resection_post_dna_damaging_treatment == "No" ~ "No",
#            TRUE ~ "Other"
#          )) %>% 
#            mutate(resect_before_treatment = case_when(
#            grepl("CVD", Resection_pre_dna_damaging_treatment) & grepl("MIBG", Resection_pre_dna_damaging_treatment) ~ "CVD/MIBG",
#            grepl("CVD", Resection_pre_dna_damaging_treatment) ~ "CVD",
#            grepl("MIBG", Resection_pre_dna_damaging_treatment) ~ "MIBG",
#            grepl("Lutate", Resection_pre_dna_damaging_treatment) ~ "Lutate",
#            grepl("TMZ", Resection_pre_dna_damaging_treatment) ~ "TMZ",
#            grepl("Uncertain", Resection_pre_dna_damaging_treatment) == "Uncertain" ~ "No Data",
#            Resection_pre_dna_damaging_treatment == "No" ~ "No",
#            TRUE ~ "Other"
#          )), 
#        aes(x="MGMT", y=log2_cpm_z, group=`Patient ID`)) +
#   geom_boxplot(mapping=aes(group=NULL), outlier.alpha = 0) +
#   geom_path(alpha=0.2, position=jp, linetype=2) +
#   geom_point(mapping=aes(colour=resect_post_treatment, shape=is_primary_or_met), position=jp) +
#   geom_text(mapping=aes(label=label), position=jp) +
#   scale_color_manual(values=c("No"="grey80", CVD="red", "CVD/MIBG"="orange", MIBG="green3", "Other"="purple", "TMZ"="gold", Lutate="blue", "No data"="black")) +
#   facet_wrap("hgnc_symbol") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=1))
# 
# 
# 
# ggplot(exp %>% mutate(Mutant=ifelse(TERT_ATRX_Mutation != "WT", "TERT/ATRX", "WT")) %>% 
#          mutate(color_code = case_when(
#            TERT_ATRX_Mutation != "WT" ~ TERT_ATRX_Mutation,
#            telhunter_log2_telcontentratio > 0.5 ~ "Long telomeres",
#            A5_ID == "E171-1" ~ "TERT over-expression",
#            TRUE ~ "WT"
#          )) %>% mutate(Ki67_Staining_numeric = as.numeric(gsub("<","",Ki67_Staining))), 
#        aes(x=paste(differential_group_sampletype, Mutant), y=Ki67_Staining_numeric, group=`Patient ID`)) +
#   geom_boxplot(mapping=aes(group=NULL)) +
#   geom_path(alpha=0.2, position=jp, linetype=2) +
#   geom_point(mapping=aes(colour=color_code, shape=differential_group_anatomy), position=jp) +
#   geom_text(mapping=aes(label=label), position=jp) +
#   facet_wrap("hgnc_symbol") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=1)) + 
#   scale_color_manual(values=c(driver_cols, `Long telomeres`=ColorPalette[["Yellow3"]], `TERT over-expression`=ColorPalette[["DarkRed1"]])) 
# 
# 
# ggplot(plot.data %>%  
#          mutate(Mutant=ifelse(TERT_ATRX_Mutation != "WT", "TERT/ATRX", "WT")) %>% 
#          filter(GeneSet %in% c("HALLMARK_GLYCOLYSIS","HALLMARK_MITOTIC_SPINDLE","HALLMARK_E2F_TARGETS",
#                                "HALLMARK_G2M_CHECKPOINT","JANSKY_cycling_neuroblast_genes","JANSKY_early_SCP_genes",
#                                "JANSKY_late_SCP_genes","ZETHOVEN_SCLCs", "ZETHOVEN_ben_met_up", 
#                                "ZETHOVEN_ben_met_down")) %>% 
#          mutate(color_code = case_when(
#            TERT_ATRX_Mutation != "WT" ~ TERT_ATRX_Mutation,
#            telhunter_log2_telcontentratio > 0.5 ~ "Long telomeres",
#            A5_ID == "E171-1" ~ "TERT over-expression",
#            TRUE ~ "WT"
#          )) %>% 
#          mutate(differential_group_sampletype = recode(differential_group_sampletype, "Other"="SFU/Local-Recurrance")), 
#        aes(x=paste(differential_group_sampletype, Mutant), y = GeneSet_score,  group = `Patient ID`)) + 
#   geom_boxplot(aes(group = NULL), outlier.alpha = 0) + 
#   geom_point(aes(color = color_code),
#              alpha = 0.7, position = pre_jitter) + 
#   geom_path(position = pre_jitter, linetype=2, alpha=0.2) +
#   facet_wrap("GeneSet", ncol = 6) + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1)) + 
#   scale_color_manual(values=c(driver_cols, `Long telomeres`=ColorPalette[["Yellow3"]], `TERT over-expression`=ColorPalette[["DarkRed1"]])) 
