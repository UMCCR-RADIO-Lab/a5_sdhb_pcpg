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

hgnc_symbols <- c("RIPK4", "RPRM", "DRG2", "USE1", "COX17","SULT4A1")


##########################
# Fetch Ensembl gene ids #
##########################

hgnc_to_ensgid <- hgnc_to_ensgid_from_biomart(hgnc_symbols = hgnc_symbols, use_cache = T, update_cache = F)

#######################
# Annotate expression #
#######################
available_ids <- intersect(rownames(wts_tcga_flynn_a5_lcpm.batch_removed), hgnc_to_ensgid$ensembl_gene_id)
exp <- wts_tcga_flynn_a5_lcpm.batch_removed[available_ids,,drop=F] %>% 
  as_tibble(rownames = "ensembl_gene_id") %>% 
  pivot_longer(cols = -ensembl_gene_id, 
               names_to = "Sample", 
               values_to = "log2_cpm") %>% 
  group_by(ensembl_gene_id) %>% 
  mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm, na.rm=T))/sd(log2_cpm, na.rm=T)) %>% 
  inner_join(wts_tcga_flynn_a5_anno) %>% 
  left_join(hgnc_to_ensgid)

exp <- exp %>%
  mutate(subtype=dplyr::recode(new_naming,
                               "A5 - Extraadrenal"="C1A1 (SDHx)",
                               "A5 - NF1"="C2A (Kinase)",
                               "A5 - Adrenal"="C1A1 (SDHx)",
                               "A5 - Head_neck"="C1A2 (SDHx-HN)",
                               "A5 - VHL"="C1B1 (VHL)",
                               "A5 - Extraadrenal_mediastinum"="C1A2 (SDHx-HN)",
                               "A5 - Unspecified"="C1A1 (SDHx)",
                               "C2B2 (MAML)"="C2B2 (MAML3)"),
         Malignancy=replace_na(Malignancy, "Unknown")) 

#exp <- exp %>% filter(subtype == "C1A1 (SDHx)")

exp <- exp %>% filter(Sample != "TCGA-SR-A6MX-05")

exp <- exp %>% mutate(Dataset = case_match(Dataset, 
                                           "TCGA" ~ "TCGA/Flynn et al.",
                                           "Flynn" ~ "TCGA/Flynn et al.",
                                           .default = Dataset))

########
# Plot #
########

ggplot(exp, aes(x=ATRX_mutation, y=log2_cpm_z)) +
  geom_boxplot() +
  #geom_jitter(mapping=aes(color=subtype, shape=Dataset), width = 0.2) +
  geom_jitter(mapping=aes(color=subtype), width = 0.2) +
  scale_color_manual(values=subtype_cols) +
  facet_grid(Dataset~hgnc_symbol, scales="free_y") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=1))  

test_results <- list()
for (gene in unique(exp$hgnc_symbol))
{
  
  for (dataset in unique(exp$Dataset))
  {
    #message("testing:", gene, " in ", dataset)
    test_results[[gene]][[dataset]] <- t.test(log2_cpm ~ ATRX_mutation, exp %>%  filter(hgnc_symbol == gene, Dataset == dataset))
    message(gene, "/", dataset, "- pval:", test_results[[gene]][[dataset]]$p.value)
  }
}