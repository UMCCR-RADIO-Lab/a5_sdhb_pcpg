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

hgnc_symbols <- c("FOLH1")
#"ENSG00000086205"

##########################
# Fetch Ensembl gene ids #
##########################

hgnc_to_ensgid <- hgnc_to_ensgid_from_biomart(hgnc_symbols = GOI, use_cache = F, update_cache = F)

#######################
# Annotate expression #
#######################
  
exp <- wts_tcga_flynn_a5_lcpm.batch_removed[hgnc_to_ensgid$ensembl_gene_id,,drop=F] %>% 
as_tibble(rownames = "ensembl_gene_id") %>% 
pivot_longer(cols = -ensembl_gene_id, 
             names_to = "Sample", 
             values_to = "log2_cpm") %>% 
inner_join(wts_tcga_flynn_a5_anno) %>% 
  left_join(hgnc_to_ensgid)

exp <- exp %>%
  mutate(subtype=dplyr::recode(new_naming,
                              "A5 - Extraadrenal"="C1A1 (SDHx)",
                              "A5 - NF1"="C2A (Kinase)",
                              "A5 - Adrenal"="C1A1 (SDHx)",
                              "A5 - Head_neck"="C1A2 (SDHx-HN)",
                              "A5 - VHL"="C1B1 (VHL)",
                              "A5 - Extraadrenal_cardiac"="C1A2 (SDHx-HN)",
                              "A5 - Unspecified"="C1A1 (SDHx)",
                              "C2B2 (MAML)"="C2B2 (MAML3)"))

########
# Plot #
########

ggplot(exp, aes(x=subtype, y=log2_cpm)) +
  geom_boxplot() +
  #geom_jitter(mapping=aes(color=subtype, shape=Dataset), width = 0.2) +
  #scale_color_manual(values=subtype_cols) +
  geom_jitter(mapping=aes(color=Dataset), width = 0.2) +
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
  inner_join(a5_anno) %>% 
  mutate(Primary_Location_Simplified=gsub("_left|_right","", Primary_Location_Simplified))
  


jp <- position_jitter(width = 0.2, seed=10)
ggplot(exp %>% arrange(Primary_Location_Simplified, `Patient ID`), aes(x=Primary_Location_Simplified, y=log2_cpm, group=`Patient ID`)) +
  geom_boxplot(mapping=aes(group=NULL)) +
  geom_path(alpha=0.2, position=jp, linetype=2) +
  geom_point(mapping=aes(color=is_primary_or_met), position=jp) +
  facet_wrap("hgnc_symbol") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5,hjust=1))
