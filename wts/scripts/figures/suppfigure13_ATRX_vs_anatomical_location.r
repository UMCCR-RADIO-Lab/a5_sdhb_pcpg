library(ggplot2)
library(patchwork)

setwd("/g/data/pq08/projects/ppgl")

#######################
# Import data loaders #
#######################

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
source("./a5/wts/scripts/data_loaders/a5_wts_dataloader.r")

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#############
# Load data #
#############

#Load A5 Samples
htseq_outs <- "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/neb/gene/"
data_loader_a5_wts_counts(htseq_outs)

a5_wts_lcpm_list <- lapply(a5_wts_dge_list, cpm, log = T)

#####################
# Prepare plot data #
#####################

hgnc_symbols <- c("ATRX")

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
  inner_join(a5_anno)


########
# Plot #
########

jp <- position_jitter(width = 0.2, seed=10)
ggplot(exp %>% arrange(primary_location_plotting, `Patient ID`), 
       aes(x=TERT_ATRX_Mutation, y=log2_cpm, group=`Patient ID`)) +
  geom_boxplot(mapping=aes(group=NULL), outlier.alpha = 0) +
  geom_path(alpha=0.2, position=jp, linetype=2) +
  geom_point(mapping=aes(colour=differential_group_sampletype_strict), position=jp) +
  ggrepel::geom_text_repel(mapping=aes(label=A5_ID), position=jp, size=4) +
  scale_color_manual(values=sampletype_strict_cols) +
  theme_bw()
