library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_gene_cn()

source("/g/data/pq08/projects/ppgl/a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
data_loader_a5_wts_counts(count_file_dir = "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/neb/gene", 
                          count_file_pattern = ".gene.counts")

log2cpm <- edgeR::cpm(a5_wts_dge_list[["SDHB"]], log=T) %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), sep="_", extra = "merge") %>% 
  pivot_longer(cols = c(-ensgid,-symbol), names_to = "A5_ID", values_to = "log2_cpm")

plot.data <- log2cpm %>% 
  filter(symbol=="TERT") %>%  
  dplyr::select(A5_ID, log2_cpm)  %>% 
  inner_join(A5_gene_cn %>% 
               mutate(A5_ID = gsub("-T0", "-", A5_ID)) %>% 
               filter(gene=="TERT") %>%  
               dplyr::select(A5_ID, maxCopyNumber)) %>% 
  mutate(maxCopyNumber=round(maxCopyNumber,0)) %>% 
  inner_join(a5_anno %>% 
               dplyr::select(A5_ID, TERT_ATRX_Mutation) %>% 
               mutate(TERT_ATRX_Mutation=ifelse(A5_ID=="E171-1", "WT", TERT_ATRX_Mutation)))

ggplot(plot.data, aes(x=maxCopyNumber, y=log2_cpm, color=TERT_ATRX_Mutation)) + geom_jitter(width = 0.1) + 
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  scale_x_continuous(expand = c(0.1,0.1)) +
  scale_color_manual(values=driver_cols) +
  xlab("TERT Copy Number") + 
  ylab("TERT Expression (log2 CPM)") 
