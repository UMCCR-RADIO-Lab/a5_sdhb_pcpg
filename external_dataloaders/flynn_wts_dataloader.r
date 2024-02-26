library(tidyverse)
library(SummarizedExperiment)
source("./public_data/data_loaders/fetch_zethoven_annotation.r")

data_loader_flynn_wts_counts <- function()
{
  message("Starting Flynn WTS data loader ...")
  
  base_dir="/g/data/pq08/projects/ppgl/public_data"
  data_dir=paste(base_dir, "wts/flynn/counts", sep="/")
  
  # load transcript counts from 2015 PPGL landscape paper
  flynn_counts <- read_tsv(paste(data_dir,"flynn_2015_ppgl_rnaseq.tsv",sep="/"), show_col_types = F)
  flynn_counts <- flynn_counts %>% dplyr::rename(ensembl_gene_id=`#Feature`) %>% dplyr::select(-Description, -Gene.Symbol)
  
 #Load zethoven annotation
  zethoven_annotation_flynn <- fetch_zethoven_annotation(datasets="Flynn")
  
  # annotation_file <- paste(base_dir,"annotation/tcga_comete_flynn_bulk_metadata.tsv", sep="/")
  # message("Loading sample annotation from ", annotation_file)
  # flynn_annotation <- read_tsv(annotation_file, show_col_types = F) %>%
  #   filter(Dataset %in% c("Flynn")) %>%
  #   dplyr::select(Barcode = "Sample.raw", Cluster,Dataset, TCGA_Cluster, Genotype, Malignancy) %>%
  #   mutate(Barcode = gsub("\\.", "-", Barcode)) %>%
  #   mutate(project_id = "PCPG") %>%
  #   left_join(zethoven_cluster_naming %>% dplyr::select(final, historical_b_variant) %>% dplyr::rename(new_naming=final), by=c("Cluster"="historical_b_variant"))
  
  #Remove samples with no counts from annotation
  zethoven_annotation_flynn <- zethoven_annotation_flynn %>% filter(Sample %in% colnames(flynn_counts))
  
  #Remove samples duplicated in the A5 cohort ("V-PH-23T"="E019")
  zethoven_annotation_flynn <- zethoven_annotation_flynn %>% filter(!(Sample %in% c("V-PH-23T")))
  
  #Remove samples with no annotation from counts
  flynn_counts <- flynn_counts %>% dplyr::select(c(ensembl_gene_id, !!zethoven_annotation_flynn$Sample))
  
  assign("flynn_wts_counts", flynn_counts, globalenv())
  assign("flynn_wts_anno", zethoven_annotation_flynn, globalenv())
  
  message("Loaded objects into global environment: 
          flynn_wts_counts (dataframe of feature counts),
          flynn_wts_anno (sample/clinical annotation)")
  
  message("Completed Flynn WTS data loader.")
}
message("Created data loader function data_loader_flynn_wts_counts()")