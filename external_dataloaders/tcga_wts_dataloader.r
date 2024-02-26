library(tidyverse)
library(SummarizedExperiment)
source("./public_data/data_loaders/fetch_zethoven_annotation.r")

#Data obtained with:
# PCPG_prj <- "TCGA-PCPG"
# # Get gene expression from all tumour TCGA samples
# query_all_tumour <- TCGAbiolinks::GDCquery(project = PCPG_prj,
#                                            data.category = "Transcriptome Profiling",
#                                            data.type = "Gene Expression Quantification",
#                                            workflow.type  = "STAR - Counts",
#                                            sample.type = c("Primary Tumor", "Metastatic", "Additional - New Primary", "Recurrent Tumor","Solid Tissue Normal"),
#                                            legacy = F)
# 
# 
# TCGAbiolinks:::GDCdownload(query_all_tumour, method = "api", files.per.chunk = 10)
# data_tumour <- TCGAbiolinks:::GDCprepare(query_all_tumour, save = T, save.filename = "tcga_counts.Rdata")

data_loader_tcga_wts_counts <- function(remove_normals=T)
{
  message("Starting TCGA WTS data loader ...")
  
  base_dir="/g/data/pq08/projects/ppgl/public_data"
  data_dir=paste(base_dir, "wts/tcga/raw_data", sep="/")
  
  #Temporary environment to stop loaded object crashing in the global environment
  tcga <- new.env()
  load(paste(data_dir,"tcga_counts.rdata", sep="/"), envir = tcga)
  tcga_rse <- tcga$data #reassign into local scope
  rm(tcga)
  
  tumour_counts <- assay(tcga_rse)
  #shorten TCGA codes from "TCGA-RW-A67X-01A-11D-A35E-05" to "TCGA-RW-A67X-01" to match annotation
  colnames(tumour_counts) <- substr(colnames(tumour_counts),1,15)
  
  #Load dataset annotation
  zethoven_annotation_tcga <- fetch_zethoven_annotation(datasets="TCGA")
  
  # annotation_file <- paste(base_dir,"annotation/tcga_comete_flynn_bulk_metadata.tsv", sep="/")
  # message("Loading sample annotation from ", annotation_file)
  # tcga_anno <- read_tsv(annotation_file, show_col_types = F) %>%
  #   filter(Dataset %in% c("TCGA")) %>%
  #   dplyr::select(Barcode = "Sample.raw", Cluster,Dataset, TCGA_Cluster, Genotype, Malignancy) %>%
  #   mutate(Barcode = gsub("\\.", "-", Barcode)) %>%
  #   mutate(project_id = "PCPG") %>%
  #   left_join(zethoven_cluster_naming %>% 
  #               dplyr::select(final, historical_b_variant) %>% 
  #               dplyr::rename(new_naming=final), 
  #             by=c("Cluster"="historical_b_variant"))
  
  
  if(remove_normals)
  {
    zethoven_annotation_tcga <- zethoven_annotation_tcga %>% filter(Genotype != "Normal")
  }
  
  #Keep only annotated samples
  tumour_counts <- tumour_counts[,colnames(tumour_counts) %in% zethoven_annotation_tcga$Sample]
  
  #move rowname to column remove .XX from ensembl gene ID
  tumour_counts <- tumour_counts %>% data.frame(check.names = F) %>% 
    tibble::rownames_to_column("ensembl_gene_id") %>% 
    mutate(ensembl_gene_id=gsub("(ENSG[0-9]+)[.][0-9]+","\\1",ensembl_gene_id))

  assign("tcga_wts_counts", tumour_counts, globalenv())
  assign("tcga_wts_anno", zethoven_annotation_tcga, globalenv())
  
  message("Loaded objects into global environment: 
          tcga_wts_counts (dataframe of feature counts),
          tcga_wts_anno (sample/clinical annotation)")

  message("Completed TCGA WTS data loader.") 
}
message("Created data loader function data_loader_tcga_wts_counts()")