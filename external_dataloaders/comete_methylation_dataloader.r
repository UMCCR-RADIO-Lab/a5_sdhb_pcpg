##Dataloader function for COMETE methylation data

library("dplyr")
library("tibble")
source("./public_data/data_loaders/fetch_zethoven_annotation.r")

download_comete_meth_geo <- function(output_dir)
{
  message("Retrieving data from GEO ...")
  library(GEOquery)
  
  meth_27k <- getGEO("GSE39198")
  meth_450k <- getGEO("GSE43293")
  
  saveRDS(meth_27, paste0(output_dir,"/comete_27.rds"))
  saveRDS(meth_450, paste0(output_dir,"comete_450.rds"))
  
  message("Data for 27k and 450k arrays downloaded to ",output_dir, " as comete_27.rds and comete_450.rds, respectively")
  
}

data_loader_comete_methylation <- function(remove_27k_availaible_on_450k = T, only_zethoven_annotated=T)
{
  message("Starting COMETE methylation data loader ...")
  
  base_dir="/g/data/pq08/projects/ppgl"
  data_dir=paste(base_dir, "public_data/methylation/comete", sep="/")
  
  message("Loading COMETE 27k GEO ExpressionSet object ...")
  meth_27k <- readRDS(paste(data_dir,"comete_27.rds", sep="/"))
  message("Loading COMETE 450k GEO ExpressionSet object ...")
  meth_450k <- readRDS(paste(data_dir,"comete_450.rds", sep="/"))
  
  message("Extracting beta-values ...")
  meth_27k_mat <- meth_27k$GSE39198_series_matrix.txt.gz@assayData$exprs
  meth_450k_mat <- meth_450k$GSE43293_series_matrix.txt.gz@assayData$exprs
  
  message("Extracting clinical annotation ...")
  meth_27k_clinical <- meth_27k$GSE39198_series_matrix.txt.gz@phenoData@data
  meth_450k_clinical <- meth_450k$GSE43293_series_matrix.txt.gz@phenoData@data
  
  #Load dataset annotation
  comete_annotation <- fetch_zethoven_annotation(datasets="E-MTAB-733") %>% mutate(Sample=Alias)
  
  # annotation_file <- paste(base_dir,"public_data/annotation/tcga_comete_flynn_bulk_metadata.tsv", sep="/")
  # message("Loading sample annotation from ", annotation_file)
  # all_annotation <- read.delim(annotation_file, sep="\t")
  # comete_annotation <- all_annotation %>% 
  #   filter(Dataset=="E-MTAB-733") %>% 
  #   mutate(Sample=Alias)  %>%
  #   left_join(zethoven_cluster_naming %>% dplyr::select(final, historical_b_variant) %>% dplyr::rename(new_naming=final), by=c("Cluster"="historical_b_variant"))
  # 
           
  zethoven_comete_annotation_27k <- comete_annotation %>% filter(Sample %in% meth_27k_clinical$title)
  zethoven_comete_annotation_450k <- comete_annotation %>% filter(Sample %in% meth_450k_clinical$title)
  
  message("Replacing GSE codes with sample names ...")
  if(!all(colnames(meth_27k_mat) == rownames(meth_27k_clinical))) {
    stop("comete_methylation_dataloade - Sample annotation and array-data column order sanity check failed")}
  colnames(meth_27k_mat) <- meth_27k_clinical$title
  colnames(meth_450k_mat) <- meth_450k_clinical$title
  
  message("Removing probes which are missing for any sample ...")
  n_all_27k <- dim(meth_27k_mat)[1]
  n_all_450k <- dim(meth_450k_mat)[1]
  meth_27k_mat <- meth_27k_mat[!is.na(rowSums(meth_27k_mat)),]
  meth_450k_mat <- meth_450k_mat[!is.na(rowSums(meth_450k_mat)),]
  n_postfilter_27k <- dim(meth_27k_mat)[1]
  n_postfilter_450k <- dim(meth_450k_mat)[1]
  message("Removed: ",n_all_27k-n_postfilter_27k, " probes from 27k; ",n_all_450k-n_postfilter_450k," probes from 450k;" )
  
  if(remove_27k_availaible_on_450k)
  {
    message("Removing samples run on both 27k and 450k from 27k beta/m-values and clinical annotation (data will remain in ExpressionSet Object)")
    available_as_450k <- gsub("Reanalyzed by: ", "", meth_27k_clinical$relation) %in% rownames(meth_450k_clinical)
    meth_27k_mat <- meth_27k_mat[,!available_as_450k]
    meth_27k_clinical <- meth_27k_clinical[!available_as_450k,]
    zethoven_comete_annotation_27k <- zethoven_comete_annotation_27k %>% filter(!(Sample %in% zethoven_comete_annotation_450k$Sample))
  }
  
  #Add batch for dataset merging
  zethoven_comete_annotation_27k <- zethoven_comete_annotation_27k %>% mutate(methylation_batch="COMETE_27k")
  zethoven_comete_annotation_450k <- zethoven_comete_annotation_450k %>% mutate(methylation_batch="COMETE_450k")
  
  #Covert to dataframe
  meth_27k_df <- meth_27k_mat %>%
    data.frame(check.names = F)%>%
    rownames_to_column("Probe")
  
  meth_450k_df <- meth_450k_mat %>%
    data.frame(check.names = F)%>%
    rownames_to_column("Probe")
  
  if(only_zethoven_annotated)
  {
    meth_27k_df <- meth_27k_df %>% dplyr::select(c("Probe", !!zethoven_comete_annotation_27k$Sample))
    meth_450k_df <- meth_450k_df %>% dplyr::select(c("Probe", !!zethoven_comete_annotation_450k$Sample))
  }
  
  message("Computing M-values ...")
  meth_27k_m_vals <- meth_27k_df %>% mutate(across(.cols = -Probe,.fns = ~log2(.x/(1-.x))))
  meth_450k_m_vals <- meth_450k_df %>% mutate(across(.cols = -Probe,.fns = ~log2(.x/(1-.x))))
  
  assign("comete_methylation_27k_exprset", meth_27k, globalenv())
  assign("comete_methylation_27k_beta", meth_27k_df, globalenv())
  assign("comete_methylation_27k_mval", meth_27k_m_vals, globalenv())
  
  assign("comete_methylation_450k_exprset", meth_450k, globalenv())
  assign("comete_methylation_450k_beta", meth_450k_df, globalenv())
  assign("comete_methylation_450k_mval", meth_450k_m_vals, globalenv())
  
  if(!only_zethoven_annotated)
  {
    #add empty rows for un-annotated samples
    zethoven_comete_annotation_27k <- bind_rows(zethoven_comete_annotation_27k, 
                                               data.frame(Sample=setdiff(meth_27k_clinical$Sample, 
                                                                         zethoven_comete_annotation_27k$Sample)))
    
    zethoven_comete_annotation_450k <- bind_rows(zethoven_comete_annotation_450k, 
                                                data.frame(Sample=setdiff(meth_450k_clinical$Sample, 
                                                                          zethoven_comete_annotation_450k$Sample)))
  }
  
  assign("comete_methylation_27k_anno", zethoven_comete_annotation_27k, globalenv())
  assign("comete_methylation_450k_anno", zethoven_comete_annotation_450k, globalenv())
  
  
  
  message("Loaded objects into global environment: 
          comete_methylation_27k_exprset (GEO ExpressionSet object for 27k arrays),
          comete_methylation_27k_beta (Beta values for 27k arrays),
          comete_methylation_27k_mval (M-values for 27k arrays),
          comete_methylation_27k_anno (Clinical annotation for 27k arrays),
          comete_methylation_450k_exprset (GEO ExpressionSet object for 450k arrays),
          comete_methylation_450k_beta (Beta values for 450k arrays),
          comete_methylation_450k_mval (M-values for 450k arrays),
          comete_methylation_450k_anno (Clinical annotation for 450k arrays)")
  
  message("Completed COMETE methylation data loader.")
}
message("Created data loader function data_loader_comete_methylation()")



