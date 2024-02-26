##Dataloader function for TCGA methylation data

library("dplyr")
library("tibble")
source("./public_data/data_loaders/fetch_zethovan_annotation.r")

# Data obtained using:
# query_all_tumour <- GDCquery(project = "TCGA-PCPG",
#                              data.category = "DNA Methylation",
#                              data.type = "Methylation Beta Value",
#                              legacy = F)
#
# GDCdownload(query_all_tumour, method = "api", files.per.chunk = 10)
#
# data_tumour <- GDCprepare(query_all_tumour, save = T, save.filename = "tcga_pcpg_beta_values.rdata")


data_loader_tcga_methylation <- function()
{
  message("Starting TCGA methylation data loader ...")
  
  base_dir="/g/data/pq08/projects/ppgl/public_data"
  data_dir=paste(base_dir, "methylation/tcga", sep="/")
 
   #Load dataset annotation
  zethovan_annotation_tcga <- fetch_zethovan_annotation(datasets="TCGA")
  
  #Add batch for dataset merging
  zethovan_annotation_tcga <- zethovan_annotation_tcga %>% mutate(methylation_batch="TCGA_450k")
  
  message("Loading RangedSummarizedExperiment object ...")
  #create a new temporary environment so the loaded object called 'data' doesn't crash into anything in global
  tcga <- new.env(parent = globalenv()) 
  load(paste(data_dir,"tcga_pcpg_beta_values.rdata", sep="/"),envir = tcga)
  meth_data <- tcga$data #reassign into local scope
  rm(tcga)
  
  message("Extracting beta-values ...")
  beta_vals <- assay(meth_data)
  beta_vals <- beta_vals[!is.na(rowSums(beta_vals)),] #Remove probes with NAs
  #shorten TCGA codes from "TCGA-RW-A67X-01A-11D-A35E-05" to "TCGA-RW-A67X-01" to match annotation
  colnames(beta_vals) <- substr(colnames(beta_vals),1,15)
  beta_vals <- data.frame(beta_vals, check.names = F) %>% rownames_to_column("Probe")
  
  message("Excluding samples without annotation from beta/m-vals ...")
  missing <- setdiff(colnames(beta_vals)[-1], zethovan_annotation_tcga$Sample)
  beta_vals <- beta_vals %>% dplyr::select(c("Probe", !!zethovan_annotation_tcga$Sample))
  message("Removed: ", toString(missing))
  
  message("Computing M-values ...")
  m_vals <- beta_vals %>% mutate(across(.cols = -Probe,.fns = ~log2(.x/(1-.x))))
  
  assign("tcga_methylation_450k_rse", meth_data, globalenv())
  assign("tcga_methylation_450k_beta", beta_vals, globalenv())
  assign("tcga_methylation_450k_mval", m_vals, globalenv())
  assign("tcga_methylation_450k_anno", zethovan_annotation_tcga, globalenv())
  
  message("Loaded objects into global environment: 
          tcga_methylation_450k_exprset (GDC RangedSummarizedExperiment Object),
          tcga_methylation_450k_beta (Beta values),
          tcga_methylation_450k_mval (M-values),
          tcga_methylation_450k_anno (Clinical annotation)")
  
  message("Completed TCGA methylation data loader.")
}
message("Created data loader function data_loader_tcga_methylation()")



