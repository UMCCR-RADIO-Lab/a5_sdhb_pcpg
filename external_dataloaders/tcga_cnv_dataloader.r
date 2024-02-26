##############################################
# Script for (down)loading TCGA CNA data set #
# Author: Aidan Flynn                        #
# Date: 01/02/2023                           #
# Languages: R                               #
##############################################

tcga_cnv_base_dir <- "/g/data/pq08/projects/ppgl"
tcga_ascat_dir=paste0(tcga_cnv_base_dir, "/public_data/copy_number_variation/tcga/ascat")
tcga_seg_dir = paste0(tcga_cnv_base_dir, "/public_data/copy_number_variation/tcga/segments")

download_tcga_ppgl_ascat <- function(output_dir) {
  library(TCGAbiolinks)
  if(dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)}
  
  query_ascat <- GDCquery(project = "TCGA-PCPG",
                          data.category = "Copy Number Variation",
                          data.type = c("Allele-specific Copy Number Segment"),
                          workflow.type  = c("ASCAT2"),
                          legacy = F)
  
  GDCdownload(query_ascat, method = "api", files.per.chunk = 5, directory=output_dir)
  
  cnv_ascat <- GDCprepare(query_ascat, 
                          save = T, 
                          save.filename = paste0(output_dir, "/tcga_ppgl_cnv_ascat.rdata"), 
                          directory=output_dir)
}

download_tcga_ppgl_cna_seg <- function(output_dir) {
  library(TCGAbiolinks)
  if(dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)}
  
  query_seg <- GDCquery(project = "TCGA-PCPG",
                        data.category = "Copy Number Variation",
                        data.type = c("Masked Copy Number Segment"),
                        legacy = F)
  
  
  GDCdownload(query_seg, method = "api", files.per.chunk = 5, directory=output_dir)
  # 
  cnv_seg <- GDCprepare(query_seg, 
                        save = T, 
                        save.filename = paste0(output_dir, "/tcga_ppgl_cnv_seg.rdata"), 
                        directory=output_dir)
}

dataloader_tcga_ppgl_cna_seg <- function(quickload=T) {
  tcga_ppgl_cna_seg_data_file <- paste0(tcga_seg_dir, "/tcga_ppgl_cnv_seg.rdata")
  if(!quickload | !file.exists(tcga_ppgl_cna_seg_data_file)){
    message("Quickload is off or pre-loaded data file not found. Attempting to download data ...")
    download_tcga_ppgl_cna_seg(output_dir = tcga_seg_dir)
  } else {
    load(tcga_ppgl_cna_seg_data_file)
    assign(x = "tcga_ppgl_cna_seg", value = data, envir = globalenv())
  }
}


dataloader_tcga_ppgl_ascat <- function(quickload=T) {
  tcga_ppgl_cna_ascat_data_file <- paste0(tcga_ascat_dir, "/tcga_ppgl_cnv_ascat.rdata")
  if(!quickload | !file.exists(tcga_ppgl_cna_ascat_data_file)){
    message("Quickload is off or pre-loaded data file not found. Attempting to download data ...")
    download_tcga_ppgl_ascat(output_dir = tcga_ascat_dir)
  } else {
    load(tcga_ppgl_cna_ascat_data_file)
    assign(x = "tcga_ppgl_cna_ascat", value = data, envir = globalenv())
  }
}