######################################################
# Script for (down)loading TCGA methylation data set #
# Author: Aidan Flynn                                #
# Date: 15/01/2023                                   #
# Languages: R                                       #
######################################################


suppressMessages(library(minfi))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(tibble))
if (!require("maxprobes")) {devtools::install_github("markgene/maxprobes")}

source("/g/data/pq08/projects/ppgl/public_data/data_loaders/fetch_zethoven_annotation.r")

###################
# Data downloader #
###################

download_tcga_ppgl_idat <- function(output_dir) {
  
  require(TCGAbiolinks)
  
  if(dir.exists(output_dir)){
    dir.create(output_dir, recursive = T)}
  
  message("Querying GDC for file information ...")
  query_idat <- GDCquery(project = "TCGA-PCPG",
                         data.category = "DNA Methylation", data.format = "IDAT",
                         legacy = F)
  
  message("Beginning download ...")
  GDCdownload(query_idat, method = "api", files.per.chunk = 5, directory=output_dir)
  
  message("Writing sample manifest ...")
  write.table(x = getResults(query = query_idat), 
              file = paste(output_dir, "sample_manifest.tsv", sep="/"), 
              sep="\t", row.names = F)
  
  # Deconstruct nested download folder into flat structure and rename files to TCGA IDs
  message("Flattening directory structure and renaming files ...")
  idat_info <- getResults(query_idat) %>%  
    rowwise() %>%  
    mutate(new_filename = gsub(pattern = ".+_noid(_(Grn|Red).idat)", 
                               replacement = paste0(cases,"\\1"), 
                               file_name))
  
  source_files <- list.files(output_dir, recursive = T, pattern = "_(Grn|Red).idat", full.names = T)
  names(source_files) <- basename(source_files)
  
  destination_files <- paste(output_dir, idat_info$new_filename, sep="/")
  names(destination_files) <- idat_info$file_name
  destination_files <- destination_files[names(source_files)]
  
  #Action file move
  purrr::walk2(.x = source_files, .y = destination_files, ~file.symlink(.x,.y))
  
  message("Download complete. IDAT files and sample manifest available in ", output_dir)
}

###############
# Data loader #
###############

data_loader_tcga_methylation <- function(data_dir="/g/data/pq08/projects/ppgl/public_data/methylation/tcga", 
                                         quickload=T, 
                                         quickload_checkpoint_dir="/g/data/pq08/projects/ppgl/public_data/methylation/tcga/quickload_checkpoints", 
                                         output_qc=F, 
                                         normalisation="Functional")
{
  det_p_threshold <- 0.01
  
  normalisation <- tolower(normalisation)
  if (!(normalisation %in% c("quantile", "functional"))) {
    stop("Normalisation method must be either 'quantile' or 'functional'")
  } 
  
  message("Starting TCGA methylation data loader ...")
  
  #Start sub-functions
  #Helper sub-function to plot mean detection p-values 
  make_mean_detection_pvals_plot <- function(targets, detP, pdf_out){
    pal <- rainbow(12)
    p_val_df <- data.frame(ID= targets$ID, p_means = colMeans(detP),group =targets$Sample_Group, col=pal[factor(targets$Sample_Group)])
    
    g1 <- ggplot(data = p_val_df, aes(x = ID, y = p_means, fill = group))+
      geom_bar(stat= "identity")+
      scale_colour_manual(values = col)+
      theme_bw(base_size = 18)+
      theme(panel.grid=element_blank(),
            strip.background = element_blank()) +
      labs(x = 'Sample', y = "Mean probe detection\n p value", fill = "Group")+
      theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size = 4))
    ggsave(plot = g1, filename = pdf_out, width = 10, height = 7, device = pdf())
  }
  
  #Helper sub-function to plot pre- and post-normalised data 
  make_norm_density_plot <- function(targets, rg_set, m_set_norm, pdf_out, normalisation){
    pdf(pdf_out, width =10)
    par(mfrow=c(1,2))
    densityPlot(rg_set, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE,
                cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    densityPlot(getBeta(m_set_norm), sampGroups=targets$Sample_Group,
                main=paste(normalisation, "Normalised"), legend=FALSE,
                cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
    dev.off()
  }
  #end sub-functions
  
  #Load zethoven dataset annotation
  zethoven_annotation_tcga <- fetch_zethoven_annotation(datasets="TCGA")
  
  #Add batch for dataset merging
  zethoven_annotation_tcga <- zethoven_annotation_tcga %>% mutate(methylation_batch="TCGA_450k")
  
  if(quickload){
    if(output_qc) { message("QC plots will not be generated while using quickload. Disable quickload to output QC") }
    
    quickload_file=paste0(quickload_checkpoint_dir, "/msetsq_", normalisation, "_filter.rds")
    if(file.exists(quickload_file)) {
      message("Quickloading ", normalisation, "-normalised quality filtered values from ", quickload_file)
      m_set_filtered <- readRDS(quickload_file)
    } else {
      stop("Quickload file not found for normalisation method ", normalisation, ". Have you run the dataloader without quickload?")
    }
    
    zethoven_annotation_tcga <- zethoven_annotation_tcga[zethoven_annotation_tcga$Sample %in% colnames(m_set_filtered),]
    zethoven_annotation_tcga <- zethoven_annotation_tcga[match(colnames(m_set_filtered), zethoven_annotation_tcga$Sample),]
    
  } else
  {

    # Read in the IDAT files
    message("Reading IDAT files...")  
    rg_set <- read.metharray.exp(base = paste(data_dir, "idat", sep="/"))
    
    #Shorten TCGA sample codes to match Zeethovan annotation
    sampleNames(rg_set) <- stringr::str_sub(sampleNames(rg_set),1,15)
    
    #remove metylation arrays not in annotation set
    unannotated_samples <- which(sampleNames(rg_set) %in% setdiff(sampleNames(rg_set), zethoven_annotation_tcga$Sample))
    message("Removing samples without annotation: ", toString(sampleNames(rg_set)[unannotated_samples]))
    rg_set <- rg_set[,-unannotated_samples]
    
    #restrict and reorder annotation
    zethoven_annotation_tcga <-  zethoven_annotation_tcga %>%  filter(Sample %in% sampleNames(rg_set)) 
    zethoven_annotation_tcga <- zethoven_annotation_tcga[match(sampleNames(rg_set), zethoven_annotation_tcga$Sample),]
    
    # P values that a probe was detected
    message("Computing detection p-values...")  
    detP <- detectionP(rg_set)
    
    qc_annotation <- zethoven_annotation_tcga %>% select(Sample, new_naming)
    colnames(qc_annotation) <- c("ID","Sample_Group")
    
    if(output_qc){
      # Plot mean detection p-values 
      message("Generating detection p-values QC plot...")  
      make_mean_detection_pvals_plot(targets=qc_annotation, detP, paste0(data_dir,"/qc/probe_expression_pvalues.pdf")) 
      
      
      # Make a QC report
      message("Generating QC report PDF...") 
      qcReport(rg_set, sampNames=qc_annotation$ID, sampGroups=qc_annotation$Sample_Group, 
               pdf= paste0(data_dir,"/qc/qc_report.pdf"))
      
    }
    
    #If during normalisation you get the following error:
    #Error in normalize.quantiles(mat[Index2, ]) : 
    #  ERROR; return code from pthread_create() is 22
    #then run: 
    # BiocManager::install("preprocessCore", configure.args="--disable-threading")
    # as per https://support.bioconductor.org/p/122925/
    if (normalisation == "quantile") {
      message("Performing quantile normalisation preprocessing...") 
      m_set_norm <- preprocessQuantile(rg_set) 
    } else if (normalisation == "functional") {
      message("Performing functional normalisation preprocessing...") 
      m_set_norm <- preprocessFunnorm(rg_set) 
    }
    
    message("Saving ", normalisation, " normalisation checkpoint data...")
    saveRDS(m_set_norm, paste0(quickload_checkpoint_dir,  "/msetsq_", normalisation, ".rds"))
    
    if(output_qc)
    {
      message("Generating normalisation QC plot...")  
      # Look at the data before and after normalisation
      make_norm_density_plot(qc_annotation, 
                             rg_set, 
                             m_set_norm, 
                             paste0(data_dir,"/qc/", normalisation, "_normalisation.pdf"),
                             normalisation)
    }
    
    # ensure probes are in the same order in the m_set_norm and detP objects
    detP <- detP[match(featureNames(m_set_norm),rownames(detP)),]
    
    # remove any probes that have failed in one or more samples
    keep <- rowSums(detP < det_p_threshold) == ncol(m_set_norm) 
    stat_summary <- table(keep)
    
    # Data to keep
    m_set_filtered <- m_set_norm[keep,]
    message("Removed failed probes with p_val > ",det_p_threshold,": n=",stat_summary["FALSE"],". Kept: n=",stat_summary["TRUE"])
    
    # Remove probes with SNPs at CpG site
    pre <- nrow(m_set_filtered)
    m_set_filtered <- dropLociWithSnps(m_set_filtered)
    post <- nrow(m_set_filtered)
    message("Removed ", pre-post ," loci with SNPs")
    
    # Get cross reactive probes from the maxprobes library
    xReactiveProbes <- unlist(maxprobes::xreactive_probes(array_type = "EPIC"))
    
    # Exclude cross reactive probes 
    keep <- !(featureNames(m_set_filtered) %in% xReactiveProbes)
    stat_summary <- table(keep)
    
    m_set_filtered <- m_set_filtered[keep,] 
    message("Removed cross-reactive probes: n=",stat_summary["FALSE"],". Kept: n=",stat_summary["TRUE"])
    
    message("Saving filtered, normalised data checkpoint ...")
    saveRDS(m_set_filtered, paste0(quickload_checkpoint_dir, "/msetsq_", normalisation, "_filter.rds"))
    
    assign("tcga_meth_rg_set", rg_set, globalenv())
    message("Loaded object tcga_meth_rg_set (un-processed red/green values) into global environment")
  }
  
  assign("Ill450k_array_annotation_hg19", getAnnotation(paste0(m_set_filtered@annotation[["array"]],"anno",".",m_set_filtered@annotation[["annotation"]])), globalenv())
  assign("tcga_methylation_450k_anno", zethoven_annotation_tcga, globalenv())
  assign("tcga_methylation_450k_filtered", m_set_filtered, globalenv())
  
  
  message("Loaded objects into global environment: 
          ill450k_array_annotation_hg19 (probe annotation),
          tcga_methylation_450k_anno (sample annotation),
          tcga_methylation_450k_filtered (dataset after ", normalisation,"  normalisation, 
                                   removing probes that are cross-reactive, 
                                   are at a SNP site, or have a detection p_val > ",det_p_threshold,")")
  
  message("Completed TCGA methylation data loader.")
}
message("Created data loader function data_loader_tcga_methylation()")



