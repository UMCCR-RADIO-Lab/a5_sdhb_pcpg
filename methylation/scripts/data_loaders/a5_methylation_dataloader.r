suppressMessages(library(missMethyl))
suppressMessages(library(minfi))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
if(!require("IlluminaHumanMethylationEPICanno.ilm10b5.hg38"))
{
  devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICmanifest") 
  devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")
}
if (!require("maxprobes")) {devtools::install_github("markgene/maxprobes")}


data_loader_a5_methylation_array <- function(quickload=T, quickload_checkpoint_dir="/g/data/pq08/projects/ppgl/a5/methylation/quickload_checkpoints", output_qc=F, remove_excluded_samples=T, normalisation="Functional")
{
  det_p_threshold <- 0.01
  
  normalisation <- tolower(normalisation)
  if (!(normalisation %in% c("quantile", "functional"))) {
    stop("Normalisation method must be either 'quantile' or 'functional'")
  } 
  
  message("Starting A5 methylation data loader ...")
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
  
  base_dir="/g/data/pq08/projects/ppgl"
  data_dir <- paste0(base_dir,"/a5/methylation/raw_data/ILMLEPIC-16614")
  
  if(!exists("a5_anno")) {
    source(paste0(base_dir,"/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r"))
    data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.org", use_cache = T)
  }
  
  excluded_samples <- a5_anno %>% filter(Exclude=="Y") %>% pull(A5_ID)
  
  # Find the sample sheet
  message("Reading targets file...")  
  targets <- read.metharray.sheet(paste0(data_dir,"/idat_symlinks"), pattern="SampleSheet.csv") %>%
    mutate (Sample_Name=gsub("E158_T01_D_(.)","E158_T0\\1_D",Sample_Name)) %>% #Fix misnamed samples
    mutate(Sample_Name = gsub("_T0", "-", Sample_Name)) %>% 
    mutate(Sample_Name = gsub("_D", "", Sample_Name))
  
  if(!("methylation_batch" %in% colnames(a5_anno))) {
    a5_anno_temp <- a5_anno %>% 
      left_join(targets %>% 
                  dplyr::select(Sample_Name, Sample_Group) %>% 
                  dplyr::rename(methylation_batch=Sample_Group), by=c("A5_ID"="Sample_Name"))
    assign("a5_anno", a5_anno_temp, globalenv())
  }
  
  if(quickload){
    if(output_qc) { message("QC plots will not be generated while using quickload. Disable quickload to output QC") }
    if(remove_excluded_samples) { message("WARNING: The 'remove_excluded_samples' flag has no effect when quickloading. " ,
                                          "If caches were generated with 'remove_excluded_samples=F' then ",
                                          "regenerate quickload caches using 'quickload=F' and ",
                                          "'remove_excluded_samples=T'") }
    
    # message("Quickloading raw methylation values from ", paste0(quickload_checkpoint_dir,"/mset_raw.rds"))
    # m_set_raw <- readRDS(paste0(quickload_checkpoint_dir,"/mset_raw.rds"))
     
    # message("Quickloading ", normalisation, "-normalised values from ", paste0(quickload_checkpoint_dir,"/msetsq_quantile.rds"))
    # m_set_norm <- readRDS(paste0(quickload_checkpoint_dir, "/msetsq_", normalisation, ".rds"))

    quickload_file <- paste0(quickload_checkpoint_dir, "/msetsq_", normalisation, "_filter.rds")
    if(file.exists(quickload_file)) {
    message("Quickloading ", normalisation, "-normalised quality filtered values from ", quickload_file)
    m_set_filtered <- readRDS(paste0(quickload_checkpoint_dir, "/msetsq_", normalisation, "_filter.rds"))
    } else {
      stop("Quickload file not found for normalisation method ", normalisation, ". Have you run the dataloader without quickload?")
    }

    message("Loading negative control probes ...")
    neg_cont_probes <- readRDS(paste0(quickload_checkpoint_dir, "/negative_control_probes.rds"))
    
    message("Loading EPIC annotation...")
    annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
    
    targets <- targets[targets$Sample_Name %in% colnames(m_set_filtered),]
    
  } else
  {
    # Following this guide modified for EPIC arrays:
    # https://bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html
    # Arrays used are infinium-methylationepic-v-1-0-b5
    
    # Load annotation
    message("Loading EPIC annotation...")  
    annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b5.hg38)
    
    if(remove_excluded_samples) {
      targets <- targets[!(targets$Sample_Name %in% excluded_samples),]
      message("Removed excluded samples:", toString(excluded_samples))
    }
    
    #IDAT symlinks created with bash commands:
    # base_dir="/g/data/pq08/projects/ppgl/a5/methylation/raw_data/ILMLEPIC-16614"
    # while read -r idat_file;
    #   do ln -s ${idat_file} ${base_dir}/idat_symlinks/$(basename ${idat_file});
    # done < <(find ${base_dir} -iname "*.idat")
    # ln -s ${base_dir}/ILMLEPIC-16614_SampleSheet.csv ${base_dir}/idat_symlinks/ILMLEPIC-16614_SampleSheet.csv
    
    # Read in the IDAT files
    message("Reading IDAT files...")  
    rg_set <- read.metharray.exp(targets=targets)
    rg_set@annotation = c(array = "IlluminaHumanMethylationEPIC", annotation = "ilm10b5.hg38")
    
    sampleNames(rg_set) <- targets$Sample_Name
    
    # P values that a probe was detected
    message("Computing detection p-values...")  
    detP <- detectionP(rg_set)
    
    if(output_qc){
      # Plot mean detection p-values 
      message("Generating detection p-values QC plot...")  
      make_mean_detection_pvals_plot(targets, detP, paste0(base_dir,"/a5/methylation/qc/probe_expression_pvalues.pdf")) 
      
      
      # Make a QC report
      message("Generating QC report PDF...") 
      qcReport(rg_set, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
               pdf= paste0(base_dir,"/a5/methylation/qc/qc_report.pdf"))
      
    }
    
    # # create a MethylSet object from the raw data for plotting
    # message("Generating raw (non-normalised) methylation values..") 
    # m_set_raw <- preprocessRaw(rg_set)
    # message("Saving raw checkpoint data...")
    # saveRDS(m_set_raw, paste0(quickload_checkpoint_dir, "/mset_raw.rds"))
    
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
      make_norm_density_plot(targets, 
                             rg_set, 
                             m_set_norm, 
                             paste0(base_dir,"/a5/methylation/qc/", normalisation, "_normalisation.pdf"),
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
    
    neg_cont_probes <- getINCs(rg_set)
    message("Saving negative control probes ...")
    saveRDS(neg_cont_probes, paste0(quickload_checkpoint_dir, "/negative_control_probes.rds"))
    
    assign("a5_meth_rg_set", rg_set, globalenv())
    message("Loaded object a5_meth_rg_set (un-processed red/green values) into global environment")
  }
  
  
  assign("epic_array_annotation_hg38", annEPIC, globalenv())
  assign("a5_methylation_targets", targets, globalenv())
  #assign("a5_methylation_raw", m_set_raw, globalenv())
  #assign("a5_methylation_quantile", m_set_norm, globalenv())
  assign("a5_methylation_filtered", m_set_filtered, globalenv())
  assign("a5_methylation_neg_cont_probes", neg_cont_probes, globalenv())
  
  
  message("Loaded objects into global environment: 
          epic_array_annotation_hg38 (probe annotation),
          a5_methylation_targets (array manifest file),
          a5_methylation_neg_cont_probes (array negative control probes from getINCs()),
          a5_methylation_filtered (dataset after ", normalisation,"  normalisation, 
                                   removing probes that are cross-reactive, 
                                   are at a SNP site, or have a detection p_val > ",det_p_threshold,")")
  
  # a5_methylation_raw (non-normalised methylation signal), 
  #         a5_methylation_quantile (quantile-normalised methylation signal), 
  
  message("Completed A5 methylation data loader.")
}
message("Created data loader function data_loader_a5_methylation_array()")