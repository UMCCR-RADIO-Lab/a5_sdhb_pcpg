##Dataloader function for COMETE small RNA data

library(dplyr)
source("./public_data/data_loaders/fetch_zethoven_annotation.r")

data_loader_comete_smallrna <- function(genome_version="hg38", only_zethoven_annotated=T, remove_outlier_samples=F)
{
  message("Starting COMETE small-RNA data loader ...")
  
  base_dir="/g/data/pq08/projects/ppgl"
  data_dir=paste(base_dir, "public_data/small_rna/comete/counts", sep="/")
  
  if(!(genome_version %in% c("grch37", "hg38")))
  {
    stop("genome_version must be 'grch37' or 'hg38'")
  }
  
  # Read in the featureCount summary
  counts_summary <- read.delim(paste0(data_dir,"/comete_subread_counts_",genome_version,"_mirbase_ref.tsv.summary"), header = T,sep = "\t", check.names = F)
  colnames(counts_summary) <- gsub("COMETE_bams[/]([A-Z0-9-]+).subread_results.bam", "\\1", colnames(counts_summary))
  
  # Read in the featureCount reads
  smallrna_read_counts <-  read.delim(paste0(data_dir,"/comete_subread_counts_",genome_version,"_mirbase_ref.tsv"), header = T,sep = "\t", comment.char = '#', check.names = F)
  colnames(smallrna_read_counts) <- gsub("COMETE_bams[/]([A-Z0-9-]+).subread_results.bam", "\\1", colnames(smallrna_read_counts))
  
  smallrna_feature_annotation <- smallrna_read_counts[,1:6]
  smallrna_read_counts <- smallrna_read_counts[,-c(2:6)]
  
  # Get the comete sample download info
  comete_annotation_file <- paste(base_dir, "public_data/small_rna/comete/annotation/comete_small_rna_refernce_e-mtab-2833.sdrf.txt", sep="/")
  comete_anno <- read.delim(comete_annotation_file, 
                            sep="\t", 
                            check.names = F) %>%
    dplyr::select(which(!duplicated(colnames(.)))) %>% 
    dplyr::rename(sample_code=`Comment[ENA_RUN]`) %>% 
    mutate(Sample = gsub("_fastq.txt.gz", "", `Comment[SUBMITTED_FILE_NAME]`))
  
  #Reorder annotation to match data
  comete_anno <- comete_anno[match(colnames(smallrna_read_counts)[-1], comete_anno$sample_code),] 
  
  #Rename with sample names
  colnames(smallrna_read_counts)[-1] <- comete_anno$Sample
  
  #Load zethoven annotation
  zethoven_annotation_comete <- fetch_zethoven_annotation(datasets="E-MTAB-733") %>% mutate(Sample=Alias)
   
  # message("Loading sample annotation from ", zethoven_annotation_file)
  # zethoven_annotation <- read.delim(zethoven_annotation_file, sep="\t")
  # 
  # zethoven_annotation_comete <- zethoven_annotation %>% 
  #   filter(Dataset=="E-MTAB-733") %>% 
  #   mutate(Sample=Alias) %>%
  #   left_join(zethoven_cluster_naming %>% dplyr::select(final, historical_b_variant) %>% dplyr::rename(new_naming=final), by=c("Cluster"="historical_b_variant"))
  # 
  
  #Filter annotation to only samples with small RNA data
  zethoven_annotation_comete <- zethoven_annotation_comete %>% filter(Sample %in% comete_anno$Sample)
  
  if(only_zethoven_annotated) {
    #remove samples without annotation
    smallrna_read_counts <- smallrna_read_counts %>% dplyr::select(c("Geneid", !!zethoven_annotation_comete$Sample))
  } else {
    #Add empty annotation rows for un-annotated samples 
    zethoven_annotation_comete <- bind_rows(zethoven_annotation_comete, 
                                            data.frame(Sample=setdiff(comete_anno$Sample, 
                                                                      zethoven_annotation_comete$Sample)))
  }
  
  
  # Remove samples that were specific to comete in the future UMAP and seem
  # to be some sort of technical effect 
  if(remove_outlier_samples) {
    comete_samples_to_remove <- readRDS("/g/data/pq08/projects/ppgl/a5/small_rna/qc/sample_blacklist/comete_batch_effect_samples.rds")
    zethoven_annotation_comete <- zethoven_annotation_comete %>% filter(!(Sample %in% comete_samples_to_remove))
    smallrna_read_counts <- smallrna_read_counts %>% dplyr::select(c("Geneid", !!zethoven_annotation_comete$Sample))
    message("Removed previously determined outlier samples:", toString(comete_samples_to_remove))
  }
  
  assign("comete_smallrna_read_counts_summary", counts_summary, globalenv())
  assign("comete_smallrna_read_counts", smallrna_read_counts, globalenv())
  assign("comete_smallrna_anno", zethoven_annotation_comete, globalenv())
  if(!exists("smallrna_feature_annotation")) {
  assign("smallrna_feature_annotation", smallrna_feature_annotation, globalenv())}
  
  message("Loaded objects into global environment: 
          comete_smallrna_read_counts_summary (featureCounts summary table),
          comete_smallrna_read_counts (raw counts table),
          comete_smallrna_anno (sample annotation),
          smallrna_feature_anno (feature annotation)")
  
  message("Completed COMETE small-RNA data loader.")
}
message("Created data loader function data_loader_comete_smallrna()")