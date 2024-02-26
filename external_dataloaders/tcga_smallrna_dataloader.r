##Dataloader function for TCGA small RNA data

library(dplyr)
source("./public_data/data_loaders/fetch_zethoven_annotation.r")

data_loader_tcga_smallrna <- function(genome_version="hg38", only_zethoven_annotated=T)
{
  message("Starting TCGA small-RNA data loader ...")
  
  base_dir="/g/data/pq08/projects/ppgl"
  data_dir=paste(base_dir, "public_data/small_rna/tcga/counts", sep="/")
  
  if(!exists("zethoven_cluster_naming")) { source(paste(base_dir,"public_data/annotation/zethoven_sn_rna_cluster_naming.r", sep="/")) }
  
  if(!(genome_version %in% c("grch37", "hg38")))
  {
    stop("genome_version must be 'grch37' or 'hg38'")
  }
  
  # Read in the featureCount summary
  counts_summary <- read.delim(paste0(data_dir,"/tcga_subread_counts_",genome_version,"_mirbase_ref.tsv.summary"), header = T,sep = "\t", check.names = F)
  colnames(counts_summary) <- gsub("TCGA_bams[/]([A-Z0-9-]+).subread_results.bam", "\\1", colnames(counts_summary))
  
  # Read in the featureCount reads
  smallrna_read_counts <-  read.delim(paste0(data_dir,"/tcga_subread_counts_",genome_version,"_mirbase_ref.tsv"), header = T,sep = "\t", comment.char = '#', check.names = F)
  #Keep only patient barcode portion of col headers 
  colnames(smallrna_read_counts) <- gsub("TCGA_bams[/]([A-Z0-9-]{15}).+subread_results.bam", "\\1", colnames(smallrna_read_counts))
  
  smallrna_feature_annotation <- smallrna_read_counts[,1:6]
  smallrna_read_counts <- smallrna_read_counts[,-c(2:6)]
  
  #Load zethoven annotation
  zethoven_annotation_tcga <- fetch_zethoven_annotation(datasets="TCGA")
  
  # #Load compendium 
  # zethoven_annotation_file <- paste(base_dir,"annotation/tcga_comete_flynn_bulk_metadata.tsv", sep="/")
  # message("Loading sample annotation from ", zethoven_annotation_file)
  # zethoven_annotation <- read.delim(zethoven_annotation_file, sep="\t")
  # zethoven_annotation_tcga <- zethoven_annotation %>% 
  #   filter(Dataset=="TCGA") %>% 
  #   mutate(Sample=substr(gsub("[.]","-", Sample.raw),1,15)) %>%
  #   left_join(zethoven_cluster_naming %>% 
  #               dplyr::select(final, historical_b_variant) %>% 
  #               dplyr::rename(new_naming=final), 
  #             by=c("Cluster"="historical_b_variant"))
  
  
  #Filter annotation to only samples with small RNA data
  zethoven_annotation_tcga <- zethoven_annotation_tcga %>% filter(Sample %in% colnames(smallrna_read_counts))
  
  if(only_zethoven_annotated) {
    #remove samples without annotation
    smallrna_read_counts <- smallrna_read_counts %>% dplyr::select(c("Geneid", !!zethoven_annotation_tcga$Sample))
  } else {
    #Add empty annotation rows for un-annotated samples 
    zethoven_annotation_tcga <- bind_rows(zethoven_annotation_tcga, 
                                            data.frame(Sample=setdiff(colnames(smallrna_read_counts)[-1], 
                                                                      zethoven_annotation_tcga$Sample)))
  }
  
  assign("tcga_smallrna_read_counts_summary", counts_summary, globalenv())
  assign("tcga_smallrna_read_counts", smallrna_read_counts, globalenv())
  assign("tcga_smallrna_anno", zethoven_annotation_tcga, globalenv())
  if(!exists("smallrna_feature_annotation")) {
    assign("smallrna_feature_annotation", smallrna_feature_annotation, globalenv())}
  
  message("Loaded objects into global environment: 
          tcga_smallrna_read_counts_summary (featureCounts summary table),
          tcga_smallrna_read_counts (raw counts table),
          tcga_smallrna_anno (sample annotation),
          smallrna_feature_anno (feature annotation)")
  
  message("Completed TCGA small-RNA data loader.")
}
message("Created data loader function data_loader_tcga_smallrna()")