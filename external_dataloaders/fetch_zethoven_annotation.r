#Helper function to read in zethoven annotation file and optionally filter for a dataset
fetch_zethoven_annotation <- function(datasets=NULL)
{
  base_dir="/g/data/pq08/projects/ppgl"
  zethoven_annotation_file <- paste(base_dir,"public_data/annotation/tcga_comete_flynn_bulk_metadata.tsv", sep="/")
 
  if(!exists("zethoven_cluster_naming")) { source(paste(base_dir,"public_data/annotation/zethoven_sn_rna_cluster_naming.r", sep="/")) }
   
  message("Loading sample annotation from ", zethoven_annotation_file)
  zethoven_annotation <- read.delim(zethoven_annotation_file)
  
  zethoven_annotation <- zethoven_annotation %>% 
    mutate(Sample=gsub("[.]","-", Sample), 
           Sample.raw=gsub("[.]","-", Sample.raw))
  
  if(!is.null(datasets)) {
    zethoven_annotation <- zethoven_annotation %>% filter(Dataset %in% datasets) 
    
    if("TCGA" %in% datasets) {
      
      #Load extended TCGA dataset annotation
      tcga_extended <- read.delim("./public_data/annotation/fishbein_2017_table_s2_tcga_annotation.tsv") %>% 
        mutate(Sample_ID=stringr::str_sub(Sample_ID,1,15))
      
      #Merge extended for dataset merging
      zethoven_annotation <- zethoven_annotation %>% 
        left_join(tcga_extended %>%  dplyr::select(-mRNA_subtype, -Gender, -Purity), 
                   by=c("Sample"="Sample_ID"), suffix = c("_zethoven", "_TCGA")) %>% 
        dplyr::select(-PPGL, -Tumor_location)
    }  
  }
  
  zethoven_annotation <- zethoven_annotation %>%
    left_join(zethoven_cluster_naming %>% 
                dplyr::select(final, historical_b_variant) %>% 
                dplyr::rename(new_naming=final), 
              by=c("Cluster"="historical_b_variant")) %>% mutate(
                new_naming=dplyr::recode(new_naming, "C2B2 (MAML)"="C2B2 (MAML3)"))
  
  return(zethoven_annotation)
}