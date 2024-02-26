library(limma)
library(edgeR)
library(dplyr)
library(purrr)

data_loader_a5_smallrna <- function(genome_version="hg38", exclude_samples=c(), minimum_RIN=6)
{
  
  #################
  # Preliminaries #
  #################
  
  message("Starting A5 small-RNA data loader ...")
  
  if(!(genome_version %in% c("grch37", "hg38")))
  {
    stop("genome_version must be 'grch37' or 'hg38'")
  }
  
  base_dir="/g/data/pq08/projects/ppgl/"
  data_dir=file.path(base_dir, "a5/small_rna/analysis/feature_counts")
 
  #####################
  # Sample annotation #
  #####################
  
  if(!exists("a5_anno")) {
    source(paste0(base_dir,"/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r"))
    data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.org", use_cache = T)
  }
   
  if(!exists("a5_qc")) {
    source(paste0(base_dir,"/a5/sample_annotation/scripts/data_loaders/a5_qc_annotation_dataloader.r"))
    data_loader_a5_qc_anno(google_account = "aidan.flynn@umccr-radio-lab.org", use_cache = T)
  }
  
  ################
  # Data Loading #
  ################
  
  # Read in the featureCount summary
  counts_summary <- read.delim(paste0(data_dir,"/a5_subread_counts_",genome_version,"_mirbase_ref.tsv.summary"), header = T,sep = "\t")
  colnames(counts_summary) <- gsub("bams[.](E[0-9]{3}).T0([0-9]).subread_results.bam", "\\1-\\2", colnames(counts_summary))
  
  # Read in the featureCount reads
  smallrna_read_counts <-  read.delim(paste0(data_dir,"/a5_subread_counts_",genome_version,"_mirbase_ref.tsv"), header = T,sep = "\t", comment.char = '#')
  colnames(smallrna_read_counts) <- gsub("bams[.](E[0-9]{3}).T0([0-9]).subread_results.bam", "\\1-\\2", colnames(smallrna_read_counts))
  
  smallrna_feature_annotation <- smallrna_read_counts[,1:6]
  smallrna_read_counts <- smallrna_read_counts[,-c(2:6)]

  smallrna_read_counts_mat <- as.matrix(smallrna_read_counts[,-1])
  rownames(smallrna_read_counts_mat) <- smallrna_read_counts$Geneid

  ###############
  # Sample Sets #
  ###############

  #######
  # Excluded samples
  #######
  
  #Samples to exclude from core set for quality/purity/genotype reasons
  exclude_global <- a5_anno %>%  filter(Exclude=="Y") %>% pull(A5_ID)
  exclude_missing_anno <- setdiff(colnames(smallrna_read_counts_mat), a5_anno$A5_ID)
  exclude_genotype <- a5_anno %>%  filter(Genotype != "SDHB") %>% pull(A5_ID)
  exclude_no_rna_data <- setdiff(a5_anno$A5_ID, colnames(smallrna_read_counts_mat))
  exclude_no_wgs_data <- c("E181-1", "E191-1")
  exclude_qc <- a5_qc %>% filter(smallrna__RIN < minimum_RIN) %>% pull(A5_ID)
  
  exclude_base <- c(exclude_no_rna_data, exclude_no_wgs_data, exclude_qc, exclude_global, exclude_missing_anno)
  
  a5_anno <- a5_anno %>% mutate(Exclude_smallRNA=ifelse(A5_ID %in% exclude_base, T, F))
  
  #######
  # Anatomy groups
  #######
  
  #Select head/neck+aortic tumours
  samples_hn <- a5_anno %>% 
    filter(differential_group_anatomy == "Head_Neck") %>% 
    pull(A5_ID)
  
  #Select abdo_thoracic tumours
  samples_abdothoracic <- a5_anno %>% 
    filter(differential_group_anatomy == "Abdominal_Thoracic") %>% 
    pull(A5_ID)
  
  #Select samples of unclear origin tumours based on UMAP
  samples_mediastinum <- a5_anno %>% 
    filter(differential_group_anatomy == "Mediastinum") %>% 
    pull(A5_ID)
  
  #Select samples of unclear origin tumours based on UMAP
  samples_ambiguous <- a5_anno %>% 
    filter(differential_group_anatomy == "Ambiguous") %>% 
    pull(A5_ID)
  
  #######
  # Create sets
  #######
  
  counts_df_list <- list()
  
  #All data
  counts_df_list[["all"]] <- smallrna_read_counts_mat
  
  #Poor QC and samples with no WGS excluded
  counts_df_list[["qc_ok"]] <- smallrna_read_counts_mat[,!(colnames(smallrna_read_counts_mat) %in% exclude_base)]
  
  #Poor QC, samples with no WGS, and NF/VHL samples excluded
  counts_df_list[["SDHB"]] <- smallrna_read_counts_mat[,!(colnames(smallrna_read_counts_mat) %in% c(exclude_base, exclude_genotype))]
  
  #Head and neck based on clinical and UMAP annotation (aortic cases excluded)
  samples_hn_keep <- setdiff(samples_hn, c(exclude_base, exclude_genotype))
  counts_df_list[["SDHB_HN"]] <- smallrna_read_counts_mat[,colnames(smallrna_read_counts_mat) %in% samples_hn_keep]
  
  #Abdominal/Thoracic based on clinical and UMAP annotation (ambiguous cases excluded)
  samples_abdothoracic_keep <- setdiff(samples_abdothoracic, c(exclude_base, exclude_genotype))
  counts_df_list[["SDHB_abdothoracic"]] <- smallrna_read_counts_mat[,colnames(smallrna_read_counts_mat) %in% samples_abdothoracic_keep]
  message("Created data subsets:", toString(names(counts_df_list)))
  
  
  ##############
  # Preprocess #
  ##############
  
  #######
  # Create DGE objects
  #######
  
  message("Converting counts to DGE objects ...")
  counts_dge_list <- lapply(counts_df_list, DGEList)
  
  # filter genes with low expression in all/nearly all samples
  message("Filtering lowly expressed genes ...")
  counts_dge_list <- map(.x = counts_dge_list, 
                         .f = function (dge_object)  { 
                           diff_groups <- a5_anno$differential_group[match(colnames(dge_object), a5_anno$A5_ID)]
                           expr_genes <- filterByExpr(dge_object, group = diff_groups)  
                           dge_object[expr_genes,, keep.lib.sizes = FALSE] })

  #######
  # TMM normalisation
  #######
  
  message("Calculating normalisation factors with edgeR ...")
  counts_dge_list <- map(counts_dge_list, calcNormFactors)
  
  
  ############################
  # Export objects to global #
  ############################
  
  assign("a5_smallrna_read_counts_summary", counts_summary, globalenv())
  assign("a5_smallrna_read_counts", smallrna_read_counts, globalenv())
  assign("a5_smallrna_dge_list", counts_dge_list, globalenv())
  if(!exists("smallrna_feature_annotation",envir = globalenv())) {
    assign("smallrna_feature_annotation", smallrna_feature_annotation, globalenv())}
  assign("a5_anno", a5_anno, globalenv())
  
  message("Loaded objects into global environment: 
          a5_smallrna_read_counts_summary (featureCounts summary table),
          a5_smallrna_read_counts (raw counts table),
          smallrna_feature_anno (feature annotation),
          a5_smallrna_dge_list  (a list of count normalised DGE objects named:
                          - all: All profiled samples
                          - qc_ok: Excluded samples removed (",toString(exclude_base),")
                          - SDHB: QC_OK with non-SDHB samples removed (",toString(exclude_genotype),")
                          - SDHB_abdothoracic: Abdominal/thoracic based on clinical and UMAP annotation (mediastinum and ambiguous cases excluded:",toString(samples_ambiguous),")
                          - SDHB_HN: Head and neck based on clinical and UMAP annotation (mediastinum and ambiguous cases excluded)"
          )
  
  message("Completed A5 small-RNA data loader.")
}
message("Created data loader function data_loader_a5_smallrna()")