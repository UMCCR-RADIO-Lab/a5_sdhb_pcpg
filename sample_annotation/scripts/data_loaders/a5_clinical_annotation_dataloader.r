library(googlesheets4)
library(dplyr)

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_add_differential_groups.r")


#Function to read the clinical data from the google sheet or offline-cache
data_loader_a5_clinical_anno <- function(google_account=NULL, use_cache=FALSE, offline_cache="/g/data/pq08/projects/ppgl/a5/offline_cache/a5_clinical_annotation.tsv")
{
  if(exists("a5_anno")) {
    message("WARNING: Object a5_anno already exists in the global environment and will be overwritten")
  }
  
  if(use_cache)
  {
    message("Loading clinical data from offline cache...")
    a5_anno <- read.delim(offline_cache, sep="\t", header=T, check.names = F)
  } else {
    if(is.null(google_account)) { stop("You must provide a google credential with access to the clinical annotation sheet")}
    gs4_auth(google_account)
    message("Fetching annotation from google sheet ...")
    a5_anno <- read_sheet("1hnXdXI29KvvuLxsaTBID1-EbE7mTSk6bG05HcSFfhgo", sheet="a5_sample_annotation", col_types = "c")
    message("Updating offline cache ...")
    write.table(a5_anno, offline_cache, sep="\t", row.names = F)
  }
  
  a5_anno <- a5_anno %>% dplyr::rename(`A5_ID`=`A5 ID`)
  a5_anno <- a5_anno %>% dplyr::filter(!is.na(`A5_ID`))
  a5_anno <- a5_anno %>% mutate(is_primary_or_met = replace(is_primary_or_met, is_primary_or_met == "Metastatic", "Metastasis"))
  
  a5_anno <- add_differential_groups(a5_anno)
 
  assign(x = "a5_anno", value = a5_anno, envir = globalenv())
  message("Loaded clinical annotation into the global environment as a5_anno")
}
message("Created data loader function data_loader_a5_clinical_anno()")

a5_to_zethoven_anno_format <- function(a5_anno) {
  a5_anno_tcga_format <- a5_anno %>% filter(Exclude == "N") %>%
    mutate(
      Dataset = "A5",
      TCGA_Cluster = "Pseudohypoxia",
      Genotype = "SDHB",
      Sex = stringr::str_to_sentence(Gender)
    ) %>%
    mutate(
      Malignancy = case_when(
        tumour_metastasised == "Yes" ~ "Malignant",
        tumour_metastasised == "No" ~ "Benign",
        TRUE ~ NA_character_
      ),
      Cluster = paste(
        "A5",
        gsub("_abdominal|_thoracic|_bladder|_[Ll]eft|_[Rr]ight", "", Primary_Location_Simplified),
        sep = " - "
      ),
      new_naming = Cluster
    ) %>%
    dplyr::rename(Sample = A5_ID) %>%
    dplyr::select(
      Sample,
      Cluster,
      Dataset,
      Sex,
      TCGA_Cluster,
      Genotype,
      Malignancy,
      new_naming
    )
  
  a5_anno_tcga_format[a5_anno_tcga_format$Sample == "E124-1", c("Genotype", "Cluster", "TCGA_Cluster", "new_naming")] <-
    c("NF1", "Kinase", "Kinase signaling", "A5 - NF1")
  a5_anno_tcga_format[a5_anno_tcga_format$Sample == "E145-1", c("Genotype", "Cluster", "TCGA_Cluster", "new_naming")] <-
    c("VHL", "VHL", "Pseudohypoxia", "A5 - VHL")
  
  return(a5_anno_tcga_format)
  
}
message("Created helper function a5_to_zethoven_anno_format()")