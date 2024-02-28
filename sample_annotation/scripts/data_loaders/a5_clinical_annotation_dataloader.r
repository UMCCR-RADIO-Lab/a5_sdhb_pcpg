library(googlesheets4)
library(dplyr)

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_add_differential_groups.r")


#Function to read the clinical data from the google sheet or offline-cache
data_loader_a5_clinical_anno <- function(google_account=NULL, 
                                         use_cache=FALSE, 
                                         offline_cache="/g/data/pq08/projects/ppgl/a5/offline_cache/a5_clinical_annotation.tsv",
                                         remove_excluded_samples=T)
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
  
  if(remove_excluded_samples == T) {
    a5_anno <- a5_anno %>% filter(Exclude == "N")
  }
  
  a5_anno <- add_differential_groups(a5_anno)
  
  a5_anno <- a5_anno %>% mutate(age_at_resection=
                       floor(lubridate::interval(lubridate::my(paste("01", `Year of birth`,sep="-")),
                                                 lubridate::dmy(`Date of resection (DD/MM/YYYY)`)) 
                             / lubridate::years(1)))
  
  a5_anno <- a5_anno %>% 
    mutate(primary_location_plotting = case_match(Primary_Location_Simplified,
      "Head_neck" ~ "Head and neck",
      "Extraadrenal_abdominal" ~ "Extraadrenal (abdominal/thoracic)",
      "Extraadrenal_thoracic" ~ "Extraadrenal (abdominal/thoracic)",
      "Extraadrenal_bladder" ~ "Extraadrenal (bladder)",
      "Adrenal" ~ "Adrenal",
      "Adrenal_right" ~ "Adrenal",
      "Adrenal_left" ~ "Adrenal" 
    ),
    primary_location_plotting = factor(primary_location_plotting, 
                                       levels = c("Adrenal", "Extraadrenal (abdominal/thoracic)", 
                                                  "Extraadrenal (bladder)", "Head and neck")),
    cell_of_origin = case_when(
      Primary_Location_Simplified == "Normal" ~ "Chromaffin",
      Primary_Location_Simplified == "Head_neck" & grepl("Non_chromaffin", Major_Cluster) ~ "Non_chromaffin",
      Primary_Location_Simplified == "Extraadrenal_thoracic" & grepl("Non_chromaffin", Major_Cluster) ~ "Non_chromaffin",
      Primary_Location_Simplified != "Head_neck" & grepl("Chromaffin", Major_Cluster) ~ "Chromaffin",
      TRUE ~ "Ambiguous"
    ),
    cell_of_origin = factor(cell_of_origin, levels = c("Chromaffin","Non_chromaffin","Ambiguous")))
  
  a5_anno <- a5_anno %>% 
    mutate(structural_variant_count = as.numeric(structural_variant_count),
           wgs_tmb = as.numeric(wgs_tmb),
           Largest_primary_dimensions_cm = as.numeric(Largest_primary_dimensions_cm),
           Age_first_dx = as.numeric(Age_first_dx),
           telhunter_log2_telcontentratio = as.numeric(telhunter_log2_telcontentratio),
           telhunter_RNA_telcontent = as.numeric(telhunter_RNA_telcontent),
           TERT_purity_adj_vaf = as.numeric(TERT_purity_adj_vaf),
           ATRX_purity_adj_vaf = as.numeric(ATRX_purity_adj_vaf),
           sample_purity = as.numeric(sample_purity),
           MKI67_log2_cpm = as.numeric(MKI67_log2_cpm)
           )
 
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
      Cluster = case_match(A5_ID, 
                           "E185-1" ~ "A5 - Head_neck",
                           "E128-1" ~ "A5 - Abdominal_Thoracic",
                           .default =paste(
        "A5",
        differential_group_anatomy,
        sep = " - "
      )),
      new_naming = dplyr::recode(Cluster,
                                         "A5 - Head_neck"="C1A2 (SDHx-HN)",
                                         "A5 - Abdominal_Thoracic"="C1A1 (SDHx)",
                                         "A5 - Thoracic_non_chromaffin"="C1A2 (SDHx-HN)",
                                         "A5 - Extraadrenal"="C1A1 (SDHx)",
                                         "A5 - NF1"="C2A (Kinase)",
                                         "A5 - Adrenal"="C1A1 (SDHx)",
                                         "A5 - VHL"="C1B1 (VHL)")
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
  
  if("E124-1" %in% a5_anno_tcga_format$Sample) {
  a5_anno_tcga_format[a5_anno_tcga_format$Sample == "E124-1", c("Genotype", "Cluster", "TCGA_Cluster", "new_naming")] <-
    c("NF1", "Kinase", "Kinase signaling", "A5 - NF1") }
  
  if("E145-1" %in% a5_anno_tcga_format$Sample) {
  a5_anno_tcga_format[a5_anno_tcga_format$Sample == "E145-1", c("Genotype", "Cluster", "TCGA_Cluster", "new_naming")] <-
    c("VHL", "VHL", "Pseudohypoxia", "A5 - VHL") }
  
  return(a5_anno_tcga_format)
  
}
message("Created helper function a5_to_zethoven_anno_format()")