#Helper script to append comparison groups for differential methylation/expression analysis
#to sample annotation dataframe

if(!exists("differential_group_levels"))
{
  differential_group_levels <- c(
  "NonMetastaticPrimary_WT_VAFOK_PurityLow",
  "NonMetastaticPrimary_WT_VAFOK_PurityOK"  ,
  "MetastaticPrimary_WT_VAFOK_PurityLow",
  "MetastaticPrimary_WT_VAFOK_PurityOK",
  "Metastasis_WT_VAFOK_PurityOK",
  "MetastaticPrimary_TERT_VAFLow_PurityOK" ,
  "MetastaticPrimary_TERT_VAFOK_PurityOK",
  "Metastasis_TERT_VAFOK_PurityOK",
  "MetastaticPrimary_ATRX_VAFOK_PurityOK"  ,
  "Metastasis_ATRX_VAFOK_PurityOK",
  "Other_WT_VAFOK_PurityLow",
  "Other_WT_VAFOK_PurityOK",
  "Other_TERT_VAFOK_PurityOK",
  "Other_TERT_VAFLow_PurityOK",
  "Other_TERT_VAFOK_PurityLow",
  "Other_ATRX_VAFOK_PurityOK",
  "Normal_WT_VAFOK_PurityOK")
  
  lockBinding("differential_group_levels", globalenv())
}

message("Created global variable 'differential_group_levels'")

add_differential_groups <- function(sample_annotation) {
  
  ###########
  # Anatomy #
  ###########
  
  sample_annotation <- sample_annotation %>%  mutate(differential_group_anatomy = case_when(
    Primary_Location_Simplified=="Head_neck" & grepl("Parasympathetic", Major_Cluster) ~ "Head_Neck",
    Primary_Location_Simplified=="Extraadrenal_thoracic_cardiac" ~ "Cardiac",
    Primary_Location_Simplified!="Head_neck" & grepl("Sympathetic", Major_Cluster) ~ "Abdominal_Thoracic",
    TRUE ~ "Ambiguous"
  ))
  
  #################
  # Sample Type  #
  #################
  
  sample_annotation <- sample_annotation %>%
    mutate(
      differential_group_sampletype = case_when(
        is_primary_or_met == "Primary" &
          tumour_metastasised == "Yes" ~ "MetastaticPrimary",
        is_primary_or_met == "Primary" &
          tumour_metastasised == "No" ~ "NonMetastaticPrimary",
        is_primary_or_met == "Metastasis" &
          tumour_metastasised == "Yes" ~ "Metastasis",
        TRUE ~ "Other")
    )
    

  ############
  # Genotype #
  ############
  
  
  sample_annotation <- sample_annotation %>%
    mutate(
      differential_group_genotype = case_when(
        differential_group_sampletype %in% c("MetastaticPrimary", "Other") & 
          TERT_ATRX_Mutation == "TERT"  ~ paste(TERT_ATRX_Mutation, 
                                                ifelse(TERT_purity_adj_vaf > 0.4, "VAFOK", "VAFLow"), sep="_"),
        differential_group_sampletype %in% c("MetastaticPrimary", "Other") & 
          TERT_ATRX_Mutation == "ATRX"  ~ paste(TERT_ATRX_Mutation, 
                                                ifelse(ATRX_purity_adj_vaf > 0.4, "VAFOK", "VAFLow"), sep="_"),
        TRUE ~ paste(TERT_ATRX_Mutation, "VAFOK", sep="_")
      )
    )
  
  
  ##########
  # Purity #
  ##########
  
  sample_annotation <- sample_annotation %>% 
    mutate(differential_group_purity = if_else(
      `sample_purity` > 0.4, 
      "PurityOK", 
      "PurityLow")) 
  
  ############
  # Combined #
  ############
  
  sample_annotation <- sample_annotation %>% 
    mutate(differential_group = paste(differential_group_sampletype,
                                      differential_group_genotype,
                                      differential_group_purity,
                                      sep="_"),
           differential_group= factor(differential_group, levels = differential_group_levels))
           
  return(sample_annotation)
}
message("Created annotation function add_differential_groups()")