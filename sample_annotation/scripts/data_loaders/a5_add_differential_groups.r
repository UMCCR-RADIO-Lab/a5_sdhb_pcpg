#Helper script to append comparison groups for differential methylation/expression analysis
#to sample annotation dataframe

if(!exists("differential_group_levels"))
{
  differential_group_levels <- c(
    "NMP_WT_VAFOK_PurityLow",
    "NMP_WT_VAFOK_PurityOK",
    
    "NMLR_WT_VAFOK_PurityOK",
    
    "SFU_ATRX_VAFOK_PurityOK",
    "SFU_ATRX_VAFLow_PurityOK",
    "SFU_TERT_VAFOK_PurityLow",
    "SFU_WT_VAFOK_PurityOK",
    
    "MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityLow",
    "MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK",
    "MetastaticPrimaryUnconfirmed_TERT_VAFLow_PurityLow", 
    "MetastaticPrimaryUnconfirmed_TERT_VAFLow_PurityOK", 
    "MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK",
    "MetastaticPrimaryUnconfirmed_WT_VAFOK_PurityLow",
    "MetastaticPrimaryUnconfirmed_WT_VAFOK_PurityOK",
    
    "MetastaticRecurranceUnconfirmed_ATRX_VAFOK_PurityOK",
    "MetastaticRecurranceUnconfirmed_TERT_VAFOK_PurityOK", 
    "MetastaticRecurranceUnconfirmed_WT_VAFOK_PurityLow",
    
    "MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK",
    "MetastaticPrimaryConfirmed_TERT_VAFLow_PurityOK",
    "MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK",
    "MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK",
    
    "MetastaticRecurranceConfirmed_TERT_VAFOK_PurityOK", 
    
    "Metastasis_ATRX_VAFOK_PurityOK",
    "Metastasis_TERT_VAFOK_PurityLow",
    "Metastasis_TERT_VAFOK_PurityOK",
    "Metastasis_WT_VAFOK_PurityOK",
    
    "Normal_WT_VAFOK_PurityOK")
  
  lockBinding("differential_group_levels", globalenv())
}

message("Created global variable 'differential_group_levels'")

add_differential_groups <- function(sample_annotation) {
  
  ###########
  # Anatomy #
  ###########
  
  sample_annotation <- sample_annotation %>%  mutate(differential_group_anatomy = case_when(
    Primary_Location_Simplified=="Normal" ~ "Normal",
    Primary_Location_Simplified=="Head_neck" & grepl("Parasympathetic", Major_Cluster) ~ "Head_Neck",
    Primary_Location_Simplified=="Extraadrenal_thoracic_mediastinum" ~ "Mediastinum",
    Primary_Location_Simplified!="Head_neck" & grepl("Sympathetic", Major_Cluster) ~ "Abdominal_Thoracic",
    TRUE ~ "Ambiguous"
  ))
  
  #################
  # Sample Type  #
  #################
  
  sample_annotation <- sample_annotation %>%
    mutate(
      differential_group_sampletype = case_when(
        is_primary_or_met == "Normal" ~ "Normal",
        is_primary_or_met %in% c("Primary", "Recurrent") &
          tumour_metastasised == "Yes" ~ "MetastaticPrimary",
        is_primary_or_met  %in% c("Primary", "Recurrent") &
          tumour_metastasised == "No" ~ "NonMetastaticPrimary",
        is_primary_or_met  %in% c("Primary", "Recurrent") &
          tumour_metastasised == "Short follow up" ~ "SFU",
        is_primary_or_met == "Metastasis" &
          tumour_metastasised == "Yes" ~ "Metastasis",
        TRUE ~ "Other")
    )
  
  #######################
  # Sample Type Strict  #
  #######################
  
  #strict = primaries and relapse only considered metastatic if 
  # primary/met relationship proven by shared genomic features
  
  sample_annotation <- sample_annotation %>%
    mutate(
      differential_group_sampletype_strict = case_when(
        #Met/primary pairs
        A5_ID %in% c("E143-3", "E146-1", "E158-1",  "E159-3", "E225-1") ~ "Metastatic primary",
        A5_ID == "E132-2" ~ "Metastatic local recurrence",
        #Additional primary
        A5_ID %in% c("E159-2") ~ "Primary (metastasis reported)",
        #normals
        is_primary_or_met == "Normal" ~ "Normal",
        #Primaries
        is_primary_or_met  == "Primary" &
          tumour_metastasised == "Yes" ~ "Primary (metastasis reported)",
        is_primary_or_met  == "Primary" &
          tumour_metastasised == "Short follow up" ~ "Primary (short follow up)",
        is_primary_or_met  == "Primary" &
          tumour_metastasised == "No" ~ "Non-metastatic primary",
        #Recurrences
        is_primary_or_met == "Recurrent" &
          tumour_metastasised == "Yes" ~ "Local recurrence (metastasis reported)",
        is_primary_or_met == "Recurrent" &
          tumour_metastasised == "No" ~ "Non-metastatic local recurrence",
        #Mets
        is_primary_or_met == "Metastasis" &
          tumour_metastasised == "Yes" ~ "Metastasis",
        
        TRUE ~ "Other"
      )) %>% 
    mutate(differential_group_sampletype_strict = 
             factor(differential_group_sampletype_strict, 
                    levels=c("Normal",
                             "Non-metastatic primary",
                             "Primary (short follow up)",
                             "Non-metastatic local recurrence",
                             "Primary (metastasis reported)",
                             "Local recurrence (metastasis reported)",
                             "Metastatic primary",
                             "Metastatic local recurrence",
                             "Metastasis"))
    )
  
  ############
  # Genotype #
  ############
  
  sample_annotation <- sample_annotation %>%  mutate(TERT_ATRX_Mutation=case_match(A5_ID, "E171-1"~"TERT", .default = TERT_ATRX_Mutation))
  
  sample_annotation <- sample_annotation %>%
    mutate(
      differential_group_genotype = case_when(
        differential_group_sampletype %in% c("MetastaticPrimary", "SFU", "Other") & 
          TERT_ATRX_Mutation == "TERT"  ~ paste(TERT_ATRX_Mutation, 
                                                ifelse(TERT_purity_adj_vaf > 0.25, "VAFOK", "VAFLow"), sep="_"),
        differential_group_sampletype %in% c("MetastaticPrimary",  "SFU", "Other") & 
          TERT_ATRX_Mutation == "ATRX"  ~ paste(TERT_ATRX_Mutation, 
                                                ifelse(ATRX_purity_adj_vaf > 0.25, "VAFOK", "VAFLow"), sep="_"),
        TRUE ~ paste(TERT_ATRX_Mutation, "VAFOK", sep="_")
      )
    )
  
  
  ##########
  # Purity #
  ##########
  
  sample_annotation <- sample_annotation %>% 
    mutate(differential_group_purity = if_else(
      `sample_purity` >= 0.5, 
      "PurityOK", 
      "PurityLow")) 
  
  ############
  # Combined #
  ############
  
  sample_annotation <- sample_annotation %>% 
    mutate(differential_group = paste(case_match(differential_group_sampletype_strict,
                                                 "Normal" ~ "Normal", 
                                                 "Non-metastatic primary" ~ "NMP",
                                                 "Non-metastatic local recurrence" ~ "NMLR",
                                                 "Primary (short follow up)" ~ "SFU", 
                                                 "Primary (metastasis reported)" ~ "MetastaticPrimaryUnconfirmed", 
                                                 "Local recurrence (metastasis reported)" ~ "MetastaticPrimaryUnconfirmed",
                                                 "Metastatic primary" ~ "MetastaticPrimaryConfirmed", 
                                                 "Metastatic local recurrence"  ~ "MetastaticPrimaryConfirmed",
                                                 "Metastasis" ~ "Metastasis"),
                                      differential_group_genotype,
                                      differential_group_purity,
                                      sep="_"))
  
  if(!all(unique(sample_annotation$differential_group) %in% differential_group_levels))
  {
    stop("Differential group sanity check failed.  Not all differential groups are in the predefined list. Missing entries:", 
         toString(setdiff(unique(sample_annotation$differential_group), differential_group_levels)))
  }
  
  sample_annotation <- sample_annotation %>% 
    mutate(differential_group= factor(differential_group, levels = differential_group_levels))
  
  return(sample_annotation)
}
message("Created annotation function add_differential_groups()")