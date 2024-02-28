######################################################################
# Script to generate contrast matrices for                           # 
# differential expression and methylation                            #
# Author: Aidan Flynn                                                #
# Date: 14/10/2022                                                   #
# Languages: R                                                       #
######################################################################

library(rlang)

#Creates contrast and design matrices for head/neck versus abdothoracic comparisons
make_hn_vs_abdominothoracic_contrasts <- function(sample_anno, exclude_samples=NULL)
{
  differential_group_anatomy <- sample_anno$differential_group_anatomy
  differential_group_driver <- sample_anno$TERT_ATRX_Mutation
  
  differential_group <- paste(differential_group_driver,differential_group_anatomy, sep="_")
  
  if(!is.null(exclude_samples))
  {
    exclude_idx <- which(as.character(sample_anno$A5_ID) %in% as.character(exclude_samples))
    differential_group[exclude_idx] <- "Exclude"
    message("Marked ", toString(exclude_samples), " for exclusion from contrasts")
  }
  
  sex <- sample_anno$Gender
  design_mat <- model.matrix(~0 + differential_group + sex)
  colnames(design_mat) <- gsub("differential_group","", colnames(design_mat))
  rownames(design_mat) <- sample_anno$A5_ID
  
  contr_matrix <- makeContrasts(
    Non_chromaffin_vs_Chromaffin = (ATRX_Abdominal_Thoracic + TERT_Abdominal_Thoracic + WT_Abdominal_Thoracic)/3 - (WT_Head_neck),
    Non_chromaffin_vs_Chromaffin_WT = (WT_Abdominal_Thoracic) - (WT_Head_neck),
    levels = colnames(design_mat))
  
  assign("contrast_matrix_hn", contr_matrix, globalenv())
  assign( "design_matrix_hn", design_mat, globalenv())
  
  message("Added contrast matrix 'contrast_matrix_hn' to the global environment")
  message("Added design matrix 'design_matrix_hn' to the global environment")
}
message("Created contrast function make_hn_vs_abdominothoracic_contrasts()")

#Creates contrast and design matrices for comparisons between ATRX/TERT/WT plus 
# metastatic/non-metastatic samples
make_genotype_sampletype_contrasts <- function(sample_anno, exclude_samples=NULL)
{
  # make a design matrix 
  if(class(sample_anno$differential_group) == "factor")
  {
    differential_group <- forcats::fct_drop(sample_anno$differential_group)
  } else {
    differential_group <- factor(sample_anno$differential_group)
  }
  
  if(!is.null(exclude_samples))
  {
    exclude_idx <- which(as.character(sample_anno$A5_ID) %in% as.character(exclude_samples))
    levels(differential_group) <- c(levels(differential_group), "Exclude")
    differential_group[exclude_idx] <- "Exclude"
    message("Marked ", toString(exclude_samples), " for exclusion from contrasts")
  }
  sex <- sample_anno$Gender
  design_mat <- model.matrix(~0 + differential_group + sex)
  colnames(design_mat) <- gsub("differential_group","", colnames(design_mat))
  rownames(design_mat) <- sample_anno$A5_ID
  
  contr_matrix <- makeContrasts(
    
    #TERT/ATRX vs Non-Metastatic Primary
    TERT_PriMet_vs_NonMetPri_WT = 
      (MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK+
          MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK+
          Metastasis_TERT_VAFOK_PurityOK)/3 -
         NMP_WT_VAFOK_PurityOK,
    
    ATRX_PriMet_vs_NonMetPri_WT = 
      (SFU_ATRX_VAFOK_PurityOK+
          MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK+
          MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK+
          Metastasis_ATRX_VAFOK_PurityOK)/4 -
         NMP_WT_VAFOK_PurityOK,
    
    #All-TERT vs Non-TERT
    TERT_All_vs_NonTERT = 
      (MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK +
         Metastasis_TERT_VAFOK_PurityOK)/3 -
      (NMP_WT_VAFOK_PurityOK +
         SFU_ATRX_VAFOK_PurityOK +
         SFU_WT_VAFOK_PurityOK +
         MetastaticPrimaryUnconfirmed_WT_VAFOK_PurityOK +
         MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK +
         Metastasis_ATRX_VAFOK_PurityOK +
         Metastasis_WT_VAFOK_PurityOK)/9,
    
    #All-ATRX vs Non-ATRX
    ATRX_All_vs_NonATRX = 
      (SFU_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK +
         Metastasis_ATRX_VAFOK_PurityOK)/4 -
      (NMP_WT_VAFOK_PurityOK +
         SFU_WT_VAFOK_PurityOK +
         MetastaticPrimaryUnconfirmed_TERT_VAFLow_PurityOK +
         MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK +
         MetastaticPrimaryUnconfirmed_WT_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_TERT_VAFLow_PurityOK +
         MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK +
         Metastasis_TERT_VAFOK_PurityOK +
         Metastasis_WT_VAFOK_PurityOK)/10,
    
    #All-TERT vs All-ATRX
    ATRX_All_vs_TERT_All = 
      (SFU_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK +
         Metastasis_ATRX_VAFOK_PurityOK)/4 -
      (MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK +
         Metastasis_TERT_VAFOK_PurityOK)/3,
    
    #Metastasis Only vs NonMetPri
    Metastasis_Only_vs_NonMetPri_WT = 
      (Metastasis_ATRX_VAFOK_PurityOK +
       Metastasis_TERT_VAFOK_PurityOK +
       Metastasis_WT_VAFOK_PurityOK)/3 -
      NMP_WT_VAFOK_PurityOK,
    
    #Metastasis and confirmed metastatic primaries vs NonMetPri
    Metastatic_Confirmed_vs_NonMetPri_WT = 
      (MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_TERT_VAFLow_PurityOK +
         MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK +
         Metastasis_ATRX_VAFOK_PurityOK +
         Metastasis_TERT_VAFOK_PurityOK +
         Metastasis_WT_VAFOK_PurityOK)/7 -
      NMP_WT_VAFOK_PurityOK,
    
    #Metastases, confirmed metastatic primaries, and unconfirmed with TERT/ATRX vs NonMetPri
    Metastatic_Likely_vs_NonMetPri_WT = 
      (MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK +
         Metastasis_ATRX_VAFOK_PurityOK +
         Metastasis_TERT_VAFOK_PurityOK +
         Metastasis_WT_VAFOK_PurityOK)/8 -
      NMP_WT_VAFOK_PurityOK,
  
    #Metastasis and confirmed or unconfirmed metastatic primaries vs NonMetPri
    Metastatic_All_vs_NonMetPri_WT = 
      (MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryUnconfirmed_TERT_VAFLow_PurityOK +
         MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK +
         MetastaticPrimaryUnconfirmed_WT_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_TERT_VAFLow_PurityOK +
         MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK +
         MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK +
         Metastasis_ATRX_VAFOK_PurityOK +
         Metastasis_TERT_VAFOK_PurityOK +
         Metastasis_WT_VAFOK_PurityOK)/11 -
      NMP_WT_VAFOK_PurityOK,
    levels = colnames(design_mat))
  
  assign("contrast_matrix_genosampletype", contr_matrix, globalenv())
  assign("design_matrix_genosampletype", design_mat, globalenv())
  
  message("Added contrast matrix 'contrast_matrix_genosampletype' to the global environment")
  message("Added design matrix 'design_matrix_genosampletype' to the global environment")
}
message("Created contrast function make_genotype_sampletype_contrasts()")


#This function produces a table of how many samples participate in each contrast group
count_contrast_members <- function(contr_matrix, design.matrix)
{
  
  #Check number of samples participating in each contrast
  nContrasts = dim(contr_matrix)[2]
  
  contrast_member_table <- data.frame(contrast=vector(mode="character"), 
                                      groupA_count=vector(mode="integer"), 
                                      groupB_count=vector(mode="integer"))
  for (i in 1:nContrasts)
  {
    groupA.idx <- which(contr_matrix[,dimnames(contr_matrix)$Contrasts[i]] > 0)
    groupB.idx <- which(contr_matrix[,dimnames(contr_matrix)$Contrasts[i]] < 0)
    groupA.levels <- dimnames(contr_matrix)$Levels[groupA.idx]
    groupB.levels <- dimnames(contr_matrix)$Levels[groupB.idx]
    
    groupA.count <- sum(rowSums(design.matrix[,groupA.levels,drop=F]) > 0)
    groupB.count <- sum(rowSums(design.matrix[,groupB.levels,drop=F]) > 0)
    
    contrast_member_table <- rbind(contrast_member_table, data.frame(contrast=dimnames(contr_matrix)$Contrasts[i], 
                                                                     groupA_count=groupA.count, 
                                                                     groupB_count=groupB.count))
  }
  
  return(contrast_member_table)
}
message("Created contrast function count_contrast_members()")


#Convert design matrix into contrast group membership list
contrastdesign_to_memberlist <- function(contrast_name, contrast_matrix, design_matrix, contrast_name_sep="_vs_")
{
  contrast_name_split=stringr::str_split_1(contrast_name,pattern = contrast_name_sep)
  
  contrast_groups <- data.frame(contrast_matrix) %>% 
    tibble::rownames_to_column("group") %>% 
    filter(group!="sexmale") %>% 
    dplyr::select(group, all_of(contrast_name)) %>% 
    mutate(side=case_when(
      !!sym(contrast_name) > 0 ~ contrast_name_split[[1]],
      !!sym(contrast_name) < 0 ~ contrast_name_split[[2]],
      !!sym(contrast_name) == 0 ~ "non_participant")) %>% 
    group_by(side) %>%
    summarise(members=list(group))
  
  contrast_groups <- setNames(contrast_groups$members,contrast_groups$side)
  
  if(length(contrast_groups[[contrast_name_split[[1]]]]) > 1) {
    groupA <- rowSums(design_matrix[,contrast_groups[[contrast_name_split[[1]]]]])
  } else {
    groupA <- design_matrix[,contrast_groups[[contrast_name_split[[1]]]]]
  }
  
  if(length(contrast_groups[[contrast_name_split[[2]]]]) > 1) {
    groupB <- rowSums(design_matrix[,contrast_groups[[contrast_name_split[[2]]]]])
  } else {
    groupB <- design_matrix[,contrast_groups[[contrast_name_split[[2]]]]]
  }
  
  if(length(contrast_groups[["non_participant"]]) > 1) {
    groupC <- rowSums(design_matrix[,contrast_groups[["non_participant"]]])
  } else {
    groupC <- design_matrix[,contrast_groups[["non_participant"]]]
  }
  
  membership <- setNames(data.frame(groupA, 
                                    groupB, 
                                    groupC), 
                         c(contrast_name_split[[1]],
                           contrast_name_split[[2]], 
                           "non_participant")) %>% 
    tibble::rownames_to_column("A5_ID") %>% 
    pivot_longer(cols = -A5_ID,names_to = "group", values_to = "participation") %>% 
    filter(participation==1) %>% dplyr::select(-participation)
  
  if(any(membership$A5_ID != dimnames(design_matrix)[[1]])) {
    warning("Return sample order does not match design matrix order")
  }
  
  return(membership)
}
message("Created contrast function contrastdesign_to_memberlist()")

simplify_design_contrast <- function(contrast, design_matrix, contrast_matrix, contrast_order = NULL){
  
  participant_list <- 
    contrastdesign_to_memberlist(contrast_name = contrast, 
                                 design_matrix = design_matrix, 
                                 contrast_matrix = contrast_matrix)
  
  if(!all(participant_list$A5_ID == rownames(design_matrix))) { stop("Participant list order sanity check failed")}
  
  sex <- ifelse(design_matrix[,"sexmale"], "male", "female")
  
  differential_group <- participant_list$group
  
  simplified_design_mat <- model.matrix(~0 + differential_group + sex)
  colnames(simplified_design_mat) <- gsub("differential_group","", colnames(simplified_design_mat))
  rownames(simplified_design_mat) <- rownames(design_matrix)
  
  if (!is.null(contrast_order)) {
    contrast_levels <- c(contrast_order, "non_participant")
  } else
  {
    contrast_levels <- c(setdiff(unique(differential_group),"non_participant"), "non_participant")
  }
  
  
  contrast_string <- setNames(paste(contrast_levels[[1]],contrast_levels[[2]], sep=" - "), contrast)
  
  simplified_contrast_matrix <- rlang::inject(makeContrasts(!!!contrast_string,  levels = colnames(simplified_design_mat)))
  
  return(list(contrast_matrix = simplified_contrast_matrix, design_matrix = simplified_design_mat))
  
}



#############
## RECODING #
#############

contrast_recode = list(
  TERT_PriMet_vs_NonMetPri_WT = 
    c(MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK = "Primary_TERT",
      MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK = "Primary_TERT",
      Metastasis_TERT_VAFOK_PurityOK = "Metastasis_TERT",
      NMP_WT_VAFOK_PurityOK = "NMP"),
  
  
  ATRX_PriMet_vs_NonMetPri_WT = 
    c(SFU_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      Metastasis_ATRX_VAFOK_PurityOK = "Metastasis_ATRX",
      NMP_WT_VAFOK_PurityOK = "NMP"),
  
  
  #All-TERT vs Non-TERT
  TERT_All_vs_NonTERT = 
    c(MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK = "Primary_TERT",
      MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK = "Primary_TERT",
      Metastasis_TERT_VAFOK_PurityOK = "Metastasis_TERT",
      NMP_WT_VAFOK_PurityOK = "Primary_WT",
      SFU_ATRX_VAFOK_PurityOK  = "Primary_ATRX",
      SFU_WT_VAFOK_PurityOK = "Primary_WT",
      MetastaticPrimaryUnconfirmed_WT_VAFOK_PurityOK = "Primary_WT",
      MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK = "Primary_WT",
      Metastasis_ATRX_VAFOK_PurityOK = "Metastasis_ATRX",
      Metastasis_WT_VAFOK_PurityOK = "Metastasis_WT"),
  
  #All-ATRX vs Non-ATRX
  ATRX_All_vs_NonATRX = 
    c(SFU_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      Metastasis_ATRX_VAFOK_PurityOK = "Metastasis_ATRX",
      NMP_WT_VAFOK_PurityOK = "Primary_WT",
      SFU_WT_VAFOK_PurityOK = "Primary_WT",
      MetastaticPrimaryUnconfirmed_TERT_VAFLow_PurityOK  = "Primary_TERT",
      MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK  = "Primary_TERT",
      MetastaticPrimaryUnconfirmed_WT_VAFOK_PurityOK  = "Primary_WT",
      MetastaticPrimaryConfirmed_TERT_VAFLow_PurityOK  = "Primary_TERT",
      MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK  = "Primary_TERT",
      MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK  = "Primary_WT",
      Metastasis_TERT_VAFOK_PurityOK = "Metastasis_TERT",
      Metastasis_WT_VAFOK_PurityOK = "Metastasis_WT"),
  
  #All-TERT vs All-ATRX
  ATRX_All_vs_TERT_All = 
    c(SFU_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK = "Primary_ATRX",
      Metastasis_ATRX_VAFOK_PurityOK = "Metastasis_ATRX",
      MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK  = "Primary_TERT",
      MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK  = "Primary_TERT",
      Metastasis_TERT_VAFOK_PurityOK = "Metastasis_TERT"),
  
  Metastatic_Confirmed_vs_NonMetPri_WT = 
    c(MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK = "MetastaticPrimary_ATRX",
      MetastaticPrimaryConfirmed_TERT_VAFLow_PurityOK  = "MetastaticPrimary_TERT",
      MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK  = "MetastaticPrimary_TERT",
      MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK  = "MetastaticPrimary_WT",
      Metastasis_ATRX_VAFOK_PurityOK = "Metastasis_ATRX",
      Metastasis_TERT_VAFOK_PurityOK = "Metastasis_TERT",
      Metastasis_WT_VAFOK_PurityOK = "Metastasis_WT",
      NMP_WT_VAFOK_PurityOK = "NMP"),
  
  #Metastases, confirmed metastatic primaries, and unconfirmed with TERT/ATRX vs NonMetPri
  Metastatic_Likely_vs_NonMetPri_WT = 
    c(MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK = "MetastaticPrimary_ATRX",
      MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK = "MetastaticPrimary_TERT",
      MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK = "MetastaticPrimary_ATRX",
      MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK = "MetastaticPrimary_TERT",
      MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK = "MetastaticPrimary_WT",
      Metastasis_ATRX_VAFOK_PurityOK = "Metastasis_ATRX",
      Metastasis_TERT_VAFOK_PurityOK = "Metastasis_TERT",
      Metastasis_WT_VAFOK_PurityOK = "Metastasis_WT",
      NMP_WT_VAFOK_PurityOK = "NMP"),
  
  #Metastasis and confirmed or unconfirmed metastatic primaries vs NonMetPri
  Metastatic_All_vs_NonMetPri_WT = 
    c(MetastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK = "MetastaticPrimary_ATRX",
      MetastaticPrimaryUnconfirmed_TERT_VAFLow_PurityOK = "MetastaticPrimary_TERT",
      MetastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK = "MetastaticPrimary_TERT",
      MetastaticPrimaryUnconfirmed_WT_VAFOK_PurityOK = "MetastaticPrimary_WT",
      MetastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK = "MetastaticPrimary_ATRX",
      MetastaticPrimaryConfirmed_TERT_VAFLow_PurityOK = "MetastaticPrimary_TERT",
      MetastaticPrimaryConfirmed_TERT_VAFOK_PurityOK = "MetastaticPrimary_TERT",
      MetastaticPrimaryConfirmed_WT_VAFOK_PurityOK = "MetastaticPrimary_WT",
      Metastasis_ATRX_VAFOK_PurityOK = "Metastasis_ATRX",
      Metastasis_TERT_VAFOK_PurityOK = "Metastasis_TERT",
      Metastasis_WT_VAFOK_PurityOK = "Metastasis_WT",
      NMP_WT_VAFOK_PurityOK = "NMP"))


recoded_contrast_strings = c(
  #All-TERT vs NonMetPri
  TERT_PriMet_vs_NonMetPri_WT = "(Primary_TERT + Metastasis_TERT)/2 - NMP",
  
  #All-ATRX vs NonMetPri
  ATRX_PriMet_vs_NonMetPri_WT = "(Primary_ATRX + Metastasis_ATRX)/2 - NMP",
  
  #All-TERT vs Non-TERT
  TERT_All_vs_NonTERT = "(Primary_TERT + Metastasis_TERT)/2 - (Primary_WT + Primary_ATRX + Metastasis_ATRX + Metastasis_WT)/4",    
  
  #All-ATRX vs Non-ATRX
  ATRX_All_vs_NonATRX = "(Primary_ATRX + Metastasis_ATRX)/2 - (Primary_WT + Primary_TERT + Metastasis_TERT + Metastasis_WT)/4",    
  
  #All-TERT vs All-ATRX
  ATRX_All_vs_TERT_All = "(Primary_ATRX + Metastasis_ATRX)/2 - (Primary_TERT + Metastasis_TERT)/2",
  
  #Metastasis Only vs NonMetPri
  Metastasis_Only_vs_NonMetPri_WT = "(Metastasis_ATRX_VAFOK_PurityOK + Metastasis_TERT_VAFOK_PurityOK + Metastasis_WT_VAFOK_PurityOK)/3 - NMP_WT_VAFOK_PurityOK",    
  
  #Metastasis and confirmed metastatic primaries vs NonMetPri
  Metastatic_Confirmed_vs_NonMetPri_WT = "(MetastaticPrimary_ATRX + MetastaticPrimary_TERT + MetastaticPrimary_WT + Metastasis_ATRX + Metastasis_TERT + Metastasis_WT)/6 - NMP",    
  
  #Metastases, confirmed metastatic primaries, and unconfirmed with TERT/ATRX vs NonMetPri
  Metastatic_Likely_vs_NonMetPri_WT = "(MetastaticPrimary_ATRX + MetastaticPrimary_TERT + MetastaticPrimary_WT + Metastasis_ATRX + Metastasis_TERT + Metastasis_WT)/6 - NMP",  
  
  #Metastasis and confirmed or unconfirmed metastatic primaries vs NonMetPri
  Metastatic_All_vs_NonMetPri_WT = "(MetastaticPrimary_ATRX + MetastaticPrimary_TERT + MetastaticPrimary_WT + Metastasis_ATRX + Metastasis_TERT + Metastasis_WT)/6 - NMP")