######################################################################
# Script to generate contrast matrices for                           # 
# differential expression and methylation                            #
# Author: Aidan Flynn                                                #
# Date: 14/10/2022                                                   #
# Languages: R                                                       #
######################################################################


#Creates contrast and design matrices for head/neck versus abdothoracic comparisons
make_hn_vs_abdominothoracic_contrasts <- function(sample_anno, exclude_samples=NULL)
{
  differential_group_anatomy <- sample_anno$differential_group_anatomy
  if(!is.null(exclude_samples))
  {
    exclude_idx <- which(as.character(sample_anno$A5_ID) %in% as.character(exclude_samples))
    differential_group_anatomy[exclude_idx] <- "Exclude"
    message("Marked ", toString(exclude_samples), " for exclusion from contrasts")
  }
  sex <- sample_anno$Gender
  design_mat <- model.matrix(~0 + differential_group_anatomy + sex)
  colnames(design_mat) <- gsub("differential_group_anatomy","", colnames(design_mat))
  rownames(design_mat) <- sample_anno$A5_ID
    
  contr_matrix <- makeContrasts(
    Parasympathetic_vs_Sympathetic = (Abdominal_Thoracic) - (Head_Neck),
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
    TERT_PriMet_vs_NonMetPri_WT = (Metastasis_TERT_VAFOK_PurityOK + MetastaticPrimary_TERT_VAFOK_PurityOK)/2 - NonMetastaticPrimary_WT_VAFOK_PurityOK,
    ATRX_PriMet_vs_NonMetPri_WT = (Metastasis_ATRX_VAFOK_PurityOK + MetastaticPrimary_ATRX_VAFOK_PurityOK)/2 - NonMetastaticPrimary_WT_VAFOK_PurityOK,
    #All-TERT vs Non-TERT
    TERT_All_vs_NonTERT = (Metastasis_TERT_VAFOK_PurityOK + MetastaticPrimary_TERT_VAFOK_PurityOK)/2 - (Metastasis_ATRX_VAFOK_PurityOK + MetastaticPrimary_ATRX_VAFOK_PurityOK + Metastasis_WT_VAFOK_PurityOK + NonMetastaticPrimary_WT_VAFOK_PurityOK)/4,
    #All-ATRX vs Non-ATRX
    ATRX_All_vs_NonATRX = (Metastasis_ATRX_VAFOK_PurityOK + MetastaticPrimary_ATRX_VAFOK_PurityOK)/2 - 
      (NonMetastaticPrimary_WT_VAFOK_PurityOK + MetastaticPrimary_WT_VAFOK_PurityOK + 
         Metastasis_WT_VAFOK_PurityOK + MetastaticPrimary_TERT_VAFLow_PurityOK + 
         MetastaticPrimary_TERT_VAFOK_PurityOK + Metastasis_TERT_VAFOK_PurityOK + 
         Other_WT_VAFOK_PurityOK)/7,
    #All-TERT vs All-ATRX
    ATRX_All_vs_TERT_All = (Metastasis_ATRX_VAFOK_PurityOK + MetastaticPrimary_ATRX_VAFOK_PurityOK)/2 - (Metastasis_TERT_VAFOK_PurityOK + MetastaticPrimary_TERT_VAFOK_PurityOK)/2,
    #All-TERT vs All-ATRX
    Metastasis_All_vs_NonMetPri_WT = (Metastasis_ATRX_VAFOK_PurityOK + Metastasis_TERT_VAFOK_PurityOK + Metastasis_WT_VAFOK_PurityOK)/3 - NonMetastaticPrimary_WT_VAFOK_PurityOK,
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