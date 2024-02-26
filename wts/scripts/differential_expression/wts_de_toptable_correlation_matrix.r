##########################################
# This script generates performs         #
# spearman correlation analysis of       #
# genes found to be DE between           #
# TERT/ATRX/Met/Non-met tumours          #
#                                        #
# Author: Aidan Flynn                    #
# Date: 01/09/2023                       #
#                                        #
##########################################

#library(ggplot2)
#library(ggrepel)
#library(patchwork)
library(furrr)

setwd("/g/data/pq08/projects/ppgl")

checkpoint_rds_cor_mat_rho <- "./a5/wts/quickload_checkpoints/wts_de_toptable_cor_mat_rho.rds"
checkpoint_rds_cor_mat_pval <- "./a5/wts/quickload_checkpoints/wts_de_toptable_cor_mat_pval.rds"

if (!exists("quickload_correlation_matrix")) { quickload_correlation_matrix <- FALSE }

if(quickload_correlation_matrix)
{
  message("Quickloading DE gene correlations")
  cor_mat_rho <- readRDS(checkpoint_rds_cor_mat_rho)
  cor_mat_pval <- readRDS(checkpoint_rds_cor_mat_pval)
} else
{
  message("Computing DE gene correlations")
#################
# Run DE script #
#################

#Performs DE
# creates globals:
# - EnsIds
# - data_loader_ensgid_to_chr()
# - wts_top_tables[]
if(!exists("wts_top_tables")) {
  source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")
}

#############################
# TopTable gene correlation #
#############################

cor_mat_rho <- list()
cor_mat_pval <- list()

adj_pval_cutoff <- 0.05
  
  if(!exists("threads")) { threads <- 12 }
  options(future.debug = FALSE)
  plan(strategy="multisession", workers=threads)
  
  for (contrast in c("TERT_PriMet_vs_NonMetPri_WT", "ATRX_PriMet_vs_NonMetPri_WT", "ATRX_All_vs_TERT_All", "Metastatic_All_vs_NonMetPri_WT"))
  {
    
    tt_sig <- wts_top_tables[["genosampletype"]][[contrast]] %>%
      filter(adj.P.Val < adj_pval_cutoff)
    
    gene_list <- unique(c(tt_sig$Gene, "ENSG00000148773.14_MKI67"))
    
    gene_list_expr <- a5_wts_lcpm_list[["SDHB_abdothoracic"]][gene_list,]

    #Create list of non-redundant comparisons
    comp_list <- 
        data.frame(
          t(combn(x = gene_list, 
                m = 2, 
                repl=FALSE)))
    colnames(comp_list) <- c("goi1","goi2")
    
    
    #Threaded iteration over gene list returning list of lists
    cor_list_i <- furrr::future_map(set_names(gene_list,gene_list),
                                    .f = \(goi1) {
                                      gene_list_j <- comp_list$goi2[comp_list$goi1 == goi1]
                                      cor_list_j <- list()

                                      for (goi2 in gene_list_j)
                                      {
                                        cor_list_j[[goi2]] <- list()
                                        corr_test_result <- cor.test(x = gene_list_expr[goi1,],y = gene_list_expr[goi2,], method=c("spearman"))

                                        cor_list_j[[goi2]]$pval <- corr_test_result$p.value
                                        cor_list_j[[goi2]]$rho <- corr_test_result$estimate
                                      }

                                      return(cor_list_j)
                                    })

    cor_mat_rho[[contrast]] <-  matrix(nrow = length(gene_list), ncol = length(gene_list), dimnames = list(gene_list,gene_list))
    cor_mat_pval[[contrast]] <- matrix(nrow = length(gene_list), ncol = length(gene_list), dimnames = list(gene_list,gene_list))
    
    #Convert list of lists into matrix
    for(i in names(cor_list_i)){
      for (j in names(cor_list_i[[i]]))
      {
        cor_mat_rho[[contrast]][i,j] <- cor_list_i[[i]][[j]][["rho"]]
        cor_mat_rho[[contrast]][j,i] <- cor_list_i[[i]][[j]][["rho"]]
        cor_mat_pval[[contrast]][i,j] <- cor_list_i[[i]][[j]][["pval"]]
        cor_mat_pval[[contrast]][j,i] <- cor_list_i[[i]][[j]][["pval"]]
      }
    }
    
    #fill in self_comps (skipped for efficiency)
    for(i in 1:ncol(cor_mat_rho[[contrast]])){
      cor_mat_rho[[contrast]][i,i] <- 1
      cor_mat_pval[[contrast]][i,i] <- 0
      }
    
  } 
  saveRDS(cor_mat_rho, checkpoint_rds_cor_mat_rho)
  saveRDS(cor_mat_pval, checkpoint_rds_cor_mat_pval)
}

