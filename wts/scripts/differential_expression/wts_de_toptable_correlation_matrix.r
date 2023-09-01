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

a5_wts_lcpm_long <- a5_wts_lcpm_list[["SDHB_abdothoracic"]] %>%
  as_tibble(rownames="ensgid_symbol") %>% 
  pivot_longer(cols = -ensgid_symbol, names_to = "A5_ID", values_to = "log2cpm") %>% 
  arrange(ensgid_symbol, A5_ID) %>% 
  mutate(ensgid_symbol = factor(ensgid_symbol),
         A5_ID = factor(A5_ID))


cor_mat_rho <- list()
cor_mat_pval <- list()

checkpoint_rds_cor_mat_rho <- "/g/data/pq08/projects/ppgl/a5/wts/quickload_checkpoints/wts_de_toptable_cor_mat_rho.rds"
checkpoint_rds_cor_mat_pval <- "/g/data/pq08/projects/ppgl/a5/wts/quickload_checkpoints/wts_de_toptable_cor_mat_pval.rds"
if(quickload_correlation_matrix)
{
  cor_mat_rho <- readRDS(checkpoint_rds_cor_mat_rho)
  cor_mat_pval <- readRDS(checkpoint_rds_cor_mat_pval)
} else
{
  
  threads=27
  options(future.debug = TRUE)
  plan(strategy="multisession", workers=threads)

  for (contrast in c("TERT_PriMet_vs_NonMetPri_WT", "ATRX_PriMet_vs_NonMetPri_WT", "Metastasis_All_vs_NonMetPri_WT"))
  {
    
    tt_sig <- wts_top_tables[["genosampletype"]][[contrast]] %>%
      filter(adj.P.Val < 0.01)
    
    gene_list <- unique(c(tt_sig$Gene, "ENSG00000148773.14_MKI67"))
    
    cor_mat_rho[[contrast]] <- cor(t(a5_wts_lcpm_list[["SDHB_abdothoracic"]][gene_list,]), use = "pairwise.complete.obs", method = "spearman")
    
    cor_mat_pval[[contrast]] <- matrix(nrow = length(gene_list), ncol = length(gene_list), dimnames = list(gene_list,gene_list))

    cor_list_i <- furrr::future_map(set_names(gene_list,gene_list),
                                    .f = \(goi1) {
                                      cor_list_j <- list()
                                      goi1_expr <- a5_wts_lcpm_long %>%
                                        filter(ensgid_symbol==goi1) %>% pull(log2cpm)

                                      for (goi2 in gene_list)
                                      {
                                        goi2_expr <- a5_wts_lcpm_long %>%
                                          filter(ensgid_symbol==goi2) %>% pull(log2cpm)

                                        corr_test_result <- cor.test(x = goi1_expr,y = goi2_expr, method=c("spearman"))

                                        cor_list_j[[goi2]] <- corr_test_result$p.value
                                      }

                                      return(cor_list_j)
                                    })


    for(i in names(cor_list_i)){
      for (j in names(cor_list_i[[i]]))
      {
        cor_mat_pval[[contrast]][i,j] <- cor_list_i[[i]][[j]]
      }
    }
    
   
    
  } 
  saveRDS(cor_mat_rho, checkpoint_rds_cor_mat_rho)
  saveRDS(cor_mat_pval, checkpoint_rds_cor_mat_pval)
}

