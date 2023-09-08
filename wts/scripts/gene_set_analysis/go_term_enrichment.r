library(clusterProfiler)
setwd("/g/data/pq08/projects/ppgl")

################
# Data Loaders #
################

#Differential Expression - WTS
if (!exists("wts_top_tables")) {
  source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")}

if !( exists("quickload_go_enrichment") { quickload_go_enrichment <- F })

checkpointfile_goenrichment <- "./a5/wts/quickload_checkpoints/goenrichment.rds"
if(quickload_go_enrichment)
{
  ego <- readRDS(checkpointfile_goenrichment)
}
else {
  ego <- enrichGO(gene          = gene,
                  universe      = names(geneList),
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  
  
  tt_p_val_cutoff <- 0.05
  tt_lfc_cutoff <- 1
  ego <- list()
  for (comp in names(wts_top_tables))
  {
    ego[[comp]] <- list()
    
    for (contr in names(wts_top_tables[[comp]])) {
      ego[[comp]][[contr]] <- list()
      
      toptable <- wts_top_tables[[comp]][[contr]] %>%
        mutate(ensgid=gsub("[.][0-9]{1,2}_.+$","",Gene))
      
      
      for (Ont in c("CC","MF","BP"))
      {
        ego[[comp]][[contr]][[Ont]] <- enrichGO(gene = toptable %>%  filter(adj.P.Val < tt_p_val_cutoff, abs(logFC) > tt_lfc_cutoff) %>%  pull(ensgid),
                                                universe      = toptable$ensgid,
                                                OrgDb         = org.Hs.eg.db,
                                                keyType       = 'ENSEMBL',
                                                ont           = Ont,
                                                pAdjustMethod = "BH",
                                                pvalueCutoff  = 0.05,
                                                qvalueCutoff  = 0.05,
                                                readable      = TRUE)
        
      }
    }
  }
  
  saveRDS(ego, checkpointfile_goenrichment)
  
}






