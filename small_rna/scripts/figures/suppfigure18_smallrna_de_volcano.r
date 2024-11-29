library(patchwork)

source("/g/data/pq08/projects/ppgl/a5/small_rna/scripts/differential_expression/a5_smallrna_differential_expression.r")


tert_vs_NMP <- plot_volcano(smallrna_top_tables[["genosampletype"]][["TERT_PriMet_vs_NonMetPri_WT"]], 20) + ggtitle("TERT Mutant Vs non-met. primary")
atrx_vs_NMP <- plot_volcano(smallrna_top_tables[["genosampletype"]][["ATRX_PriMet_vs_NonMetPri_WT"]], 20) + ggtitle("ATRX Mutant Vs non-met. primary")
tert_vs_atrx <- plot_volcano(smallrna_top_tables[["genosampletype"]][["ATRX_All_vs_TERT_All"]], 20) + ggtitle("TERT Vs ATRX")
met_vs_nonmet <- plot_volcano(smallrna_top_tables[["genosampletype"]][["Metastatic_All_vs_NonMetPri_WT"]], 20) + ggtitle("All metastatic (pri/met) vs non-metastatic primaries")

tert_vs_NMP + 
  atrx_vs_NMP + 
  tert_vs_atrx + 
  met_vs_nonmet + 
  plot_layout(nrow=2, guides = "collect")
