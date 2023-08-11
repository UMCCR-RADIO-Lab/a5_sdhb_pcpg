library(ggplots)
library(patchwork)

setwd("/g/data/pq08/projects/ppgl/")

#Subplot scripts
source("./a5/wts/scripts/gene_set_analysis/mpas_activity_score_E167_ras_mut.r")
source("./a5/tertiary/scripts/figures/E167_E169_mgmt_meth_expr.r")
source("./a5/wts/scripts/figures/E167_mlh1_expr.r")

E167_panel <- gg_mgmt_meth$E167 +
  gg_mgmt_expr$E167 +
  gg_mlh1_expr +
  gg_mpas_score + scale_color_manual(values = c("E124-1" = "black", 
                                              "E167-1" = "black",
                                              "E167-2" = "black", 
                                              "Other" = "grey")) +
  xlab("") + guides(color=F)


ggsave(plot = E167_panel, filename = "./a5/tertiary/results/figures/E167_vignette_meth_expr_panel.pdf", device = pdf(), width = 10, height = 24, units = "cm")
