setwd("/g/data/pq08/projects/ppgl")


################
# Data Loaders #
################

quickload_diff_meth = T; quickload_gsea <- T; quickload_dmr <- T;
source("./a5/methylation/scripts/a5_methylation_analysis_v2.r")

source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")

######################
# Plotting functions #
######################

plot_volcano_expr<- function(tt, n_label) {
  
  genes_to_label <- tt %>%
    slice_min(adj.P.Val, n = n_label) %>% 
    pull(Gene)
  
  tt_plot <- tt %>% 
    mutate(significant = if_else(adj.P.Val < 0.05, "Adj. P Value < 0.05", "Not significant")) %>% 
    mutate(label = if_else(Gene %in% genes_to_label, TRUE, FALSE)) %>% 
    mutate(gene_symbol = ifelse(grepl("^A[CL][0-9]", gene_symbol), 
                                stringr::str_extract(string = Gene, pattern = "ENSG[0-9]+"), 
                                gene_symbol))
  
  volcano <- ggplot(tt_plot, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col=significant)) +
    scale_color_manual(values=c("red", "black")) +
    geom_text_repel(
      data = . %>% filter(label == TRUE),
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.1, "lines"),
      segment.color = "grey",
      max.overlaps = Inf,
      size = 4.5,
      aes(label=gene_symbol))
  
  return(volcano)
  
}

plot_volcano_dmr <- function(tt, n_label, ranking_stat="Fisher", meandiff_cutoff=0.25) {
  
  tt_plot <- tt %>% 
    mutate(significant = if_else(#!!sym(ranking_stat) < 10^-16 | 
      (!!sym(ranking_stat) < 0.05 & abs(meandiff) > meandiff_cutoff), 
      "Adj. P Value < 0.05", 
      "Not significant")) %>% 
    arrange(!!sym(ranking_stat)) %>% 
    mutate(label = if_else(row_number() < n_label & significant == "Adj. P Value < 0.05", TRUE, FALSE))
  
  volcano <- ggplot(tt_plot, aes(meandiff, -log10(!!sym(ranking_stat)))) +
    geom_point(aes(col=significant)) +
    scale_color_manual(values=c("red", "black")) +
    geom_text_repel(
      data = . %>% filter(label == TRUE),
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.1, "lines"),
      segment.color = "grey",
      max.overlaps = Inf,
      size = 4.5,
      aes(label=overlapping.genes))
  
  return(volcano)
  
}

gg_expr <- plot_volcano_expr(wts_top_tables[["Parasympathetic_vs_Sympathetic"]][["Parasympathetic_vs_Sympathetic"]], 30) + ggtitle("Head and Neck vs Abdo-thoracic")

gg_dmr <- plot_volcano_dmr(GenomicRanges::as.data.frame(dmr_lists[["hn"]][["Parasympathetic_vs_Sympathetic"]]), 30) + ggtitle("Head and Neck vs Abdo-thoracic")


gg_composed <- gg_expr + gg_dmr + plot_layout(guides = "collect")  & theme(aspect.ratio = 1)

ggsave(plot = gg_composed, 
       filename = "./a5/tertiary/results/figures/top_dex_dmr_volcano_plot.pdf",
       width =180, height = 90, units = "mm",device = pdf(),  scale = 2
)
