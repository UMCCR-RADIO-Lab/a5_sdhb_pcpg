source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")
source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#####################
# Genes of Interest #
#####################

adj_pval_cutoff = 0.05

GOI.met <- 
  wts_top_tables[["genosampletype"]][["Metastatic_All_vs_NonMetPri_WT"]] %>% 
  filter(adj.P.Val < adj_pval_cutoff) %>%  
  mutate(source ="Metastatic_All_vs_NonMetPri_WT")


GOI.tert <- wts_top_tables[["genosampletype"]][["TERT_All_vs_NonTERT"]] %>% 
  filter(adj.P.Val < adj_pval_cutoff) %>% 
  mutate(source="TERT_All_vs_NonTERT")


GOI.atrx <-
  wts_top_tables[["genosampletype"]][["ATRX_All_vs_NonATRX"]] %>% 
  filter(adj.P.Val < adj_pval_cutoff) %>%  
  mutate(source="ATRX_All_vs_NonATRX")


GOI <- bind_rows(GOI.met, GOI.tert, GOI.atrx)

#GOI <- GOI %>% group_by(Gene) %>% mutate(source = ifelse(n() == 3, "Common", source)) 

plot_data <- GOI %>% 
  mutate(gene_biotype=ifelse(grepl("pseudogene", gene_biotype), "pseudogene", gene_biotype)) %>% 
  mutate(direction = ifelse(logFC > 0, "Up-regulated", "Down-regulated")) %>% 
  dplyr::select(Gene, gene_biotype, gene_symbol, ensembl_gene_id, direction, source) %>% 
  distinct() %>%
  group_by(source) %>% 
  mutate(source_total=n())  %>%
  group_by(source, gene_biotype, direction, source_total) %>%
  dplyr::count() %>%
  mutate(prop=n/source_total) %>% 
  filter(gene_biotype %in% c("lncRNA","protein_coding", "pseudogene")) 



plot_data <- plot_data %>% 
  mutate(source = factor(as.character(source), 
                         levels = c("Metastatic_All_vs_NonMetPri_WT", 
                                    "TERT_All_vs_NonTERT", "ATRX_All_vs_NonATRX")),
         source=forcats::fct_recode(source, 
                                    "All Metastasis vs NMP"="Metastatic_All_vs_NonMetPri_WT", 
                                    "TERT vs Rest"="TERT_All_vs_NonTERT", 
                                    "ATRX vs Rest"="ATRX_All_vs_NonATRX")) 

plot_data <- plot_data 

fill_pal <- c("All Metastasis vs NMP"=ColorPalette[["DarkRed2"]],
              "TERT vs Rest"=driver_cols[["TERT"]],
              "ATRX vs Rest"=driver_cols[["ATRX"]])

col_pal <- c("Up-regulated"="red",
              "Down-regulated"="blue")

ggplot(plot_data, 
       aes(x=gene_biotype, y=prop, fill=source)) + 
  geom_col(position = position_dodge(), linewidth=1) + 
  scale_fill_manual(values = fill_pal) +
  scale_color_manual(values = col_pal) +
  theme_bw() + 
  theme(axis.text.x =element_text(angle = 90, hjust=1, vjust=0.5)) + 
  facet_wrap("direction") + ylab("Proportion of all DE genes")

ggplot(plot_data, 
       aes(x=direction, y=prop, fill=source)) + 
  geom_col(position = position_dodge(), linewidth=1) + 
  scale_fill_manual(values = fill_pal) +
  scale_color_manual(values = col_pal) +
  theme_bw() + 
  theme(axis.text.x =element_text(angle = 90, hjust=1, vjust=0.5)) + 
  facet_wrap("gene_biotype") + ylab("Proportion of all DE genes")
