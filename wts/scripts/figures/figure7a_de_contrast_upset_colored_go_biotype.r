library(ggupset)

if(!exists("padj_cutoff")) { padj_cutoff = 0.05 }
if(!exists("lfc_cutoff")) { lfc_cutoff = 1 }


source("/g/data/pq08/projects/ppgl/a5/./wts/scripts/figures/a5_wts_de_annotated_volcano_plots.r")

contrasts_to_use <- c("TERT_PriMet_vs_NonMetPri_WT","ATRX_PriMet_vs_NonMetPri_WT","Metastatic_All_vs_NonMetPri_WT", "ATRX_All_vs_TERT_All")

de_genes <- purrr::map2(.x = wts_top_tables$genosampletype[contrasts_to_use], 
                        .y=names(wts_top_tables$genosampletype[contrasts_to_use]), 
                        .f = \(tt, source) { tt %>% 
                            filter(adj.P.Val < padj_cutoff, abs(logFC) > lfc_cutoff) %>% 
                            dplyr::select(Gene, gene_biotype,  ensembl_gene_id) %>% mutate(source=source) %>% left_join(goterms) }) %>% 
  bind_rows()


upset_data <- de_genes %>% 
  group_by(Gene) %>% 
  summarise(source=list(source))

gg_upset <- ggplot(upset_data, aes(x=source)) +
  geom_bar() +
  #geom_text(mapping=aes(y=1, label=after_stat(count))) +
  scale_x_upset() + ylab("Genes (n)")

# go_biotype_cols <- c("Cell adhesion" = c25[[1]],
#                      "Cell cycle/proliferation" = c25[[2]],
#                      "Cell migration" = c25[[3]],
#                      "Cellular developmental process" = c25[[4]],
#                      "Cellular response to DNA damage stimulus" = c25[[5]],
#                      "Chromatin/chromosome organization" = c25[[6]],
#                      "Nervous system development/neurogenesis" = c25[[7]],
#                      "None" = "gray30",
#                      "Other" = "gray70",
#                      "Programmed cell death" = c25[[10]],
#                      "Pseudogene" = c25[[11]],
#                      "Regulation of DNA-templated transcription" = c25[[23]],
#                      "lncRNA" = c25[[18]],
#                      "Angiogenesis" = c25[[14]],
#                      "DNA replication"  = c25[[15]])

plot_data <- upset_data %>% 
  rowwise() %>% 
  mutate(source_string = toString(source))

combo_order <- upset_data %>% 
  rowwise() %>% 
  mutate(source_string = toString(source)) %>%  
  group_by(source_string) %>% 
  dplyr::count() %>% 
  arrange(desc(n)) %>% 
  pull(source_string)

plot_data <- plot_data %>% mutate(source_string = factor(source_string, levels=combo_order)) 

plot_data <- plot_data %>%  
  ungroup() %>% 
  left_join(de_genes %>%  
              dplyr::select(Gene, gene_biotype, parent_go_term) %>%  unique()) %>% 
  mutate(GO_biotype = case_when(gene_biotype == "protein_coding" ~ parent_go_term,
                                grepl("pseudogene",gene_biotype) ~ "Pseudogene",
                                gene_biotype == "lncRNA" ~ "lncRNA", 
                                TRUE ~ "Other")) %>% 
  group_by(source_string, GO_biotype) %>% 
  dplyr::count()

gg_go_biotype <- ggplot(plot_data, aes(x=source_string, y=n, fill=GO_biotype)) + 
  geom_col() +
  scale_fill_manual(values=go_biotype_cols)



gg_upset_gtable <- ggplotGrob(gg_upset)
gg_go_biotype_gtable <- ggplotGrob(gg_go_biotype)

gg_go_biotype_gtable_panel <-  gtable::gtable_filter(gg_go_biotype_gtable, "panel")
gg_go_biotype_gtable_axis_l <-  gtable::gtable_filter(gg_go_biotype_gtable, "axis-l")

gg_upset_gtable$grobs[[6]] <- gg_go_biotype_gtable_panel$grobs[[1]]
gg_upset_gtable$grobs[[3]] <- gg_go_biotype_gtable_axis_l$grobs[[1]]

gg_go_biotype_gtable_guide_box <-  gtable::gtable_filter(gg_go_biotype_gtable, "guide-box")
gg_upset_gtable <- gtable::gtable_add_cols(gg_upset_gtable, gg_go_biotype_gtable$widths[[9]], 9)
gg_upset_gtable <- gtable::gtable_add_rows(gg_upset_gtable, heights = gg_go_biotype_gtable$heights[[9]], 15)

gg_upset_gtable <- gtable::gtable_add_grob(x = gg_upset_gtable,
                                           grobs = gg_go_biotype_gtable_guide_box,
                                           t = 7,
                                           l = 9,
                                           b = 9,
                                           r = 10,
                                           z = 14,
                                           name="guide-box")

grid.newpage()
grid.draw(gg_upset_gtable)

