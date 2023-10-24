library(GO.db)
library(DBI)

#####################
# Functions Loaders #
#####################

get_go_children <- function (go_parent_id) 
{
  con <- GO_dbconn()
  sql_query = paste0("select goparent.term as parent_go_term, 
                           goparent.go_id as parent_go_id, 
                           gochild.go_id as child_go_id,
                           gochild.term 
                           from 
                           go_term as goparent 
                           inner join go_bp_offspring on go_bp_offspring._id=goparent._id 
                           inner join go_term as gochild on go_bp_offspring._offspring_id=gochild._id where parent_go_id IN (",
                     paste0("'", go_parent_id,"'", collapse = ","), ");")
  go_children <- dbGetQuery(con, sql_query)
  return(go_children)
}

#####################
# Umbrella GO Terms #
#####################

go_parents <- matrix(
  data = c(
    1, "GO:0007049","cell cycle",
    2, "GO:0008283","cell population proliferation",
    3, "GO:0012501","programmed cell death",
    4, "GO:0016477","cell migration",
    5, "GO:0098602","cell adhesion",
    6, "GO:0007155","cell adhesion",
    7, "GO:0006974","DNA damage response",
    8, "GO:0006281","DNA repair",
    9, "GO:0006260", "DNA replication",
    10, "GO:0051276","chromosome organization",
    11, "GO:0006325","chromatin organization",
    12, "GO:0016570","histone modification",
    13, "GO:0032774", "RNA biosynthetic process",
    14, "GO:0048869","cellular developmental process",
    15, "GO:0022008","neurogenesis",
    16, "GO:0001525", "angiogenesis",
    17, "GO:0007399","nervous system development",
    18, "GO:0016192", "vesicle-mediated transport"),
  byrow = T, ncol=3)

dimnames(go_parents) = list(c(go_parents[,1]), c("rank", "go_id", "go_term"))

go_parents <- data.frame(go_parents)

go_parent_child = get_go_children(go_parents$go_id)

################
# Data Loaders #
################

#Differential Expression - WTS
source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")

#Differential Expression - top de gene correlations
quickload_correlation_matrix = T
source("./a5/wts/scripts/differential_expression/wts_de_toptable_correlation_matrix.r")
mki67_cor <- c(cor_mat_rho$TERT_All_vs_NonTERT["ENSG00000148773.14_MKI67",],
               cor_mat_rho$ATRX_All_vs_NonATRX["ENSG00000148773.14_MKI67",],
               cor_mat_rho$Metastasis_All_vs_NonMetPri_WT["ENSG00000148773.14_MKI67",])
mki67_cor <- mki67_cor[!duplicated(names(mki67_cor))]

#####################
# Genes of Interest #
#####################

gene_select <- "topbottom_logfc" #"top_abslogfc"
adj_pval_cutoff <- 0.05

GOI_lfc_pipe <- . %>% filter(adj.P.Val < adj_pval_cutoff) %>% 
  arrange(desc(abs(logFC))) %>% 
  mutate(Symbol = gsub("ENSG[0-9]+([.][0-9]+)?_(.+)","\\2", Gene)) %>% 
  dplyr::select(Gene, Symbol, adj.P.Val, logFC) %>% 
  {
    if(gene_select == "top_abslogfc")
    { . }
    else if(gene_select == "topbottom_logfc")
    { arrange(., logFC) }
  } %>% 
  mutate(mki67_correlation=mki67_cor[Gene])

GOI.met <- 
  wts_top_tables[["genosampletype"]][["Metastasis_All_vs_NonMetPri_WT"]] %>% 
  GOI_lfc_pipe %>%  
  mutate(source ="Metastasis_All_vs_NonMetPri_WT")


GOI.tert <- wts_top_tables[["genosampletype"]][["TERT_All_vs_NonTERT"]] %>% 
  GOI_lfc_pipe %>% 
  mutate(source="TERT_All_vs_NonTERT")


GOI.atrx <-
  wts_top_tables[["genosampletype"]][["ATRX_All_vs_NonATRX"]] %>% 
  GOI_lfc_pipe %>%  
  mutate(source="ATRX_All_vs_NonATRX")


GOI <- bind_rows(GOI.met, GOI.tert, GOI.atrx)

#GOI <- GOI %>% group_by(Gene) %>% mutate(source = ifelse(n() == 3, "Common", source)) 

GOI <- GOI %>% 
  mutate(source = factor(as.character(source), 
                         levels = c("Metastasis_All_vs_NonMetPri_WT", 
                         "TERT_All_vs_NonTERT", "ATRX_All_vs_NonATRX")),
         source=forcats::fct_recode(source, 
                                    "All Metastasis vs NMP"="Metastasis_All_vs_NonMetPri_WT", 
                                    "TERT vs Rest"="TERT_All_vs_NonTERT", 
                                    "ATRX vs Rest"="ATRX_All_vs_NonATRX")) 

goterms <- ensgid_to_goterm_from_biomart(ens_gids = gsub("[.][0-9]{1,2}_.+$","",GOI$Gene), 
                                         permitted_domains = "biological_process")
goterms_test <- goterms

#Count number of evidence codes for each go_term and keep the highest (ties will be kept)
goterms_test <- goterms_test %>% 
  group_by(ensembl_gene_id, go_id, name_1006, definition_1006, go_linkage_type, namespace_1003) %>% 
  mutate(n_evidence = n()) %>% 
  slice_max(n=1, order_by = n_evidence)

  #mutate(n_evidence = ifelse(parent_go_term %in% c("Other", "None"), 0, n()))

#Annotate GO terms with any matching parent terms
goterms_test <- 
  goterms_test %>% 
  left_join(go_parent_child %>%  
              dplyr::select(child_go_id, parent_go_id, parent_go_term), 
            by = c("go_id" = "child_go_id")) %>% 
  mutate(parent_go_id = ifelse(go_id %in% go_parents$go_id, go_id, parent_go_id),
         parent_go_term = ifelse(go_id %in% go_parents$go_id, name_1006, parent_go_term))

#Add parent ranks and fill in missing values
goterms_test <-  
  goterms_test %>% 
  left_join(go_parents %>% 
              dplyr::select(go_id, parent_go_rank=rank), 
            by = c("parent_go_id"="go_id")) %>% 
  mutate(parent_go_term = case_when( 
    is.na(go_id) ~ "None",
    !is.na(go_id) & is.na(parent_go_id) ~ "Other",
    TRUE ~ parent_go_term),
    parent_go_rank = case_when( 
      is.na(go_id) ~ 99,
      !is.na(go_id) & is.na(parent_go_id) ~ 100,
      TRUE ~ as.numeric(parent_go_rank)))


#Sort by parent GO terms by assigned rank and keep the top
goterms_test <- goterms_test %>% 
  ungroup() %>% 
dplyr::select(ensembl_gene_id, parent_go_id, parent_go_term, parent_go_rank) %>% 
  distinct() %>% 
  group_by(ensembl_gene_id) %>% 
  arrange(parent_go_rank) %>% 
  slice_head(n=1)

#Recode some GO terms to reduce categories
goterms_test <- goterms_test %>% 
  ungroup() %>% 
  mutate(parent_go_term = dplyr::recode(as.character(parent_go_term),
                                        "cell cycle" = "cell cycle/proliferation",
                                        "cell population proliferation" = "cell cycle/proliferation",
                                        "nervous system development" = "nervous system development/neurogenesis",
                                        "neurogenesis" = "nervous system development/neurogenesis",
                                        "chromosome organization" = "chromatin/chromosome organization",
                                        "chromatin organization" ="chromatin/chromosome organization",
                                        "vesicle-mediated transport" = "Other",
                                        "histone modification" = "Other"
  ),
  parent_go_term = stringr:::str_to_sentence(parent_go_term),
  parent_go_term = gsub("[Dd]na", "DNA", parent_go_term),
  parent_go_term = gsub("[Rr]na", "RNA", parent_go_term)) %>% 
  mutate(parent_go_term = factor(as.character(parent_go_term), levels=sort(unique(.$parent_go_term)))) %>% 
  distinct() %>% arrange(desc(parent_go_rank))


#Annotate genes with go info
plot_data <- GOI  %>% 
  mutate(ensembl_gene_id=gsub("[.][0-9]{1,2}_.+$","",Gene)) %>%  left_join(goterms_test) %>% 
  group_by(Gene, Symbol, mki67_correlation, source, ensembl_gene_id, parent_go_id, parent_go_term, parent_go_rank) %>% 
  summarize(adj.P.Val= median(adj.P.Val), logFC=median(logFC))


go_levels <- 
  c("Angiogenesis",
    "Cell adhesion",
    "Cell cycle/proliferation",
    "Cell migration",
    "Cellular developmental process",
    "Cellular response to DNA damage stimulus",
    "Chromatin/chromosome organization",
    "DNA repair",
    "DNA replication",
    "Nervous system development/neurogenesis  ",
    "Programmed cell death",
    "RNA biosynthetic process",
    "Other",
    "None")

plot_data <- plot_data %>%  
  mutate(parent_go_term = factor(as.character(parent_go_term), levels=go_levels)) %>%  
  arrange(source, desc(parent_go_term))

c26 = c("dodgerblue2",
        "green4", 
        "#E31A1C", #red
        "#6A3D9A", # purple
        "#FF7F00", # orange
        "gold1",
        "skyblue2",
        "#FB9A99", # lt pink
        "palegreen2",
        "#CAB2D6", # lt purple
        "#FDBF6F", # lt orange
        "khaki2",
        "grey40",
        "grey70",
        "maroon",
        "orchid1",
        "deeppink1",
        "blue1",
        "steelblue4",
        "darkturquoise",
        "green1",
        "yellow4",
        "yellow3",
        "darkorange4",
        "brown",
        "black",
        "aquamarine1",
        "darkslategrey")

col_scale <- setNames(c26[1:length(go_levels)], go_levels)

ggplot(plot_data, aes(x= source, y=mki67_correlation, color=parent_go_term)) + geom_jitter(width=0.2, size=1.2) + 
  scale_color_manual(values = col_scale) + 
  guides(color=guide_legend(override.aes = list(size=2))) +
  theme_bw()
