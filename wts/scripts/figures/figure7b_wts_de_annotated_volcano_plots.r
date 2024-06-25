library(tidyverse)
library(patchwork)
library(GO.db)
library(DBI)

setwd("/g/data/pq08/projects/ppgl/")

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
    13, "GO:0048869","cellular developmental process",
    14, "GO:0022008","neurogenesis",
    15, "GO:0001525", "angiogenesis",
    16, "GO:0007399","nervous system development",
    17, "GO:0016192", "vesicle-mediated transport",
    #18, "GO:0032774", "RNA biosynthetic process"
    18, "GO:0006355", "Regulation of DNA-templated transcription"
    ),
  byrow = T, ncol=3)

dimnames(go_parents) = list(c(go_parents[,1]), c("rank", "go_id", "go_term"))

go_parents <- data.frame(go_parents)

go_parent_child = get_go_children(go_parents$go_id)

#############
# Gene Sets #
#############

Hs.H <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
Hs.C8 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C8")

# put the genesets into list form
Hs.H.list <- lapply(split(Hs.H, f = Hs.H$gs_name), function(w) { unique(w$gene_symbol) })

#Janskey

signature_genes_jansky <- read_csv("./a5/sn_rna_seq/reference/published_datasets/jansky_2021/jansky_sympathoadrenal_genes.csv", skip = 1)

jansky.genesets.list <- list(
  JANSKY_late_chromaffin_genes = signature_genes_jansky %>% dplyr::filter(cluster == "late Chromaffin cells") %>% pull(gene),
  JANSKY_early_chromaffin_genes = signature_genes_jansky %>% dplyr::filter(cluster == "Chromaffin cells") %>% pull(gene),
  JANSKY_connecting_chromaffin_genes = signature_genes_jansky %>% dplyr::filter(cluster == "connecting Progenitor cells")  %>% pull(gene),
  JANSKY_bridge_genes = signature_genes_jansky %>% dplyr::filter(cluster == "Bridge") %>% pull(gene),
  JANSKY_early_neuroblast_genes = signature_genes_jansky %>% dplyr::filter(cluster == "Neuroblasts") %>% pull(gene),
  JANSKY_late_neuroblast_genes = signature_genes_jansky %>% dplyr::filter(cluster == "late Neuroblasts") %>% pull(gene),
  JANSKY_cycling_SCP_genes = signature_genes_jansky %>% dplyr::filter(cluster == "cycling SCPs") %>% pull(gene),
  JANSKY_cycling_neuroblast_genes = signature_genes_jansky %>% dplyr::filter(cluster == "cycling Neuroblasts") %>% pull(gene),
  JANSKY_early_SCP_genes = signature_genes_jansky %>% dplyr::filter(cluster == "SCPs") %>% pull(gene),
  JANSKY_late_SCP_genes = signature_genes_jansky %>% dplyr::filter(cluster == "late SCPs") %>% pull(gene))

zethovan_pseudobulk_de <- read_csv("./a5/sn_rna_seq/reference/published_datasets/zethoven_2022/zethoven_suppdata_6.csv")

# filter out genes below significance threshold cutoff of 0.05 adj p value
# filter out genes with log fold change below 3 
zethovan_pseudobulk_de <- zethovan_pseudobulk_de %>% 
  dplyr::filter(grepl("_normal" , Contrast)) %>%
  mutate(Contrast = Contrast %>% str_replace_all("\\.", "_") %>%
           str_replace("Sustentacular_cells", "SCLCs") %>% 
           str_remove("_normal")) %>% 
  dplyr::filter(adj.P.Val < 0.05) %>% 
  dplyr::filter(logFC > 3)

# chose logFC cutoff to have a good number of genes in each set (200-400 genes)
# this can be tweaked further if necessary
zethoven.celltypes.list <- lapply(split(zethovan_pseudobulk_de,
                                        f = zethovan_pseudobulk_de$Contrast),
                                  function(y) y$Gene)
names(zethoven.celltypes.list) <- paste0("ZETHOVEN_", names(zethoven.celltypes.list))

all_genesets <- c(Hs.H.list, jansky.genesets.list, zethoven.celltypes.list)

gene_geneset_lookup <- data.frame(GeneSet=rep(names(all_genesets), times= purrr::map(all_genesets, length)), gene_symbol=purrr::reduce(all_genesets, c))
gene_geneset_lookup <- gene_geneset_lookup[!duplicated(gene_geneset_lookup$gene_symbol),]

genesets_to_use <- c(
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_MITOTIC_SPINDLE",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "JANSKY_cycling_neuroblast",
  "ZETHOVEN_SCLCs",
  "JANSKY_late_SCPs",
  "JANSKY_early_SCPs"
)

gene_geneset_lookup <- gene_geneset_lookup %>% filter(GeneSet %in% genesets_to_use)

################
# Data Loaders #
################

if (!exists("wts_top_tables")) {
output_qc <- FALSE; output_tables <- FALSE; output_plots <- FALSE
source("./a5/wts/scripts/differential_expression/a5_wts_differential_expression.r")
}

############
# GO Terms #
############

GOI <- purrr::map(.x = c(wts_top_tables$Chromaffin_vs_Non_chromaffin, wts_top_tables$genosampletype), .f = \(tt) { tt %>% filter(adj.P.Val < 0.05) %>% pull(ensembl_gene_id) }) %>% 
  purrr::reduce(.f = c) %>% unique()

goterms <- ensgid_to_goterm_from_biomart(ens_gids = gsub("[.][0-9]{1,2}_.+$","",GOI), 
                                         permitted_domains = "biological_process")


#Count number of evidence codes for each go_term and keep the highest (ties will be kept)
goterms <- goterms %>% 
  group_by(ensembl_gene_id, go_id, name_1006, definition_1006, go_linkage_type, namespace_1003) %>% 
  mutate(n_evidence = n()) %>% 
  slice_max(n=1, order_by = n_evidence)

#Annotate GO terms with any matching parent terms
goterms <- 
  goterms %>% 
  left_join(go_parent_child %>%  
              dplyr::select(child_go_id, parent_go_id, parent_go_term), 
            by = c("go_id" = "child_go_id")) %>% 
  mutate(parent_go_id = ifelse(go_id %in% go_parents$go_id, go_id, parent_go_id),
         parent_go_term = ifelse(go_id %in% go_parents$go_id, name_1006, parent_go_term))

#Add parent ranks and fill in missing values
goterms <-  
  goterms %>% 
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
goterms <- goterms %>% 
  ungroup() %>% 
  dplyr::select(ensembl_gene_id, parent_go_id, parent_go_term, parent_go_rank) %>% 
  distinct() %>% 
  group_by(ensembl_gene_id) %>% 
  arrange(parent_go_rank) %>% 
  slice_head(n=1)

#Recode some GO terms to reduce categories
goterms <- goterms %>% 
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
    "Nervous system development/neurogenesis",
    "Programmed cell death",
    #"RNA biosynthetic process",
    "Regulation of DNA-templated transcription",
    "Other",
    "None")


############
# Plotting #
############



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



plot_volcano <- function(tt, n_label=30, genes_to_label=NULL, lfc_cutoff=1.5, padj_cutoff=0.05, color_feature, color_scale){
  
  if(is.null(genes_to_label)) {
  genes_to_label <- tt %>% filter(adj.P.Val < padj_cutoff & abs(logFC) > lfc_cutoff) %>% 
    slice_min(adj.P.Val, n = n_label) %>% 
    pull(Gene)
  } else {
    genes_to_label <- tt %>%  filter(gene_symbol %in% genes_to_label) %>%  pull(Gene)
  }
  
  if(!is.factor(tt[[color_feature]])) { stop("Color feature columns must be a factor")}
  
  tt_plot <- tt %>% #flag the most significant genes to label
    mutate(color_feature = if_else(adj.P.Val < padj_cutoff & abs(logFC) > lfc_cutoff , !!sym(color_feature), "Not significant"),
           color_feature = factor(as.character(color_feature), levels=levels(!!sym(color_feature)))) %>% 
    mutate(label = if_else(Gene %in% genes_to_label, TRUE, FALSE)) 
  
  volcano <- ggplot(tt_plot, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(col=color_feature)) +
    scale_color_manual(values=color_scale) +
    geom_text_repel(
      data = . %>% filter(label == TRUE),
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.1, "lines"),
      #segment.color = "grey",
      max.overlaps = Inf,
      size = 4,
      nudge_y = 0.2,
      aes(label=gene_symbol)) +
    geom_hline(yintercept=-log10(padj_cutoff), linetype=2,color="red")+
    geom_vline(xintercept=c(lfc_cutoff, -1 * lfc_cutoff), linetype=2,color="red")+
    theme_bw() + theme(aspect.ratio = 1)
  
  return(volcano)
  
}


if(!exists("lfc_cutoff")) {
  lfc_cutoff=1
}

if(!exists("padj_cutoff")) {
  padj_cutoff=0.05
}

gg_volcano_go_benmet <- list()
gg_volcano_geneset_benmet <- list()
genes_to_label <- list(Metastatic_All_vs_NonMetPri_WT = c("TERT","ATRX"),
                       TERT_PriMet_vs_NonMetPri_WT = c("TERT"), 
                       ATRX_PriMet_vs_NonMetPri_WT = c("DRG2", "RPRM"))
for (contrast in c("Metastatic_All_vs_NonMetPri_WT", "TERT_PriMet_vs_NonMetPri_WT", "ATRX_PriMet_vs_NonMetPri_WT"))
{
  
  plot_data_base <- wts_top_tables[["genosampletype"]][[contrast]] %>% 
    left_join(gene_geneset_lookup) %>% 
    left_join(goterms)  %>% 
    mutate(gene_symbol = ifelse(grepl("^A[CL][0-9]", gene_symbol), 
                                stringr::str_extract(string = Gene, pattern = "ENSG[0-9]+"), 
                                gene_symbol))
  
  plot_data <- plot_data_base %>% mutate(parent_go_term = factor(as.character(parent_go_term), levels=c(go_levels, "Not significant")))
  col_set = setNames(c26[1:length(go_levels)], go_levels)
  
  gg_volcano_go_benmet[[contrast]] <- 
    plot_volcano(tt = plot_data,
                 padj_cutoff = padj_cutoff,
                 lfc_cutoff = lfc_cutoff,
                 genes_to_label = genes_to_label[[contrast]],
                 color_feature = "parent_go_term",
                 color_scale = col_set) + ggtitle(contrast)
  
  plot_data <- plot_data_base %>% mutate(GeneSet = replace_na(GeneSet, "None"),
                                         GeneSet = factor(as.character(GeneSet), levels=c(unique(GeneSet), "Not significant")))
  col_set = setNames(c26[1:length(unique(plot_data$GeneSet))], unique(plot_data$GeneSet))
  col_set[["None"]] = "grey50"
  
  gg_volcano_geneset_benmet[[contrast]] <- 
    plot_volcano(tt = plot_data, 
                 padj_cutoff = padj_cutoff,
                 lfc_cutoff = lfc_cutoff,
                 genes_to_label = genes_to_label[contrast],
                 color_feature = "GeneSet",
                 color_scale = col_set) + ggtitle(contrast)
  
}


gg_volcano_go_nonvschromaffin <- list()
gg_volcano_geneset_nonvschromaffin <- list()
for (contrast in c("Chromaffin_vs_Non_chromaffin"))
{
  
  plot_data_base <- wts_top_tables[["Chromaffin_vs_Non_chromaffin"]][[contrast]] %>% 
    left_join(gene_geneset_lookup) %>% 
    left_join(goterms) %>% 
    mutate(gene_symbol = ifelse(grepl("^A[CL][0-9]", gene_symbol), 
                                stringr::str_extract(string = Gene, pattern = "ENSG[0-9]+"), 
                                gene_symbol))
  
  plot_data <- plot_data_base %>% mutate(parent_go_term = factor(as.character(parent_go_term), levels=c(go_levels, "Not significant")))
  col_set = setNames(c26[1:length(go_levels)], go_levels)
  gg_volcano_go_nonvschromaffin[[contrast]] <- 
    plot_volcano(tt = plot_data, 
                 padj_cutoff = padj_cutoff,
                 lfc_cutoff = lfc_cutoff,
                 n_label = 30, 
                 color_feature = "parent_go_term",
                 color_scale = col_set) + ggtitle(contrast)
  
  plot_data <- plot_data_base %>% mutate(GeneSet = replace_na(GeneSet, "None"),
                                         GeneSet = factor(as.character(GeneSet), levels=c(unique(GeneSet), "Not significant")))
  col_set = setNames(c26[1:length(unique(plot_data$GeneSet))], unique(plot_data$GeneSet))
  col_set[["None"]] = "grey50"
  gg_volcano_geneset_nonvschromaffin[[contrast]] <- 
    plot_volcano(tt = plot_data, 
                 n_label = 30, 
                 padj_cutoff = padj_cutoff,
                 lfc_cutoff = lfc_cutoff,
                 color_feature = "GeneSet",
                 color_scale = col_set) + ggtitle(contrast)
  
}


#####################################
# ATRX/TERT/Met_vs_NMP intersection #
#####################################

contrast = "Metastatic_All_vs_NonMetPri_WT"
plot_data_base <- wts_top_tables[["genosampletype"]][[contrast]] %>% 
  left_join(goterms)  %>% 
  mutate(gene_symbol = ifelse(grepl("^A[CL][0-9]", gene_symbol), 
                              stringr::str_extract(string = Gene, pattern = "ENSG[0-9]+"), 
                              gene_symbol))

plot_data <- plot_data_base %>% mutate(parent_go_term = factor(as.character(parent_go_term), levels=c(go_levels, "Not significant")))

GOI <- purrr::map(.x = c(wts_top_tables$genosampletype[c("TERT_PriMet_vs_NonMetPri_WT","ATRX_PriMet_vs_NonMetPri_WT","Metastatic_All_vs_NonMetPri_WT")]), 
                  .f = \(tt) { tt %>% filter(adj.P.Val < padj_cutoff, abs(logFC) > lfc_cutoff) %>% pull(ensembl_gene_id) }) %>% 
  purrr::reduce(.f = intersect) %>% unique()

hgnc_lookup <- wts_top_tables$genosampletype$Metastatic_All_vs_NonMetPri_WT %>% dplyr::select(ensembl_gene_id, gene_symbol)

plot_data <- plot_data %>%  
  mutate(color_field=ifelse(ensembl_gene_id %in% GOI, as.character(parent_go_term), "Not intersecting"),
         color_field=factor(as.character(color_field), levels=c(go_levels, "Not significant", "Not intersecting")))
plot_data <- plot_data %>%  arrange(desc(color_field))

col_set = go_biotype_cols
col_set[["Other"]] <- "black"
col_set[["None"]] <- "black"
col_set[["Not intersecting"]] <- "grey70"
col_set[["Not significant"]] <- "grey50"

cellcycle_genes <- unlist(purrr:::map(as.list(Revelio::revelioTestData_cyclicGenes), .f = \(x) as.character(x[!is.na(x)])))
cellcycle_genes.use <- intersect(cellcycle_genes, hgnc_lookup$gene_symbol[hgnc_lookup$ensembl_gene_id %in% GOI])
gg_volcano_go_benmet[[contrast]] <- 
  plot_volcano(tt = plot_data,
               padj_cutoff = padj_cutoff,
               lfc_cutoff = lfc_cutoff,
               genes_to_label = c("CDK1", "BUB1", "TOP2A", "AURKB", "EZH2", "MKI67", "FOXM1",
                                  "LHX9","CDH19","PENK","FOXD3","DCHS2","SSTR1","FZD9"),
               color_feature = "color_field",
               color_scale = col_set) + ggtitle(contrast)
