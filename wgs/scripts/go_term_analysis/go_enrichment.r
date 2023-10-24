library(clusterProfiler)
setwd("/g/data/pq08/projects/ppgl")

################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/wts/scripts/data_loaders/a5_wts_dataloader.r")

if (!exists("a5_somatic_variants_keep")) {
  source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
  data_loader_somatic_variants(quickload = T)
}


gene_universe <- AnnotationDbi::select(org.Hs.eg.db, columns = "ENSEMBL", keys=keys(org.Hs.eg.db)) 
gene_universe <- gene_universe %>%  filter(!is.na(ENSEMBL)) %>% dplyr::rename("ensembl_gene_id"="ENSEMBL")
data_loader_get_biotypes(gene_universe$ensembl_gene_id, use_cache = T)
gene_universe <- ensid_to_biotype %>%  filter(gene_biotype =="protein_coding") %>%  pull(ensembl_gene_id)

if (!exists("quickload_go_enrichment")) { quickload_go_enrichment <- T }
message("quickload_go_enrichment:", quickload_go_enrichment)

checkpointfile_goenrichment <- "./a5/wgs/quickload_checkpoints/goenrichment.rds"
if(quickload_go_enrichment) {
  ego <- readRDS(checkpointfile_goenrichment)
} else {
  ego <- list()
  
  all_genes <- unique(a5_somatic_variants_keep %>% filter(A5_ID != "E167-1") %>% pull(Gene))
  
  ego[["All"]] <- list()
  for (Ont in c("CC","MF","BP"))
  {
    ego[["All"]][[Ont]] <- enrichGO(gene = all_genes,
                           universe      = gene_universe,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = 'ENSEMBL',
                           ont           = Ont,
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.05,
                           qvalueCutoff  = 0.05,
                           readable      = TRUE)
    ego[["All"]][[Ont]] <- data.frame(ego[["All"]][[Ont]]) %>% mutate(source="All", ontology=Ont)
    
  }
  
  met_genes <- a5_somatic_variants_keep %>%  
    inner_join(a5_anno %>%  
                 dplyr::select(A5_ID, differential_group_sampletype_strict)) %>% 
    filter(A5_ID != "E167-1", 
           differential_group_sampletype_strict %in% c("Metastasis", "Metastatic local recurrence", "Metastatic primary")) %>% 
    pull(Gene) %>%  unique()
  
  ego[["Met"]] <- list()
  for (Ont in c("CC","MF","BP"))
  {
    ego[["Met"]][[Ont]] <- enrichGO(gene = met_genes,
                                    universe      = gene_universe,
                                    OrgDb         = org.Hs.eg.db,
                                    keyType       = 'ENSEMBL',
                                    ont           = Ont,
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.05,
                                    readable      = TRUE)
    
    ego[["Met"]][[Ont]] <- data.frame(ego[["Met"]][[Ont]]) %>% mutate(source="Met", ontology=Ont)
    
  }
  
  nonmetpri_genes <- a5_somatic_variants_keep %>%  
    inner_join(a5_anno %>%  
                 dplyr::select(A5_ID, differential_group_sampletype_strict)) %>% 
    filter(A5_ID != "E167-1", 
           differential_group_sampletype_strict %in% c("Non-metastatic primary")) %>% 
    pull(Gene) %>%  unique()
  
  ego[["Nonmetpri"]] <- list()
  for (Ont in c("CC","MF","BP"))
  {
    ego[["Nonmetpri"]][[Ont]] <- enrichGO(gene = nonmetpri_genes,
                                    universe      = gene_universe,
                                    OrgDb         = org.Hs.eg.db,
                                    keyType       = 'ENSEMBL',
                                    ont           = Ont,
                                    pAdjustMethod = "BH",
                                    pvalueCutoff  = 0.05,
                                    qvalueCutoff  = 0.05,
                                    readable      = TRUE)
    ego[["Nonmetpri"]][[Ont]] <- data.frame(ego[["Nonmetpri"]][[Ont]]) %>% mutate(source="Nonmetpri", ontology=Ont)
    
  }
  
  saveRDS(ego, checkpointfile_goenrichment)
  
}


ego <- bind_rows(bind_rows(ego$All),
bind_rows(ego$Met),
bind_rows(ego$Nonmetpri))



