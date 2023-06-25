library(biomaRt)
library(topGO)
library(org.Hs.eg.db)

#############
# IMPORTANT #
#############

# This script requires the 'A5_gridss' and 'A5_vcfs.keep' objects loaded by "load_all_data.R" to be in the environment

##############
# /IMPORTANT #
##############

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))#,host = "http://asia.ensembl.org")) 


#Feed gene list to biomart requesting a table with ensembl_gene_id and external_gene_name (aka. gene symbol)
hgnc2ens <- getBM(filters= "external_gene_name", 
                  attributes= c("ensembl_gene_id", "external_gene_name", "chromosome_name"), 
                  values=unique(c(A5_gridss$GeneStartName, A5_gridss$GeneEndName)),
                  mart= mart)
hgnc2ens <- hgnc2ens %>%  filter(chromosome_name %in% c(1:22,"X","Y")) %>%  dplyr::select(-chromosome_name)

hgnc2ens <- bind_rows(hgnc2ens,
  A5_vcfs.keep %>%  ungroup() %>%  filter(Normal_AF==0, Tumour_AF > 0.15) %>% dplyr::select(Gene, PCGR_SYMBOL) %>% dplyr::rename(ensembl_gene_id=Gene, external_gene_name=PCGR_SYMBOL)) %>%  distinct()

GoTerms <- getBM(filters= "ensembl_gene_id",
                 attributes= c("ensembl_gene_id",
                               "go_id",
                               "name_1006",
                               "definition_1006",
                               "go_linkage_type",
                               "namespace_1003",
                               "external_gene_name"),
                 values=hgnc2ens$ensembl_gene_id,mart= mart)


#######################
# sample_gene_mapping #
#######################


sample_gene_mapping <- 
  bind_rows(A5_vcfs.keep %>%  ungroup() %>% dplyr::select(A5_ID, Gene) %>% dplyr::rename(ensembl_gene_id=Gene)
            ,
            A5_gridss %>% mutate(A5_ID=gsub("T0(.)","\\1",A5_ID)) %>% dplyr::select(A5_ID, GeneStartName, GeneEndName) %>% 
              pivot_longer(cols=c(GeneStartName, GeneEndName), names_to = "discard", values_to="external_gene_name") %>% 
              dplyr::select(-discard) %>%  distinct() %>% 
              inner_join(hgnc2ens) %>%  dplyr::select(-external_gene_name)) %>% 
  left_join(Anno %>%  dplyr::select(A5_ID, is_primary_or_met, tumour_metastasised))

sample_gene_mapping <- sample_gene_mapping %>%  filter(A5_ID != "E167-1")

genes_in_benign <- sample_gene_mapping %>% filter(is_primary_or_met=="Primary", tumour_metastasised!="Yes") %>% pull(ensembl_gene_id) %>%  unique()
genes_in_metastatic_only <- sample_gene_mapping %>% filter(!(ensembl_gene_id %in% genes_in_benign), tumour_metastasised=="Yes") %>% pull(ensembl_gene_id) %>%  unique()

#############
# All Genes #
#############

AllGenes <- select(org.Hs.eg.db,columns = c("ENSEMBL","SYMBOL","REFSEQ"),
                   keys=keys(org.Hs.eg.db)) %>%  
  filter(!is.na(ENSEMBL),!is.na(SYMBOL),!is.na(REFSEQ)) %>% 
  pull(ENSEMBL) %>% unique()

geneID2GO <- getBM(filters= "ensembl_gene_id", 
                   attributes= c("ensembl_gene_id", "external_gene_name", "chromosome_name", "go_id"), 
                   values=AllGenes,
                   mart= mart)

geneID2GO <- geneID2GO %>%  
  filter(go_id != "",
         ensembl_gene_id != "", 
         external_gene_name != "", 
         chromosome_name %in% c(1:22,"X","Y")) %>% 
  dplyr::select(-external_gene_name, -chromosome_name)

geneID2GO <- split(x = geneID2GO, f = geneID2GO$ensembl_gene_id)
geneID2GO <- lapply(geneID2GO, function (golist) { golist$go_id })

geneList <- factor(as.integer(AllGenes %in% hgnc2ens$ensembl_gene_id), levels=c(0,1))
names(geneList) <- AllGenes


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO,  nodeSize = 20)

GOdata <- updateGenes(object = GOdata, geneList = geneList[GOdata@feasible])

#GOdata@allScores <- as.numeric(GOdata@allGenes %in% metastatic_only)

#GOdata@geneSelectionFun <- function(x) { return (x==1) }

ts <- termStat(GOdata)

ts <- ts %>%  tibble::rownames_to_column(var = "GO_ID") %>% mutate(ratio=Significant/Expected) %>%  arrange(desc(ratio))

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                    orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 1000)

allRes <- allRes %>% filter(Annotated < 50)
allRes$classicFisher_FDR <- p.adjust(allRes$classicFisher, method = p.adjust.methods, n = length(allRes$classicFisher))

allRes.sig <- allRes %>%  filter(classicFisher_FDR <0.05)

########################################
# Met Only Genes - All Gene Background #
########################################

geneList <- factor(as.integer(AllGenes %in% genes_in_metastatic_only), levels=c(0,1))
names(geneList) <- AllGenes


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO,  nodeSize = 20)

GOdata <- updateGenes(object = GOdata, geneList = geneList[GOdata@feasible])

#GOdata@allScores <- as.numeric(GOdata@allGenes %in% metastatic_only)

#GOdata@geneSelectionFun <- function(x) { return (x==1) }

ts <- termStat(GOdata)

ts <- ts %>%  tibble::rownames_to_column(var = "GO_ID") %>% mutate(ratio=Significant/Expected) %>%  arrange(desc(ratio))

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 1000)

allRes <- allRes %>% filter(Annotated < 50)
allRes$classicFisher_FDR <- p.adjust(allRes$classicFisher, method = p.adjust.methods, n = length(allRes$classicFisher))

allRes.sig <- allRes %>%  filter(classicFisher_FDR <0.05)


#############################################
# Met Only Genes - Mutated Genes Background #
#############################################

geneList <- factor(as.integer(hgnc2ens$ensembl_gene_id %in% genes_in_metastatic_only), levels=c(0,1))
names(geneList) <- hgnc2ens$ensembl_gene_id


GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = geneID2GO,  nodeSize = 20)

GOdata <- updateGenes(object = GOdata, geneList = geneList[GOdata@feasible])

#GOdata@allScores <- as.numeric(GOdata@allGenes %in% metastatic_only)

#GOdata@geneSelectionFun <- function(x) { return (x==1) }

ts <- termStat(GOdata)

ts <- ts %>%  tibble::rownames_to_column(var = "GO_ID") %>% mutate(ratio=Significant/Expected) %>%  arrange(desc(ratio))

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher,
                   orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 1000)

allRes <- allRes %>% filter(Annotated < 50)
allRes$classicFisher_FDR <- p.adjust(allRes$classicFisher, method = p.adjust.methods, n = length(allRes$classicFisher))

allRes.sig <- allRes %>%  filter(classicFisher_FDR <0.05)




#################
# Curated Genes #
#################


c(A5_vcfs.keep.curated$Gene), "ENSG00000054267"

setdiff(shiva_curated_genes.sv, A5_vcfs.keep.curated$SYMBOL)

sample_gene_mapping <- 
  bind_rows(A5_vcfs.keep.curated %>%  ungroup() %>% dplyr::select(A5_ID, PCGR_SYMBOL)
            ,
            A5_gridss.curated %>% mutate(A5_ID=gsub("T0(.)","\\1",A5_ID)) %>% dplyr::select(A5_ID, GeneStartName, GeneEndName) %>% 
              pivot_longer(cols=c(GeneStartName, GeneEndName), names_to = "discard", values_to="PCGR_SYMBOL") %>% 
              dplyr::select(-discard) %>% 
              filter(PCGR_SYMBOL %in% shiva_curated_genes.sv) %>%  distinct()) %>% 
  left_join(Anno %>%  dplyr::select(A5_ID, is_primary_or_met))
  
sample_gene_mapping <- sample_gene_mapping %>%  filter(A5_ID != "E167-1")

sample_gene_mapping.summary <- sample_gene_mapping %>% group_by(PCGR_SYMBOL) %>% 
  summarise(PriMet=paste(unique(is_primary_or_met), collapse =","), 
            A5_ID=paste(unique(A5_ID), collapse =","))



mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = "http://asia.ensembl.org")) 

#Feed gene list to biomart requesting a table with ensembl_gene_id and external_gene_name (aka. gene symbol)
hgnc2ens <- getBM(filters= "external_gene_name", attributes= c("ensembl_gene_id", "external_gene_name", "chromosome_name"), values=c(shiva_curated_genes.variant, shiva_curated_genes.sv),mart= mart)
hgnc2ens <- hgnc2ens %>%  filter(chromosome_name %in% c(1:22,"X","Y"))

GoTerms <- getBM(filters= "ensembl_gene_id",
                                 attributes= c("ensembl_gene_id",
                                               "go_id",
                                               "name_1006",
                                               "definition_1006",
                                               "go_linkage_type",
                                               "namespace_1003",
                                               "external_gene_name"),
                                 values=hgnc2ens$ensembl_gene_id,mart= mart)


View(GoTerms %>% filter(namespace_1003 != "cellular_component", go_linkage_type != "IEA"))

GoTerms %>% 
  filter(namespace_1003 != "cellular_component", go_linkage_type != "IEA") %>% 
  left_join(sample_gene_mapping.summary, by=c("external_gene_name"="PCGR_SYMBOL")) %>% 
group_by(name_1006,definition_1006, go_linkage_type, namespace_1003) %>% 
  summarise(Genes=paste(unique(external_gene_name), collapse =","), 
            count=length(unique(external_gene_name)),
            PriMet=paste(unique(PriMet), collapse =","),
            A5_ID=paste(unique(A5_ID), collapse =",")) %>% 
  filter(name_1006 != "") %>% 
  write.table("C:/ResearchData/RADIO/A5/secondary_analysis/mutated_genes/Shiva_curated_GOTerms_by_term.tsv", sep="\t", row.names = F)

