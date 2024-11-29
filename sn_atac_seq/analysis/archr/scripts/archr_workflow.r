#######################
# Script to run motif #
# enrichment analysis #
# on snATAC using     # 
# ArchR               #
#                     #
# Author: Aidan Flynn #
# Date: 20/10/2024    #
#                     #
#######################

#########################################################
# Note: This script requires a large amount of          #
# RAM to load both snATAC and snRNA data simultaneously #
#########################################################

setwd("/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/archr")

################
# Dependencies #
################

library(ArchR)
library(tidyr)
library(dplyr)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)


#########
# Setup #
#########

set.seed(1)

addArchRThreads(threads = 6) 

system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))

addArchRGenome("hg38")

################
# Data loaders #
################

#Clinical data
source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(use_cache = T)

#snRNA data loader
source("/g/data/pq08/projects/ppgl/a5/sn_rna_seq/scripts/data_loaders/a5_snrna_dataloader.r")
data_loader_a5_snrna(quickload = T)
a5_snrna <- snrna_annotate_cell_types(snrna_object = a5_snrna, output_qc = F, 
                                      qc_out_dir = "/g/data/pq08/projects/ppgl/a5/sn_rna_seq/qc")

###########################
# Convert counts to Arrow #
###########################

input_dirs <- list.dirs("/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/cellranger/counts_hg38", 
                        full.names = T, 
                        recursive = F)
input_dirs <- input_dirs[grepl("E[0-9]{3}|NAM", input_dirs, perl = T)]
sample_names <- basename(input_dirs)

input_files <- file.path(input_dirs, "outs/fragments.tsv.gz")
names(input_files) <- sample_names

ArrowFiles <- createArrowFiles(
  inputFiles = input_files,
  sampleNames = names(input_files),
  outputNames = file.path("./raw_arrow_files",names(input_files)),
  minTSS = 4, 
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

###################
# Filter Doublets #
###################

#Custom modification to use amulet output in ArchR doublet filtering
filterDoublets_amulet <- function (ArchRProj = NULL, 
                                   cutEnrich = 1, 
                                   cutScore = -Inf, 
                                   filterRatio = 1, 
                                   amulet_doublets=NULL) 
{
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = cutEnrich, name = "cutEnrich", valid = c("numeric"))
  ArchR:::.validInput(input = cutScore, name = "cutScore", valid = c("numeric"))
  ArchR:::.validInput(input = filterRatio, name = "filterRatio", valid = c("numeric"))
  if (any(grepl("filterDoublets", names(ArchRProj@projectSummary)))) {
    stop("Already ran filterDoublets on ArchRProject! Cannot be re-ran on an ArchRProject!")
  }
  df <- ArchR::getCellColData(ArchRProj, c("Sample", "DoubletEnrichment", 
                                           "DoubletScore"))
  
  df$AmuletDoublet<- FALSE
  if (!is.null(amulet_doublets)) {
    df$AmuletDoublet[rownames(df) %in% amulet_doublets] <- TRUE
  }
  
  splitDF <- split(seq_len(nrow(df)), as.character(df$Sample))
  cellsFilter <- lapply(splitDF, function(y) {
    x <- df[y, , drop = FALSE]
    n <- nrow(x)
    x <- x[order(x$DoubletEnrichment, decreasing = TRUE),]
    if (!is.null(cutEnrich)) {
      x <- x[which(x$DoubletEnrichment >= cutEnrich | x$AmuletDoublet), ]
    }
    if (!is.null(cutScore)) {
      x <- x[which(x$DoubletScore >= cutScore | x$AmuletDoublet), ]
    }
    if (nrow(x) > 0) {
      head(rownames(x), filterRatio * n * (n/1e+05))
    }
    else {
      NULL
    }
    
  }) %>% unlist(use.names = FALSE)
  message("Filtering ", length(cellsFilter), " cells from ArchRProject!")
  tabRemove <- table(df[cellsFilter, ]$Sample)
  tabAll <- table(df$Sample)
  samples <- unique(df$Sample)
  for (i in seq_along(samples)) {
    if (!is.na(tabRemove[samples[i]])) {
      message("\t", samples[i], " : ", tabRemove[samples[i]], 
              " of ", tabAll[samples[i]], " (", round(100 * 
                                                        tabRemove[samples[i]]/tabAll[samples[i]], 1), 
              "%)")
    }
    else {
      message("\t", samples[i], " : ", 0, " of ", tabAll[samples[i]], 
              " (0%)")
    }
  }
  if (length(cellsFilter) > 0) {
    ArchRProj@cellColData <- ArchRProj@cellColData[rownames(ArchRProj@cellColData) %ni% 
                                                     cellsFilter, , drop = FALSE]
  }
  ArchRProj <- addProjectSummary(ArchRProj = ArchRProj, name = "filterDoublets", 
                                 summary = c(cutEnrich = cutEnrich, cutScore = cutScore, 
                                             filterRatio = filterRatio))
  ArchRProj
}

addArchRThreads(threads = 1) 
ArrowFiles <- list.files("./raw_arrow_files", full.names = T, pattern =  "*.arrow", recursive = F)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "project",
  copyArrows = TRUE
)

proj <- addDoubletScores(
  input = proj,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  #knnMethod = "LSI",
  #force = TRUE
)

#Load amulet output
amulet_output <- list.files("../amulet", 
                            recursive = T, 
                            pattern = "MultipletCellIds_01.txt", 
                            full.names = T)
names(amulet_output) <- stringr::str_extract(pattern = "(E|NAM)[0-9]{3}(-.)?", 
                                             string = amulet_output)

#read doublet lists and conform naming to ArchR style
amulet_doublets <- purrr::map2(.x = amulet_output, .y=names(amulet_output), 
                               .f = \(x,y) { paste(y, readr::read_lines(x), sep='#') })
amulet_doublets <- purrr::reduce(.f = c, .x = amulet_doublets)

#Tag doublets
proj <- filterDoublets_amulet(ArchRProj = proj, 
                              amulet_doublets = amulet_doublets)

saveArchRProject(proj, load=F)


################################
# scRNA cell identity transfer #
################################


#Unsupervised snRNA integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = a5_snrna,
  addToArrow = FALSE,
  groupRNA = "cell_type",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Un))
rnaClust <- colnames(cM)[apply(cM, 1 , which.max)]
lift_over <- data.frame(rnaClust, atacClust=rownames(cM)) #Assignments

#Supervised groups
groupList <- SimpleList(
  Tumor = SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% lift_over$atacClust[lift_over$rnaClust=="Tumor"]],
    RNA = rownames(a5_snrna@meta.data)[a5_snrna@meta.data$cell_type=="Tumor"]
  ),
  Chromaffin = SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% lift_over$atacClust[lift_over$rnaClust=="Chromaffin cells"]],
    RNA = rownames(a5_snrna@meta.data)[a5_snrna@meta.data$cell_type=="Chromaffin cells"]
  ),
  Adrenocortical = SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% lift_over$atacClust[lift_over$rnaClust=="Adrenocortical cells"]],
    RNA = rownames(a5_snrna@meta.data)[a5_snrna@meta.data$cell_type=="Adrenocortical cells"]
  ),
  Lymphocytes = SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% lift_over$atacClust[lift_over$rnaClust=="Lymphocytes"]],
    RNA = rownames(a5_snrna@meta.data)[a5_snrna@meta.data$cell_type=="Lymphocytes"]
  ),
  Myeloid = SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% lift_over$atacClust[lift_over$rnaClust=="Myeloid cells"]],
    RNA = rownames(a5_snrna@meta.data)[a5_snrna@meta.data$cell_type=="Myeloid cells"]
  ),
  Endothelial = SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% lift_over$atacClust[lift_over$rnaClust=="Endothelial cells"]],
    RNA = rownames(a5_snrna@meta.data)[a5_snrna@meta.data$cell_type=="Endothelial cells"]
  ),
  Fibro_SCLC = SimpleList(
    ATAC = proj$cellNames[proj$Clusters %in% lift_over$atacClust[lift_over$rnaClust %in% c("Fibroblasts","SCLCs")]],
    RNA = rownames(a5_snrna@meta.data)[a5_snrna@meta.data$cell_type %in% c("Fibroblasts","SCLCs")]
  )  
)

addArchRThreads(threads = 6) 

#Supervised snRNA integration
proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = a5_snrna,
  addToArrow = TRUE,
  force=T,
  groupRNA = "cell_type",
  groupList = groupList, 
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)

pal <- paletteDiscrete(values = a5_snrna$cell_type)

gg_predictedGroup_Un <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "predictedGroup_Un", 
  pal = pal
)

gg_predictedGroup_predictedGroup <- plotEmbedding(
  proj, 
  colorBy = "cellColData", 
  name = "predictedGroup", 
  pal = pal
)

#Cell type marker genes
markerGenes  <- c(
  "FAP", "PDGFRB", "PDGFRA", "ACTA2", "COL1A1", #fibroblast_markers
  "MSR1", "CD163", "CCL3", #mono_macro_markers
  "PNMT", "TH", "DBH","CHGA", "CHGB", #chromaffin_markers
  "STAR", "CYP11B1", "CYP11A1", #adrenocortical_markers
  "FLT1", "EPAS1", #endothelial_markers
  "CDH19", "SOX10", "S100B", "VIM", #sustentacular_markers
  "CD2", "CD3E" , "MS4A1" #lymphocyte_markers
)

p1 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneIntegrationMatrix", 
  name = markerGenes, 
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

p1c <- lapply(p1, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

p2 <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  continuousSet = "horizonExtra",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

p2c <- lapply(p2, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

pdf(file = "/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/archr/plots/snrna_transfer_clusters.pdf", height = 10, width = 10, onefile = T)
print(gg_predictedGroup_Un)
print(gg_predictedGroup)
do.call(cowplot::plot_grid, c(list(ncol = 5), p1c))
do.call(cowplot::plot_grid, c(list(ncol = 5), p2c))
dev.off()

cM <- confusionMatrix(proj$Clusters, proj$predictedGroup)
labelOld <- rownames(cM)

labelNew <- colnames(cM)[apply(cM, 1, which.max)]

#Assign labels and override SCLCs
proj$ClustersRefined <- mapLabels(proj$Clusters, newLabels = labelNew, oldLabels = labelOld)
proj$ClustersRefined[proj$predictedGroup == "SCLCs"] <- "SCLCs"

saveArchRProject(proj, load=F)

################
# Peak finding #
################

addArchRThreads(threads = 1) 

setwd("/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/archr")

a5_anno.use <- bind_rows(a5_anno %>%  
                           dplyr::select(A5_ID, 
                                         differential_group_sampletype_strict, 
                                         TERT_ATRX_Mutation) %>% 
                           filter(A5_ID %in% proj@cellColData$Sample),
                         data.frame(A5_ID="NAM018",
                                    differential_group_sampletype_strict="Normal", 
                                    TERT_ATRX_Mutation="WT"))

anno_holder <- data.frame(A5_ID=as.character(proj@cellColData$Sample)) %>% left_join(a5_anno.use)

proj@cellColData$differential_group_sampletype_strict <- anno_holder$differential_group_sampletype_strict
proj@cellColData$TERT_ATRX_Mutation <- anno_holder$TERT_ATRX_Mutation
proj@cellColData$differential_group_sampletype_strict[proj@cellColData$predictedGroup != "Tumor"] <- "Normal"
proj@cellColData$TERT_ATRX_Mutation[proj@cellColData$predictedGroup != "Tumor"] <- "WT"
#Override a few cells in the normal chromaffin that get classified as tumour
proj@cellColData$predictedGroup[as.character(proj@cellColData$Sample) == "NAM018" & proj@cellColData$predictedGroup == "Tumor"] <- "Chromaffin cells"

saveArchRProject(proj, load=F)

addArchRThreads(threads = 4) 

proj_chromaffin_only <- subsetArchRProject(
  ArchRProj = proj,
  cells = rownames(proj@cellColData)[proj@cellColData$predictedGroup=="Tumor" | proj@cellColData$predictedGroup=="Chromaffin cells"],
  outputDirectory = "proj_chromaffin_only",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

proj_chromaffin_only@cellColData$GenotypeContrastGroup <- paste(proj_chromaffin_only@cellColData$predictedGroup, 
                                                                proj_chromaffin_only@cellColData$TERT_ATRX_Mutation, 
                                                                sep="-")


proj_chromaffin_only <- addGroupCoverages(ArchRProj = proj_chromaffin_only, 
                                          groupBy = "GenotypeContrastGroup",
                                          minRep = 3, 
                                          maxRep = 10, 
                                          minCells = 100, 
                                          maxCells = 500, 
                                          sampleRatio = 0.5,
                                          force = T)


proj_chromaffin_only <- addReproduciblePeakSet(ArchRProj = proj_chromaffin_only,
                                               groupBy = "GenotypeContrastGroup",
                                               reproducibility = "(n+1)/2",
                                               pathToMacs2="/g/data/pq08/projects/ppgl/a5/software/macs2/bin/macs2",
                                               plot = F,
                                               force = T)

proj_chromaffin_only <- addPeakMatrix(proj_chromaffin_only, force = T)

saveArchRProject(ArchRProj = proj_chromaffin_only, load = F)

##################
# Motif analysis #
##################

proj_chromaffin_only <- addMotifAnnotations(ArchRProj = proj_chromaffin_only, 
                                            motifSet = "cisbp", 
                                            name = "Motif", 
                                            force = T)

contrasts <- list(tert_vs_chrom = c("Tumor-TERT", "Chromaffin cells-WT"),
                  atrx_vs_chrom = c("Tumor-ATRX", "Chromaffin cells-WT"),
                  tert_vs_atrx = c("Tumor-TERT", "Tumor-ATRX"))

markers <- list()
for (contrast in names(contrasts))
{
  markers[[contrast]] <- getMarkerFeatures(
    ArchRProj = proj_chromaffin_only, 
    useMatrix = "PeakMatrix",
    groupBy = "GenotypeContrastGroup",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = contrasts[[contrast]][[1]],
    bgdGroups = contrasts[[contrast]][[2]]
  )
}

enriched_motifs <- list()
depleted_motifs <- list()

for (contrast in names(contrasts))
{
  enriched_motifs[[contrast]] <- peakAnnoEnrichment(
    seMarker = markers[[contrast]],
    ArchRProj = proj_chromaffin_only,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC >= 1"
  )
}

for (contrast in names(contrasts))
{
  depleted_motifs[[contrast]] <- peakAnnoEnrichment(
    seMarker = markers[[contrast]],
    ArchRProj = proj_chromaffin_only,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC <= -1"
  )
}

gg_motif <- list()
for (contrast in names(enriched_motifs))
{
  gg_motif[[contrast]] <- list()
  for (direction in c("enriched", "depleted")){
    if (direction == "enriched") {
      df <- data.frame(TF = rownames(enriched_motifs[[contrast]]), mlog10Padj = assay(enriched_motifs[[contrast]])[,1])
    } else
    {
      df <- data.frame(TF = rownames(depleted_motifs[[contrast]]), mlog10Padj = assay(depleted_motifs[[contrast]])[,1])
    }
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))
    
    contrast_members <- str_split_1(contrast, "_vs_")
    i = ifelse(direction == "enriched", 1, 2)
    plot_title <- paste(contrast, paste0("enriched ", contrast_members[[i]]), sep =" - ")
    
    gg_motif[[contrast]][[direction]] <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 3,
        nudge_x = 2,
        color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
      ggtitle(plot_title)
  }
}

gg_motif$tert_vs_chrom$enriched +
  gg_motif$tert_vs_chrom$depleted + 
  gg_motif$atrx_vs_chrom$enriched +
  gg_motif$atrx_vs_chrom$depleted +
  gg_motif$tert_vs_atrx$enriched +
  gg_motif$tert_vs_atrx$depleted +
  plot_layout(ncol = 2)


top_tables <- data.frame(TF=vector(mode = "character"),
                         mlog10Padj=vector(mode = "numeric"),
                         Enrichment=vector(mode = "numeric"),
                         contrast=vector(mode = "character"),
                         direction=vector(mode = "character"))
for (contrast in names(enriched_motifs))
{
  contrast_members <- str_split_1(contrast, "_vs_")
  df_enr <- data.frame(TF = rownames(enriched_motifs[[contrast]]), 
                       mlog10Padj = assay(enriched_motifs[[contrast]],1)[,1],
                       Enrichment = assay(enriched_motifs[[contrast]],3)[,1],
                       contrast=contrast,
                       direction=paste0("enriched ", contrast_members[[1]]))
  df_dep <- data.frame(TF = rownames(depleted_motifs[[contrast]]), 
                       mlog10Padj = assay(depleted_motifs[[contrast]],1)[,1],
                       Enrichment = assay(depleted_motifs[[contrast]],3)[,1],
                       contrast=contrast,
                       direction=paste0("enriched ", contrast_members[[2]]))
  
  top_tables <- bind_rows(top_tables, df_enr, df_dep)
}

motifs <- top_tables %>% group_by(contrast, direction) %>% 
  arrange(desc(mlog10Padj)) %>%  slice_head(n=10) 
motifs$TF <- gsub("_[0-9]+$", "", motifs$TF)

