################
# Data loading #
################

#This script relies on an archr project folder created by 
# /g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/archr/scripts/archr_workflow.r

setwd("/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/archr")
proj_chromaffin_only <- loadArchRProject(path = "proj_chromaffin_only")

##################################
# Find enriched/depleted motifs #
##################################

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

################################
# Motif enrichment score plots #
################################

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

gg_enrich <- gg_motif$tert_vs_chrom$enriched +
  gg_motif$tert_vs_chrom$depleted + 
  gg_motif$atrx_vs_chrom$enriched +
  gg_motif$atrx_vs_chrom$depleted +
  gg_motif$tert_vs_atrx$enriched +
  gg_motif$tert_vs_atrx$depleted +
  plot_layout(nrow = 1)


##########################
# Motif enrichment UMAPS #
##########################

motifs_of_interest <- c("AR","CEBPZ","CTCFL","CTCF","EGR1","FOSL1",
            "JUNB","JUND","KLF6","NFYB","NR3C1","PBX3",
            "PGR","RFX2","SP1","SP2","TFAP2A","TFAP2B",
            "TFAP2C","WT1","ZNF148")
motif_select <- paste(motifs_of_interest, collapse="|")
markerMotifs <- getFeatures(proj_chromaffin_only, select = motif_select, useMatrix = "MotifMatrix")
motif_select <- paste0("z:",paste(motifs_of_interest, collapse="_|z:"), "_")
markerMotifs <- markerMotifs[grepl(motif_select, markerMotifs)]

p <- plotEmbedding(
  ArchRProj = proj_chromaffin_only, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj_chromaffin_only)
)

p3 <- plotEmbedding(
  ArchRProj = proj_chromaffin_only, 
  colorBy = "cellColData", 
  name = "Sample", 
  embedding = "UMAP"
) 
p4 <- plotEmbedding(
  ArchRProj = proj_chromaffin_only, 
  colorBy = "cellColData", 
  name = "GenotypeContrastGroup", 
  embedding = "UMAP"
)

gg_motif_umap <- lapply(p, function(x){
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

gg_motif_umap$Samples <- p3 + ggtitle("") + xlab("") + ylab("")
gg_motif_umap$GenotypeContrastGroup <- p4 + ggtitle("") + xlab("") + ylab("")
purrr::reduce(.x = gg_motif_umap, .f = `+`) 
