library(Seurat)
library(Signac)
#library(rtracklayer)
library(tidyverse)
library(patchwork)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomeInfoDb)


generate_snatac_qc <- function(snatac_object, qc_out_dir) {
  

  tss_group_plot <- TSSPlot(snatac_object, group.by = "high.tss", assay = "cellranger") + NoLegend()
  
  frag_hist_nuc_group <- FragmentHistogram(object = snatac_object, group.by = "nucleosome_group", assay = "cellranger")
  
  # QC plot
  multiqc_violin <- VlnPlot(
    object = snatac_object,
    features = c("pct_reads_in_peaks", "peak_region_fragments","passed_filters",
                 "TSS.enrichment", "blacklist_ratio", "nucleosome_signal"),
    pt.size = 0,
    ncol = 3,
    group.by = "sampleID")
  
  # vln plot the qc cutoffs 
  vln_1 <- VlnPlot(
    object = snatac_object,
    group.by = "sampleID", 
    features = "passed_filters", 
    pt.size = 0) +
    geom_hline(yintercept=5000, linetype = "dashed", colour = "red") +
    NoLegend() +
    labs(title = "fragments", x = NULL)#+ 
  #scale_fill_manual(values = sample_colours)
  
  # NOTE: Blacklist ratio is not supported by cellranger v2.0
  # there is a bug in cellranger atac aggr that is making blacklist_region_fragments the same as TSS_fragments 
  
  vln_2 <- VlnPlot(
    object = snatac_object,
    group.by = "sampleID", 
    features = "blacklist_ratio", 
    pt.size = 0) +
    geom_hline(yintercept=0.005, linetype = "dashed", colour = "red") +
    NoLegend() +
    labs(title = "blacklist ratio", x = NULL) 
  
  vln_3 <- VlnPlot(
    object = snatac_object,
    group.by = "sampleID", 
    features = "TSS.enrichment", 
    pt.size = 0) +
    geom_hline(yintercept=2, linetype = "dashed", colour = "red") +
    NoLegend() +
    labs(title = "TSS enrichment score", x = NULL) 
  
  # individual fragment histograms
  Frag_hist <- FragmentHistogram(snatac_object, assay = "cellranger") 
  
  # Make a QC plot for fragment length of all samples on the same plot 
  # using geom_smooth
  
  reads <- Signac:::MultiGetReadsInRegion(snatac_object, assay = "cellranger", region = "chr1-1-2000000")
  
  # add sample names to the fragments dataframe
  fragments_df <- snatac_object@meta.data %>%
    dplyr::select(sampleID, sequencing_library) %>% 
    rownames_to_column(var = "cell") %>% 
    left_join(reads, by = "cell")
  
  fragment_dist <- ggplot(fragments_df) +
    geom_density(aes(x = length, colour = sampleID)) +
    theme_classic()
  
  pdf(file = file.path(qc_out_dir,"sn_atac_qc.pdf"), width = 8, height=7, onefile = T,)
  print(tss_group_plot)
  print(frag_hist_nuc_group)
  print(multiqc_violin)
  print(vln_1)
  print(vln_2)
  print(vln_3)
  print(fragment_dist)
  dev.off()
  
}

perform_cell_type_transfer <- function(snatac_object, snrna_object,  output_qc, qc_out_dir)
{
  
  message("Starting snRNA-snATAC cell-type label transfer")
  
  # Latent semantic indexing
  message("Performing latent semantic indexing ...")
  snatac_object <- FindTopFeatures(snatac_object, min.cutoff = 10)
  snatac_object <- RunTFIDF(snatac_object)
  snatac_object <- RunSVD(snatac_object, reduction.key = "LSIcellranger_", reduction.name = "lsi.cellranger")
  
  message("Finding clusters ...")
  # the first component is correlated to sequencing depth, therefore it is removed from downstream analysis 
  snatac_object <- RunUMAP(snatac_object, reduction = "lsi.cellranger", dims = 2:30, reduction.key = "UMAPcellranger_",reduction.name = "umap.cellranger")
  snatac_object <- FindNeighbors(snatac_object, reduction = "lsi.cellranger", dims = 2:50)
  snatac_object <- FindClusters(snatac_object, algorithm = 3, resolution = 1.5)
  
  #----
  # generate gene-activity matrix 
  #----
  message("Generating gene-activity matrix ...")
  #compute gene activity
  gene.activities <- GeneActivity(snatac_object)
  # add the gene activity matrix to the Seurat object as a new assay and normalize it
  snatac_object[["RNA"]] <- CreateAssayObject(counts = gene.activities)
  snatac_object <- NormalizeData(
    object = snatac_object,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = median(snatac_object$nCount_RNA)
  )
  DefaultAssay(snatac_object) <- "RNA"
  
  #----
  # annotate atac nuclei with snRNA-seq reference 
  #----
  
  Idents(snrna_object) <- snrna_object$cell_type # make sure correct cell identities are set
  # perform label transfer 
  message("Performing cell-type label transfer ...")
  transfer.anchors <- FindTransferAnchors(
    reference = snrna_object,
    query = snatac_object,
    reduction = "cca"
  )
  
  predicted.labels <- TransferData(
    anchorset = transfer.anchors,
    refdata = snrna_object$cell_type,
    weight.reduction = snatac_object[["lsi.cellranger"]],
    dims = 2:30
  )
  
  snatac_object <- AddMetaData(object = snatac_object, metadata = predicted.labels)
  # create a metadata column for describing individual tumour samples and predicted cell types 
  snatac_object@meta.data <- snatac_object@meta.data %>% 
    mutate(neoplastic_cell_sample = case_when(
      snatac_object$sampleID == "NAM018" ~ "Normal",
      predicted.id == "Tumor"~ sampleID,
      .default = "Normal"))
  
  # make a column for the normal cell type or tumor sample ID
  # this annotation is used for peak calling, so that it is performed
  # individually for each cell type and tumor sample
  
  snatac_object@meta.data <- snatac_object@meta.data %>% 
    mutate(cluster=case_when(
      sampleID == "NAM018" & predicted.id == "Tumor" ~ "Chromaffin cells", # correct any normal cells classified as tumor
      sampleID == "NAM018" ~ predicted.id,
      predicted.id == "Tumor"~ sampleID,
      TRUE ~ predicted.id))
  
  Idents(snatac_object) <- snatac_object$predicted.id
  
  
  if(output_qc)
  { 
    
    pdf(file = file.path(qc_out_dir,"cell_type_transfer_qc.pdf"), onefile = T, width = 7, height = 8)
    DepthCor(snatac_object, reduction = 'lsi.cellranger')
    
    DimPlot(snatac_object)
    DimPlot(snatac_object, group.by="cluster", cols = c25)
    DimPlot(snatac_object,group.by = "predicted.id", cols = c25)
    
    #----
    # UMAP plots 
    #----
    
    # plot samples
    atac_umap_sample <- DimPlot(snatac_object, group.by = "sampleID", pt.size = 0.5, label = T) +
      ggtitle("atac sample") + NoLegend()
    
    # plot clusters
    atac_umap_cluster <- DimPlot(snatac_object, group.by = "seurat_clusters", pt.size = 0.5, label = T) +
      ggtitle("atac cluster") + NoLegend()
    
    # plot cell types 
    atac_umap_predictedid <- DimPlot(snatac_object, group.by = "predicted.id", pt.size = 0.5) +
      ggtitle("atac predicted cell type")
    
    atac_umap_celltype <- DimPlot(snatac_object, group.by = "neoplastic_cell_sample", pt.size = 0.5, label = F) +
      ggtitle("Tumour samples (neoplastic cells)")
    
    atac_umap_sample + atac_umap_cluster + atac_umap_celltype
    
    FeaturePlot(snatac_object, features = c("prediction.score.max"))  + atac_umap_celltype
    
    VlnPlot(snatac_object, features = c("prediction.score.max"), group.by = "seurat_clusters")
    
    hist(snatac_object$prediction.score.max)
    
    abline(v = 0.5, col = "red")
    
    # Marker gene expression DotPlot
    fibroblast_markers <- c("FAP", "PDGFRB", "PDGFRA", "ACTA2", "COL1A1")
    mono_macro_markers <-c("MSR1", "CD163", "CCL3")
    chromaffin_markers <- c("PNMT", "TH", "DBH","CHGA", "CHGB")
    adrenocortical_markers <- c("STAR", "CYP11B1", "CYP11A1")
    endothelial_markers <- c("EPAS1", "FLT1")
    sclc_markers <- c("CDH19", "SOX10", "S100B", "VIM")
    lymphocyte_markers <- c("CD2", "CD3E" , "MS4A1")
    
    marker_genes <- c(
      adrenocortical_markers,
      chromaffin_markers,
      endothelial_markers,
      fibroblast_markers,
      lymphocyte_markers,
      mono_macro_markers,
      sclc_markers)
    
    # make a dotplot for marker gene activity in the atac cell types  
    DotPlot(
      snatac_object,
      features = rev(marker_genes),
      group.by = "predicted.id",
      cols = c("grey", "red")) +
      scale_colour_distiller(palette = "RdYlBu")+
      ggtitle("snATAC-seq gene activity") + theme(axis.text.x = element_text(angle = 90, hjust=0.9, vjust=0.5))
    
    # make a dotplot for marker gene activity in the atac clusters
    DotPlot(
      snatac_object,
      features = rev(marker_genes),
      group.by = "seurat_clusters",
      cols = c("grey", "red"))+
      scale_colour_distiller(palette = "RdYlBu")+
      ggtitle("snATAC-seq gene activity") + theme(axis.text.x = element_text(angle = 90, hjust=0.9, vjust=0.5))
    
    FeaturePlot(snatac_object, features = c("TERT"))
    
    
    dev.off() 
  }
  
  # filter out the low-scoring cells 
  snatac_object <- subset(snatac_object,
                      subset = prediction.score.max > 0.5) 
  
  return(snatac_object)
  
}

data_loader_a5_snatac <- function(quickload=T,
                                  quickload_checkpoint_dir="/g/data/pq08/projects/ppgl/a5/sn_atac_seq/quickload_checkpoints/", 
                                  aggregated_count_h5="/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/cellranger/aggregated_hg38/outs/filtered_peak_bc_matrix.h5",
                                  barcode_meta_data="/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/cellranger/aggregated_hg38/outs/singlecell.csv",
                                  fragments_file="/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/cellranger/aggregated_hg38/outs/fragments.tsv.gz",
                                  aggregation_meta_data="/g/data/pq08/projects/ppgl/a5/sn_atac_seq/analysis/cellranger/aggregated_hg38/outs/aggregation_csv.csv",
                                  remove_samples=c("E140-1","E146-1","E171-1","E225-1"),
                                  genome_version = "hg38",
                                  seurat_min_cells=3,
                                  seurat_min_features=200,
                                  output_qc=F,
                                  snrna_cell_type_transfer = T,
                                  qc_out_dir="/g/data/pq08/projects/ppgl/a5/sn_atac_seq/qc",
                                  qc_thresholds=list(frags_passed_filters_min = 5000, 
                                                     TSS_enrichment = 2,
                                                     nucleosome_signal = 4)
)
{
  
  message("Starting A5 sn-ATAC data loader...")
  
  base_dir="/g/data/pq08/projects/ppgl"
  
  if(!exists("a5_anno")) {
    source(file.path(base_dir,"a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r"))
    data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.org", use_cache = T, remove_excluded_samples = T)
  }
  
  source(file.path(base_dir,"./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r"))
  
  if(quickload)
  {
    message("Read sn_ATAC data from quickload checkpoint ... ")
    a5_snatac <- readRDS(file.path(quickload_checkpoint_dir, "sn_atac_hg38_filtered.rds"))
    
  } else {
    
    message("Starting aggregate data read-in...")
    set.seed(1234)
    
    #Load count data
    counts <- Read10X_h5(filename = aggregated_count_h5)
    metadata <- read.csv(
      file = barcode_meta_data,
      header = TRUE,
      row.names = 1
    )
    
    #Create assay
    chrom_assay <- CreateChromatinAssay(
      counts = counts,
      sep = c(":", "-"),
      genome = genome_version,
      fragments = fragments_file,
      min.cells = seurat_min_cells,
      min.features = seurat_min_features
    )
    
    #Create Seurat object
    a5_snatac <- CreateSeuratObject(
      counts = chrom_assay,
      assay = "cellranger",
      meta.data = metadata
    )
    

    
    # extract gene annotations from EnsDb
    # get gene annotations for hg38
    annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotation) <- "UCSC"
    genome(annotation) <- "hg38"
    Annotation(a5_snatac) <- annotation
    
    #----
    # label barcodes according to sample of origin 
    #----
    
    #Aggregation tables contain the original sample identities 
    # (number appended to each cell barcode correspond to the row number of the aggregation table)
    aggregation_meta <- read.csv(aggregation_meta_data)
    colnames(aggregation_meta)[1] <- "Sample_Index"
    aggregation_meta$Sample_Index <- as.character(aggregation_meta$Sample_Index)
    aggregation_meta$library_id <- gsub("_", "-", aggregation_meta$library_id)
    
    #Join original sample IDs to table of barcodes
    Barcodes <- data.frame(Barcode=colnames(a5_snatac)) %>% 
      separate(Barcode, into = c("Barcode","Sample_Index"), sep="-") %>% 
      left_join(aggregation_meta %>% dplyr::select(Sample_Index, library_id))
    
    #store sample index in metadata
    a5_snatac[["sequencing_library"]] <- Barcodes$Sample_Index
    
    #store original sample ID info in the metadata
    a5_snatac[["sampleID"]] <- Barcodes$library_id
    
    Idents(a5_snatac) <- a5_snatac$sampleID
    
    # Annotate clinical data 
    a5_anno.use <- a5_anno %>%
      dplyr::filter(A5_ID %in% aggregation_meta$library_id)
    
    # add the clinical data to the metadata
    a5_snatac_md <- a5_snatac@meta.data
    a5_snatac_md <- a5_snatac_md %>%
      rownames_to_column(var = "barcode") %>%
      left_join(a5_anno.use, by = c("sampleID"="A5_ID"))
    rownames(a5_snatac_md) <- a5_snatac_md$barcode
    a5_snatac@meta.data <- a5_snatac_md
    
    
    #----
    # QC metric computation 
    #----
    
    # compute nucleosome signal score per cell
    a5_snatac <- NucleosomeSignal(object = a5_snatac)
    # compute TSS enrichment score per cell
    a5_snatac <- TSSEnrichment(a5_snatac, fast = F, assay = "cellranger")
    # add blacklist ratio and fraction of reads in peaks
    a5_snatac$pct_reads_in_peaks <- a5_snatac$peak_region_fragments / a5_snatac$passed_filters * 100
    # blacklist_ratio calculation was changed from peak_region_fragments to passed_filters as I want to avoid bias towards the more common cell types
    a5_snatac$blacklist_ratio <- a5_snatac$blacklist_region_fragments / a5_snatac$passed_filters 
    # TSS Plot
    a5_snatac$high.tss <- ifelse(a5_snatac$TSS.enrichment > qc_thresholds$TSS_enrichment, "High", "Low")
    
    # Nucleosome plot 
    a5_snatac$nucleosome_group <- ifelse(a5_snatac$nucleosome_signal > qc_thresholds$nucleosome_signal, 
                                         paste("NS > ", qc_thresholds$nucleosome_signal), 
                                         paste("NS < ", qc_thresholds$nucleosome_signal))
    
    #############
    # QC Output #
    #############
    
    if(output_qc)
    {
      if(!exists("ColorPalette")) {
        source(paste0(base_dir,"/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r"))
      }
      generate_snatac_qc(a5_snatac, qc_out_dir)
      
    }
    
    #----
    # Filter out poor quality Samples
    #----
    a5_snatac$sample_qc <- ifelse(a5_snatac$sampleID %in% c(remove_samples), "poor","ok")
    a5_snatac <- subset(
      a5_snatac,
      sample_qc == "ok")
    
    #----
    # Filter out poor quality cells 
    #----
    message("Removing low quality cells")
    #table(a5_snatac$passed_filters > qc_thresholds$frags_passed_filters_min & a5_snatac$TSS.enrichment > qc_thresholds$TSS_enrichment) # passed filters is the fragment number
    
    a5_snatac$cell_qc <-  if_else(a5_snatac$passed_filters < qc_thresholds$frags_passed_filters_min &
                                    a5_snatac$TSS.enrichment < qc_thresholds$TSS_enrichment,
                                  "poor","ok")
    # compare the cells i'm keeping vs removing
    # FragmentHistogram(object = a5_snatac, group.by = "cell_qc")+
    #   TSSPlot(a5_snatac, group.by = "cell_qc") + NoLegend()
    
    # I want to filter based on total reads (passed_filters) and TSS enrichment
    # avoiding filtering based on the cellranger peaks because I'm calling peaks with MACS2 later on
    # reads within peaks could bias the counts of under-represented cell types
    a5_snatac <- subset(
      a5_snatac,
      cell_qc == "ok")
    
    if(snrna_cell_type_transfer){
      snrna_loaded_locally <- F
      if(!exists("a5_snrna")) {
        source(file.path(base_dir,"a5/sn_rna_seq/scripts/data_loaders/a5_snrna_dataloader.r"))
        data_loader_a5_snrna(quickload=T)
        assign(x = "a5_snrna", value = snrna_annotate_cell_types(a5_snrna), envir = globalenv()) 
        snrna_loaded_locally <- T
      }
      
      a5_snatac <- perform_cell_type_transfer(snatac_object = a5_snatac,snrna_object = a5_snrna, output_qc, qc_out_dir)
      
      if (snrna_loaded_locally) {rm(a5_snrna, envir = globalenv())}
    }
    
    #reset default assay
    DefaultAssay(a5_snatac) <- "cellranger"
    
    message("Writing quickload checkpoint ... ")
    saveRDS(a5_snatac, file.path(quickload_checkpoint_dir, "sn_atac_hg38_filtered.rds"))
    
    
  }
  
  assign("a5_snatac", a5_snatac, globalenv())
  message("Loaded objects into global environment: 
          a5_snatac - Seurat object with snATAC counts and metadata")
  
  message("Completed A5 snATAC data loader.")
}
message("Created data loader function data_loader_a5_snatac()")