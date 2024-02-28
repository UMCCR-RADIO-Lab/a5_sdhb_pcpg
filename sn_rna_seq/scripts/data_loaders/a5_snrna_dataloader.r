library(Seurat)
library(ggplot2)
library(patchwork)

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

dev_off_safe <- function() { if(names(dev.cur()) != "null device") {dev.off()} }

#Function to generate quality control plots and tables
generate_snrna_qc <- function(sn_data,
                              qc_out_dir,
                              nCount_RNA_upper,
                              nCount_RNA_lower,
                              nFeature_RNA_upper,
                              nFeature_RNA_lower,
                              percent_mt_cutoff)
{
  
  
  metadata <- sn_data@meta.data
  n_samples <- length(unique(metadata$orig.ident))
  
  set.seed(1234)
  sample_colours = setNames(sample(x = ColorPalette,size = n_samples,replace = F), unique(metadata$orig.ident))
  
  ############
  # QC Plots #
  ############
  
  # Visualize the number of cell counts per sample
  ncells_bar <- metadata %>% 
    ggplot(aes(x=A5_ID, fill=A5_ID)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5, face="bold")) +
    ggtitle("NCells") + 
    scale_fill_manual(values = sample_colours)
  ggsave(plot = ncells_bar, filename = file.path(qc_out_dir, "ncells_per_sample_barplot.pdf"))
  dev_off_safe()
  
  # Visualize the number UMIs/transcripts per cell
  nCount_density <- metadata %>% 
    ggplot(aes(color=A5_ID, x=nCount_RNA, fill= A5_ID)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    geom_vline(xintercept = nCount_RNA_lower, linetype = "dashed", colour = "red") + 
    scale_colour_manual(values = sample_colours) + 
    scale_fill_manual(values = sample_colours) +
    theme(legend.position = "none")
  
  # Visualize the distribution of genes detected per cell via histogram
  nfeatures_density <- metadata %>% 
    ggplot(aes(color=A5_ID, x=nFeature_RNA, fill= A5_ID)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() +
    ylab("Cell density") +
    geom_vline(xintercept = nFeature_RNA_lower, linetype = "dashed", colour = "red") + 
    scale_colour_manual(values = sample_colours) + 
    scale_fill_manual(values = sample_colours) +
    theme(legend.position = "none")
  
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
  complexity_density <- metadata %>%
    ggplot(aes(x=log10.genes.per.umi, color = A5_ID, fill=A5_ID)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    scale_colour_manual(values = sample_colours) + 
    scale_fill_manual(values = sample_colours)
  
  ggsave(plot = (nfeatures_density + nCount_density + complexity_density), 
         filename = file.path(qc_out_dir, "nfeatures_ncount_complexity_density.pdf"),
         width = 15,
         height= 6)
  dev_off_safe()
  
  # visualise percentage mitochondrial, number of counts and number of genes as violin plots 
  nfeature_vln <- VlnPlot(sn_data, features = "nFeature_RNA", group.by = "A5_ID", pt.size = 0) +
    geom_hline(yintercept = c(nFeature_RNA_lower, nFeature_RNA_upper), linetype = "dashed", colour = "red") +
    scale_fill_manual(values = sample_colours) + 
    labs(title = "nFeature_RNA", x = NULL) 
  ncount_vln <- VlnPlot(sn_data, features = "nCount_RNA", group.by = "A5_ID", pt.size = 0) +
    geom_hline(yintercept = nCount_RNA_lower, linetype = "dashed", colour = "red") +
    scale_fill_manual(values = sample_colours) + 
    coord_cartesian(ylim = c(0 , 25000)) +
    labs(title = "nCount_RNA", x = NULL) 
  mt_vln <- VlnPlot(sn_data, features = "percent.mt", group.by = "A5_ID", pt.size = 0) +
    geom_hline(yintercept = percent_mt_cutoff, linetype = "dashed", colour = "red") +
    scale_fill_manual(values = sample_colours) + 
    coord_cartesian(ylim = c(0 , 40)) +
    labs(title = "percent.mt", x = NULL)
  
  ggsave(plot = (nfeature_vln + 
                   ncount_vln + 
                   mt_vln +
                   plot_layout(ncol = 3) &
                   theme(legend.position = "none")), 
         filename = file.path(qc_out_dir, "nfeatures_ncount_violin.pdf"),
         width=15,
         height = 6)
  dev_off_safe()
  
  featurescatter1 <- FeatureScatter(sn_data, feature1 = "nCount_RNA", feature2 = "percent.mt", group.by = 'predicted_doublet')
  featurescatter2 <- FeatureScatter(sn_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = 'predicted_doublet')
  
  ggsave(plot = (featurescatter1 + featurescatter2) & scale_color_discrete(name='predicted_doublet'), 
         filename = file.path(qc_out_dir, "feature_scatter.png"), device = png(), dpi = 300,
         width=16,
         height=10)
  dev_off_safe()
  
  ############
  # QC Table #
  ############
  
  qc_data <- bind_cols(orig.ident=metadata$orig.ident, sn_data[[c("nCount_RNA","nFeature_RNA","percent.mt")]], )
  
  
  compute_percentiles <- function(qc)
  {
    qc %>% 
      dplyr::select(-orig.ident) %>% 
      as.list %>% 
      purrr::map(quantile) %>% 
      as.data.frame() %>% as_tibble(rownames = "quantile") %>% 
      tidyr::pivot_longer(cols = -quantile, names_to="metric") %>% 
      arrange(metric) %>% mutate(metric=paste(metric, gsub("%","th percentile",quantile))) %>% 
      dplyr::select(-quantile) 
  }
  
  
  #Quantile data
  
  quantile_data <- qc_data %>% 
    group_by(orig.ident) %>% 
    {
      rlang::set_names(dplyr::group_split(.), pull(group_keys(.), orig.ident))
    } %>% 
    purrr::map(.f = compute_percentiles) %>% bind_rows(.id="orig.ident") %>% 
    bind_rows(compute_percentiles(qc_data) %>% mutate(orig.ident="All"))
  
  
  ## nCount_RNA
  nCount_RNA_per_sample <- 
    qc_data %>%
    mutate(nCount_RNA_category = cut(nCount_RNA,
                                     breaks = c(0,
                                                nCount_RNA_lower,
                                                nCount_RNA_upper,
                                                10 ^ 10),
                                     labels = c(
                                       paste(0, nCount_RNA_lower, sep = "-"),
                                       paste(nCount_RNA_lower, nCount_RNA_upper, sep ="-"),
                                       paste(">", nCount_RNA_upper)
                                     )
    )) %>%
    group_by(orig.ident, nCount_RNA_category) %>% 
    dplyr::count() %>%  
    group_by(orig.ident) %>% 
    mutate(proportion = round(n / sum(n),4))
  
  nCount_RNA_all <- nCount_RNA_per_sample %>% 
    group_by(nCount_RNA_category) %>% 
    summarise(n=sum(n)) %>% 
    ungroup() %>% 
    mutate(proportion = round(n / sum(n),4)) %>% mutate(orig.ident="All Samples")
  
  nCount_RNA <- bind_rows(nCount_RNA_per_sample, nCount_RNA_all) %>% 
    dplyr::rename(nCount_RNA = nCount_RNA_category)
  
  ## nFeature_RNA
  nFeature_RNA_per_sample <- qc_data %>%
    mutate(nFeature_RNA_category = cut(nFeature_RNA,
                                       breaks = c(0,
                                                  nFeature_RNA_lower,
                                                  nFeature_RNA_upper,
                                                  10 ^ 10),
                                       labels = c(paste(0, nFeature_RNA_lower, sep = "-"),
                                                  paste(nFeature_RNA_lower, nFeature_RNA_upper, sep = "-"),
                                                  paste(">", nFeature_RNA_upper)
                                       )
    )) %>%
    group_by(orig.ident, nFeature_RNA_category) %>% 
    dplyr::count() %>%  
    group_by(orig.ident) %>% mutate(proportion = round(n / sum(n),4))
  
  nFeature_RNA_all <- nFeature_RNA_per_sample %>% 
    group_by(nFeature_RNA_category) %>% 
    summarise(n=sum(n)) %>% 
    ungroup() %>% 
    mutate(proportion = round(n / sum(n),4)) %>% mutate(orig.ident="All Samples")
  
  nFeature_RNA <- bind_rows(nFeature_RNA_per_sample, nFeature_RNA_all)  %>% 
    dplyr::rename(nFeature_RNA = nFeature_RNA_category)
  
  ## percent.mt
  percent_mt_per_sample <- 
    qc_data %>%
    mutate(percent_mt_category = cut(percent.mt,
                                     breaks = c(0,
                                                percent_mt_cutoff,
                                                100),
                                     labels = c(paste(0, percent_mt_cutoff, sep = "-"),
                                                paste(percent_mt_cutoff, 100, sep = "-")), 
                                     include.lowest=T)) %>%
    group_by(orig.ident, percent_mt_category) %>% 
    dplyr::count() %>%  
    group_by(orig.ident) %>% 
    mutate(proportion = round(n / sum(n),4)) 
  
  percent_mt_all <- percent_mt_per_sample %>% 
    group_by(percent_mt_category) %>% 
    summarise(n=sum(n)) %>% 
    ungroup() %>% 
    mutate(proportion = round(n / sum(n),4)) %>% 
    mutate(orig.ident="All Samples")
  
  percent_mt <- bind_rows(percent_mt_per_sample, percent_mt_all) %>% 
    dplyr::rename(percent_mt = percent_mt_category)
  
  # write to file
  qc_file_path <- file.path(qc_out_dir,"qc_metrics_quantiles.tsv")
  readr::write_delim(x = quantile_data, file = qc_file_path, col_names = T)
  
  qc_file_path <- file.path(qc_out_dir,"qc_metrics_ncount.tsv")
  readr::write_delim(x = nCount_RNA, file = qc_file_path, col_names = T)
  
  qc_file_path <- file.path(qc_out_dir,"qc_metrics_nfeature.tsv")
  readr::write_delim(x = nFeature_RNA, file = qc_file_path, col_names = T)
  
  qc_file_path <- file.path(qc_out_dir,"qc_metrics_percent_mt.tsv")
  readr::write_delim(x = percent_mt, file = qc_file_path, col_names = T)
  
}

data_loader_a5_snrna <- function(quickload=T,
                                 cell_ranger_count_dir="/g/data/pq08/projects/ppgl/a5/sn_rna_seq/analysis/cellranger_hg38/intronic_counted/",
                                 scrublet_out_dir="/g/data/pq08/projects/ppgl/a5/sn_rna_seq/analysis/scrublet/hg38/intronic_counted/",
                                 quickload_checkpoint_dir="/g/data/pq08/projects/ppgl/a5/sn_rna_seq/quickload_checkpoints", 
                                 seurat_min_cells=3,
                                 seurat_min_features=200,
                                 exclude_samples=c(),
                                 output_qc=F,
                                 qc_out_dir="/g/data/pq08/projects/ppgl/a5/sn_rna_seq/qc",
                                 sample_rename = c(),
                                 perform_feature_regression=T,            
                                 qc_thresholds=list(nCount_RNA_upper = 16000, 
                                                    nCount_RNA_lower = 750, 
                                                    nFeature_RNA_upper = 6000, 
                                                    nFeature_RNA_lower = 500,
                                                    percent_mt_cutoff=10),
                                 regression_vars=c("S.Score", "G2M.Score", "nCount_RNA", "percent.mt")
)
{
  
  message("Starting A5 snRNA data loader.")
  
  base_dir="/g/data/pq08/projects/ppgl"
  
  if(!exists("a5_anno")) {
    source(paste0(base_dir,"/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r"))
    data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.org", use_cache = T, remove_excluded_samples = T)
  }
  
  if(!("differential_group" %in% colnames("a5_anno"))) {
    a5_anno <- add_differential_groups(a5_anno)
  }
  
  if(quickload)
  {
    
    
    if(perform_feature_regression)
    {
      quick_load_file <- file.path(quickload_checkpoint_dir,"a5_scrna_quality_filtered_regressed.rds") 
    } else 
    {
      quick_load_file <- file.path(quickload_checkpoint_dir,"a5_scrna_quality_filtered.rds")
    }
    
    if(file.exists(quick_load_file)) {
      message("Performing sn_RNA quickload from:", quick_load_file)
      message("--------------------")
      message("Warning: During quickload parameters such as regression_vars, exclude_samples, sample_rename, and
            any qc thresholds will reflect the values at the time of the checkpoint not those currently specified")
      message("--------------------")
      
      a5_snrna <- readRDS(quick_load_file)
      
    } else {
      message("Quickload file (", quick_load_file,') not found. Ensure you have run the data loader without quickload first')
    }
    
  } else {
    
    ####################
    # Helper functions #
    ####################
    
    #Helper to read counts into seurat objects #
    make_seurat <- function(sample_name, matrix_dir, scrublet_csv, min_cells, min_features){
      #' read seurat object from 10x output and add scrublet predictions to metadata 
      sample.data <- Read10X(data.dir = matrix_dir)
      sample.seurat <- CreateSeuratObject(counts = sample.data, 
                                          project = sample_name, 
                                          min.cells = min_cells, 
                                          min.features = min_features)
      # add scrublet results to metadata
      sample.scrublet <- read.csv(scrublet_csv, sep = ",")
      rownames(sample.scrublet) <- colnames(sample.data)
      sample.seurat <- AddMetaData(sample.seurat, metadata = sample.scrublet)
      return(sample.seurat)
    }
    
    ########################
    # Populate input files #
    ########################
    
    # make a list with all the sample names 
    sample_names <- list.dirs(path = cell_ranger_count_dir, full.names = F, recursive = F)
    
    sample_names <- setdiff(sample_names, exclude_samples)
    
    matrix_dirs <- file.path(cell_ranger_count_dir, 
                             sample_names, 
                             "outs/filtered_feature_bc_matrix/")
    
    scrublet_csvs <- file.path(scrublet_out_dir, 
                               paste0(sample_names, 
                                      "_scrublet_output_table.csv"))
    
    #Rename samples as required
    if(!is.null(sample_rename)) {
      sample_names <- recode(sample_names, !!!sample_rename)
    }
    
    #Set vector names
    names(sample_names) <- sample_names
    names(matrix_dirs) <- sample_names
    names(scrublet_csvs) <- sample_names
    
    #########################
    # Create Seurat Objects #
    #########################
    
    # loop through the samples, load them into seurat objects
    seurat_objects <- purrr::pmap(.l = list(sample_name=sample_names, 
                                            matrix_dir=matrix_dirs, 
                                            scrublet_csv=scrublet_csvs),
                                  .f = \(sample_name,matrix_dir,scrublet_csv) {
                                    make_seurat (sample_name,
                                                 matrix_dir,
                                                 scrublet_csv, 
                                                 min_cells=seurat_min_cells, 
                                                 min_features = seurat_min_features)})
    
    
    #########################
    # Merge Sample Datasets #
    #########################
    
    # merge all of the samples into one object
    a5_snrna <- merge(x=seurat_objects[[1]], 
                      y=seurat_objects[-1],
                      add.cell.ids = names(seurat_objects),
                      project = "A5_single_nuclei")
    
    rm(seurat_objects)
    
    
    ########################
    # Calculate QC metrics #
    ########################
    
    # calculate % mitochondrial gene expression
    a5_snrna[["percent.mt"]] <- PercentageFeatureSet(a5_snrna, pattern = "^MT-")
    # calculate log10 genes detected per umi 
    a5_snrna[["log10.genes.per.umi"]] <- log10(a5_snrna$nFeature_RNA) / log10(a5_snrna$nCount_RNA)
    
    ##################################
    # Sample annotation and metadata #
    ##################################
    
    # filter to retain only the snRNA-seq samples
    a5_anno_snrna <- a5_anno %>%
      dplyr::filter(A5_ID %in% sample_names) 
    
    # get available metadata for the pheo-atlas paper samples
    zethoven_metadata <- readr::read_tsv(file.path(base_dir, "a5/sample_annotation/metadata_sc_samples.tsv") ,show_col_types = F) 
    
    # harmonize the important single cell metadata fields with A5 metadata  
    zethoven_samples <- c("P018-PGL1", "P018-PGL3", "E240", "E243", "E018")
    
    zethoven_metadata <- zethoven_metadata %>%
      dplyr::rename(A5_ID = Sample) %>% 
      filter(A5_ID %in% c(zethoven_samples)) %>% 
      dplyr::select(A5_ID, 
                    Gender, 
                    tumour_metastasised=Malignancy, 
                    Primary_Location_Simplified=Location) %>% 
      mutate(#A5_ID = ifelse(!grepl("-", A5_ID), paste(A5_ID,"1", sep="-"), A5_ID),
        Gender = recode(Gender,"M" = "male","F" = "female"),
        Genotype = case_when(
          A5_ID %in% c("E240", "E243") ~ "WT", 
          A5_ID %in% c("E018") ~ "RET",
          A5_ID %in% c("P018-PGL1", "P018-PGL3") ~ "SDHB"),
        is_primary_or_met = case_when(
          A5_ID %in% c("E240", "E243") ~ "Normal", 
          A5_ID %in% c("P018-PGL1", "P018-PGL3", "E018") ~ "Primary"),
        tumour_metastasised = case_when(
          A5_ID %in% c("E240", "E243") ~ "Normal", 
          A5_ID %in% c("P018-PGL1", "P018-PGL3") ~ "No",
          A5_ID %in% c("E018") ~ "Short follow up"),
        is_head_and_neck = case_when(
          A5_ID %in% c("E240", "E243") ~ "Normal", 
          A5_ID %in% c("P018-PGL1", "P018-PGL3", "E018") ~ "PC_or_PGL"),
        Major_Cluster = case_when(
          A5_ID %in% c("E240", "E243") ~ "Normal", 
          A5_ID %in% c("P018-PGL1", "P018-PGL3", "E018") ~ "Chromaffin"),
        `Patient ID` = gsub("-.+", "", A5_ID),
        PublicationID = A5_ID,
        Primary_Location_Simplified = gsub("AT-PGL",
                                           "Extraadrenal_abdominal", 
                                           Primary_Location_Simplified),
        TERT_ATRX_Mutation = case_when(
          A5_ID %in% c("E240", "E243", "P018-PGL1", "P018-PGL3") ~ "WT", 
          A5_ID %in% c("E018") ~ "ATRX"),
        TERT_purity_adj_vaf = case_when(
          A5_ID %in% c("E018") ~ 0.28,
          TRUE ~ 0),
        ATRX_purity_adj_vaf = 0,
        sample_purity=0.51 #arbitrary value above 0.5 for use by the differential groups function
      ) %>% add_differential_groups()
    
    # join the A5 and singlecell metadata
    a5_anno_snrna <- bind_rows(a5_anno_snrna, zethoven_metadata)
    
    # add the clinical data to the metadata, retaining the rownames 
    metadata <- a5_snrna@meta.data
    metadata <- metadata %>%
      # add cell barcode column 
      tibble::rownames_to_column(var = "barcode") %>%
      # add a sample column 
      mutate(A5_ID = orig.ident) %>% 
      # join the clinical data to the metadata 
      left_join(a5_anno_snrna, by = "A5_ID") %>% 
      # add a column describing the library preparation chemistry that was used 
      mutate(chemistry = if_else(A5_ID %in% c("E140-1",
                                              "E143-1",
                                              "E171-1",
                                              "E225-1"),
                                 "3prime", "5prime"))
    
    # make the metadata colnames compatible with seurat commands
    compatiblenames <- setNames(colnames(metadata), make.names(colnames(metadata)))
    metadata <- metadata %>%
      rename(!!!compatiblenames)
    
    rownames(metadata) <- metadata$barcode
    a5_snrna@meta.data <- metadata[,-1]
    
    ######################
    # scMatch annotation #
    ######################
    
    # TODO: re-run this on the new data?
    # #read in the top scMatch annotations for each cell, add annotations to metadata
    # scMatch_annos <- read.csv("Data/snRNA-seq/scMatch/results/annotation_result_keep_all_genes/human_Spearman_top_ann.csv")
    # #make barcodes the same
    # scMatch_annos$cell <- gsub('.', '-',scMatch_annos$cell, 
    #                            fixed = TRUE) 
    # #check that cell order in the scMatch is still correct
    # table(colnames(rna)==(scMatch_annos$cell))
    # names(scMatch_annos$cell.type) <- (scMatch_annos$cell)
    # rna <- AddMetaData(rna, metadata = scMatch_annos$cell.type, col.name = 'scMatch_predicted_cell_type')
    
    #############
    # QC Output #
    #############
    
    if(output_qc)
    {
      if(!exists("ColorPalette")) {
        source(paste0(base_dir,"/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r"))
      }
      rlang:::inject(generate_snrna_qc(sn_data=a5_snrna, qc_out_dir, !!!qc_thresholds))
      
    }
    
    #####################
    # Quality Filtering #
    #####################
    
    a5_snrna <- subset(
      a5_snrna,
      subset = nFeature_RNA > qc_thresholds[["nFeature_RNA_lower"]] &
        nCount_RNA > qc_thresholds[["nCount_RNA_lower"]] &
        percent.mt < qc_thresholds[["percent_mt_cutoff"]] & 
        predicted_doublet == "False")
    
    
    saveRDS(a5_snrna, file.path(quickload_checkpoint_dir,"a5_scrna_quality_filtered.rds"))
    
    #########################
    # Cell Cycle Regression #
    #########################
    
    if (perform_feature_regression)
    {
      
      # Normalisation and Scaling- required for cell cycle scoring
      a5_snrna <- NormalizeData(a5_snrna)
      a5_snrna <- FindVariableFeatures(a5_snrna, selection.method = "vst", nfeatures = 2000)
      a5_snrna <- ScaleData(a5_snrna, features = rownames(a5_snrna))
      
      # make list of genes for s and g2 phases
      s.genes <- cc.genes$s.genes
      g2m.genes <- cc.genes$g2m.genes
      
      #assign phase based on phase-specific gene expression
      a5_snrna <- CellCycleScoring(a5_snrna, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
      
      
      #To visualise the cell-cycle-related variance, run PCA using the s and g2 gene lists
      if(output_qc)
      {
        a5_snrna <- RunPCA(a5_snrna, features = c(s.genes, g2m.genes), approx=F)
        
        plot1 <- DimPlot(object = a5_snrna, group.by = "Phase") + ggtitle("Phase - before regression - cell cycle genes only")
        
        a5_snrna <- RunPCA(a5_snrna)
        plot2 <- DimPlot(object = a5_snrna, group.by = "Phase") + ggtitle("Phase - before regression - all genes") #Clearly the G2M cells are on the right side of PC1
        
        plot3 <- FeaturePlot(a5_snrna, features = "nCount_RNA", order = TRUE) + 
          scale_color_gradientn(values=scales::rescale(c(0,500,1000,5000,10000,20000,100000)),
                                colors = c("grey","red", "orange", "purple", "green", "blue"))
        plot4 <- FeaturePlot(a5_snrna, features = "percent.mt", order = TRUE)
      }
      
      #regress out cell cyle phase, percentage mitochrondrial rna, sequencing depth 
      #Log normalise and scale: this will create a normalised and regressed dataset that I can use for downstream analysis
      a5_snrna <- ScaleData(a5_snrna, vars.to.regress = c(regression_vars), features = rownames(a5_snrna))
      saveRDS(a5_snrna, file.path(quickload_checkpoint_dir,"a5_scrna_quality_filtered_regressed.rds"))
      
      if(output_qc)
      {
        a5_snrna <- RunPCA(a5_snrna)
        plot5 <- DimPlot(object = a5_snrna, group.by = "Phase") + ggtitle("Phase - after regression - all genes")
        
        plot6 <- FeaturePlot(a5_snrna, features = "nCount_RNA", order = TRUE) + 
          scale_color_gradientn(values=scales::rescale(c(0,500,1000,5000,10000,20000,100000)),
                                colors = c("grey","red", "orange", "purple", "green", "blue"))
        plot7 <- FeaturePlot(a5_snrna, features = "percent.mt", order = TRUE)
        
        pdf(file = file.path(qc_out_dir,"pca_and_cell_cycle_regression.pdf"), onefile = T)
        purrr::walk(.x=list(plot1,plot2,plot3,plot4,plot5,plot6,plot7), .f = print)
        dev_off_safe()
      }
      
    }
    
  }
  
  ################################
  #Export to  global environment #
  ################################
  
  assign("a5_snrna", a5_snrna, globalenv())
  
  message("Loaded objects into global environment: 
          a5_snrna - Seurat object with snRNA counts and metadata")
  
  message("Completed A5 snRNA data loader.")
}
message("Created data loader function data_loader_a5_snrna()")

snrna_annotate_cell_types <- function(snrna_object, output_qc=FALSE, qc_out_dir=NULL)
{
  
  ################
  # Marker Genes #
  ################
  
  fibroblast_markers <- c("FAP", "PDGFRB", "PDGFRA", "ACTA2", "COL1A1")
  mono_macro_markers <-c("MSR1", "CD163", "CCL3")
  chromaffin_markers <- c("PNMT", "TH", "DBH","CHGA", "CHGB")
  adrenocortical_markers <- c("STAR", "CYP11B1", "CYP11A1")
  endothelial_markers <- c("FLT1", "EPAS1")
  sustentacular_markers <- c("CDH19", "SOX10", "S100B", "VIM")
  lymphocyte_markers <- c("CD2", "CD3E" , "MS4A1")
  
  #################
  # Find Clusters #
  #################
  
  
  if (output_qc)
  {
    pdf(file.path(qc_out_dir,"pca_elbow_plot.pdf"))
    ElbowPlot(snrna_object, ndims = 50)
    dev.off()
  }
  
  seed = 42
  set.seed(seed)
  snrna_object <- RunUMAP(snrna_object, dims = 1:25, seed.use = seed)
  snrna_object <- FindNeighbors(snrna_object, dims = 1:25)
  snrna_object <- FindClusters(snrna_object, resolution = 1) 
  
  if (output_qc)
  {
    pdf(file.path(qc_out_dir,"cluster_qc.pdf"), width = 15, height = 10) #dev.off is after all plots
    
    clusters_dimplot_original <- DimPlot(snrna_object, group.by = 'seurat_clusters', label = TRUE) + NoLegend() +ggtitle('Seurat Cluster')
    samples_dimplot <- DimPlot(snrna_object, group.by = 'A5_ID', label = T) +ggtitle('Sample')
    chemistry_dimplot <- DimPlot(snrna_object, group.by = 'chemistry') + ggtitle('Chemistry')
    
    print(samples_dimplot +
            clusters_dimplot_original +
            chemistry_dimplot +
            plot_layout(ncol=2,nrow = 2))
    
    
    # check quality control metrics across the clusters 
    qc_violin <-  VlnPlot(snrna_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",  'doublet_score'), ncol = 4, pt.size=0, group.by = 'seurat_clusters')
    qc_featureplot <- FeaturePlot(snrna_object, features = c('percent.mt', 'nCount_RNA', 'doublet_score'))
    qc_featureplot2 <- FeaturePlot(snrna_object, features = c('nFeature_RNA')) + 
      scale_color_gradientn(values=scales::rescale(c(0,500,1000,5000,10000,20000,100000)),
                            colors = c("grey","red", "orange", "purple", "green", "blue"))
    
    cc_dimplot <-  DimPlot(snrna_object, group.by = 'Phase')
    
    print(qc_violin)
    print(qc_featureplot)
    print(qc_featureplot2)
    print(cc_dimplot)
    
    marker_genes_dp <- DotPlot(
      snrna_object,
      group.by = "seurat_clusters",
      features = rev(c(
        chromaffin_markers,
        adrenocortical_markers,
        sustentacular_markers,
        mono_macro_markers,
        endothelial_markers,
        lymphocyte_markers,
        fibroblast_markers))
    ) +
      RotatedAxis()
    
    
    print(marker_genes_dp)
    
    
  }
  
  #----
  # DGE for all clusters
  #----
  
  # to help annotate, identify top DEGs for each cluster 
  # cluster_markers <- FindAllMarkers(snrna_object,
  #                                   only.pos = TRUE,
  #                                   min.pct = 0.25,
  #                                   logfc.threshold = 0.25)
  
  #write_tsv(cluster_markers, file = "results/differential_gene_expression/cluster_markers.tsv")
  #cluster_markers <- read_tsv("results/differential_gene_expression/cluster_markers.tsv")
  # get top 10 significant genes for each cluster 
  # top_5 <- cluster_markers %>% 
  #   filter(p_val_adj < 0.05) %>% 
  #   group_by(cluster) %>% 
  #   slice_max(n = 5, order_by = avg_log2FC)
  # 
  # top_5 %>% filter(cluster == 27)
  
  #----
  # Annotate clusters
  #----
  
  # TODO: change this to just annotate tumour clusters as "Tumour"
  # can ten make another metadata column where they are annotated per sample of origin
  
  #annotate the clusters according to cell type:
  #based these annotations on the expression of marker genes (dotPlot) as well as scMatch predictions
  
  # all cell types -  c("Tumour", "Chromaffin cells", "Adrenocortical cells", "Endothelial cells", "Fibroblasts", "SCLCs", "Myeloid cells", "Lymphocytes")
  
  new.cluster.ids <- c(
    "Tumor", #0
    "Tumor", #1
    "Tumor", #2
    "Tumor", #3
    "Tumor", #4
    "Tumor", #5 
    "Tumor", #6
    "Tumor", #7
    "Tumor", #8
    "Tumor", #9
    "Myeloid cells", #10
    "Tumor", #11
    "Chromaffin cells", #12 
    "Tumor", #13
    "Tumor", #14
    "Adrenocortical cells", #15
    "SCLCs", #16
    "Tumor",#17
    "Fibroblasts", #18
    "Endothelial cells", #19
    "Tumor", #20
    "Adrenocortical cells", #21
    "Lymphocytes", #22
    "Tumor", #23
    "Tumor", #24
    "Endothelial cells", #25
    "Tumor", #26 
    "Fibroblasts", #27
    "Tumor", #28
    "SCLCs", #29
    "Lymphocytes", #30
    "Tumor", #31
    "Tumor", #32
    "Tumor" #33
    
  )
  
  if (digest::digest(snrna_object$seurat_clusters, algo="md5") != "36b6e32adda0b888a695546ad4e48cfb"){
    if (names(dev.cur()) != "null device") { dev.off() }
    stop("The checksum for cell to cluster assignment does not match the stored checksum. This may be due to a change in 
    package versions or input data. The cluster to cell-type assignments hardcoded in this function are only valid when 
    the checksum is true.  You will need to manually re-assess the cluster assignments and update the checksum")
  }
  
  Idents(snrna_object) <- snrna_object$seurat_clusters
  names(x = new.cluster.ids) <- levels(x = snrna_object)
  snrna_object <- RenameIdents(object = snrna_object, new.cluster.ids)
  snrna_object$cell_type <- Idents(snrna_object)
  
  if(output_qc){
    # clusters labelled according to scMatch results and these marker genes
    
    marker_genes_dp <- DotPlot(
      snrna_object,
      features = rev(c(
        chromaffin_markers,
        adrenocortical_markers,
        sustentacular_markers,
        mono_macro_markers,
        endothelial_markers,
        lymphocyte_markers,
        fibroblast_markers))
    ) +
      RotatedAxis()
    
    annotated_dimplot <- DimPlot(snrna_object, cols = cell_type_colours, label=T) #+ NoLegend()
    
    print((annotated_dimplot + samples_dimplot) /
            marker_genes_dp)
    dev.off() #close PDF opened further up
  }
  
  return(snrna_object)
  
}
