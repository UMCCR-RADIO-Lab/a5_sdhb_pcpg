python_style_zip <- function(...) {
  mapply(list, ..., SIMPLIFY = FALSE)
}

get_MANE_refseq_ids <- function()
{
  refseq_ids <- readr::read_delim(file = "/g/data/pq08/reference/GRCh38/mane/MANE.GRCh38.v1.0.refseq_genomic.gff.gz", 
                                  delim="\t", 
                                  comment = "#", 
                                  col_names =  c("seqname","source","feature","start",
                                                 "end","score","strand","frame","attribute"),
                                  col_types = "ccccccccc") %>%
    filter(feature=="mRNA") %>% 
    separate_rows(attribute,sep = ";") %>% 
    separate(col = attribute, 
             into=c("attribute_name", "attribute_value"),
             sep="=") %>% 
    filter(attribute_name == "transcript_id") %>%
    mutate(attribute_value=gsub("[.].+$","",attribute_value)) %>% 
    pull(attribute_value)
  return(refseq_ids)
}

get_region_probes <- function(gene_symbol=NULL, 
                              region= NULL,
                              probe_annotation,
                              available_probes){
  
  if (!is.null(gene_symbol) )
  {
    # entrez_id <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, 
    #                                    keys = gene_symbol, 
    #                                    columns="ENTREZID", 
    #                                    keytype = "SYMBOL")[["ENTREZID"]]
    # 
    # if(length(entrez_id) > 1 ) { stop ("multimapper") }
    # 
    # tx_annotation <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, keys=entrez_id,   
    #                                        columns=c("GENEID", "TXCHROM", "TXSTART", "TXEND", "TXNAME", "CDSNAME"), keytype="GENEID")
    # region_chr <- tx_annotation$TXCHROM[[1]]
    # region_min_coord <- min(tx_annotation$TXSTART)
    # region_max_coord <- max(tx_annotation$TXEND)
    
    region_probes <- probe_annotation %>% GenomicRanges::as.data.frame() %>% 
      filter(grepl(paste0(gene_symbol,";|",gene_symbol,"$"), UCSC_RefGene_Name))
    
    if (nrow(region_probes) == 0)
    {
      warning("No probe annotation found with gene symbol ", gene_symbol,". Returning NULL")
      return(NULL)
    }
    
  }
  
  if (!is.null(region) )
  {
    region_array <- stringr::str_split_1(region, pattern = "[:-]")
    region_chr <- ifelse(grepl("^chr", region_array[[1]]), region_array[[1]], paste0("chr", region_array[[1]]))
    region_min_coord <- as.numeric(min(region_array[2:3]))
    region_max_coord <- as.numeric(max(region_array[2:3]))
    
    region_probes <- probe_annotation %>% GenomicRanges::as.data.frame() %>% 
      mutate(pos=as.numeric(pos)) %>% 
      filter(chr==region_chr,
             pos >= region_min_coord,
             pos <= region_max_coord)
  }
  
  region_probes <- region_probes %>% 
    filter(Name %in% available_probes) 
  
  if (nrow(region_probes) == 0)
  {
    warning("No probes available for gene symbol ", gene_symbol,". Returning NULL")
    return(NULL)
  }
  
  if( !exists("mane_refseqid", envir = globalenv()))
  {
    mane_refseqid <- get_MANE_refseq_ids()
    assign("mane_refseqid", mane_refseqid, envir = globalenv())
  }
  
  if(all(region_probes$UCSC_RefGene_Name == ""))
  {
    region_probes$UCSC_RefGene_Name <- "Unannotated"
  }
  
  if(all(region_probes$UCSC_RefGene_Accession == ""))
  {
    region_probes$UCSC_RefGene_Accession <- "Unannotated"
  }
  
  if(all(region_probes$UCSC_RefGene_Group == ""))
  {
    region_probes$UCSC_RefGene_Group <- "Unannotated"
  }
  
  region_probes <- region_probes %>% 
    rowwise() %>% 
    mutate(UCSC_annotation = paste(
      lapply(
        python_style_zip(str_split_1(UCSC_RefGene_Name, pattern = ";" ), 
                         str_split_1(UCSC_RefGene_Accession, pattern = ";" ), 
                         str_split_1(UCSC_RefGene_Group, pattern = ";" )), 
        paste, collapse="|"), 
      collapse=";")) %>% 
    ungroup() %>% 
    separate_rows(UCSC_annotation,sep = ";") %>% 
    separate(UCSC_annotation, 
             into=c("UCSC_RefGene_Name","UCSC_RefGene_Accession","UCSC_RefGene_Group"), 
             sep = "[|]") 
  
  #Keep MANE where available else pick single transcript
  region_accessions <- unique(region_probes[,c("UCSC_RefGene_Name","UCSC_RefGene_Accession")]) %>% 
    arrange(UCSC_RefGene_Name, UCSC_RefGene_Accession) %>% 
    mutate(mane=UCSC_RefGene_Accession %in% mane_refseqid) %>% 
    group_by(UCSC_RefGene_Name) %>% mutate(any_mane=any(mane)) %>% 
    filter(mane | (!any_mane & row_number()==1)) %>% 
    filter(UCSC_RefGene_Name != "", UCSC_RefGene_Accession  != "")
  
  if(any(region_probes$UCSC_RefGene_Accession %in% mane_refseqid))
  {
    if (!is.null(gene_symbol)){
      region_probes <- region_probes %>% filter(UCSC_RefGene_Name == gene_symbol, UCSC_RefGene_Accession %in% UCSC_RefGene_Accession) 
    } else {
      region_probes <- region_probes %>% filter(UCSC_RefGene_Accession %in% region_accessions$UCSC_RefGene_Accession) 
    }
  }
  
  region_probes <- region_probes %>% dplyr::select(Name, chr, pos, UCSC_RefGene_Name, UCSC_RefGene_Accession, UCSC_RefGene_Group) %>% 
    group_by(Name, chr, pos, UCSC_RefGene_Name) %>%  
    summarise(UCSC_RefGene_Group=paste(unique(UCSC_RefGene_Group), collapse = "/"))
  
  return(region_probes)
}

summarize_probes <- function(plot.data, grouping_column, colour_column, shape_column) {
  
  #Group probes into promoter/body and annotate probe count
  plot.data.summarised <- plot.data %>% 
    mutate(UCSC_RefGene_Group = case_when(
      grepl("TSS|UTR|1stExon", UCSC_RefGene_Group) & grepl("Body", UCSC_RefGene_Group) ~ "Promoter/Body",
      grepl("TSS|UTR|1stExon", UCSC_RefGene_Group) ~ "Promoter",
      grepl("Body", UCSC_RefGene_Group) ~ "Body",
      TRUE ~ "Other/Unannotated"
    )) %>% 
    group_by(Sample_ID, UCSC_RefGene_Name, UCSC_RefGene_Group) %>% 
    mutate(UCSC_RefGene_Group_with_count = paste0(UCSC_RefGene_Group, "(n=", n(),")"),
           UCSC_RefGene_Group_count = paste0("n=", n()))
  
  #group probes by sample, probe annotation columns, and any plot feature columns (injected) 
  injectables <- c(sym(grouping_column))
  if(!is.null(colour_column)) { injectables <- c(injectables, sym(colour_column)) }
  if (!is.null(shape_column)){ injectables <- c(injectables, sym(shape_column)) }
  injectables <- unique(injectables)
  plot.data.summarised <- 
    rlang::inject(
      group_by(plot.data.summarised, 
               Sample_ID, 
               !!!injectables, 
               UCSC_RefGene_Name, UCSC_RefGene_Group, 
               UCSC_RefGene_Group_with_count, 
               UCSC_RefGene_Group_count)) 
  
  #summarize probe values into mean and sd
  plot.data.summarised <- plot.data.summarised %>% 
    summarise(probe_id="mean_all_probes", 
              pos=10^9, #Large number to ensure summaries sort to end
              b_sd=sd(b_val, na.rm = T),
              m_sd=sd(m_val, na.rm = T),
              b_val=mean(b_val, na.rm = T), 
              m_val=mean(m_val, na.rm = T))
  
  plot.data.summarised <- plot.data.summarised %>% 
    group_by(!!sym(grouping_column)) %>% 
    mutate(m_z=(m_val-mean(m_val, na.rm = T))/sd(m_val)) 
  
  return(plot.data.summarised)
}

plot_methylation <- function(gene_symbol=NULL, 
                             region=NULL, 
                             sample_annotation, 
                             grouping_column="group",
                             plot_title = NULL, 
                             show_mean_summarised_only=T, 
                             label_outliers=T,
                             outlier_z_threshold=2,
                             label_small_sample_sets=T,
                             small_set_cutoff = 5,
                             samples_to_label=c(), 
                             b_vals, 
                             m_vals,
                             array_type=c("EPIC","450K"),
                             plot_mode=c("beta","m"),
                             colour_column=NULL,
                             shape_column=NULL,
                             colour_scale=NULL,
                             shape_scale=NULL)
{
  
  
  
  
  
  if (is.null(gene_symbol) & is.null(region)) { stop("Gene symbol or region must be supplied") }
  
  if (!is.null(gene_symbol) & !is.null(region)) { stop("One of Gene symbol or region must be supplied, not both") }
  
  if (is.null(plot_title)) { plot_title = ifelse(is.null(region),  gene_symbol, region) }
  
  if (!(grouping_column %in% colnames(sample_annotation))) { 
    stop("No column named ", grouping_column, " in supplied sample annotation. Please specify an appropriate column name.")}
  
  if (!is.null(colour_column) && !(colour_column %in% colnames(sample_annotation))) { 
    stop("No column named ", colour_column, " in supplied sample annotation. Please specify an appropriate colour column name.")}
  
  if (!is.null(shape_column) && !(shape_column %in% colnames(sample_annotation))) { 
    stop("No column named ", shape_column, " in supplied sample annotation. Please specify an appropriate shape column name.")}
  
  if (!("Sample_ID" %in% colnames(sample_annotation))) { 
    stop("No column named Sample_ID in supplied sample annotation.")}
  
  array_type <- match.arg(toupper(array_type), c("450K","EPIC"))  
  plot_mode <- match.arg(tolower(plot_mode), c("beta","m"))
  
  if(array_type=="450K"){
    probe_annotation <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } else {
    probe_annotation <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }
  
  region_probes <- get_region_probes(gene_symbol=gene_symbol, 
                                     region= region,
                                     probe_annotation = probe_annotation,
                                     available_probes = rownames(m_vals))
  
  if(is.null(region_probes)) { return(NULL) }
  
  probe_data_b <- b_vals[region_probes$Name,] %>% data.frame(check.names = F) %>% 
    tibble::rownames_to_column("probe_id") %>% 
    pivot_longer(cols=-probe_id, names_to = "Sample_ID", values_to = "b_val")
  
  probe_data_m <- m_vals[region_probes$Name,] %>% data.frame(check.names = F) %>% 
    tibble::rownames_to_column("probe_id") %>% 
    pivot_longer(cols=-probe_id, names_to = "Sample_ID", values_to = "m_val")
  
  probe_data <-  probe_data_b %>% inner_join(probe_data_m)
  
  plot.data <- probe_data %>%  
    inner_join(region_probes %>% dplyr::select(probe_id=Name, chr, pos, UCSC_RefGene_Name, UCSC_RefGene_Group)) %>%  
    inner_join(sample_annotation)
  
  plot.data.summarised <- summarize_probes(plot.data, grouping_column, colour_column, shape_column)
  
  #Annotate samples for labelling
  plot.data.summarised <- plot.data.summarised %>% 
    group_by(!!sym(grouping_column)) %>% 
    mutate(label=ifelse((label_outliers & abs(m_z) > outlier_z_threshold) | 
                          (label_small_sample_sets & n() <= small_set_cutoff) | 
                          Sample_ID %in% samples_to_label , Sample_ID, NA))
  
  plot.data <- plot.data %>% dplyr::bind_rows(plot.data, plot.data.summarised)
  
  
  if(plot_mode == "beta") { 
    y_source = "b_val" 
  } else { 
    y_source = "m_val" 
  } 
  
  injectables <- list()
  if(!is.null(colour_column)) { injectables <- c(injectables, color=sym(colour_column)) }
  if (!is.null(shape_column)){ injectables <- c(injectables, shape=sym(shape_column)) }
  
  if(show_mean_summarised_only)
  {
    
    plot.data <- plot.data %>% filter(probe_id=="mean_all_probes")
    
    gg_meth <- ggplot(plot.data, 
                      aes(x=!!sym(grouping_column), 
                          y=!!sym(y_source))) + 
      geom_boxplot(outlier.alpha = 0)
    
    gg_meth <- gg_meth + 
      geom_jitter(mapping=rlang::inject(aes(!!!injectables)), 
                  width = 0.2)
    
    gg_meth <- gg_meth + 
      geom_text(x=0.7,y=0.1, aes(label=UCSC_RefGene_Group_count), size=3) +
      geom_text(aes(label=label), nudge_y = (max(plot.data[[y_source]])-min(plot.data[[y_source]]))*0.02, size=3) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1)) + 
      ggtitle(plot_title)
    
    if(length(unique(plot.data$UCSC_RefGene_Name)) > 1 ) {
      gg_meth <- gg_meth + facet_grid(UCSC_RefGene_Name~UCSC_RefGene_Group) 
    } else 
    {
      gg_meth <- gg_meth + facet_wrap("UCSC_RefGene_Group", ncol=3) 
    }
    
  } else {
    plot.data <- plot.data %>%  
      mutate(probe_id=paste0(probe_id,"(",UCSC_RefGene_Group,")")) %>% 
      arrange(pos) %>% mutate(probe_id = factor(probe_id, levels=unique(.$probe_id)))
    
    gg_meth <- ggplot(plot.data, 
                      aes(x=!!sym(grouping_column), 
                          y=!!sym(y_source))) + 
      geom_boxplot(outlier.alpha = 0)    
    
    gg_meth <- gg_meth + 
      geom_jitter(mapping=rlang::inject(aes(!!!injectables)), 
                  width = 0.2)
    
    gg_meth <- gg_meth + 
      geom_text(aes(label=label), nudge_y = (max(plot.data[[y_source]])-min(plot.data[[y_source]]))*0.02, size=3) +
      facet_wrap("probe_id") +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1)) + 
      ggtitle(plot_title)
  }
  
  if(plot_mode == "beta"){
    gg_meth <- gg_meth + coord_cartesian(ylim = c(0,1))
  }
  
  if(!is.null(colour_column) & !is.null(colour_scale)) { 
    gg_meth <- gg_meth + scale_color_manual(values = colour_scale)
  }
  
  if (!is.null(shape_column) & !is.null(shape_scale)){ 
    gg_meth <- gg_meth + scale_shape_manual(values = shape_scale)
  }
  
  return(gg_meth)
  
}


plot_methylation_vs_expr <- function(gene_symbol=NULL, 
                                     probes=NULL,
                                     sample_annotation,
                                     summarise_probes_by_region=T,
                                     grouping_column="group",
                                     plot_title = NULL, 
                                     samples_to_label=c(), 
                                     b_vals, 
                                     m_vals,
                                     log2_cpm,
                                     label_outliers=T,
                                     outlier_z_threshold=2,
                                     array_type=c("EPIC","450K"),
                                     plot_mode=c("beta","m"),
                                     colour_column=NULL,
                                     shape_column=NULL,
                                     colour_scale=NULL,
                                     shape_scale=NULL,
                                     export_plot_data = NULL)
{
  

  if (is.null(gene_symbol) & is.null(probes)) { stop("Gene symbol or probe list must be supplied") }

  if (is.null(plot_title)) { plot_title = gene_symbol }
    
  if (!(grouping_column %in% colnames(sample_annotation))) { 
    stop("No column named ", grouping_column, " in supplied sample annotation. Please specify an appropriate column name.")}
 
  if (!is.null(colour_column) && !(colour_column %in% colnames(sample_annotation))) { 
    stop("No column named ", colour_column, " in supplied sample annotation. Please specify an appropriate colour column name.")}
  
  if (!is.null(shape_column) && !(shape_column %in% colnames(sample_annotation))) { 
    stop("No column named ", shape_column, " in supplied sample annotation. Please specify an appropriate shape column name.")}
  
  if (!("Sample_ID" %in% colnames(sample_annotation))) { 
    stop("No column named Sample_ID in supplied sample annotation.")}
  
  if(stringr::str_detect(string = rownames(log2_cpm)[[1]], pattern="ENSG[0-9]+[.][0-9]+_.+"))
  {
    log2_cpm <- log2_cpm %>% data.frame(check.names = F) %>% tibble::rownames_to_column("ensgid_symbol") %>% 
      separate(ensgid_symbol, into=c("ensgid", "symbol"), extra = "merge", sep="_") %>%  dplyr::select(-ensgid)
  } else {
    if (rownames(log2_cpm)[[1]] != "1") {
    log2_cpm <- log2_cpm %>% data.frame(check.names = F) %>% tibble::rownames_to_column("symbol") 
    } else {
      stop("log2_cpm must be supplied as a matrix or dataframe with gene symbols as rownames")
    }
  }
  
  array_type <- match.arg(toupper(array_type), c("450K","EPIC"))  
  plot_mode <- match.arg(tolower(plot_mode), c("beta","m"))
  
  if(array_type=="450K"){
    probe_annotation <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } else {
    probe_annotation <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }
  
  if (is.null(probes))
  {
  region_probes <- get_region_probes(gene_symbol=gene_symbol, 
                                     probe_annotation = probe_annotation,
                                     available_probes = rownames(m_vals))
  } else {
    region_probes = data.frame(probe_annotation[probe_annotation$Name %in% probes, c("Name", "chr", "pos", "UCSC_RefGene_Name", "UCSC_RefGene_Group")])
    
    region_probes = region_probes %>% 
      tidyr::separate_rows(UCSC_RefGene_Name, UCSC_RefGene_Group, sep=";") %>% 
      distinct()
    
    if(length(unique(region_probes$UCSC_RefGene_Name)) > 1 & is.null(gene_symbol)) {
      stop("Supplied probes map to more than one gene. Unable to fetch gene expression. Please supply gene symbol.")
    }
    
    if(all(unique(region_probes$UCSC_RefGene_Name) == "")) {
      stop("Supplied probes do not map to any gene symbol. Unable to fetch gene expression.")
    }
    
    if(is.null(gene_symbol)){ 
      gene_symbol <- unique(region_probes$UCSC_RefGene_Name)
    }
    
  }
  
  if(is.null(region_probes)) { return(NULL) }
  
  probe_data_b <- b_vals[region_probes$Name, ,drop=F] %>% data.frame(check.names = F) %>% 
    tibble::rownames_to_column("probe_id") %>% 
    pivot_longer(cols=-probe_id, names_to = "Sample_ID", values_to = "b_val")
  
  probe_data_m <- m_vals[region_probes$Name,,drop=F] %>% data.frame(check.names = F) %>% 
    tibble::rownames_to_column("probe_id") %>% 
    pivot_longer(cols=-probe_id, names_to = "Sample_ID", values_to = "m_val")
  
  probe_data <-  probe_data_b %>% inner_join(probe_data_m)
  
  if (gene_symbol %in% log2_cpm$symbol)
  {
    expr_data <- log2_cpm %>%  filter(symbol==gene_symbol) %>% 
      pivot_longer(cols=-symbol, names_to = "Sample_ID", values_to = "log2cpm")
    
  } else {
    warning("Gene symbol ", gene_symbol, " not found in log2cpm object. Returning NULL")
    return(NULL)
  }
  
  plot.data <- probe_data %>%  
  inner_join(region_probes %>% dplyr::select(probe_id=Name, chr, pos, UCSC_RefGene_Name, UCSC_RefGene_Group)) %>%  
  inner_join(sample_annotation)

  if(summarise_probes_by_region)
  {
    
    plot.data <- summarize_probes(plot.data, grouping_column, colour_column, shape_column)
    
    plot.data <- plot.data %>% inner_join(expr_data, by=c("UCSC_RefGene_Name"="symbol","Sample_ID"="Sample_ID"))
    
    #group samples for labeling
    plot.data <- plot.data %>% group_by(UCSC_RefGene_Name, !!sym(grouping_column))
      
  } else {
    
    plot.data <- plot.data %>% inner_join(expr_data, by=c("UCSC_RefGene_Name"="symbol","Sample_ID"="Sample_ID"))
    
    #group samples for labeling
    plot.data <- plot.data %>% group_by(probe_id, !!sym(grouping_column)) 
    
  }
  
  #Annotated samples for outlier labeling based on applied groups
  plot.data <- plot.data %>% mutate(m_z=(m_val-mean(m_val, na.rm = T))/sd(m_val),
         log2cpm_z=(log2cpm-mean(log2cpm))/sd(log2cpm),
         label=ifelse((label_outliers & (abs(m_z) > outlier_z_threshold | abs(log2cpm_z) > outlier_z_threshold)) | 
                        Sample_ID %in% samples_to_label , Sample_ID, NA))
  
  if(!is.null(export_plot_data)){
    assign(x = export_plot_data, value = plot.data, envir = globalenv())
    message("Exported plot data to:", export_plot_data)
  }
  
  if(plot_mode == "beta") { 
    y_source = "b_val" 
  } else { 
    y_source = "m_val" 
  } 
  
  injectables <- list()
  if(!is.null(colour_column)) { injectables <- c(injectables, color=sym(colour_column)) }
  if (!is.null(shape_column)){ injectables <- c(injectables, shape=sym(shape_column)) }
  
  if(summarise_probes_by_region)
  {
    
    gg_meth <- ggplot(plot.data, aes(y=log2cpm, x=!!sym(y_source))) 
    
    gg_meth <- gg_meth + 
      geom_point(mapping=rlang::inject(aes(!!!injectables)))
    
    gg_meth <- gg_meth + 
      geom_text(y=max(plot.data$log2cpm)-(0.02*max(plot.data$log2cpm)),
                x=0.5, 
                aes(label=UCSC_RefGene_Group_count)) +
      geom_text(aes(label=label), nudge_y = (max(plot.data[["log2cpm"]])-min(plot.data[["log2cpm"]]))*0.02, size=3) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1)) + 
      ggtitle(plot_title)  +
      facet_wrap("UCSC_RefGene_Group", ncol=3) 
    
  } else {
    
    plot.data <- plot.data %>%  
      mutate(probe_id=paste0(probe_id,"(",UCSC_RefGene_Group,")")) %>% 
      arrange(pos) %>% mutate(probe_id = factor(probe_id, levels=unique(.$probe_id)))
    
    gg_meth <- ggplot(plot.data, aes(y=log2cpm, x=!!sym(y_source))) 
    
    gg_meth <- gg_meth + 
      geom_point(mapping=rlang::inject(aes(!!!injectables)))
    
    gg_meth <- gg_meth + 
      geom_text(aes(label=label), nudge_y = (max(plot.data[["log2cpm"]])-min(plot.data[["log2cpm"]]))*0.02, size=3) +
      theme_bw() + 
      theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1)) + 
      ggtitle(plot_title)  +
      facet_wrap("probe_id", ncol=5) 
    
  }
  
  
  if(!is.null(colour_column) & !is.null(colour_scale)) { 
    gg_meth <- gg_meth + scale_color_manual(values = colour_scale)
  }
  
  if (!is.null(shape_column) & !is.null(shape_scale)){ 
    gg_meth <- gg_meth + scale_shape_manual(values = shape_scale)
  }
  
  if(plot_mode == "beta"){
    gg_meth <- gg_meth + coord_cartesian(xlim = c(0,1))
  }
  
  return(gg_meth)
  
}

