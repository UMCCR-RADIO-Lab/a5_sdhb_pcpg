library(tidyverse)
library(limma)
library(edgeR)
library(googlesheets4)
library(biomaRt)
library(org.Hs.eg.db)
library(SummarizedExperiment)

#Function to read and merge transcript count files from htseq 
data_loader_a5_wts_counts <- function(count_file_dir=NULL, count_file_pattern=".gene.counts")
{  
  message("Starting A5 WTS data loader ...")
  
  #####################
  # Sample annotation #
  #####################
  
  if(!exists("a5_anno")) {
    source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
    data_loader_a5_clinical_anno(use_cache = T) }
  
  if(!("differential_group" %in% names(a5_anno))) {
    a5_anno <- add_differential_groups(a5_anno) }
  
  ###############
  # Load Counts #
  ###############
  
  # read in the counts for each sample
  message("Reading in count files ...")
  counts_df <- list.files(count_file_dir, full.names = T, pattern = count_file_pattern) %>% 
    map(function(x){
      sample_name <- sub(pattern = count_file_pattern, replacement = "", basename(x))
      message("Processing: ", sample_name)
      read_tsv(x,col_types = "ccd") %>%
        dplyr::rename(!!sample_name := count)
    }) %>% 
    purrr::reduce(left_join, by = c("gene_id", "gene_symbol")) %>% 
    mutate(gene = paste0(gene_id, "_", gene_symbol)) %>% 
    # remove the rows without gene symbols (QC rows from htseq)
    dplyr::filter(!is.na(gene_symbol)) %>% 
    dplyr::select(-gene_id, -gene_symbol) %>% 
    column_to_rownames(var = "gene")
  
  ###############
  # Sample Sets #
  ###############
  
  #######
  # Excluded samples
  #######
  
  #Samples to exclude from core set for quality/purity/genotype reasons
  exclude_global <- a5_anno %>%  filter(Exclude=="Y") %>% pull(A5_ID)
  exclude_missing_anno <- setdiff(colnames(counts_df), a5_anno$A5_ID)
  exclude_genotype <- a5_anno %>%  filter(Genotype != "SDHB") %>% pull(A5_ID)
  exclude_no_rna_data <- setdiff(a5_anno$A5_ID, colnames(counts_df))
  exclude_no_wgs_data <- c("E181-1", "E191-1")
  exclude_qc <- c("E144-1", "E159-3")

  exclude_base <- c(exclude_no_rna_data, exclude_no_wgs_data, exclude_qc, exclude_global, exclude_missing_anno)
  

  a5_anno <- a5_anno %>% mutate(Exclude_RNA=ifelse(A5_ID %in% exclude_base, T, F))
  
  #######
  # Anatomy groups
  #######
  
  #Select head/neck
  samples_hn <- a5_anno %>% 
    filter(differential_group_anatomy == "Head_neck") %>% 
    pull(A5_ID)
  
  #Select abdo_thoracic tumours
  samples_abdothoracic <- a5_anno %>% 
    filter(differential_group_anatomy == "Abdominal_Thoracic") %>% 
    pull(A5_ID)
  
  #Select samples of unclear origin tumours based on UMAP
  samples_thoracic_non_chromaffin <- a5_anno %>% 
    filter(differential_group_anatomy == "Thoracic_non_chromaffin") %>% 
    pull(A5_ID)
  
  #Select samples of unclear origin tumours based on UMAP
  samples_ambiguous <- a5_anno %>% 
    filter(differential_group_anatomy == "Ambiguous") %>% 
    pull(A5_ID)
  
  #######
  # Create sets
  #######
  
  counts_df_list <- list()
  
  #All data
  counts_df_list[["all"]] <- counts_df
  
  #Poor QC and samples with no WGS excluded
  counts_df_list[["qc_ok"]] <- counts_df[,!(colnames(counts_df) %in% exclude_base)]
  
  #Poor QC, samples with no WGS, and NF/VHL samples excluded
  counts_df_list[["SDHB"]] <- counts_df[,!(colnames(counts_df) %in% c(exclude_base, exclude_genotype))]
  
  #Head and neck based on clinical and UMAP annotation (aortic cases excluded)
  samples_hn_keep <- setdiff(samples_hn, c(exclude_base, exclude_genotype))
  counts_df_list[["SDHB_HN"]] <- counts_df[,colnames(counts_df) %in% samples_hn_keep]
  
  #Abdominal/Thoracic based on clinical and UMAP annotation (ambiguous cases excluded)
  samples_abdothoracic_keep <- setdiff(samples_abdothoracic, c(exclude_base, exclude_genotype))
  counts_df_list[["SDHB_abdothoracic"]] <- counts_df[,colnames(counts_df) %in% samples_abdothoracic_keep]
  message("Created data subsets:", toString(names(counts_df_list)))
  
  ##############
  # Preprocess #
  ##############
  
  #######
  # Create DGE objects
  #######
  
  message("Converting counts to DGE objects ...")
  counts_dge_list <- lapply(counts_df_list, DGEList)
  
  # filter genes with low expression in all/nearly all samples
  message("Filtering lowly expressed genes ...")
  counts_dge_list <- map(.x = counts_dge_list, 
                          .f = function (dge_object)  { 
                            diff_groups <- a5_anno$differential_group[match(colnames(dge_object), a5_anno$A5_ID)]
                            expr_genes <- filterByExpr(dge_object, group = diff_groups)  
                            dge_object[expr_genes,, keep.lib.sizes = FALSE] })
  
  #######
  # TMM normalisation
  #######
  
  message("Calculating normalisation factors with edgeR ...")
  counts_dge_list <- map(counts_dge_list, calcNormFactors)
  
  ####################
  # Export to Global #
  ####################
  
  assign(x = "a5_wts_counts", value = counts_df, envir = globalenv())
  assign(x = "a5_wts_dge_list", value = counts_dge_list, envir = globalenv())
  assign(x = "a5_anno", value = a5_anno, envir = globalenv())
  
  message("Loaded objects into global environment: 
          a5_wts_counts (dataframe of raw feature counts),
          a5_wts_dge_list (a list of count normalised DGE objects named:
                          - all: All profiled samples
                          - qc_ok: Excluded samples removed (",toString(exclude_base),")
                          - SDHB: QC_OK with non-SDHB samples removed (",toString(exclude_genotype),")
                          - SDHB_abdothoracic: Abdominal/thoracic based on clinical and UMAP annotation (thoracic_non_chromaffin and ambiguous cases excluded:",toString(samples_ambiguous),toString(samples_thoracic_non_chromaffin),")
                          - SDHB_HN: Head and neck based on clinical and UMAP annotation (thoracic_non_chromaffin and ambiguous cases excluded)"
          )
  
  message("Completed A5 WTS data loader.")
}
message("Created data loader function data_loader_a5_wts_counts()")

#Function to fetch the biotype of a gene from ensembl using biomart
biotype_from_biomart <- function(EnsIds, 
                                 use_cache=FALSE, 
                                 update_cache=TRUE, 
                                 ensembl_mirror="https://asia.ensembl.org",
                                 offline_cache="/g/data/pq08/projects/ppgl/a5/offline_cache/biotype_from_biomart.tsv")
{
  if(use_cache)
  {
    message("Loading biotype data from offline cache...")
    return_data <- read.delim(offline_cache, sep="\t", header=T, check.names = F)
    if (any(!(as.character(na.omit(EnsIds)) %in% as.character(na.omit(return_data$ensembl_gene_id))))) {
      missing_genes <- setdiff(EnsIds, return_data$ensembl_gene_id)
      message_genes <- missing_genes
      if (length(missing_genes) > 200) { message_genes <- message_genes[1:200]; message_genes[201] <- paste("+", length(missing_genes)-200, "genes")}
      message("WARNING: Not all requested gene IDs are present in the offline cache. This may be because the cache was built with a different 
              query set or because not all IDs return a result. Rerun in online mode to attempt to update. Missing genes:", toString(message_genes))
    }
    return_data <- data.frame(ensembl_gene_id=EnsIds) %>% left_join(return_data)
    return(return_data)
    
  } else {
  message("Fetching biotypes from biomart")
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = ensembl_mirror))
  
  biotype_type <- getBM(filters= "ensembl_gene_id",
                        attributes= c("ensembl_gene_id",
                                      "gene_biotype"),
                        values=EnsIds,mart= mart)
  
  
  
  return_data <- data.frame(ensembl_gene_id=EnsIds) %>% left_join(biotype_type)
  
  if(update_cache) {
    write.table(return_data, offline_cache, sep="\t", row.names = F)
  }
  
  }
  
  return(return_data)
}
message("Created helper function biotype_from_biomart()")

#Function to fetch the biotype of a gene from org.Hs.eg.db
biotype_from_orghs <- function(EnsIds)
{
  message("Fetching biotypes from org.Hs.eg")
  ens_org_list <- as.list(org.Hs.egENSEMBL[mappedkeys(org.Hs.egENSEMBL)])
  ens_org_df <- data.frame(entrez_id=names(ens_org_list), ensg_id=unlist(lapply(ens_org_list,toString)))
  ens_org_df <- ens_org_df %>% separate_rows(ensg_id, sep=", ")
  
  genetype_org_list <- as.list(org.Hs.egGENETYPE[mappedkeys(org.Hs.egGENETYPE)]); 
  genetype_org_df <- data.frame(entrez_id=names(genetype_org_list), bio_type=unlist(lapply(genetype_org_list,toString)))
  
  ens_genetype <- data.frame(ensembl_gene_id=EnsIds)  %>% 
    left_join(ens_org_df, by=c("ensembl_gene_id"="ensg_id")) %>%
    left_join(genetype_org_df) %>% dplyr::select(-entrez_id)
  
  return(ens_genetype)
}
message("Created helper function biotype_from_orghs()")

#Function to fetch the biotype of a gene from ensembl, falling back to org.Hs.eg.db if using biomart fails
data_loader_get_biotypes <- function(EnsIds, use_cache=FALSE)
{
  ensid_to_biotype <- NULL  
tryCatch(
    expr = 
      {
        ensid_to_biotype <- biotype_from_biomart(EnsIds, use_cache=use_cache)
      },
    error = function(cond){ 
      message("Problem with biomart, failing over to org.Hs.eg - ", cond)
      }
  )
  if(is.null(ensid_to_biotype))
  {
    ensid_to_biotype <- biotype_from_orghs(EnsIds) %>% rename(gene_biotype=bio_type)
  }
  ensid_to_biotype <- ensid_to_biotype %>% distinct() %>% mutate(gene_biotype=replace_na(gene_biotype,"Not available"))

assign(x = "ensid_to_biotype", value = ensid_to_biotype, envir = globalenv())
message("Loaded EnsemblID to gene-biotype mapping table into the global environment as ensid_to_biotype")
}
message("Created helper function data_loader_get_biotypes()")

data_loader_ensgid_to_chr <- function(EnsIds, 
                                      use_cache=FALSE, 
                                      update_cache=TRUE, 
                                      ensembl_mirror="https://asia.ensembl.org",
                                      offline_cache="/g/data/pq08/projects/ppgl/a5/offline_cache/ensgid_to_chr_from_biomart.tsv")
{
  if(use_cache)
  {
    message("Loading ensembl_id to chromosome conversion table from offline cache...")
    ensgid_to_chr <- read.delim(offline_cache, sep="\t", header=T, check.names = F)
    if (any(!(as.character(na.omit(EnsIds)) %in% as.character(na.omit(ensgid_to_chr$ensembl_gene_id))))) {
      missing_genes <- setdiff(EnsIds, ensgid_to_chr$ensembl_gene_id)
      message_genes <- missing_genes
      if (length(missing_genes) > 200) { message_genes <- message_genes[1:200]; message_genes[201] <- paste("+", length(missing_genes)-200, "genes")}
      message("WARNING: Not all requested gene IDs are present in the offline cache. This may be because the cache was built with a different 
              query set or because not all IDs return a result. Rerun in online mode to attempt to update. Missing genes:", toString(message_genes))
    }
    ensgid_to_chr <- data.frame(ensembl_gene_id=EnsIds) %>% left_join(ensgid_to_chr)
  } else {
    message("Fetching chromosomes from biomart")
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = ensembl_mirror, ))
    
    ensgid_to_chr <- getBM(filters= "ensembl_gene_id",
                          attributes= c("ensembl_gene_id",
                                        "chromosome_name",
                                        "start_position",
                                        "end_position"),
                          values=EnsIds,mart= mart)
    
    if(nrow(ensgid_to_chr) == 0) { 
      message("Biomart query returned no results. Check Ensembl Gene ID formatting.")
      return(NULL)
      }
    
    ensgid_to_chr <- data.frame(ensembl_gene_id=EnsIds) %>% left_join(ensgid_to_chr)
    
    if(update_cache) {
    write.table(ensgid_to_chr, offline_cache, sep="\t", row.names = F)
    }
    
  }
  assign(x = "ensgid_to_chr", value = ensgid_to_chr, envir = globalenv())
  message("Loaded EnsemblID to chr mapping table into the global environment as ensgid_to_chr")
}
message("Created helper function data_loader_ensgid_to_chr()")

#Function to fetch the ensembl gene ID from hgnc symbol using biomart
hgnc_to_ensgid_from_biomart <- function(hgnc_symbols, 
                                        ensembl_mirror="https://asia.ensembl.org",
                                        use_cache=FALSE, 
                                        update_cache=TRUE,
                                        offline_cache="/g/data/pq08/projects/ppgl/a5/offline_cache/hgnc_to_ensgid_from_biomart.tsv")
{
  if(use_cache)
  {
    message("Loading HGNC to ensembl_id conversion table from offline cache...")
    hgnc_to_ensgid <- read.delim(offline_cache, sep="\t", header=T, check.names = F)
    if (any(!(as.character(na.omit(hgnc_symbols)) %in% as.character(na.omit(hgnc_to_ensgid$hgnc_symbol))))) {
      missing_genes <- setdiff(hgnc_symbols, hgnc_to_ensgid$hgnc_symbol)
      message_genes <- missing_genes
      if (length(missing_genes) > 200) { message_genes <- message_genes[1:200]; message_genes[201] <- paste("+", length(missing_genes)-200, "genes")}
      message("WARNING: Not all requested gene IDs are present in the offline cache. This may be because the cache was built with a different 
              query set or because not all IDs return a result. Rerun in online mode to attempt to update. Missing genes:", toString(message_genes))
    }
    hgnc_to_ensgid <- hgnc_to_ensgid %>% filter(hgnc_symbol %in% hgnc_symbols)
  } else {
    message("Fetching ensembl_ids from biomart")
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = ensembl_mirror))
    
    hgnc_to_ensgid <- getBM(filters= "external_gene_name",
                            attributes= c("ensembl_gene_id",
                                          "external_gene_source",
                                          "external_gene_name"),
                            values=hgnc_symbols,mart= mart)
    
    hgnc_to_ensgid <- hgnc_to_ensgid %>% 
      filter(external_gene_source=="HGNC Symbol") %>% 
      dplyr::select(-external_gene_source) %>% 
      dplyr::rename("hgnc_symbol"="external_gene_name")
    
    if(update_cache) {
      write.table(hgnc_to_ensgid, offline_cache, sep="\t", row.names = F)
    }
  }
  return_data <- data.frame(hgnc_symbol=hgnc_symbols) %>% 
    left_join(hgnc_to_ensgid, by=c("hgnc_symbol")) 
  return(return_data)
}
message("Created helper function hgnc_to_ensgid_from_biomart()")

#Function to fetch the ensembl gene ID from hgnc symbol using biomart
ensgid_to_hgnc_from_biomart <- function(ens_gids, 
                                        ensembl_mirror="https://asia.ensembl.org",
                                        use_cache=FALSE, 
                                        update_cache=TRUE,
                                        offline_cache="/g/data/pq08/projects/ppgl/a5/offline_cache/ensgid_to_hgnc_from_biomart.tsv")
{
  if(use_cache)
  {
    message("Loading ensembl_id to HGNC conversion table from offline cache...")
    ensgid_to_hgnc <- read.delim(offline_cache, sep="\t", header=T, check.names = F)
    if (any(!(as.character(na.omit(ens_gids)) %in% as.character(na.omit(ensgid_to_hgnc$ensembl_gene_id))))) {
      missing_genes <- setdiff(ens_gids, ensgid_to_hgnc$ensembl_gene_id)
      message_genes <- missing_genes
      if (length(missing_genes) > 200) { message_genes <- message_genes[1:200]; message_genes[201] <- paste("+", length(missing_genes)-200, "genes")}
      message("WARNING: Not all requested gene IDs are present in the offline cache. This may be because the cache was built with a different 
              query set or because not all IDs return a result. Rerun in online mode to attempt to update. Missing genes:", toString(message_genes))
    }
    ensgid_to_hgnc <- ensgid_to_hgnc %>% filter(ensembl_gene_id %in% ens_gids)
  } else {
    message("Fetching ensembl_ids from biomart")
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = ensembl_mirror))
    
    ensgid_to_hgnc <- getBM(filters= "ensembl_gene_id",
                            attributes= c("ensembl_gene_id",
                                          "external_gene_name",
                                          "external_gene_source"),
                            values=ens_gids,mart= mart)
    
    ensgid_to_hgnc <- ensgid_to_hgnc %>% 
      filter(external_gene_source=="HGNC Symbol") %>% 
      dplyr::select(-external_gene_source) %>% 
      dplyr::rename("hgnc_symbol"="external_gene_name")
    
    if(update_cache) {
      write.table(ensgid_to_hgnc, offline_cache, sep="\t", row.names = F)
    }
  }
  return_data <- data.frame(ensembl_gene_id=ens_gids) %>% 
    left_join(ensgid_to_hgnc) 
  return(return_data)
}
message("Created helper function ensgid_to_hgnc_from_biomart()")


#Function to fetch the gene ontology terms associated with and ensemble gene ID using biomart
ensgid_to_goterm_from_biomart <- function(ens_gids, 
                                        ensembl_mirror="https://asia.ensembl.org",
                                        use_cache=TRUE, 
                                        update_cache=FALSE,
                                        offline_cache="/g/data/pq08/projects/ppgl/a5/offline_cache/ensgid_to_gene_ontology_from_biomart.tsv",
                                        permitted_evidence_codes=c("EXP","IDA","IPI","IMP","IGI","IEP",
                                                                   "HTP","HDA","HMP","HGI","HEP",
                                                                   "IBA","IBD","IKR","IRD","ISS","ISO","ISA","ISM","IGC",
                                                                   "RCA","TAS","NAS"), #see https://geneontology.org/docs/guide-go-evidence-codes/
                                        permitted_domains=c("cellular_component","biological_process","molecular_function")
                                        )
{
  if(use_cache)
  {
    message("Loading gene ontology table from offline cache...")
    ensgid_to_go <- read.delim(offline_cache, sep="\t", header=T, check.names = F)
    if (any(!(as.character(na.omit(ens_gids)) %in% as.character(na.omit(ensgid_to_go$ensembl_gene_id))))) {
      missing_genes <- setdiff(ens_gids, ensgid_to_go$ensembl_gene_id)
      message_genes <- missing_genes
      if (length(missing_genes) > 200) { message_genes <- message_genes[1:200]; message_genes[201] <- paste("+", length(missing_genes)-200, "genes")}
      message("WARNING: Not all requested gene IDs are present in the offline cache. This may be because the cache was built with a different 
              query set or because not all IDs return a result. Rerun in online mode to attempt to update. Missing genes:", toString(message_genes))
    }
    
    ensgid_to_go <- ensgid_to_go %>% filter(ensembl_gene_id %in% ens_gids)
    
  } else {
    message("Fetching gene ontology terms from biomart")
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl",host = ensembl_mirror))
    
    mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "asia")
    
    ensgid_to_go <- biomaRt::getBM(filters= "ensembl_gene_id",
                               attributes= c("ensembl_gene_id",
                                             "go_id", 
                                             "name_1006",
                                             "definition_1006",
                                             "go_linkage_type",
                                             "namespace_1003"),
                               values=ens_gids, mart = mart)
    

    if(update_cache) {
      write.table(ensgid_to_go, offline_cache, sep="\t", row.names = F)
    }
  }
  
  ensgid_to_go <- ensgid_to_go %>% 
    filter(go_linkage_type %in% permitted_evidence_codes, 
           namespace_1003 %in% permitted_domains)
  
  return_data <- data.frame(ensembl_gene_id=ens_gids) %>% 
    left_join(ensgid_to_go) 
  
  return(return_data)
}
message("Created helper function ensgid_to_goterm_from_biomart()")

###############
# RNA Fusions #
###############

data_loader_arriba <- function(arriba_out_dir=NULL, sample_to_exclude=vector(mode="character"))
{
  arriba_file_pattern = "_fusions.tsv"
  pcawg_drivers <- read.delim("/g/data/pq08/reference/gene_lists/pcawg_mutational_drivers.txt")
  
  message("Loading Arriba RNA Fusion Data...")
  
  arriba_files <- list.files(arriba_out_dir, pattern = arriba_file_pattern, recursive = T, full.names = T)
  names(arriba_files) <- gsub(".+/(E[0-9]{3}-[0-9])_fusions.tsv","\\1", arriba_files) 
  a5_arriba <- purrr:::map(arriba_files, ~readr:::read_delim(.x, show_col_types = F)) %>% bind_rows(.id="A5_ID") %>% dplyr::rename(`gene1`=`#gene1`)
  
  a5_arriba <- a5_arriba %>% filter(!(A5_ID %in% sample_to_exclude))
  a5_arriba_keep_pcawg <- a5_arriba %>% filter((gene1 %in% pcawg_drivers$Gene | gene2 %in% pcawg_drivers$Gene), confidence=="high", (split_reads1 >2 | split_reads2 >2 | discordant_mates >2))
  a5_arriba_keep_pcawg_brief <-  
    bind_rows(
      a5_arriba_keep_pcawg  %>%
        mutate(Annotation=paste0("RNA-Fusion:",gene2,"|",type ,"|", confidence, " conf","|", reading_frame)) %>% 
        dplyr::select(A5_ID, gene1, Annotation) %>% 
        distinct() %>% 
        dplyr::rename(Gene=gene1), 
      a5_arriba_keep_pcawg  %>%
        mutate(Annotation=paste0("RNA-Fusion:",gene1,"|",type ,"|", confidence, " conf","|", reading_frame)) %>% 
        dplyr::select(A5_ID, gene2, Annotation) %>% distinct() %>% dplyr::rename(Gene=gene2)
    ) %>% 
    separate_rows(Gene, sep=",") %>% 
    mutate(Gene=gsub("\\([0-9]+\\)","",Gene)) %>% distinct()
  
  assign("a5_arriba", a5_arriba, envir = globalenv()) 
  assign("a5_arriba_keep_pcawg", a5_arriba_keep_pcawg, envir = globalenv()) 
  assign("a5_arriba_keep_pcawg_brief", a5_arriba_keep_pcawg_brief, envir = globalenv()) 
  
  message("Created objects a5_arriba, a5_arriba_keep_pcawg (filtered to PCAWG CG list), a5_arriba_keep_pcawg_brief (summarised version)")
  
}
message("Created data loader function data_loader_arriba()")
