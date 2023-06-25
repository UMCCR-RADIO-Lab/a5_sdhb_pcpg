library(googlesheets4)
library(biomaRt)
library(parallel)
library(tidyr)
library(dplyr)

basedir="/g/data/pq08/projects/ppgl/a5"

source(paste0(basedir, "/wgs/scripts/data_loaders/small_variant_dataloader.r"))
source(paste0(basedir, "/wgs/scripts/data_loaders/manta_vcf_dataloader.r"))
source(paste0(basedir, "/wgs/scripts/data_loaders/purple_gene_cn_dataloader.r"))
source(paste0(basedir, "/wgs/scripts/data_loaders/linx_output_dataloader.r"))
source(paste0(basedir, "/wgs/scripts/data_loaders/cna_event_labeller.r"))
source(paste0(basedir, "/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r"))

if (!exists("pcawg_drivers_hgnc_syn"))
{
  pcawg_drivers <- read.delim("/g/data/pq08/reference/gene_lists/pcawg_mutational_drivers.txt")
  
  #message("Loading Gene Symbol Conversions from Ensembl...")
  
  # mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  # pcawg_drivers_hgnc_syn <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol", "external_synonym","chromosome_name"),values=gsub("[.][0-9][0-9]?$","", pcawg_drivers$Ensembl),mart= mart)
  # pcawg_drivers_hgnc_syn <- pcawg_drivers_hgnc_syn %>% pivot_longer(cols=c(hgnc_symbol, external_synonym), names_to="SymbolType",values_to="Symbol") %>% distinct() %>% mutate(chromosome_name=paste0("chr",chromosome_name))
}

if (!exists("a5_anno"))
{
  data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = TRUE)
}

####################
# Excluded Samples #
####################
exclude_samples <- a5_anno %>% dplyr::filter(Exclude=="Y") %>% pull(A5_ID) #c("E181-1")
#"E124-1", "E145-1",

#####################
# Curated Gene List #
#####################

shiva_curated_genes.variant <- c("AKAP12","ATRX","ANXA6", "ARID4B", "AXL","BCR","BNC2","CBX7","CDK10","CDKN2A","CEP152","CFTR","CHST9",
                                 "CRYBG1","CXCL14","DACT2","DICER1","DNMT3A","DPH1","ELF3","EP300", "EPAS1","EYA4","FAT1","FAT4",
                                 "FOXP1","FOXQ1","GJB2","GLI2","GLI3","GNAS","GPRC5A","GRB7","HDAC6","IGF1R","KDR","KLF4",
                                 "KLK10","KMT2A","KMT2B","KMT2C","KMT2D", "KMT2E", "KRAS","LRP1B","LSM1","MALT1","MAP2K4","MCM7","MDM2","MECP2", "MYCN",
                                 "MMP7","MYO18B","NF1","NFATC1","NKX2-5","NOTCH2","NRAS","PCDH17","PIM3","PML","PPP1R12B",
                                 "PPP1R13L","PRRT2","PRUNE2","PTCH1","PTK7","PTPRG","RASA1","RBM45","RET","SASH1","SATB2","SETDB1",
                                 "SMAD4","SMC1A","SMO","SOX3","SOX6","SOX7","TERT","TP53","TP63","TTF1","UNC5D","VGLL4","VWA5A")

shiva_curated_genes.sv <- c("SP140","BARD1","NFE2L2","TAF1D","AKAP9","ARID4B","ASXL2","BCL2L11","BRCA2","CACNA1D","CDKN2A",
                            "CREBBP","DNMT3B","E2F3","EFTUD2","EPAS1","EPHA5","EPHB1","ERCC3","FBXO11","FLT1","FOXP1","HGF",
                            "KIF3C","KMT2C ","LRP1B","MAP2K4","NSD1","PM20D1","PREX2","PRUNE1","PTPRD","RNF14","SLC6A1","SLX4",
                            "SMYD3","SPRY4","SPTBN2","TTF1","YES1","ZBTB20","ZNF33A")


##########################
# Germline PCGR Variants #
##########################

data_loader_germline_variants <- function(threads = 1)
{
  message("Loading Germline variant data...")
  
  #germline_cpsr_vcf_dir=paste0(basedir, "/wgs/symlinks/germline_cpsr_vcf")
  germline_cpsr_tsv_dir=paste0(basedir, "/wgs/symlinks/germline_cpsr_tsv")
  
  #A5_germline_cpsr <- load_germline_vcfs(germline_cpsr_vcf_dir, combine = T, threads = threads) # 
  A5_germline_cpsr <- load_germline_cpsr_tsv(germline_cpsr_tsv_dir)
  
  A5_germline_cpsr_keep <- A5_germline_cpsr %>% 
    filter(!(A5_ID %in% exclude_samples), 
           !(CLINVAR_CLASSIFICATION %in% c("Benign","Likely_Benign")),
           !(grepl("Q80", CDS_CHANGE) & SYMBOL=="AR"))
  
  A5_germline_cpsr_keep_pcawg <- A5_germline_cpsr_keep %>%  filter(SYMBOL %in% pcawg_drivers$Gene)
  
  A5_germline_cpsr_keep_pcawg_brief <- A5_germline_cpsr_keep_pcawg %>% mutate(Annotation=paste0("Germline:",CONSEQUENCE, "(",CLINVAR_CLASSIFICATION,")")) %>%
    dplyr::select(A5_ID, SYMBOL, Annotation) %>%  dplyr::rename(Gene=SYMBOL) %>% 
    group_by(A5_ID, Gene) %>% summarise(Annotation=paste(Annotation, sep = ",", collapse = ","))
  
  assign("A5_germline_cpsr", A5_germline_cpsr, envir = globalenv())
  assign("A5_germline_cpsr_keep_pcawg", A5_germline_cpsr_keep_pcawg, envir = globalenv())
  assign("A5_germline_cpsr_keep_pcawg_brief", A5_germline_cpsr_keep_pcawg_brief, envir = globalenv())
}
####################
# Somatic Variants #
####################
snv_quickloadRDS <- "/g/data/pq08/projects/ppgl/a5/wgs/quickload_checkpoints/somatic_vcfs.rds"

data_loader_somatic_variants <- function(somatic_vcf_dir=NULL, quickload=TRUE, 
                                         blacklist="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/blacklists/blacklist_readsupport_gteq3_samplesupport_gteq3.tsv")
{
 if(quickload==FALSE & is.null(somatic_vcf_dir)) {
   stop("You must provide a somatic VCF directory when quickload==FALSE")
 }
   
  AnnoColNames <- c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature",
                    "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
                    "Protein_position","Amino_acids","Codons","Existing_variation", "ALLELE_NUM",
                    "DISTANCE","STRAND","FLAGS","PICK","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID",
                    "CANONICAL","MANE","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","RefSeq",
                    "DOMAINS","HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF",
                    "gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF",
                    "gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF",
                    "CLIN_SIG","SOMATIC","PHENO","CHECK_REF","NearestExonJB") 
  
  message("Loading Somatic Mutation Data...")
  
  CGI <- read.delim("/g/data/pq08/reference/gene_lists/cgi/catalog_of_validated_oncogenic_mutations.tsv")
  CGI <- CGI %>%
    separate_rows(gdna, sep="__") %>%
    mutate(gtemp=gsub("(chr[0-9XY]{1,2}):g.([0-9]+)([ATCG]+)>([ATCG]+)", "\\1-\\2-\\2-\\3-\\4",gdna),
           gtemp=gsub("(chr[0-9XY]{1,2}):g.([0-9]+)_([0-9]+)(.+)", "\\1-\\2-\\3--\\4",gtemp),
           gtemp=gsub("(chr[0-9XY]{1,2}):g.([0-9]+)(.+)", "\\1-\\2---\\3",gtemp)) %>%
    separate(gtemp, into=c("seqnames","start", "end", "REF", "ALT"), sep="-", remove=T) %>%
    mutate(start=as.numeric(start), end=as.numeric(end)) %>%
    dplyr::rename(CGI_context=context, CGI_source=source, CGI_info=info)
  
  
  keepregx.all <- "inframe_insertion|missense_variant|inframe_deletion|splice_acceptor_variant|stop_lost|stop_gained|splice_donor_variant|frameshift_variant|start_lost"
  keepregx.special.tert <- "upstream_gene_variant"
  keepregx.special.atrx <- "splice_region_variant|regulatory_region_variant"
  
  if (!quickload)
  {
    a5_somatic_variants <- load_somatic_vcfs(somatic_vcf_dir=somatic_vcf_dir,combine = T, threads = 3)
    a5_somatic_variants$start <- as.numeric(a5_somatic_variants$start)
    a5_somatic_variants$end <- as.numeric(a5_somatic_variants$end)
    a5_somatic_variants$TUMOR_AF <- as.numeric(a5_somatic_variants$TUMOR_AF)
    
    #a5_somatic_variants <- a5_somatic_variants[,!unlist(lapply(a5_somatic_variants, function (x) { all(is.na(x))}))]
    
    a5_somatic_variants <- a5_somatic_variants %>% mutate(A5_ID=gsub("T0","", A5_ID))
    a5_somatic_variants <- a5_somatic_variants %>% 
      left_join(CGI %>% dplyr::select(seqnames, start, end, ALT, CGI_context, CGI_source, CGI_info), by=c("seqnames", "start", "end", "ALT"))
    
    a5_somatic_variants <- a5_somatic_variants %>% 
      separate_rows(CSQ,sep=",") %>% 
      separate(col=CSQ, into=AnnoColNames, sep="[|]") %>%  
      group_by(A5_ID, seqnames, start,  ALT) %>% 
      mutate(IMPACT=factor(IMPACT, levels=c("MODIFIER", "LOW", "MODERATE", "HIGH" ))) %>% 
      arrange(desc(IMPACT), desc(CANONICAL), desc(MANE)) %>% 
      slice_head(n=1)
    
    saveRDS(a5_somatic_variants, snv_quickloadRDS)
  }
  else {
    a5_somatic_variants <- readRDS(snv_quickloadRDS)
  }
  
  a5_somatic_variants <- a5_somatic_variants %>% filter(!(A5_ID %in% exclude_samples))
  
  #Append TERT mutation incorrectly filtered by UMMCRise bug due to high mutation load
  missing_tert <- a5_somatic_variants %>% 
    filter(A5_ID=="E167-2", PCGR_SYMBOL=="TERT") %>% 
    mutate(A5_ID="E167-1", Tumour_AF="0.32", Tumour_ADJAF="0.3", Tumour_RD="45", Tumour_AD="21", Tumour_DP="66") 
  
  a5_somatic_variants <- bind_rows(a5_somatic_variants, missing_tert)
  
  a5_somatic_variants_keep <- a5_somatic_variants %>% filter((grepl(keepregx.all,PCGR_CONSEQUENCE) & PCGR_SYMBOL != "TERT") | 
                                                               (grepl(keepregx.special.atrx,PCGR_CONSEQUENCE) &  PCGR_SYMBOL == "ATRX") | 
                                                               (grepl(keepregx.special.tert,PCGR_CONSEQUENCE) &  PCGR_SYMBOL == "TERT"))
  
  if (!is.null(blacklist))
  {
    blacklisted_variants <- read.delim(blacklist)
    a5_somatic_variants <- a5_somatic_variants %>% 
      left_join(blacklisted_variants %>% 
                  dplyr::rename("blacklist_n_normal_observed"=n_above_threshold) %>% 
                  mutate(blacklist=TRUE),
                by = c("seqnames", "start", "end", "ALT"))
    
    a5_somatic_variants$blacklist_n_normal_observed[is.na(a5_somatic_variants$blacklist_n_normal_observed)] <- 0
    a5_somatic_variants$blacklist[is.na(a5_somatic_variants$blacklist)] <- F
    
  }
  
  a5_somatic_variants_keep_curated <- a5_somatic_variants_keep %>%  
    filter(PCGR_SYMBOL %in% shiva_curated_genes.variant)
  
  #failed validation
  a5_somatic_variants_keep_curated <- a5_somatic_variants_keep_curated %>% 
    filter(!(A5_ID=="E227-T01" && PCGR_SYMBOL=="KRAS"), 
           !(A5_ID=="E122-T01" && PCGR_SYMBOL=="NF1"))
  
  #Remove all E167-T01 variants except KRAS
  a5_somatic_variants_keep_curated <- a5_somatic_variants_keep_curated %>% 
    filter(A5_ID != "E167-T01" || 
             (A5_ID != "E167-T01" & PCGR_SYMBOL=="KRAS") || 
             (A5_ID != "E167-T01" & PCGR_SYMBOL=="TERT" & PCGR_SYMBOL=="TERT" & start == 1295113))
  
  a5_somatic_variants_keep_pcawg <- a5_somatic_variants_keep %>% filter(PCGR_SYMBOL %in% pcawg_drivers$Gene)
  
  a5_somatic_variants_keep_pcawg_brief <- 
    a5_somatic_variants_keep_pcawg %>%  
    distinct() %>% 
    mutate(PCGR_ClinSig=ifelse(!is.na(PCGR_CLINVAR_CLNSIG) & PCGR_CLINVAR_CLNSIG !="uncertain_significance",
                               paste0(",PCGR_ClinSig:",PCGR_CLINVAR_CLNSIG),
                               ""),
           CGI=ifelse(!is.na(CGI_source),
                      paste0(",CGI:",CGI_source,"-",CGI_context),
                      ""),
           Annotation=paste0(PCGR_CONSEQUENCE,"(", round(TUMOR_AF,2),")", 
                             CGI, 
                             PCGR_ClinSig)) %>%
    dplyr::select(A5_ID, PCGR_SYMBOL, Annotation) %>%  
    dplyr::rename(Gene=PCGR_SYMBOL) %>% 
    group_by(A5_ID, Gene) %>% 
    summarise(Annotation=paste(Annotation, sep = ",", 
                               collapse = ","))
  
  assign("a5_somatic_variants", a5_somatic_variants, envir = globalenv())
  assign("a5_somatic_variants_keep", a5_somatic_variants_keep, envir = globalenv()) 
  assign("a5_somatic_variants_keep_pcawg", a5_somatic_variants_keep_pcawg, envir = globalenv()) 
  assign("a5_somatic_variants_keep_pcawg_brief", a5_somatic_variants_keep_pcawg_brief, envir = globalenv()) 
  assign("a5_somatic_variants_keep_curated", a5_somatic_variants_keep_curated, envir = globalenv()) 
  
  message("Created objects", 
          " a5_somatic_variants(all variants), ",
          "a5_somatic_variants_keep (filtered for consequence), ",
          "a5_somatic_variants_keep_curated (filtered for consequence and Shiva's curated list), ", 
          "a5_somatic_variants_keep_pcawg (filtered for consequence and PCAWG CG list), ",
          "a5_somatic_variants_keep_pcawg_brief (summarised version)")
  
}

######################
# Structral Variants #
######################

## MANTA

data_loader_sv_manta <- function()
{
  
  message("Loading Manta SV Data...")
  
  #taskset 0xffff R
  system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))
  cl <- makeCluster(3, outfile="")
  clusterEvalQ(cl, library("dplyr"))
  clusterEvalQ(cl, library("tidyr"))
  clusterEvalQ(cl, source("C:/Users/AFFLY/OneDrive - The University of Melbourne/ResearchData/RADIO/A5/Scripts/dataloaders/load_mantavcfs.R"))
  
  A5_manta.files <- fetchPurpleMantaFileNames(basedir)
  #A5_manta.files <- A5_manta.files[names(A5_manta.files) != "E138"]
  A5_manta.vcfs <- readPurpleMantaVCFs(A5_manta.files)
  A5_manta.df <- parLapply(cl = cl, X = A5_manta.vcfs,fun = PurpleMantaVCFtoDataFrame)
  stopCluster(cl)
  
  #saveRDS(A5_manta.df,"C:/ResearchData/RADIO/A5/Data/WGS/Tool outputs/Manta/A5_manta.df.rds")
  #A5_manta.df <- readRDS("C:/ResearchData/RADIO/A5/Data/WGS/Tool outputs/Manta/A5_manta.df.rds")
  FocalDelSize=1000000
  A5_manta.df_keep <- lapply(A5_manta.df, function (mantasvs, GenesOfInterest) {
    mantasvs <- mantasvs[!(mantasvs$SVTYPE %in% c("DEL","DUP")) |  
                           (mantasvs$SVTYPE %in% c("DEL","DUP") & apply(mantasvs[, c("BPI_START","BPI_END")], 1, max)-apply(mantasvs[, c("BPI_START","BPI_END")], 1, min) < FocalDelSize),]
    
    mantasvs <- mantasvs[sapply(stringr::str_split(mantasvs$PURPLE_AF,pattern = ","), max) > 0.05,]
    
    boolean_store <- vector(length = nrow(mantasvs),mode = "logical")
    for (i in 1:length(GenesOfInterest))
    {
      boolean_store <- boolean_store | grepl(pattern = paste0("&",GenesOfInterest[i],"&|^",GenesOfInterest[i],"&|&",GenesOfInterest[i],"$|^",GenesOfInterest[i],"$"),x = mantasvs$Gene_Name)
    }
    return(mantasvs[boolean_store,])
  }, GenesOfInterest=pcawg_drivers$Gene) 
  A5_manta.df_keep <- bind_rows(A5_manta.df_keep, .id = "A5_ID")
  
  A5_manta.df_keep_pcawg_brief <- A5_manta.df_keep_pcawg %>%  mutate(Coord=paste0(seqnames,":",start,"-",end)) %>% select(A5_ID, Gene, SVTYPE, Coord) %>%  distinct() %>%  group_by(A5_ID, Gene, SVTYPE) %>% summarise(Coord=paste(Coord, sep="", collapse = ";")) 
  ##Needs work ####
  #A5_manta.df_keep_pcawg %>%  mutate(Size=ifelse(SVTYPE!="BND", BPI_END-BPI_START,0)) %>%  select(A5_ID, Gene, seqnames,  start, BPI_START, BPI_END, REF,   ALT, SVTYPE, Size) %>% distinct() %>% group_by(A5_ID, Gene) %>% mutate(nEvents=n()) %>% arrange(A5_ID, Gene) %>% filter(!(SVTYPE=="DEL" & nEvents==1)))
  
  # mantasv_in_pcawg_gene <-
  #   unlist(
  #     lapply(
  #       lapply(
  #         lapply(A5_manta.df_keep$`GENE(s)`, #Split genes on &
  #                stringr::str_split_fixed,
  #                pattern="&",
  #                n=Inf),
  #         "%in%",
  #         pcawg_drivers$Gene), #check if genes match between lists
  #       any)) #return true if any match
  
  assign("A5_manta.df_keep", A5_manta.df_keep, envir = globalenv()) 
  assign("A5_manta.df_keep_pcawg", A5_manta.df_keep_pcawg, envir = globalenv()) 
  assign("A5_manta.df_keep_pcawg_brief", A5_manta.df_keep_pcawg_brief, envir = globalenv()) 
  
  message("Created objects A5_manta.df_keep (filtered to BND/INV and focal DEL/DUP <1MB), A5_manta.df_keep_pcawg (further filtered to PCAWG CG list), A5_manta.df_keep_pcawg_brief  (summarised version)")
  
}

## GRIDSS
data_loader_sv_gridsslinx <- function(quickload=TRUE, threads=1)
{
  linx_quickloadRDS=paste0(basedir, "/wgs/quickload_checkpoints/linx_sv_quickload.rdata")
  
  message("Loading GRIDSS SV Data...")
  
  if(quickload & file.exists(linx_quickloadRDS))
  {
    A5_gridss <- readRDS(linx_quickloadRDS)
  } else {
    pre_processed_linx_data <- paste0(basedir, "/wgs/quickload_checkpoints/all_linx_svs.rdata")
    if(file.exists(pre_processed_linx_data)) {
      message("Found preprocessed linx data in:", pre_processed_linx_data, ". To completely reprocess the data, delete preprocessed data file and rerun the data loader")
      A5_gridss.svs <- readRDS(pre_processed_linx_data)
    } else {
      A5_gridss.svs <- load_linx_all_samples(linx_base_dir = paste0(basedir, "/wgs/analysis/linx"),
                                             purple_base_dir = paste0(basedir, "/wgs/analysis/purple"),
                                             threads=threads)
      saveRDS(A5_gridss.svs, pre_processed_linx_data)
    }
    
    A5_gridss <- lapply(A5_gridss.svs, function (sample_svs) { bind_rows(lapply(sample_svs, summarise_sv_as_row)) })
    A5_gridss <- bind_rows(A5_gridss, .id = "A5_ID")
    
    saveRDS(A5_gridss, linx_quickloadRDS)
    
  }
  
  gridss_blacklist <- readr::read_delim(paste0(basedir, 
                                               "/wgs/analysis/gridss/blacklist/blacklist_gteq2.tsv"),
                                        col_types = "cdccddccdcdd")
  
  #Extract second chromosome info from ALT format
  gridss_blacklist <- gridss_blacklist %>% 
    mutate(CHROM2_POS2=ifelse(
      stringr::str_detect(ALT,"chr"), 
      gsub(".*(chr[0-9XYM]{1,2}:[0-9]+).*","\\1",ALT),
      "-:0")) %>% 
    separate(CHROM2_POS2, into = c("chrom2", "start2"), sep=":") %>% 
    mutate(start2 = as.numeric(start2)) %>% 
    dplyr::rename(chrom1=CHROM, start1=POS)
  
  
  A5_gridss <- A5_gridss %>% mutate(across(.cols = matches("^start|^end|SR|RP"), .fns = as.numeric))
  
  A5_gridss <- A5_gridss %>%
    mutate(chrom2 = replace_na(chrom2, "-"),
           start2 = replace_na(start2, 0)) %>%
    left_join(gridss_blacklist)
  
  
  gridss_recurrent_cutoff <- 2
  
  #Remove SVs called purely on copy number by purple 
  A5_gridss_keep <- A5_gridss %>%  filter(!grepl("purple", vcfId))
  
  #Remove SVs called in more than a threshold number of unrelated patients
  A5_gridss_keep <- A5_gridss_keep %>%
    filter(is.na(observed_patient_count) | observed_patient_count <= gridss_recurrent_cutoff)
  
  #Apply support filters
  A5_gridss_keep <- A5_gridss_keep %>% filter(SR>3 | RP >3 | (SR >1 & RP >1))
  
  
  A5_gridss_keep_pcawg <- A5_gridss_keep %>%  filter(GeneStartName %in% pcawg_drivers$Gene | GeneEndName %in% pcawg_drivers$Gene)
  
  StartCols <- colnames(A5_gridss_keep_pcawg)[grep("Start",colnames(A5_gridss_keep_pcawg))]
  EndCols <- colnames(A5_gridss_keep_pcawg)[grep("End",colnames(A5_gridss_keep_pcawg))]
  
  A5_gridss_keep_pcawg.long <- bind_rows(A5_gridss_keep_pcawg %>% 
                                           mutate(fusion_pair=ifelse(isFusion, paste(GeneStartName,GeneEndName,sep="_"), NA)) %>% 
                                           dplyr::select(!contains("GeneEnd")) %>% rename_with(.fn = ~gsub("Start" , "", .x)),
                                         A5_gridss_keep_pcawg  %>% 
                                           mutate(fusion_pair=ifelse(isFusion, paste(GeneStartName,GeneEndName,sep="_"), NA)) %>% 
                                           dplyr::select(!contains("GeneStart")) %>% rename_with(.fn = ~gsub("End" , "", .x))) %>%  distinct()
  
  
  A5_gridss_keep_curated <- A5_gridss_keep %>% filter((GeneStartName %in% shiva_curated_genes.sv & GeneStartDisrupted) | (GeneEndName %in% shiva_curated_genes.sv & GeneEndDisrupted))
  
  assign("A5_gridss", A5_gridss, envir = globalenv()) 
  assign("A5_gridss_keep", A5_gridss_keep, envir = globalenv()) 
  assign("A5_gridss_keep_pcawg", A5_gridss_keep_pcawg, envir = globalenv()) 
  assign("A5_gridss_keep_pcawg.long", A5_gridss_keep_pcawg.long, envir = globalenv()) 
  assign("A5_gridss_keep_curated", A5_gridss_keep_curated, envir = globalenv()) 
  
  message("Loaded objects into global environment: 
          A5_gridss (dataframe of all SVs in gridss/purple/linx output),
          A5_gridss_keep (filtered to keep only SVs seen in <=",gridss_recurrent_cutoff," patients 
                          with 'SR>3 | RP>1 | (SR >1 & RP >1)')
          A5_gridss_keep_pcawg (keep list filtered to PCAWG cancer gene list),
          A5_gridss_keep_pcawg.long  (start/end gene stacked),
          A5_gridss_keep_curated (keep filtered to Shiva's curation list)")
  
}

####################
# Gene Copy Number #
####################

data_loader_gene_cn <- function()
{
  
  message("Loading Gene CN Data...")
  
  FocalDelSize=250000
  
  purple_cn_dir <- paste0(basedir, "/wgs/symlinks/purple_cn")
  
  A5_gene_cn <- readPurpleGeneCN(fetchPurpleGeneCNFileNames(purple_cn_dir))
  A5_gene_cn$A5_ID <- gsub(".purple.cnv.gene.tsv","",A5_gene_cn$A5_ID)
  #A5_gene_cn <- A5_gene_cn %>% mutate(A5_ID=gsub("(E[0-9]{3})$","\\1_1", A5_ID), A5_ID=gsub("_","-", A5_ID))
  A5_gene_cn <- A5_gene_cn %>% filter(!(A5_ID %in% exclude_samples))
  ploidyadjust <- A5_gene_cn %>% group_by(A5_ID) %>% summarise(MedCN=median(maxCopyNumber)) %>% mutate(ploidy=round(MedCN,0),ploidyadj=ploidy-2) %>% dplyr::select(-MedCN)
  A5_gene_cn <- A5_gene_cn %>% left_join(ploidyadjust) %>% mutate(minCopyNumber_ploidyAdj=minCopyNumber-ploidyadj,maxCopyNumber_ploidyAdj=maxCopyNumber-ploidyadj)
  A5_gene_cn <- A5_gene_cn %>% left_join(a5_anno %>% dplyr::select(A5_ID, Gender))
  A5_gene_cn_keep_pcawg <- A5_gene_cn %>% filter(minCopyNumber_ploidyAdj < 1.75 | maxCopyNumber_ploidyAdj > 2.25, gene %in% pcawg_drivers$Gene)
  A5_gene_cn_keep_pcawg_brief <- A5_gene_cn_keep_pcawg %>%
    mutate(
      Event=case_when(
        (minCopyNumber > 1.8 & minCopyNumber < 2.25) & minMinorAlleleCopyNumber < 0.3 ~ "CNLOH",
        minCopyNumber > 2.8 & minMinorAlleleCopyNumber < 0.3 ~ "GainLOH",
        minCopyNumber < 0.25 ~ "CNHomoDel",
        minCopyNumber > 0.25 & minCopyNumber < 0.85  ~ "CNHomoDelSubClonal",
        (minCopyNumber != maxCopyNumber) & (minCopyNumber/maxCopyNumber < 0.75) ~ "CNgeneBreak",
        minCopyNumber > 0.75 & minCopyNumber < 1.3 & (minRegionEnd-minRegionStart) < FocalDelSize  ~ "CNDelFocal",
        minCopyNumber > 1.3 & minCopyNumber < 1.75  & (minRegionEnd-minRegionStart) < FocalDelSize  ~ "CNDelFocalSubClonal",    
        minCopyNumber > 0.75 & minCopyNumber < 1.3 & !(chromosome %in% c("X","Y") & Gender=="male") ~ "CNDel",
        minCopyNumber > 1.3 & minCopyNumber < 1.75 & !(chromosome %in% c("X","Y") & Gender=="male") ~ "CNDelSubClonal",
        minCopyNumber > 0.75 & minCopyNumber < 1.3 & (chromosome %in% c("X","Y") & Gender=="male") ~ "X/Y-Male",
        minCopyNumber > (2*ploidy) & maxCopyNumber > (2*ploidy) ~ "CNGain",
        TRUE ~ "Exclude"), 
      Annotation=ifelse(minCopyNumber != maxCopyNumber, paste0(Event,"(",round(minCopyNumber,2),"/",round(maxCopyNumber,2),")"), paste0(Event,"(",round(minCopyNumber,2),")"))
    ) %>% 
    filter(Event != "Exclude") %>% 
    dplyr::select(A5_ID, gene, Annotation) 
  
  assign("A5_gene_cn", A5_gene_cn, envir = globalenv()) 
  assign("A5_gene_cn_keep_pcawg", A5_gene_cn_keep_pcawg, envir = globalenv()) 
  assign("A5_gene_cn_keep_pcawg_brief", A5_gene_cn_keep_pcawg_brief, envir = globalenv()) 
  
  message("Created objects A5_gene_cn, A5_gene_cn_keep_pcawg (filtered to PCAWG CG list and non-diploid), A5_gene_cn_keep_pcawg_brief (summarised version)")
  
}


################
# CNA segments #
################

chromothripsis_regions <- data.frame(A5_ID=c("E123-T01","E128-T02","E138-T01","E180-T01","E180-T01","E231-T01","E198-T01","E198-T01"),
                                     chromosome=c("chr5","chr1","chr2","chr3","chr11","chr3","chr17","chr22"),
                                     start=c(0,124000000,0,0,0,50000000,0,17000000),
                                     end=c(182000000,248658000,242100000,197500000,135066000,198500000,10000000,50000000))

data_loader_cna_segs <- function()
{
  
  purple_seg_dir <- paste0(basedir, "/wgs/symlinks/purple_seg")
  
  A5_seg.files <- list.files(purple_seg_dir, pattern = ".purple.segment.tsv", full.names = T ,recursive = T)
  names(A5_seg.files) <- gsub(".purple.segment.tsv", "", basename(A5_seg.files))
  A5_seg <- lapply(A5_seg.files, read.delim)
  
  A5_seg <- bind_rows(A5_seg, .id="A5_ID")
  A5_seg <- A5_seg %>% mutate(A5_ID=gsub("T0","",A5_ID))
  A5_seg <- A5_seg %>%  dplyr::filter(chromosome!="chrY")
  
  A5_seg_keep <- A5_seg  %>%  filter(germlineStatus != "NOISE" & 
                                       (
                                         (end-start > 1000000 & (bafCount > 20 | chromosome=="chrX") & depthWindowCount > 100) | 
                                           support %in% c("BND","DEL","SGL","DUP","INV","MULTIPLE","INF")
                                       )
  )
  
  A5_seg_keep <- A5_seg_keep %>% left_join(a5_anno %>%  dplyr::select(A5_ID, Gender))
  A5_seg_keep <- A5_seg_keep %>% group_by(A5_ID) %>% mutate(mean_tumorCopyNumber=mean(tumorCopyNumber)) %>%  ungroup()
  
  chr_offsets <- A5_seg_keep %>% 
    mutate(chromosome=factor(as.character(chromosome), levels=paste0("chr",c(1:22,"X")))) %>% 
    arrange(chromosome) %>% group_by(A5_ID, chromosome) %>% 
    summarise(max_end=max(end)) %>% 
    ungroup() %>% group_by(chromosome) %>% 
    summarise(chr_size=median(max_end)) %>% mutate(offset = cumsum(as.numeric(chr_size))-chr_size)
  
  A5_seg_keep <- A5_seg_keep %>% left_join(chr_offsets) %>% mutate(start_offset=start+offset, end_offset=end+offset)
  
  A5_seg_keep <- A5_seg_keep %>% mutate(
    Class = classify_cna_event(minorAlleleCopyNumber,majorAlleleCopyNumber, Gender,chromosome, tumorCopyNumber, mean_tumorCopyNumber) #function sourced from cna_event_labeller.R
  ) %>%
    mutate(
      Class = ifelse((
        Class == "Minor Subclonal Loss" &
          (
            !is.na(lead(Class)) & lead(Class) == "Subclonal Loss"
          ) |
          (!is.na(lag(Class))  &
             lag(Class) == "Subclonal Loss")
      ),
      "Subclonal Loss",
      Class),
      Class = ifelse((
        Class == "Minor Subclonal Loss" &
          (
            !is.na(lead(Class)) & lead(Class) == "Diploid/Haploid-X"
          ) |
          (!is.na(lag(Class))  &
             lag(Class) == "Diploid/Haploid-X")
      ),
      "Diploid/Haploid-X",
      Class
      )
    ) %>%
    mutate(Class = factor(Class, levels = cn_event_types))
  
  
  
  #Manual overide for chromothripsis regions
  
  for (cr in 1:nrow(chromothripsis_regions))
  {
    A5_seg_keep <- A5_seg_keep %>% mutate(Class=ifelse(
      (A5_ID==chromothripsis_regions$A5_ID[[cr]] | A5_ID==gsub("T0", "", chromothripsis_regions$A5_ID[[cr]])) & 
        (chromosome==chromothripsis_regions$chromosome[[cr]] & 
           ((start >= chromothripsis_regions$start[[cr]]) & (end <= chromothripsis_regions$end[[cr]]))), "Chromothripsis", as.character(Class)))
  }
  
  # A5_seg_keep <- A5_seg_keep %>% mutate(Class=case_when(
  #   ((A5_ID=="E123-1" | A5_ID=="E123-T01") & chromosome=="chr5") ~ "Chromothripsis",
  #   ((A5_ID=="E128-2" | A5_ID=="E128-T02") & (chromosome=="chr1" & start > 124*10^6)) ~ "Chromothripsis",
  #   ((A5_ID=="E138-1" | A5_ID=="E138-T01") & chromosome=="chr2") ~ "Chromothripsis",
  #   ((A5_ID=="E180-1" | A5_ID=="E180-T01") & (chromosome=="chr3" | chromosome=="chr11")) ~ "Chromothripsis",
  #   ((A5_ID=="E231-1" | A5_ID=="E231-T01") & (chromosome=="chr3" & start > 50*10^6)) ~ "Chromothripsis",
  #   ((A5_ID=="E198-1" | A5_ID=="E198-T01") & (chromosome=="chr17" & start < 10*10^6)) ~ "Chromothripsis",
  #   ((A5_ID=="E198-1" | A5_ID=="E198-T01") & (chromosome=="chr22" & ((start > 17*10^6) & (end < 50*10^6)))) ~ "Chromothripsis",
  #   TRUE ~ as.character(Class)
  # ))
  
  
  min_gap <-5000000
  A5_seg_keep.merged <- A5_seg_keep %>% mutate(Class=gsub("Minor Subclonal Loss","Subclonal Loss", Class)) %>% arrange(A5_ID, chromosome,start,end) %>% group_by(A5_ID, chromosome) %>% 
    mutate(
      # lead_class=lead(Class),     
      # lag_class=lag(Class),
      #      lag_start=lag(start),
      #      lag_end=lag(end),
      #      diff_start_lend=start-lag(end),
      #      diff_lstart_end=lead(start) - end,
      seg_boundary=case_when(
        (Class != lead(Class) & Class != lag(Class)) ~ "Standalone",
        ((start-lag(end) > min_gap | is.na(lag(end))) & (is.na(lead(start)) | lead(start) - end > min_gap)) ~ "Standalone",
        (Class != lead(Class) & (start-lag(end) > min_gap | is.na(lag(end)))) ~ "Standalone",
        (Class != lag(Class) & (is.na(lead(start)) | lead(start) - end > min_gap)) ~ "Standalone",
        ((Class != lag(Class)) | (start-lag(end) > min_gap) | is.na(lag(start))) ~ "Begin",
        Class != lead(Class) | lead(start) - end > min_gap | is.na(lead(start)) ~ "End",
        TRUE ~ NA_character_
      )
    ) %>% 
    filter(!is.na(seg_boundary)) %>% 
    mutate(new_seg_start=ifelse(seg_boundary %in% c("Standalone", "Begin"), start, lag(start)),
           new_seg_end=ifelse(seg_boundary %in% c("Standalone", "End"), end, lead(end))) %>% 
    dplyr::select(-start, -end, -seg_boundary) %>% 
    group_by(A5_ID,chromosome,Gender,Class, new_seg_start, new_seg_end, offset) %>% slice_head(n=1) %>% ungroup() %>% 
    arrange(A5_ID,chromosome, new_seg_start, new_seg_end)
  
  #Annotate merged segments with median of values from constituant segments
  A5_seg_keep.merged$observedBAF <- NA
  A5_seg_keep.merged$minorAlleleCopyNumber <- NA
  A5_seg_keep.merged$majorAlleleCopyNumber <- NA
  A5_seg_keep.merged$tumorCopyNumber <- NA
  for (i in 1:nrow(A5_seg_keep.merged))  
  {
    SOI <- A5_seg_keep[A5_seg_keep$A5_ID==A5_seg_keep.merged$A5_ID[i] &
                         A5_seg_keep$chromosome==A5_seg_keep.merged$chromosome[i] &
                         A5_seg_keep$start >= A5_seg_keep.merged$new_seg_start[i] &
                         A5_seg_keep$end <= A5_seg_keep.merged$new_seg_end[i],]
    
    A5_seg_keep.merged$observedBAF[i] <- median(SOI$observedBAF) 
    A5_seg_keep.merged$minorAlleleCopyNumber[i] <- median(SOI$minorAlleleCopyNumber, na.rm = T) 
    A5_seg_keep.merged$majorAlleleCopyNumber[i] <- median(SOI$majorAlleleCopyNumber, na.rm = T)
    A5_seg_keep.merged$tumorCopyNumber[i] <- median(SOI$tumorCopyNumber, na.rm = T)
  }
  
  
  A5_seg_keep.merged <- A5_seg_keep.merged %>% 
    dplyr::rename(start=new_seg_start, end=new_seg_end) %>%  
    mutate(start_offset=start+offset, end_offset=end+offset) %>% 
    mutate(Class=factor(Class,levels=cn_event_types)) %>% 
    arrange(A5_ID, chromosome,start,end) %>% 
    mutate(Class=case_when(
      Class=="Subclonal Gain" & majorAlleleCopyNumber < 1.15 ~ factor("Diploid/Haploid-X", levels=cn_event_types),
      Class=="Minor Subclonal Loss" & minorAlleleCopyNumber > 0.95 ~ factor("Diploid/Haploid-X", levels=cn_event_types),
      TRUE ~ Class
    ))
  
  assign("A5_seg", A5_seg, envir = globalenv())
  assign("A5_seg_keep", A5_seg_keep, envir = globalenv())
  assign("A5_seg_keep.merged", A5_seg_keep.merged, envir = globalenv())
  assign("chr_offsets", chr_offsets, envir = globalenv())
  assign("chromothripsis_regions", chromothripsis_regions, envir = globalenv())
  
  message("Created objects A5_seg, A5_seg_keep (size/quality filtered and annotated), A5_seg_keep.merged (fused similar segments), chr_offsets (linear chr pos for plotting")
  
  
}

###########
# Cleanup #
###########

# rm(list = c("A5_mafs","A5_zexp","A5_gene_cn", "A5_manta.df_keep", "A5_manta.vcfs"))
# gc()

