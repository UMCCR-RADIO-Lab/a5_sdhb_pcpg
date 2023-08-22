library(dplyr)
library(tidyr)

setwd("/g/data/pq08/projects/ppgl/")

################
# Data Loaders #
################

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
if(!exists("a5_anno")) {
  data_loader_a5_clinical_anno(use_cache = T) }

source("./a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_somatic_variants(quickload = T)

source("./a5/wgs/scripts/paired_sample_analysis/process_paired_pileups.r")

blacklist <- read.delim("./a5/wgs/analysis/mpileup/blacklist/blacklists/blacklist_readsupport_gteq3_samplesupport_gteq3.tsv")

###############################
# Filter deleterious variants #
###############################

deleterious_vars <- a5_somatic_variants_keep %>%  
  ungroup() %>% 
  dplyr::select(CHROM=seqnames, POS=start, REF, ALT) %>% 
  distinct() %>% 
  mutate(Deleterious="Deleterious")


############################
# Remove unrelated primary #
############################

#pileup_vafs$E159 <- pileup_vafs$E159 %>%  dplyr::select(-`E159-T02_pileup_vaf`)

######################
# Generate summaries #
######################

summaries <- list()
min_read_support=3

for (patient in names(pileup_vafs)) {
  
  #Attach alt counts to VAFs and pivot to have one sample per row
  vaf_data <- pileup_vafs[[patient]] %>% 
    left_join(
      pileup_summaries[[patient]] %>% 
        dplyr::select(c("CHROM", "POS", "REF", "ALT", ends_with("altCount")))) %>% 
    rename_with(.fn = \(x) gsub("([0-9])_","\\1__", x)) %>% 
    pivot_longer(cols = -c(CHROM,POS,REF,ALT),  
                 names_to = c("Sample", "Metric"), 
                 names_sep = "__", values_to = "value") %>% 
    pivot_wider(id_cols = c(CHROM, POS, REF, ALT, Sample), 
                names_from = "Metric", 
                values_from = "value")
  
    
    #set VAF of variants with support below threshold in tumour to zero
    vaf_data <- vaf_data %>% 
      mutate(pileup_vaf=ifelse(!grepl("-B01",Sample) & altCount < min_read_support, 0, pileup_vaf)) 
    
    #pivot to one row per variant
    vaf_data <- vaf_data %>% 
      mutate(Sample=paste0(Sample, "_pileup_vaf")) %>% 
    pivot_wider(id_cols = c(CHROM,POS,REF,ALT), 
                names_from = Sample, 
                values_from = pileup_vaf)
    
  
    #Annotate and filter blacklist
    vaf_data <- vaf_data %>%  
    left_join(blacklist) %>% 
    filter(is.na(n_above_threshold)) 
  
    #Filter 
      vaf_data <- vaf_data %>%   
      filter(!!sym(paste0(patient, "-B01_pileup_vaf")) == 0)
    
  temp_colname <- gsub("T0(.)_pileup_vaf","\\1",colnames(vaf_data))
  pub_ids <- a5_anno %>%  filter(A5_ID %in% intersect(temp_colname,a5_anno$A5_ID)) %>% dplyr::select(A5_ID, PublicationID)
  pub_ids <- set_names(pub_ids$PublicationID, pub_ids$A5_ID)
  
  temp_colname <- recode(temp_colname, !!!pub_ids)
  colnames(vaf_data) <- temp_colname
  vaf_data <- vaf_data %>%  dplyr::select(-matches("-B01_pileup_vaf"))
  
  vaf_data <- vaf_data %>%  
    pivot_longer(cols = starts_with(patient), 
                 names_to = "Sample", values_to = "VAF") %>% 
    filter(VAF > 0) %>% 
    group_by(CHROM, POS, REF, ALT) %>% 
    summarise(Privacy=paste(Sample, sep="/", collapse="/"))
  
  vaf_data <- vaf_data %>% 
    left_join(deleterious_vars) 
  vaf_data$Deleterious[is.na(vaf_data$Deleterious)] <- "Non-deleterious"
  
  vaf_data <- vaf_data %>% 
    ungroup() %>% 
    mutate(nVar_Total=n(), 
           nVar_Total_Deleterious=(vaf_data %>% filter(Deleterious=="Deleterious") %>% nrow())) %>%  
    group_by(Privacy, nVar_Total, nVar_Total_Deleterious, Deleterious) %>%  
    dplyr::count() 
  
  walk(.x = pub_ids, .f =  function(pub_id) {
    temp_all <- vaf_data %>%  ungroup() %>% filter(grepl(pub_id, Privacy)) %>%  summarise(nMut=sum(n)) %>% pull(nMut)
    temp_del <- vaf_data %>%  ungroup() %>% filter(grepl(pub_id, Privacy), Deleterious=="Deleterious") %>%  summarise(nMut=sum(n)) %>% pull(nMut)
    temp <- vaf_data %>%  mutate("{pub_id}_nTotal" := temp_all, "{pub_id}_nDel" := temp_del)
    assign(x = "vaf_data", value = temp,envir = globalenv())
  })
  
  
  vaf_data <- vaf_data %>% 
    pivot_wider(id_cols = c(Privacy, nVar_Total, nVar_Total_Deleterious, paste0(pub_ids, "_nTotal"), paste0(pub_ids, "_nDel")), 
                names_from = Deleterious, 
                values_from = n) %>% 
    dplyr::rename("nVar_Deleterious_PrivacyClass"=Deleterious, "nVar_NonDeleterious_PrivacyClass"=`Non-deleterious`) %>% 
    mutate(pcnt_GlobalMut_Deleterious=(nVar_Deleterious_PrivacyClass/nVar_Total_Deleterious)*100, 
           pcnt_GlobalMut_AllMut=((nVar_Deleterious_PrivacyClass+nVar_NonDeleterious_PrivacyClass)/nVar_Total)*100)
  
  walk(.x = pub_ids, .f =  function(pub_id) {
    temp <- vaf_data %>%  
      mutate("{pub_id}_pcntTotal" := ifelse(grepl(pub_id, Privacy), 
                                            ((nVar_Deleterious_PrivacyClass+nVar_NonDeleterious_PrivacyClass)/!!sym(paste0(pub_id, "_nTotal")))*100,
                                            NA),
             "{pub_id}_pcntDel" := ifelse(grepl(pub_id, Privacy),
                                          (nVar_Deleterious_PrivacyClass/!!sym(paste0(pub_id, "_nDel")))*100,
                                          NA))
    assign(x = "vaf_data", value = temp,envir = globalenv())
  })
  
  colnames(vaf_data) <- gsub(paste0(patient,"-"),"", colnames(vaf_data))
  
  summaries[[patient]] <- vaf_data
}

summaries <- bind_rows(summaries)


write.table(x = summaries, 
            file = "./a5/wgs/results/paired_sample_analysis/variant_overlap_counts.tsv", 
            sep="\t", 
            row.names = F, 
            quote=F)

#LEGACY
# 
# 
# 
# source(paste0(basedir,"/wgs/scripts/utility/compare_vcf_calls.r"))
# source(paste0(basedir,"/wgs/scripts/data_loaders/wgs_dataloaders.r"))
# source(paste0(basedir,"/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r"))
# library(patchwork)
# library(ggplot2)
# library(googlesheets4)
# library(furrr)
# 
# data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)
# 
# 
# #List of A5 patients with paired samples
# Samples <- list(
#   E122=c("E122-T01","E122-T02"),
#   E128=c("E128-T01","E128-T02"),
#   E132=c("E132-T01","E132-T02"),
#   E143=c("E143-T01","E143-T02","E143-T03"),
#   E146=c("E146-T01","E146-T02"),
#   E158=c("E158-T01","E158-T02"),
#   E159=c("E159-T01","E159-T02","E159-T03","E159-T04"),
#   E166=c("E166-T01","E166-T02"),
#   E167=c("E167-T01","E167-T02"),
#   E169=c("E169-T01","E169-T02"),
#   E225=c("E225-T01","E225-T02"),
#   E229=c("E229-T01","E229-T02")
# )
# 
# 
# with_primary <- a5_anno %>% filter(`Patient ID` %in% names(Samples)) %>%  filter(is_primary_or_met=="Primary") %>% pull(`Patient ID`) %>% unique()
# with_met <- a5_anno %>% filter(`Patient ID` %in% names(Samples)) %>%  filter(is_primary_or_met=="Metastasis") %>% pull(`Patient ID`) %>% unique()
# 
# Samples <- Samples[intersect(with_primary,with_met)]
# 
# 
# #Possible consequence values deleterious variants
# coding_consequences <- c(                  
#   "splice_acceptor_variant",
#   "splice_region_variant",
#   "splice_donor_variant",
#   "stop_gained",
#   "stop_lost",
#   "start_lost",
#   "missense_variant",
#   "frameshift_variant",
#   "inframe_insertion",
#   "inframe_deletion",
#   "synonymous_coding")
# 
# coding_consequences.regex <- paste(coding_consequences, collapse ="|")
# 
# threads=3
# plan(strategy = "multisession", workers=threads)
# 
# #Pairwise matching of VCFs - returns a list of dataframes with variant calls marked present/absent for each sample from each patient
# 
# 
# matched_variant_tables <- 
#   furrr::future_map(.x = Samples, 
#                     .f = function(patient_samples) {
#                       patient = stringr::str_extract(patient_samples[[1]], "E...")
#                       message("Started variant matching for: ", patient)
#                       
#                       vcf_files <- paste0(basedir, "/wgs/analysis/bcbio/",
#                                           patient,
#                                           "/umccrised/",
#                                           gsub("T0","",patient_samples),
#                                           "__",
#                                           patient_samples,
#                                           "/small_variants/",
#                                           gsub("T0","",patient_samples),
#                                           "__",
#                                           patient_samples,"-somatic-PASS.vcf.gz")
#                       names(vcf_files) <- patient_samples
#                       
#                       return_table <-
#                         multi_compare_small_variant_vcf(
#                           vcfs = vcf_files,
#                           normal_names = gsub("T.+", "B01", patient_samples),
#                           tumour_names = patient_samples,
#                           genome = "BSgenome.Hsapiens.UCSC.hg38",
#                           Keep_Only_Highest_Consequence = T,
#                           collapse_annotation = T)
#                       
#                       
#                       message("Completed variant matching for: ", patient)
#                       
#                       return(return_table)
#                     })
# 
# #Fix and samples names that have dashes replaced with dots
# matched_variant_tables <-
#   lapply(matched_variant_tables, function (table) {
#     return(rename_with(
#       table,
#       .fn = ~ gsub("(E[0-9]{3}).([TB][0-9]{2})", "\\1-\\2", .x),
#       .cols = matches("_Called$|_VAF$")
#     ))
#   })
# 
# re_id <- a5_anno %>% dplyr::select(A5_ID, `Patient ID`, PublicationID) %>% mutate(A5_ID=gsub("-","-T0", A5_ID)) %>% 
#   filter(`Patient ID` %in% names(Samples))
# 
# for(p in names(Samples))
# {
#   
#   for(i in 1:ncol(matched_variant_tables[[p]]))
#   {
#     
#     for(r in 1:nrow(re_id))
#     {
#       colnames(matched_variant_tables[[p]])[[i]] <- gsub(re_id$A5_ID[[r]], 
#                                                          re_id$PublicationID[[r]], 
#                                                          colnames(matched_variant_tables[[p]])[[i]])
#     }
#   } 
# }
# 
# 
# for(p in names(Samples))
# {
#   primary_col <- paste0(p, "-P1_Called")
#   met_cols <- colnames(matched_variant_tables[[p]])[grepl("-M._Called",colnames(matched_variant_tables[[p]]))]
#   nmets <- length(met_cols)
#   
#   if(nmets==1)
#   {
#     matched_variant_tables[[p]] %>% ungroup() %>% 
#       mutate(any_met=paste(!!sym(met_cols))) %>% 
#       dplyr::select(all_of(c(primary_col,  met_cols, "any_met"))) %>% 
#       group_by(!!!rlang::syms(c(primary_col, "any_met"))) %>%  dplyr::count()
#   }
#   
#   if(nmets==2)
#   {
#     matched_variant_tables[[p]] %>% ungroup() %>% 
#       mutate(any_met=paste(!!!rlang::syms(met_cols))) %>% 
#       dplyr::select(all_of(c(primary_col,  met_cols, "any_met"))) %>% 
#       group_by(!!!rlang::syms(c(primary_col, "any_met"))) %>%  dplyr::count()
#   }
#   
#   if(nmets==1)
#   {
#     matched_variant_tables[[p]] %>% ungroup() %>% 
#       mutate(any_met=paste(!!rlang::syms(test_cols))) %>% 
#       dplyr::select(all_of(c(primary_col,  "E225-M1_Called", "any_met"))) %>% 
#       filter(any_met)
#   }
# }
