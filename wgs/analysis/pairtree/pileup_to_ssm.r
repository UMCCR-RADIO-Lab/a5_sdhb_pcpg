library(tidyr)
source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_somatic_variants(quickload = T)
data_loader_cna_segs()

Samples <- list(
  E122=c("E122-T01","E122-T02"),
  E128=c("E128-T01","E128-T02"),
  E132=c("E132-T01","E132-T02"),
  E136=c("E136-T01","E136-T02"),
  E143=c("E143-T01","E143-T02","E143-T03"),
  E146=c("E146-T01","E146-T02"),
  E158=c("E158-T01","E158-T02"),
  E159=c("E159-T01","E159-T03","E159-T04"), #"E159-T02"
  #E166=c("E166-T01","E166-T02"),
  E167=c("E167-T01","E167-T02"),
  E169=c("E169-T01","E169-T02"),
  E225=c("E225-T01","E225-T02"),
  E229=c("E229-T01","E229-T02")
)


#################
# Read mpileups #
#################

mpileup_dir <- "/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/paired_samples/count_summaries"

mpileup_files <- list.files(mpileup_dir, pattern = "*.txt", full.names = T)
names(mpileup_files) <- gsub(".txt","", basename(mpileup_files))

mpileup_counts <- purrr::map(mpileup_files, readr::read_delim)

mpileup_counts <- mpileup_counts[intersect(names(mpileup_counts), names(Samples))]

#############
# Blacklist #
#############

variant_blacklist <- readr::read_tsv("/g/data/pq08/projects/ppgl/a5/wgs/analysis/blacklist_from_vcf/output/blacklist_reject_pr_ratio_gt1.tsv")
variant_blacklist <- variant_blacklist %>% 
  dplyr::rename(CHROM=Chr, POS=Pos, ALT=Alt)


blacklisted_variants <- variant_blacklist %>%  
  #filter(mutect2_REJECT > 0 | strelka2_REJECT > 0 | vardict_REJECT > 0) %>%  
  filter(pass_reject_ratio > 1) %>%  
  ungroup() %>% 
  dplyr::select(CHROM,POS, ALT) %>% 
  distinct() #%>% 
  #dplyr::rename(CHROM=seqnames, POS=start) 

######################
# Generate SSM files #
######################

pileup_to_ssm <- function (Patient_ID, counts) {
  
  A5_IDs <- gsub(pattern = "-T0", replacement = "-", Samples[[Patient_ID]])
  
  normal_alt_col <- which(grepl("B01_altCount", colnames(counts)))
  counts <- counts[counts[,normal_alt_col] == 0,]
  
  #remove variants using high stringency blacklist
  counts <- counts %>% anti_join(blacklisted_variants)
  
  #Pivot to long form and filter variants with any normal support and unused columns
  counts <- counts %>% 
    pivot_longer(cols = -c(CHROM, POS, REF, ALT)) %>% 
    filter(!grepl("refCount",name),
           grepl(paste(Samples[[Patient_ID]],  collapse="|"), name), #remove blood and extra primaries
           CHROM != "chrM")
    
  
  #Remove variants supported by <3 but greater than 0 reads in any sample as they are ambiguous
  counts <- counts %>% group_by(CHROM,POS, REF, ALT) %>% mutate(Garbage=ifelse(any(value <= 3 & value > 0), "Y", "N"))
  
  #Remove variants from excluded unrelated primaries
  zero_count_vars <- counts %>% group_by(CHROM,POS, REF, ALT) %>% 
    filter(grepl("_altCount", name)) %>% 
    filter(all(value == 0)) %>% 
    dplyr::select(CHROM, POS, REF, ALT) %>% 
    distinct()
  counts <- counts %>% anti_join(zero_count_vars)
  
  #correct A5_IDs, pivot wider as one sample per line, and append publicationIDs
  counts <- counts %>% mutate(name=gsub("-T0","-", name)) %>% 
    separate_wider_delim(cols = name, delim = "_", names = c("A5_ID", "metric")) %>% 
    left_join(a5_anno %>% dplyr::select(A5_ID,PublicationID, sample_purity)) 
  
  #Pivot wider and annotate VAF
  counts <- counts %>% 
    pivot_wider(id_cols = c(CHROM, POS, REF, ALT, A5_ID, PublicationID, sample_purity, Garbage), names_from = metric, values_from = value) %>% 
    mutate(VAF=altCount/TotalCount)
  
  ##################
  # Preprocess CNA #
  ##################
  
  #annotate segment CN 
  #counts$.gr <- makeGRangesFromDataFrame(counts, seqnames.field = "CHROM", start.field = "POS", end.field = "POS")
  patient_segs <- a5_seg_keep %>%  filter(A5_ID %in% A5_IDs) %>% 
    group_by(A5_ID) %>% 
    { 
      nm <- group_keys(.)[["A5_ID"]]
      . %>% group_split() %>% set_names(nm)
    }()
  
  
  patient_segs_gr <- map(patient_segs, .f = \(x) GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = T))
  
  
  # Annotate variants with CNA
  patient_mut_seg_annotated <- list()
  for (s in names(patient_segs_gr)) {
    
    temp <- counts %>%  
      filter(A5_ID == s)
    temp <- GenomicRanges::makeGRangesFromDataFrame(df = temp, 
                                                    seqnames.field = "CHROM", 
                                                    start.field = "POS", 
                                                    end.field = "POS",
                                                    keep.extra.columns = T)
    
    
    hits <- findOverlaps(subject = temp, query = patient_segs_gr[[s]])
    
    mcols(temp)["minorAlleleCopyNumber"] <- NA
    mcols(temp)["majorAlleleCopyNumber"] <- NA
    mcols(temp)[["minorAlleleCopyNumber"]][subjectHits(hits)] <- mcols(patient_segs_gr[[s]])[["minorAlleleCopyNumber"]][queryHits(hits)]
    mcols(temp)[["majorAlleleCopyNumber"]][subjectHits(hits)] <- mcols(patient_segs_gr[[s]])[["majorAlleleCopyNumber"]][queryHits(hits)]
    
    patient_mut_seg_annotated[[s]] <- GenomicRanges::as.data.frame(temp)
    
    
  }
  
  patient_mut_seg_annotated <- bind_rows(patient_mut_seg_annotated)
  patient_mut_seg_annotated <- patient_mut_seg_annotated %>%  arrange(PublicationID, seqnames, start)
  
  
  # Compute variant probability priors
  patient_mut_seg_annotated <- patient_mut_seg_annotated %>% 
    mutate(minorAlleleCopyNumber = round(minorAlleleCopyNumber),
           majorAlleleCopyNumber = round(majorAlleleCopyNumber)) %>% 
    mutate(var_read_prob=case_when(
      seqnames == "chrY" ~ 1,
      is.na(minorAlleleCopyNumber) | is.na(majorAlleleCopyNumber)  | majorAlleleCopyNumber == 0  ~ 0.5, #Missing CN data
      minorAlleleCopyNumber == 1 & majorAlleleCopyNumber == 1 ~ 0.5, #Diploid region
      minorAlleleCopyNumber == 0 & majorAlleleCopyNumber >= 2 & VAF > (sample_purity/2) ~ 1, #mutation before CNLOH
      minorAlleleCopyNumber == 0 & majorAlleleCopyNumber >= 2 & VAF < (sample_purity/2) ~ 1/majorAlleleCopyNumber, #mutation after CNLOH
      minorAlleleCopyNumber == 0 & majorAlleleCopyNumber == 1 ~ 1, #Loss
      .default= 1/(majorAlleleCopyNumber + minorAlleleCopyNumber)
    ))
  
  patient_mut_seg_annotated <- patient_mut_seg_annotated %>% dplyr::select(CHROM=seqnames,POS=start, REF, ALT, PublicationID, altCount,TotalCount, var_read_prob, Garbage) %>% 
    pivot_longer(cols = c(altCount, TotalCount, var_read_prob), names_to = "metric", values_to = "value")
  
  #Sort then collapse counts for each variant into a comma separated list of counts per sample
  patient_mut_seg_annotated <- patient_mut_seg_annotated %>%  
    arrange(CHROM, POS, metric, PublicationID) %>% 
    group_by(CHROM, POS, REF, ALT, metric, Garbage) %>% 
    summarise(sample_order = paste(PublicationID, collapse = ","), collapsed_values = paste(value, collapse = ",")) %>% 
    pivot_wider(id_cols = c(CHROM, POS, REF, ALT, sample_order, Garbage), names_from = metric, values_from = collapsed_values) 
  
  #Update column names to match SSM format and annotate gene names and row-number ids 
  ssm <- patient_mut_seg_annotated %>% 
    dplyr::rename(total_reads=TotalCount, var_reads=altCount) %>% 
    inner_join(a5_somatic_variants %>% ungroup()  %>%  dplyr::select(CHROM=seqnames, POS=start, PCGR_SYMBOL)  %>%  distinct()) %>% 
    ungroup() %>% 
    mutate(name=case_when(
      is.na(PCGR_SYMBOL) ~ paste(CHROM, POS, sep="_"),
      duplicated(PCGR_SYMBOL) | duplicated(PCGR_SYMBOL, fromLast=T) ~ paste(PCGR_SYMBOL,POS, sep="_"),
      .default = PCGR_SYMBOL
    )) %>% 
    mutate(id=paste0("s", row_number())) %>% 
    arrange(CHROM, POS) %>% 
    dplyr::select(id, name, var_reads, total_reads, var_read_prob, sample_order, Garbage)
  
  
  return(ssm)
}

#Process pilups to SSM
ssm <- purrr::map2(names(mpileup_counts), mpileup_counts, pileup_to_ssm)
names(ssm) <- names(mpileup_counts)

#sub samples E167 to save compute
set.seed(5)
E167_nmut <- nrow(ssm$E167)
keep <- ((1:E167_nmut %in% sample(1:E167_nmut, 2500)) | !grepl("^0", ssm$E167$var_reads))
ssm$E167 <- ssm$E167[keep,]

#############
# Write SSM/config #
#############
pairtree_input_path <- "/g/data/pq08/projects/ppgl/a5/wgs/analysis/pairtree/input"
purrr::walk2(.x = names(ssm), 
             .y = ssm, 
             .f = function(patientID, ssm)
             {
               
               config_dir <- file.path(pairtree_input_path,patientID)
               
               if(!dir.exists(config_dir)) { dir.create(config_dir)}
               
               ssm_path <- file.path(config_dir, paste(patientID, "ssm", sep="."))
               config_path <- file.path(config_dir, paste(patientID, "params.json", sep="."))
               
               #config
               sample_array <- gsub(',', '","', ssm$sample_order[1])
               garbage_array <- paste(ssm %>% filter(Garbage == "Y") %>%  pull(id), collapse='","')
               writeLines(paste0('{ "samples":["', sample_array, '"], "clusters":[], "garbage":["',garbage_array,'"]}'), con = config_path)
               
               #SSM
               readr::write_delim(ssm %>%  dplyr::select(-sample_order, -Garbage) , file = ssm_path, delim ="\t")
             })
