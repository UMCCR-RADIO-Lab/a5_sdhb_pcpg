basedir <- "/g/data/pq08/projects/ppgl/a5"
setwd(basedir)
renv::activate("./")

source(paste0(basedir,"/wgs/scripts/utility/compare_vcf_calls.r"))
source(paste0(basedir,"/wgs/scripts/data_loaders/wgs_dataloaders.r"))
source(paste0(basedir,"/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r"))
library(patchwork)
library(ggplot2)
library(googlesheets4)
library(furrr)

data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)


#List of A5 patients with paired samples
Samples <- list(
  E122=c("E122-T01","E122-T02"),
  E128=c("E128-T01","E128-T02"),
  E132=c("E132-T01","E132-T02"),
  E143=c("E143-T01","E143-T02","E143-T03"),
  E146=c("E146-T01","E146-T02"),
  E158=c("E158-T01","E158-T02"),
  E159=c("E159-T01","E159-T02","E159-T03","E159-T04"),
  E166=c("E166-T01","E166-T02"),
  E167=c("E167-T01","E167-T02"),
  E169=c("E169-T01","E169-T02"),
  E225=c("E225-T01","E225-T02"),
  E229=c("E229-T01","E229-T02")
)


with_primary <- a5_anno %>% filter(`Patient ID` %in% names(Samples)) %>%  filter(is_primary_or_met=="Primary") %>% pull(`Patient ID`) %>% unique()
with_met <- a5_anno %>% filter(`Patient ID` %in% names(Samples)) %>%  filter(is_primary_or_met=="Metastasis") %>% pull(`Patient ID`) %>% unique()

Samples <- Samples[intersect(with_primary,with_met)]


#Possible consequence values deleterious variants
coding_consequences <- c(                  
  "splice_acceptor_variant",
  "splice_region_variant",
  "splice_donor_variant",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "missense_variant",
  "frameshift_variant",
  "inframe_insertion",
  "inframe_deletion",
  "synonymous_coding")

coding_consequences.regex <- paste(coding_consequences, collapse ="|")

threads=3
plan(strategy = "multisession", workers=threads)

#Pairwise matching of VCFs - returns a list of dataframes with variant calls marked present/absent for each sample from each patient


matched_variant_tables <- 
  furrr::future_map(.x = Samples, 
                    .f = function(patient_samples) {
                      patient = stringr::str_extract(patient_samples[[1]], "E...")
                      message("Started variant matching for: ", patient)
                      
                      vcf_files <- paste0(basedir, "/wgs/analysis/bcbio/",
                                          patient,
                                          "/umccrised/",
                                          gsub("T0","",patient_samples),
                                          "__",
                                          patient_samples,
                                          "/small_variants/",
                                          gsub("T0","",patient_samples),
                                          "__",
                                          patient_samples,"-somatic-PASS.vcf.gz")
                      names(vcf_files) <- patient_samples
                      
                      return_table <-
                        multi_compare_small_variant_vcf(
                          vcfs = vcf_files,
                          normal_names = gsub("T.+", "B01", patient_samples),
                          tumour_names = patient_samples,
                          genome = "BSgenome.Hsapiens.UCSC.hg38",
                          Keep_Only_Highest_Consequence = T,
                          collapse_annotation = T)
                      
                      
                      message("Completed variant matching for: ", patient)
                      
                      return(return_table)
                    })

#Fix and samples names that have dashes replaced with dots
matched_variant_tables <-
  lapply(matched_variant_tables, function (table) {
    return(rename_with(
      table,
      .fn = ~ gsub("(E[0-9]{3}).([TB][0-9]{2})", "\\1-\\2", .x),
      .cols = matches("_Called$|_VAF$")
    ))
  })

re_id <- a5_anno %>% dplyr::select(A5_ID, `Patient ID`, PublicationID) %>% mutate(A5_ID=gsub("-","-T0", A5_ID)) %>% 
  filter(`Patient ID` %in% names(Samples))

for(p in names(Samples))
{
  
  for(i in 1:ncol(matched_variant_tables[[p]]))
  {
    
    for(r in 1:nrow(re_id))
    {
      colnames(matched_variant_tables[[p]])[[i]] <- gsub(re_id$A5_ID[[r]], 
                                                         re_id$PublicationID[[r]], 
                                                         colnames(matched_variant_tables[[p]])[[i]])
    }
  } 
}


for(p in names(Samples))
{
  primary_col <- paste0(p, "-P1_Called")
  met_cols <- colnames(matched_variant_tables[[p]])[grepl("-M._Called",colnames(matched_variant_tables[[p]]))]
  nmets <- length(met_cols)
  
  if(nmets==1)
  {
    matched_variant_tables[[p]] %>% ungroup() %>% 
      mutate(any_met=paste(!!sym(met_cols))) %>% 
      dplyr::select(all_of(c(primary_col,  met_cols, "any_met"))) %>% 
      group_by(!!!rlang::syms(c(primary_col, "any_met"))) %>%  dplyr::count()
  }
  
  if(nmets==2)
  {
    matched_variant_tables[[p]] %>% ungroup() %>% 
      mutate(any_met=paste(!!!rlang::syms(met_cols))) %>% 
      dplyr::select(all_of(c(primary_col,  met_cols, "any_met"))) %>% 
      group_by(!!!rlang::syms(c(primary_col, "any_met"))) %>%  dplyr::count()
  }
  
  if(nmets==1)
  {
    matched_variant_tables[[p]] %>% ungroup() %>% 
      mutate(any_met=paste(!!rlang::syms(test_cols))) %>% 
      dplyr::select(all_of(c(primary_col,  "E225-M1_Called", "any_met"))) %>% 
      filter(any_met)
  }
}
