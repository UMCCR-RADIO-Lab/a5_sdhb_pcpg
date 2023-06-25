setwd("/g/data/pq08/projects/ppgl")

source("./a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
library(patchwork)
library(ggplot2)
library(googlesheets4)
library(VariantAnnotation)
library(furrr)

threads=6
options(future.debug = FALSE)
options(future.globals.maxSize=10^10) 
plan(strategy="multisession", workers=threads) 

data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

#Fetch VCF names
# vcf_files <- list.files("./a5/wgs/symlinks/somatic_vcf/bcbio_ensemble",full.names = T)
# names(vcf_files) <- stringr::str_extract(vcf_files, "E...-.")
# 
# #Swap UMCCRised VCF (incomplete due to high TMB in E167-T01) with BCBio VCF
# vcf_files[["E167-1"]] <- "./a5/wgs/analysis/bcbio/E167/final/2021-10-28_E167/E167-1-ensemble-annotated.vcf.gz"
# 
# vcfs <- furrr::future_map(.x = vcf_files, .f = readVcf)

data_loader_somatic_variants(quickload = T)

rm(list=c("a5_somatic_variants_keep", "a5_somatic_variants_keep_curated", 
          "a5_somatic_variants_keep_pcawg", "a5_somatic_variants_keep_pcawg_brief"))

sample_combo <- crossing(Sample1=unique(a5_somatic_variants$A5_ID),
                         Sample2=unique(a5_somatic_variants$A5_ID)) %>% 
  filter(Sample1 != Sample2)

a5_somatic_variants_list <- a5_somatic_variants %>% 
  group_by(A5_ID) %>% 
  group_split() %>% 
  purrr::set_names(purrr::map_chr(., ~.x$A5_ID[1])) %>% 
  map(.f = function (vcf) {
    tumour <- gsub("-(.)","-T0\\1",vcf$A5_ID[1])
    normal <- gsub("-(.)","-B0\\1",vcf$A5_ID[1])
    return(vcf %>% ungroup() %>% 
             dplyr::select(seqnames, start, end, REF, ALT, PCGR_SYMBOL, Tumour_AF, Normal_AF, blacklist) %>% 
             dplyr::rename("{tumour}_AF":=Tumour_AF, "{normal}_AF":=Normal_AF) %>% 
             mutate(across(.cols = ends_with("AF"),
                           .fns =  as.numeric)))
  })



quick_compare <- function(s1, s2, a5_somatic_variants_list) {
  
  s1t_af_colname <- paste0(gsub("-(.)","-T0\\1",s1),"_AF")
  s2t_af_colname <- paste0(gsub("-(.)","-T0\\1",s2),"_AF")
  
  
  current_comp <- a5_somatic_variants_list[[s1]] %>% 
    full_join(a5_somatic_variants_list[[s2]], by=c("seqnames", "start", "end", "REF", "ALT", "PCGR_SYMBOL", "blacklist")) %>% 
    mutate(across(.cols = ends_with("AF"),
                  .fns =  \(x) replace_na(x,0))) %>% 
    mutate(called=case_when(
      !!sym(s1t_af_colname) > 0 & !!sym(s2t_af_colname) > 0 ~ "Both",
      !!sym(s1t_af_colname) > 0 ~ "s1",
      !!sym(s2t_af_colname) > 0 ~ "s2"
    )) %>% 
    filter(!!sym(s1t_af_colname) > 0 | !!sym(s2t_af_colname) > 0) #Weird zero VAF variants
  
  return_comp <- current_comp %>% mutate(pair=paste(s1,s2,sep="_")) %>%  
    group_by(pair,called, blacklist) %>% 
    dplyr::count() %>% 
    mutate(blacklist=ifelse(blacklist, "OnBlackList","NotOnBlackList")) %>% 
    pivot_wider(id_cols = pair, names_from = c(called,blacklist), values_from = n)
  
  return_comp$VAF <- list(current_comp %>% dplyr::select(all_of(c(s1t_af_colname, s2t_af_colname, "called"))))
  return(return_comp)
}

#all_comps <- readRDS("~/all_comps_temp.rds")

all_comps <- furrr::future_map2(sample_combo$Sample1,sample_combo$Sample2, .f = \(s1,s2) quick_compare(s1,s2,a5_somatic_variants_list))
all_comps <- bind_rows(all_comps)


all_comps_novaf <- all_comps %>% dplyr::select(-VAF)
all_comps_novaf[all_comps_novaf$pair=="E167-1_E167-2",c(2,3,4,5,6,7)] <- list(2416,42,1919103,147,5434,416)
all_comps_novaf$Both_NotOnBlackList  = replace_na(all_comps_novaf$Both_NotOnBlackList,0)
all_comps_novaf$Both_OnBlackList  = replace_na(all_comps_novaf$Both_OnBlackList,0)
all_comps_novaf <- all_comps_novaf %>% mutate(
  prop_overlap=(Both_NotOnBlackList+
                  Both_OnBlackList)/
    (Both_NotOnBlackList+
       Both_OnBlackList+
       s1_NotOnBlackList+
       s1_OnBlackList+
       s2_NotOnBlackList+
       s2_OnBlackList))

all_comps_novaf <- all_comps_novaf %>% separate(pair, into=c("s1name","s2name"), sep="_", remove = F)

#remove redundant A/B B/A comparisons
all_comps_novaf <- all_comps_novaf %>% 
  rowwise() %>% 
  mutate(key=paste(sort(x = c(s1name,s2name), decreasing = T), collapse="_")) %>% 
  ungroup() %>% filter(!duplicated(key))

all_comps_novaf <- all_comps_novaf %>% 
  mutate(Patient=ifelse(gsub("-.","",s1name)==gsub("-.","",s2name),gsub("-.","",s1name), "Unmatched"),
    Matched=ifelse(gsub("-.","",s1name)==gsub("-.","",s2name), "Matched", "Unmatched"))



all_comps_novaf <- all_comps_novaf %>% 
  inner_join(a5_anno %>% 
               dplyr::select(A5_ID, s1_type=is_primary_or_met, s1publabel=PublicationID),
             by=c("s1name"="A5_ID")) %>% 
  inner_join(a5_anno %>% 
               dplyr::select(A5_ID, s2_type=is_primary_or_met, s2publabel=PublicationID),
             by=c("s2name"="A5_ID")) %>% 
  mutate(comptype=ifelse(Matched=="Matched", 
                     paste(s1_type,s2_type, sep="_"), "Unmatched"))

all_comps_novaf <- all_comps_novaf %>% mutate(label=ifelse(Matched=="Matched", paste(s1publabel,s2publabel, sep="_"), NA))



ggplot(all_comps_novaf, aes(x=Patient,y=prop_overlap*100)) + geom_jitter(width=0.2)

gg_scatter_1v1 <- ggplot(all_comps_novaf,
       aes(x=Matched,y=log10(Both+0.5), label=label, color=Patient)) + 
  geom_jitter(width=0.2, height=0.05, alpha=0.1, size=0.01) + 
  ggrepel::geom_text_repel(size=3, max.overlaps = 100) +
  ylab("nVar found Both - log10") +theme_bw() +
  ggtitle("Number variants overlapping in 1v1 comparisons")

gg_scatter_1v1 <- ggplot(mapping=aes(x=Matched,y=log10(Both+0.5))) + 
  geom_violin(data = all_comps_novaf %>% filter(Matched=="Unmatched")) + 
  geom_point(data = all_comps_novaf %>% filter(Matched=="Matched"),
             mapping=aes(color=Patient), size=3) +
  ggrepel::geom_text_repel(data=all_comps_novaf %>% filter(Matched=="Matched"), 
                           mapping=aes(label=label), size=3, max.overlaps = 100) +
  ylab("nVar found Both - log10") + theme_bw() +
  ggtitle("Number variants overlapping in 1v1 comparisons")

gg_dist <- ggplot(all_comps_novaf %>% filter(Matched=="Unmatched") %>% group_by(Both) %>% dplyr::count() %>% ungroup() %>%  mutate(prop=n/sum(n)),
       aes(x=factor(Both),y=prop*100)) + geom_col() + 
  ylab("percentage of comparisons") +
  xlab("Number of overlapping variants") +
  ggtitle("Distribution of overlap counts in unmatched comparisons") + theme_bw()

gg_scatter_1v1 + gg_dist

ggplot(all_comps_novaf %>% mutate(label=ifelse(Patient=="Matched", pair, NA)), 
       aes(x=Both)) + 
  geom_histogram() + 
  scale_y_log10() + 
  xlab("nVar in common") +
  ylab("nComparisons") +
  ggtitle("histogram of frequency of overlap count across comparisons")
  
all_comps_novaf %>% mutate(prop_overlap_onblacklist=Both_OnBlackList/(Both_OnBlackList+Both_NotOnBlackList)) %>% 
ggplot(mapping = aes(x=prop_overlap_onblacklist, color=Matched)) + geom_histogram(stat = "count", binwidth = 0.1)
all_comps_novaf %>% mutate(prop_overlap_onblacklist=Both_OnBlackList/(Both_OnBlackList+Both_NotOnBlackList)) %>% 
  ggplot(mapping = aes(x=prop_overlap_onblacklist, color=Matched)) + geom_density(adjust=0.05)

ggplot(data= all_comps_novaf  %>%  filter(Matched=="Matched"), 
       mapping = aes(x=Both_OnBlackList, 
                     y=Both_NotOnBlackList, 
                     color=comptype)) + 
  geom_point() +
  ggrepel::geom_text_repel(mapping = aes(label=pair), size=2) +
  scale_color_manual(values=RColorBrewer::brewer.pal(n = 8,name = "Set1")[c(1:5,7)]) +
  xlab("n_On_Blacklist")+ 
  ylab("n_Not_On_Blacklist") + 
  theme_bw() + 
  ggtitle("Matched")


ggplot(data= all_comps_novaf  %>%  filter(Matched=="Unmatched"), 
       mapping = aes(x=Both_OnBlackList, 
                     y=Both_NotOnBlackList)) + 
  geom_jitter(size=0.5, alpha=0.5) +
  xlab("n_On_Blacklist")+ 
  ylab("n_Not_On_Blacklist") + theme_bw() + ggtitle("Unmatched")


source(paste0(basedir,"/wgs/scripts/utility/compare_vcf_calls.r"))

Samples <- list(
  #E122=c("E122-T01","E122-T02"),
  #E128=c("E128-T01","E128-T02"),
  #E132=c("E132-T01","E132-T02"),
  #E143=c("E143-T01","E143-T02","E143-T03"),
  #E146=c("E146-T01","E146-T02"),
  #E158=c("E158-T01","E158-T02"),
  E159=c("E159-T01","E159-T02","E159-T03","E159-T04")
  #E166=c("E166-T01","E166-T02"),
  #E167=c("E167-T01","E167-T02"),
  #E169=c("E169-T01","E169-T02"),
  #E225=c("E225-T01","E225-T02"),
  #E229=c("E229-T01","E229-T02")
)

#Pairwise matching of VCFs - returns a list of dataframes with variant calls marked present/absent for each sample from each patient

deleterious_coding_consequences <- c(                  
  "splice_acceptor_variant",
  "splice_region_variant",
  "splice_donor_variant",
  "stop_gained",
  "stop_lost",
  "start_lost",
  "missense_variant",
  "frameshift_variant",
  "inframe_insertion",
  "inframe_deletion")

deleterious_coding_consequences.regex <- paste(deleterious_coding_consequences, collapse ="|")

matched_variant_tables <- 
  furrr::future_map(.x = Samples, 
                    .f = function(patient_samples, delregex) {
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
                      
                      #Swap UMCCRised VCF (incomplete due to high TMB in E167-T01) with BCBio VCF
                      if(patient == "E167")
                      {
                        vcf_files[["E167-T01"]] <-
                          paste0(basedir,"/wgs/analysis/bcbio/E167/final/2021-10-28_E167/E167-1-ensemble-annotated.vcf.gz")
                        vcf_files[["E167-T02"]] <-
                          paste0(basedir,"/wgs/analysis/bcbio/E167/final/2021-10-28_E167/E167-2-ensemble-annotated.vcf.gz")
                      }
                      
                      
                      return_table <-
                        multi_compare_small_variant_vcf(
                          vcfs = vcf_files,
                          normal_names = gsub("T.+", "B01", patient_samples),
                          tumour_names = patient_samples,
                          genome = "BSgenome.Hsapiens.UCSC.hg38",
                          Keep_Only_Highest_Consequence = T,
                          collapse_annotation = T)
                      
                      
                      
                      
                      # #Override missed TERT mutation due to UMCCRISE filtering it out
                      # if (patient == "E167") {
                      #   temp <- return_table
                      #   temp[temp$seqnames == "chr5" &
                      #          temp$start == 1295113, "E167.T01_Called"] <- 'Y'
                      #   return_table <- temp
                      #   rm(temp)
                      # }
                      message("Completed variant matching for: ", patient)
                      
                      return(return_table)
                    }, delregex=deleterious_coding_consequences.regex)

#Fix and samples names that have dashes replaced with dots
matched_variant_tables <-
  lapply(matched_variant_tables, function (table) {
    return(rename_with(
      table,
      .fn = ~ gsub("(E[0-9]{3}).([TB][0-9]{2})", "\\1-\\2", .x),
      .cols = matches("_Called$|_VAF$")
    ))
  })


#Annotate black listed variants
blacklist="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/blacklists/blacklist_readsupport_gt3_samplesupport_gt3.tsv"
blacklisted_variants <- read.delim(blacklist)

matched_variant_tables <-
  lapply(matched_variant_tables, function (table, black_list = blacklisted_variants) {
    table <- table %>% 
      left_join(black_list %>% 
                  dplyr::rename("blacklist_n_normal_observed"=n_above_threshold) %>% 
                  mutate(blacklist=TRUE),
                by = c("seqnames", "start", "end", "ALT"))
    
    table$blacklist_n_normal_observed[is.na(table$blacklist_n_normal_observed)] <- 0
    table$blacklist[is.na(table$blacklist)] <- F
    return(table)
  })

E159_venn <- list(`E159-T01`=matched_variant_tables$E159 %>% 
       filter(`E159-T01_Called`=="Y") %>% 
       transmute(id=paste(seqnames, start,end, ALT, sep="-")) %>% 
       pull(id),
     `E159-T02`=matched_variant_tables$E159 %>% 
       filter(`E159-T02_Called`=="Y") %>% 
       transmute(id=paste(seqnames, start,end, ALT, sep="-")) %>% 
       pull(id),
     `E159-T03`=matched_variant_tables$E159 %>% 
       filter(`E159-T03_Called`=="Y") %>% 
       transmute(id=paste(seqnames, start,end, ALT, sep="-")) %>% 
       pull(id),
     `E159-T04`=matched_variant_tables$E159 %>% 
       filter(`E159-T04_Called`=="Y") %>% 
       transmute(id=paste(seqnames, start,end, ALT, sep="-")) %>% 
       pull(id))

ggVennDiagram::ggVennDiagram(E159_venn, category.names = names(E159_venn))


all_vafs <- tibble(s1=vector(mode="character"),
                       s2=vector(mode="character"),
                       patient=vector(mode="character"),
                       variant=vector(mode="character"),
                       s1_VAF=vector(mode="numeric"),
                   s2_VAF=vector(mode="numeric"),
                   blacklist=vector(mode="logical"), 
                   blacklist_n_normal_observed=vector(mode="numeric"))
for (patient in names(Samples))
{
  pairs <- make_pairs(Samples[[patient]]) #make unique combinations
  
  for (i in 1:nrow(pairs))
  {
    s1_vaf_col <- paste0(pairs[i,"Sample1"],"_VAF")
    s2_vaf_col <- paste0(pairs[i,"Sample2"],"_VAF")
    
    temp <- matched_variant_tables[[patient]] %>% 
      ungroup() %>% 
      dplyr::rename("s1_VAF"=!!sym(s1_vaf_col),"s2_VAF"=!!sym(s2_vaf_col)) %>% 
      mutate(variant=paste(seqnames,start, ALT)) %>% 
      dplyr::select(variant, s1_VAF, s2_VAF, blacklist, blacklist_n_normal_observed) %>% 
      mutate(s1=pairs[i,"Sample1"], 
             s2=pairs[i,"Sample2"],
             patient=patient)
    all_vafs <- bind_rows(all_vafs, temp)
  }
  
  
}

all_vafs <- all_vafs %>% mutate(s1= gsub("-T0","-", s1), s2= gsub("-T0","-", s2)) %>% 
  inner_join(a5_anno %>% 
               dplyr::select(A5_ID, s1_type=is_primary_or_met, s1publabel=PublicationID),
             by=c("s1"="A5_ID")) %>% 
  inner_join(a5_anno %>% 
               dplyr::select(A5_ID, s2_type=is_primary_or_met, s2publabel=PublicationID),
             by=c("s2"="A5_ID")) %>% mutate(comptype=paste(s1_type,s2_type, sep="_")) %>% 
  mutate(pair=paste(s1publabel,s2publabel, sep="_"))
  

repeat_variants <- a5_somatic_variants %>% 
  ungroup() %>% 
  mutate(variant=paste(seqnames,start, ALT)) %>% 
  mutate(patient=gsub("-.","",A5_ID)) %>% 
  dplyr::select(patient, variant) %>% 
  distinct() %>% 
  group_by(variant) %>% 
  dplyr::count() %>% 
  filter(n>1)

all_vafs <- all_vafs %>% left_join(repeat_variants, by="variant") %>% 
  mutate(n=replace_na(n,1))



vaf_summary_plots <- list()
for (current_patient in c("E128","E159","E229"))
{
  pairs <- all_vafs %>% 
    filter(comptype=="Primary_Primary", patient==current_patient) %>% 
    dplyr::select(s1,s2) %>% distinct()
  
  for (i in 1:nrow(pairs))
  {
    current_s1 = pairs$s1[i]
    current_s2 = pairs$s2[i]
    
    s1_purity <- a5_anno %>% filter(A5_ID==current_s1) %>% pull(sample_purity)
    s2_purity <- a5_anno %>% filter(A5_ID==current_s2) %>% pull(sample_purity)
    
    s1_pubid <- a5_anno %>% filter(A5_ID==current_s1) %>% pull(PublicationID)
    s2_pubid <- a5_anno %>% filter(A5_ID==current_s2) %>% pull(PublicationID)

    
    gg_s1_dens <- ggplot(all_vafs %>% filter(s1==current_s1) %>% 
             filter(comptype=="Primary_Primary", s1_VAF > 0), 
           aes(x=s1_VAF)) + 
      geom_density(adjust=0.5) + theme_bw()  + 
      coord_cartesian(xlim = c(0,1))
    
    gg_s2_dens <- ggplot(all_vafs %>% filter(s2==current_s2) %>% 
             filter(comptype=="Primary_Primary", s2_VAF > 0), 
           aes(y=s2_VAF)) + 
      geom_density(adjust=0.5) + theme_bw() + scale_x_reverse() +
      coord_cartesian(ylim = c(0,1))
    
    gg_s1_s2_vaf <- ggplot(all_vafs %>% filter(comptype=="Primary_Primary", patient==current_patient), 
           aes(x=s1_VAF, y=s2_VAF,  color=factor(n), shape=blacklist)) + 
      geom_point() + 
      geom_vline(xintercept=s1_purity, linetype=2) +
      geom_hline(yintercept=s2_purity, linetype=2) +
      scale_color_discrete(name="nPatients observed") +
      scale_shape_manual(name="Blacklist", values=c(`TRUE`=8, `FALSE`=16)) +
      coord_cartesian(xlim = c(0,1),ylim = c(0,1))
    
    vaf_summary_plots[[paste(current_s1, current_s2, sep="_")]] <- gg_s2_dens + 
      gg_s1_s2_vaf + 
      plot_spacer() + 
      gg_s1_dens  + 
      plot_layout(nrow = 2, ncol = 2, widths=c(1,3), heights = c(3,1)) + 
      plot_annotation(title = paste(s1_pubid, "vs", s2_pubid))
    

  }
  
  
}


matched_variant_tables$E229 %>% filter(`E229-T01_Called`=="Y",`E229-T02_Called`=="Y")


all_vafs %>% 
  filter(comptype=="Primary_Primary", patient==current_patient, s1_VAF>0.05, s2_VAF>0.05, blacklist==F) %>% 
  separate(variant, into = c("CHROM","POS","ALT")) %>% 
  dplyr::select("CHROM","POS","ALT") %>% 
  mutate(POS=as.numeric(POS)) %>% 
  inner_join(pileup_summary_alt_only) %>% View()
