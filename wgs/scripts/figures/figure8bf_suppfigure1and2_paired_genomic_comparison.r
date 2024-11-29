library(patchwork)
library(ggplot2)
library(googlesheets4)
library(furrr)

setwd("/g/data/pq08/projects/ppgl/a5")


##################
# Helper scripts #
##################

source("./utility_scripts/compare_vcf_calls.r")

###############
# Data loader #
###############

source("./wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_cna_segs()
data_loader_sv_gridsslinx(quickload = T)

source("./sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

##########
# Colors #
##########

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

################
# Dictionaries #
################

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

#Classes of CNA events
cn_event_types <- c("Diploid/Haploid-X", "Loss",  "Subclonal Loss", "Minor Subclonal Loss",  
                    "Loss + Subclonal CNLOH", "CNLOH", "Gain", "Subclonal Gain", "Gain+LOH", "WGD", "WGD+Gain",
                    "Hom. Del.", "Other", "Chromothripsis")

#Genes of interest for annotation
GOI <- read.delim("/g/data/pq08/reference/gene_lists/pcawg_mutational_drivers.txt") %>% pull(Gene) %>% as.character()
GOI.E167 <- data.frame(gene=c("NRAS", "TERT", "MLH1"),
                       position=c(114716126, 1295113, 37050485))

#List of A5 patients with paired samples
Samples <- list(
  E122=c("E122-T01","E122-T02"),
  E128=c("E128-T01","E128-T02"),
  E132=c("E132-T01","E132-T02"),
  E136=c("E136-T01","E136-T02"),
  E143=c("E143-T01","E143-T02","E143-T03"),
  E146=c("E146-T01","E146-T02"),
  E158=c("E158-T01","E158-T02"),
  E159=c("E159-T01","E159-T02","E159-T03","E159-T04"),
  #E166=c("E166-T01","E166-T02"),
  E167=c("E167-T01","E167-T02"),
  E169=c("E169-T01","E169-T02"),
  E225=c("E225-T01","E225-T02"),
  E229=c("E229-T01","E229-T02")
)



## Match Variants between VCFs

threads=4
plan(strategy = "multisession", workers=threads)

#Pairwise matching of VCFs - returns a list of dataframes with variant calls marked present/absent for each sample from each patient

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
blacklist_normal_file <- "./wgs/analysis/mpileup/blacklist/blacklists/blacklist_readsupport_gteq3_samplesupport_gteq3.tsv"
blacklist_normal <- read.delim(blacklist_normal_file)
blacklist_tumour_file <- "./wgs/analysis/blacklist_from_vcf/output/blacklist_reject_pr_ratio_gt1.tsv"
blacklist_tumour <- read.delim(blacklist_tumour_file)

blacklisted_variants <-blacklist_tumour %>%  
  filter(pass_reject_ratio > 1) %>% 
  left_join(blacklist_normal, 
            by=c("Chr"="CHROM", "Pos"="POS", "Alt"="ALT")) %>% 
  dplyr::select(Chr, Pos, Alt) %>%  mutate(blacklist=TRUE)


matched_variant_tables <-
  lapply(matched_variant_tables, function (table, black_list = blacklisted_variants) {
    table <- table %>% 
      anti_join(black_list,
                by = c("seqnames"="Chr", "start"="Pos", "ALT"="Alt"))
    
    return(table)
  })


plan(strategy = "sequential")

#####
# Merge the annotation from the UMCCRised vcf (cancer genes only)
#####

#UMCCRised VCFs
vcf_files <- paste0(basedir, "/wgs/analysis/bcbio/",
                    "E167",
                    "/umccrised/",
                    gsub("T0","",Samples[["E167"]]),
                    "__",
                    Samples[["E167"]],
                    "/small_variants/",
                    gsub("T0","",Samples[["E167"]]),
                    "__",
                    Samples[["E167"]],"-somatic-PASS.vcf.gz")
names(vcf_files) <- Samples[["E167"]]

#Read in VCFs and extract annotation information
E167_umccrised_data <- purrr::map(vcf_files, readVcf) %>% 
  purrr::map(function(vcf_data) {
    data.frame(seqnames=seqnames(vcf_data), 
               start=start(vcf_data),
               end=end(vcf_data),
               info(vcf_data)[c("PCGR_SYMBOL","PCGR_TIER","PCGR_CONSEQUENCE")])
  }) %>% 
  bind_rows() %>% 
  distinct()

#Annotate matched data with UMCCRised subset of VCFs (cancer genes only due to high TMB)
matched_variant_tables$E167 <- 
  matched_variant_tables$E167 %>% 
  left_join(E167_umccrised_data)


##Guess at annotation for variants not in UMCCRise
matched_variant_tables$E167 <- matched_variant_tables$E167 %>% 
  distinct() %>% 
  group_by(seqnames) %>%  
  group_split() %>% 
  furrr::future_map(.f = function(chr_data) {
    chr_data %>% rowwise() %>% mutate(PCGR_TIER=ifelse(
      any(stringr::str_detect(pattern = deleterious_coding_consequences, string=Annotation)), 
      "TIER_4", 
      "NONCODING")) }) %>% bind_rows()

#########################################
## Preprocess data for input to ggplots #
#########################################

plan(strategy = "multisession", workers=threads)

#Result lists
comparison_grids <- list()
dotplot_data <- list()

for (patient in names(Samples))
{
  #####################
  ## Comparison Grid ##
  #####################
  #nSample by nSample - intersection show count of common variants
  
  pairs <- make_pairs(Samples[[patient]]) #make unique combinations
  pairs <- rbind(pairs, cbind(Sample1=Samples[[patient]],Sample2=Samples[[patient]])) #add self pairings
  pairs$nShared <- 0
  
  #count events present in specified pairings
  for(i in 1:nrow(pairs))
  {
    pairs$nShared[[i]] <- matched_variant_tables[[patient]] %>% ungroup() %>% 
      filter(!!sym(paste(pairs[i,1],"Called", sep="_")) == "Y", 
             !!sym(paste(pairs[i,2],"Called", sep="_")) == "Y") %>% 
      dplyr::count() %>% pull(n)
  }
  
  #Make reciprocal pairs to fill in grid and sort/factorise to ensure order
  pairs <- unique(rbind(pairs, data.frame(Sample1=pairs$Sample2,Sample2=pairs$Sample1, nShared=pairs$nShared)))
  pairs$Sample1 <- as.character(pairs$Sample1)
  pairs$Sample2 <- as.character(pairs$Sample2)
  pairs$Sample1 <- factor(pairs$Sample1, levels = sort(unique(pairs$Sample1)))
  pairs$Sample2 <- factor(pairs$Sample2, levels = sort(unique(pairs$Sample2)))
  
  comparison_grids[[patient]] <- pairs
  
  ######################################
  ## Variant Membership Dot Plot Data ##
  ######################################
  
  
  matched_variant.longform <- matched_variant_tables[[patient]] %>% 
    ungroup() %>% 
    dplyr::select(PCGR_SYMBOL, seqnames, start, end, REF, ALT,PCGR_TIER, matches(c("E...-T0._Called", "E...-T0._VAF"))) %>% 
    mutate(across(.cols = ends_with("_Called"), .fns=~ifelse(.x=="Y", 1,0))) %>% #convert called to numeric so it rbinds with VAF
    #Pivot into long form separating  sample name and data-type into columns 
    pivot_longer(cols = matches(c("_Called$", "T[0-9]{2}_VAF$")), 
                 names_sep = "_",
                 names_to = c("A5_ID", "DataType"), 
                 values_to = "Value") %>% 
    #Pivot wider to unmix datatypes but keep samples in longform
    pivot_wider(names_from = DataType, values_from = Value)  %>% 
    #Filter to only call events
    filter(Called==1) %>% 
    group_by(seqnames) %>% group_split() %>% 
    furrr::future_map(.f=function(chr_data){
      chr_data %>% group_by(seqnames, start, end, REF, ALT) %>% 
        mutate(Privacy=paste(A5_ID,sep = "/", collapse = "/"),
               nSamplesCalled=n()) %>%  ungroup()}) %>% bind_rows()
  
  
  dotplot_data[[patient]] <- matched_variant.longform  %>% 
    mutate(PCGR_SYMBOL=replace_na(PCGR_SYMBOL, "")) %>% 
    ## Create annotation groups based on variant tiers and genes of interest list
    ## E167 is treated separately due to extreme mutation load
    mutate(Annotation=case_when((
      is.na(PCGR_TIER) | PCGR_TIER == "NONCODING") ~ "Non-Coding/Syn.",
      PCGR_TIER %in% c("TIER_1","TIER_2","TIER_3") ~ PCGR_SYMBOL,
      PCGR_TIER == "TIER_4" & PCGR_SYMBOL %in% GOI ~ PCGR_SYMBOL,
      PCGR_TIER == "TIER_4" ~ "TIER_4"),
      #handle E167-1 high TMB special case
      Annotation=case_when(
        A5_ID %in% c("E167-T01","E167-T02") & (PCGR_SYMBOL %in% GOI.E167$gene & start %in% GOI.E167$position) ~ PCGR_SYMBOL,
        #Keep annotation for variants also present in E167-2
        A5_ID == "E167-T01" & nSamplesCalled==2 ~ Annotation,
        A5_ID == "E167-T01" & Annotation=="Non-Coding/Syn." ~ Annotation,
        A5_ID == "E167-T01" ~ PCGR_TIER,
        TRUE ~ Annotation),
      AnnotationClass=ifelse(Annotation %in% c("TIER_1", "TIER_2", 
                                               "TIER_3", "TIER_4",
                                               "Non-Coding/Syn."),"Summary","Specific"),
      AnnotationClass=factor(as.character(AnnotationClass), levels = c("Specific","Summary")),
      VAF=ifelse(AnnotationClass=="Specific",VAF,NA),
      Position=ifelse(AnnotationClass=="Specific",start,""))
  
  #Reduce annotation Specificity to just coding/noncoding - optional
  dotplot_data[[patient]] <- dotplot_data[[patient]] %>%  mutate(Annotation=gsub("TIER_.", "Coding", Annotation))
  
  dotplot_data[[patient]] <- dotplot_data[[patient]] %>%
    #Keep only one representative entry for each annotation group and annotate with variant count
    group_by(Annotation, Position, A5_ID, Privacy) %>% mutate(variant_group_size=n()) %>% slice_head(n=1) %>%
    arrange(desc(AnnotationClass), desc(nSamplesCalled), Privacy) %>%
    mutate(Annotation=ifelse(AnnotationClass == "Summary", paste0(Annotation," (n=",variant_group_size,")"), Annotation)) %>%
    mutate(Annotation=gsub("TIER_", "Tier ", Annotation)) %>%
    ungroup() %>% 
    #Handle edge cases where the same gene is mutated in multiple samples at different positions
    mutate(dupSymbol=(duplicated(PCGR_SYMBOL) | duplicated(PCGR_SYMBOL, fromLast=T)),
           dupSymPos=(duplicated(paste(PCGR_SYMBOL,Position)) | duplicated(paste(PCGR_SYMBOL,Position), fromLast=T)),
           Annotation=ifelse(AnnotationClass=="Specific" & dupSymbol & !dupSymPos, 
                             paste0(Annotation, "(g.",Position,")"),
                             Annotation)) %>% 
    mutate(Annotation=factor(Annotation, levels = unique(.$Annotation))) %>%
    mutate(LabelA=paste(Annotation, Privacy, sep="_"), LabelA=factor(LabelA, levels = unique(LabelA)))
  
}

################
# VariantPlots #
################

variant_dot_plots <- list()
variant_grid_plots <- list()
composed_plots <- list()
for (patient in names(Samples))
{
  #Recoding A5_IDs to publication IDs 
  name_recode <- 
    a5_anno %>% 
    filter(`Patient ID`==patient) %>% 
    dplyr::select(A5_ID,PublicationID) %>% 
    { set_names(x = .$PublicationID,nm = .$A5_ID) }
  name_recode_levels <- c(sort(name_recode[grepl("P",name_recode)]),
                          sort(name_recode[grepl("R",name_recode)]),
                          sort(name_recode[grepl("M",name_recode)]))
  
  comparison_grid_current <- comparison_grids[[patient]] %>% 
    mutate(across(.cols = c(Sample1,Sample2), 
                  .fns = \(x){ factor(dplyr::recode(gsub("T0","",x),
                                                    !!!name_recode), 
                                      levels=name_recode_levels)}))
  
  dotplot_current <- dotplot_data[[patient]] %>% 
    mutate(A5_ID= factor(dplyr::recode(gsub("T0","",A5_ID),!!!name_recode), 
                         levels=name_recode_levels))
  
  #Make plots
  variant_grid_plots[[patient]] <-
    ggplot(comparison_grid_current, aes(x = Sample1, y = Sample2, label =
                                          nShared)) +
    geom_text(size = 3) +
    geom_hline(yintercept = seq(1.5, length(unique(
      comparison_grid_current$Sample1
    ))), color = "grey") +
    geom_vline(xintercept = seq(1.5, length(unique(
      comparison_grid_current$Sample1
    ))), color = "grey") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_blank(),
      #axis.text.x=element_blank(),
      axis.text.x = element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1
      ),
      #axis.ticks.x=element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_y_discrete(limits = rev(levels(comparison_grid_current$Sample2))) #+
  #scale_x_discrete(position = "top")
  
  
  variant_dot_plots[[patient]] <-
    ggplot(dotplot_current, aes(x = A5_ID, y = LabelA)) +
    geom_point() +
    theme_bw() +
    theme(
      #axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    )  + 
    #This is to handle when two different privacy states have the same Annotation, 
    #using a combined privacy/annotation label prevents them colliding,
    #the combined label is then overwritten by the annotation alone
    scale_y_discrete(labels=dotplot_current %>% ungroup() %>% dplyr::select(LabelA, Annotation) %>% distinct() %>% pull(Annotation))
 
  composed_plots[[patient]] <-
    variant_dot_plots[[patient]] + theme(plot.margin = margin(6, 6, 0, 6)) +
    variant_grid_plots[[patient]] +
    plot_layout(nrow = 2, heights = c(1, 0.3))
}

#############
# CNA plots #
#############


zip <- function(...) {
  mapply(list, ..., SIMPLIFY = FALSE)
}
alternating_labels <- c(unlist(zip(paste0("chr", seq(1,20,2)),"")),"","","chrX")

cna_plots_vertical <- list()
for (patient in names(Samples))
{
  
  #Recoding A5_IDs to publication IDs 
  name_recode <- 
    a5_anno %>% 
    filter(`Patient ID`==patient) %>% 
    dplyr::select(A5_ID,PublicationID) %>% 
    { set_names(x = .$PublicationID,nm = .$A5_ID) }
  name_recode_levels <- c(sort(name_recode[grepl("P",name_recode)]),
                          sort(name_recode[grepl("R",name_recode)]),
                          sort(name_recode[grepl("M",name_recode)]))
  
  cna_plots_vertical[[patient]] <- ggplot(a5_seg_keep %>% 
                                            mutate(Class=forcats:::fct_recode(.f=Class,  
                                                                              Loss="Loss + Subclonal CNLOH",
                                                                              Gain="WGD+Gain",
                                                                              Gain="Gain+LOH",
                                                                              #None="Minor Subclonal Loss", 
                                                                              None="Diploid/Haploid-X")) %>% 
                                            filter(A5_ID %in% gsub("T0","",Samples[[patient]])) %>% 
                                            mutate(A5_ID=dplyr::recode(gsub("T0","",A5_ID),!!!name_recode),
                                                   A5_ID=factor(A5_ID, levels = name_recode_levels)), 
                                          aes(x=A5_ID, xend=A5_ID, y=start_offset, yend=end_offset, color=Class)) + 
    geom_segment(linewidth=2.5) +
    geom_hline(data = chr_offsets, mapping=aes(yintercept=offset), linewidth=0.1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
          axis.ticks.y = element_blank(),
          panel.grid = element_blank()) +
    scale_y_reverse(breaks = chr_offsets$offset, labels = alternating_labels, expand=c(0,0)) +
    scale_color_manual(values=cna_palette) + 
    ylab("Copy Number")
}

############################
# Structural variant plots #
############################

sv_comparison <- list()
sv_dotplot_data <- list()
for (patient in names(Samples))
{
  pairs <- make_pairs(Samples[[patient]])
  
  pair_comps <- list()
  
  for (i in 1:nrow(pairs))
  {
    
    pair_svs <- A5_gridss_keep %>%  filter(A5_ID %in% pairs[i, ])
    
    pair_comps[[paste(pairs[i, 1], pairs[i, 2], sep = "_")]] <-
      pair_svs %>% 
      dplyr::select(A5_ID, chrom1, start1, end1, A5_ID, chrom2, start2, end2) %>% 
      distinct() %>%
      mutate(end2=replace_na(end2,0)) %>% 
      pivot_wider(id_cols = c(chrom1, start1, end1, chrom2, start2, end2), 
                  names_from = A5_ID, values_from = A5_ID) %>% 
      mutate("{pairs[i, 1]}":=ifelse( is.na(!!sym(pairs[i, 1])), F, T),
             "{pairs[i, 2]}":=ifelse( is.na(!!sym(pairs[i, 2])), F, T)) %>% 
      dplyr::rename("{paste(pairs[i, 1],'Called', sep='_')}" := pairs[i, 1]) %>% 
      dplyr::rename("{paste(pairs[i, 2],'Called', sep='_')}" := pairs[i, 2]) 
    
  }

  #Merge 1v1 comparisons together where more than 2 tumours
  sv_comparison[[patient]] <-
    purrr::reduce(
      .x = pair_comps,
      .f = function(comp1, comp2) {
        comp1 %>%
          full_join(comp2) %>%
          mutate(across(matches("_Called"),
                        ~replace_na(.x, FALSE)))
      }
    )
  
  #Format dot plot data
  sv_dotplot_data[[patient]] <- sv_comparison[[patient]] %>% 
    group_by(across(matches("_Called"))) %>% 
    dplyr::count() %>%
    unite(col="Combo", matches("_Called"), sep="/", remove=F) %>% 
    mutate(Annotation=paste0("SV (n=",n,")")) %>% 
    dplyr::select(-n) %>% 
    pivot_longer(cols = c(-Annotation,-Combo), names_to="A5_ID", values_to="Present") %>% filter(Present) %>% 
    mutate(LabelA=paste(Annotation, Combo, sep="_")) %>% 
    group_by(LabelA) %>% 
    mutate(nSamples=n()) %>% arrange(desc(nSamples)) %>% 
    mutate(LabelA=factor(LabelA, levels=unique(.$LabelA))) %>% 
    mutate(A5_ID=factor(gsub("_Called", "", A5_ID), levels=sort(Samples[[patient]])))
  
}

sv_dot_plots <- list()
for (patient in names(Samples))
{
  
  #Recoding A5_IDs to publication IDs 
  name_recode <- 
    a5_anno %>% 
    filter(`Patient ID`==patient) %>% 
    dplyr::select(A5_ID,PublicationID) %>% 
    { set_names(x = .$PublicationID,nm = .$A5_ID) }
  
  name_recode_levels <- c(sort(name_recode[grepl("P",name_recode)]),
                          sort(name_recode[grepl("R",name_recode)]),
                          sort(name_recode[grepl("M",name_recode)]))
  
  sv_dotplot_current <- sv_dotplot_data[[patient]] %>% 
    mutate(A5_ID=factor(dplyr::recode(gsub("T0","",A5_ID), !!!name_recode), 
                        levels = name_recode_levels))
  
  sv_dot_plots[[patient]] <-
    ggplot(sv_dotplot_current, aes(x = A5_ID, y = LabelA)) +
    geom_point() +
    theme_bw() +
    theme(
      #axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.x = element_blank(),
      axis.title.y = element_blank(),
      axis.title.x = element_blank()
    ) + 
    scale_x_discrete(drop=F) +
    #This is to handle when two different privacy states have the same Annotation, 
    #using a combined privacy/annotation label prevents them colliding,
    #the combined label is then overwritten by the annotation alone
    scale_y_discrete(labels=sv_dotplot_data[[patient]] %>% ungroup() %>% dplyr::select(LabelA, Annotation) %>% distinct() %>% pull(Annotation))
}

##################
# Composed Plots #
##################

composed_plots <- list()
for (patient in names(Samples))
{
  composed_plots[[patient]] <- 
    variant_dot_plots[[patient]] + 
    theme(plot.margin = margin(6,6,0,6), axis.ticks.length.x = unit(0.1,"cm")) + 
    variant_grid_plots[[patient]] + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.length.x = unit(0.25,"cm"), 
          plot.margin = margin(0,6,0,6)) + 
    sv_dot_plots[[patient]] + 
    theme(axis.text.x.top = element_blank(), 
          axis.ticks.length.x = unit(0.25,"cm"), 
          plot.margin = margin(0,6,0,6)) + 
    cna_plots_vertical[[patient]] + 
    theme(axis.ticks.y = element_blank(), 
          plot.margin = margin(0,6,6,6)) +
    plot_layout(nrow=4, heights = c(1.3, 0.3, 0.5, 2))
  print(composed_plots[[patient]])
}

# Write Plots
plot_config <- 
  as_tibble(
    matrix(
      c(
        "E122", 7,5,
        "E128", 6,5,
        "E132", 5,5,
        "E136", 5,5,
        "E143", 8,5.5,
        "E146", 5,5,
        "E158", 6,5,
        "E159", 8,5.5,
        "E166", 7,5,
        "E167", 6,5,
        "E169", 6.5,5,
        "E225", 6.5,5,
        "E229", 5,5
      ),
      byrow = T, 
      ncol = 3, 
      dimnames = list(c(),c("PatientID","height","width"))))

purrr::pwalk(.l = plot_config, 
             .f = \(PatientID, height, width) {
               ggsave(composed_plots[[PatientID]], 
                      filename = file.path(basedir,"wgs/results/paired_sample_analysis/",
                                           paste0(PatientID,
                                                  "_pairplot_composed.pdf")
                      ),
                      height = as.numeric(height),
                      width =  as.numeric(width),
                      units = "in")
             }
)
