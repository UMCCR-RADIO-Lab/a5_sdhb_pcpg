library(googlesheets4)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(patchwork)

################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_somatic_variants(quickload = T)
data_loader_cna_segs()

source("/g/data/pq08/projects/ppgl/a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
data_loader_a5_wts_counts(count_file_dir = "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/neb/gene", 
                          count_file_pattern = ".gene.counts")

##########
# Colors #
##########

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")


##################
# Read in pileup #
##################

E143_pileup <- readr::read_delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/paired_samples/count_summaries/E143.txt") %>% 
  pivot_longer(cols = c(-CHROM,-POS,-REF,-ALT), names_to = c("A5_ID","metric"), names_pattern = "(E143-...)_(.+)", values_to = "count") %>% 
  pivot_wider(names_from = metric, values_from = count) %>% mutate(A5_ID = gsub("-T0" , "-", A5_ID)) %>% 
  dplyr::rename(pileup_altCount = altCount, pileup_refCount = refCount, pileup_TotalCount = TotalCount) %>% 
  mutate(pileup_VAF = pileup_altCount/pileup_TotalCount)

E143_pileup <- E143_pileup %>% inner_join(
  E143_pileup %>%  
    dplyr::select(A5_ID, CHROM, POS, REF, ALT, pileup_VAF) %>% 
    pivot_wider(names_from = A5_ID, values_from = pileup_VAF) %>%  
    mutate(across(.cols=c(`E143-B01`, `E143-1`, `E143-2`, `E143-3`), .fns=~replace_na(.x,0))) %>% 
    mutate(Privacy= case_when(
      `E143-1` > 0 &  `E143-2` > 0 & `E143-3` > 0.3  ~ "Shared HVP",
      `E143-1` > 0 &  `E143-2` > 0 & `E143-3` > 0  ~ "Shared LVP",
      `E143-1` > 0 &  `E143-2` == 0 & `E143-3` > 0  ~ "P/M1",
      `E143-1` == 0 &  `E143-2` > 0 & `E143-3` > 0  ~ "P/M2",
      `E143-1` > 0 &  `E143-2` > 0 & `E143-3` == 0  ~ "M1/M2",
      `E143-1` == 0 &  `E143-2` == 0 & `E143-3` > 0  ~ "P",
      `E143-1` > 0 &  `E143-2` == 0 & `E143-3` == 0  ~ "M1",
      `E143-1` == 0 &  `E143-2` > 0 & `E143-3` == 0  ~ "M2"
    )) %>% 
    filter(!is.na(Privacy)) %>% dplyr::select(CHROM, POS, REF, ALT, Privacy))


E143_mut <- a5_somatic_variants %>% filter(A5_ID %in% c("E143-1","E143-2","E143-3")) %>% 
  full_join(E143_pileup, by=c("A5_ID", "seqnames"="CHROM", "start"="POS", "ALT")) %>% 
  mutate(end=ifelse(is.na(end), start, end)) 

#####################
# Annotate variants #
#####################

E143_mut <- E143_mut %>% 
  inner_join(a5_anno %>% dplyr::select(A5_ID, PublicationID, sample_purity)) %>% 
  mutate(Tumour_AF=as.numeric(Tumour_AF),
         PublicationID=factor(PublicationID, levels=c("E143-P1","E143-M1","E143-M2")),
         color=ifelse(!is.na(PCGR_SYMBOL) & PCGR_SYMBOL=="TERT", "red", "black"),
         linetype=ifelse(!is.na(PCGR_SYMBOL) & PCGR_SYMBOL=="TERT", 1, 2),
         linewidth=ifelse(!is.na(PCGR_SYMBOL) & PCGR_SYMBOL=="TERT", 2, 1),
         alpha=ifelse(!is.na(PCGR_SYMBOL) & PCGR_SYMBOL=="TERT", 1, 0.1),
         varid=paste(seqnames, start, ALT)) %>% 
  arrange(PublicationID, Tumour_AF) %>% ungroup() %>% 
  filter(pileup_VAF > 0)

##################
# Preprocess CNA #
##################

E143_segs <- a5_seg_keep %>% 
  filter(A5_ID %in% c("E143-1","E143-2","E143-3")) %>% 
  group_by(A5_ID) %>% 
  { 
    nm <- group_keys(.)[["A5_ID"]]
    . %>% group_split() %>% set_names(nm)
  }()


E143_segs_gr <- map(E143_segs, .f = \(x) GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = T))

##################
# Preprocess CNA #
##################

E143_mut_seg_annotated <- list()
for (s in names(E143_segs)) {
  
  temp <- E143_mut %>%  
    filter(A5_ID == s) %>% 
    dplyr::select(PCGR_SYMBOL, seqnames, start, end, A5_ID, PublicationID, TUMOR_AF, 
                  pileup_VAF, varid, color, 
                  linetype, linewidth, alpha, Privacy)
  temp <- GenomicRanges::makeGRangesFromDataFrame(df = temp, 
                                                  keep.extra.columns = T)
  
  
  
  hits <- findOverlaps(subject = temp, query = E143_segs_gr[[s]])
  
  mcols(temp)["cna_type"] <- NA
  mcols(temp)[["cna_type"]][subjectHits(hits)] <- mcols(E143_segs_gr[[s]])[["Class"]][queryHits(hits)]
  
  E143_mut_seg_annotated[[s]] <- GenomicRanges::as.data.frame(temp)
  
  
}
E143_mut_seg_annotated <- bind_rows(E143_mut_seg_annotated)
E143_mut_seg_annotated <- E143_mut_seg_annotated %>%  arrange(PublicationID, seqnames, start)
E143_mut_seg_annotated <- E143_mut_seg_annotated %>% 
  mutate(is_tert=ifelse(!is.na(PCGR_SYMBOL) & PCGR_SYMBOL=="TERT", "TERT", "Other"))

E143_mut_seg_annotated_keep <- E143_mut_seg_annotated %>% 
  inner_join(a5_somatic_variants_keep %>% 
               ungroup() %>% 
               filter(A5_ID %in% c("E143-1","E143-2","E143-3")) %>% 
               dplyr::select(seqnames, start, end, ALT) %>% distinct())


############
# Plotting #
############


#VAF Panel
pj <- position_jitter(width = 0.2, seed = 42)
gg_vaf <- ggplot(data = E143_mut_seg_annotated, 
                 mapping=aes(x=PublicationID, 
                             y=pileup_VAF, 
                             group=varid, 
                             linetype=linetype,
                             linewidth=linewidth)) + 
  geom_point(mapping = aes(color=cna_type),position = pj, size=0.5) +
  geom_path(position = pj, aes(color=is_tert, alpha=alpha)) + 
  geom_segment(data= E143_mut %>% dplyr::select(PublicationID, sample_purity) %>% distinct(), 
               mapping=aes(x=as.numeric(PublicationID) - 0.5,
                           xend=as.numeric(PublicationID) + 0.5,
                           y=sample_purity,
                           yend=sample_purity,
                           group=NULL, 
                           color=NULL,
                           linetype=NULL,
                           linewidth=NULL)) +
  scale_linetype_identity() +
  scale_linewidth_continuous(range = c(0.5,1.5)) +
  theme_bw() +
  scale_alpha_identity() +
  #scale_color_identity() +
  scale_color_manual(values=c(cna_palette, "TERT"="red")) +
  coord_cartesian(ylim = c(0,1)) + ylab("VAF") +
  guides(linewidth="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))

#C-circle panel
gg_tcr <- ggplot(a5_anno %>% 
                   filter(A5_ID %in% c("E143-1","E143-2","E143-3")) %>% 
                   mutate(telhunter_log2_telcontentratio=as.numeric(telhunter_log2_telcontentratio),
                          PublicationID=factor(PublicationID, 
                                               levels=c("E143-P1","E143-M1","E143-M2"))), 
                 aes(x=PublicationID, y= telhunter_log2_telcontentratio, fill=c_circle_result)) + 
  geom_col() + 
  scale_fill_manual(name="C-circle Status",values = c("No data"="darkgrey", 
                                                      "Negative"=ColorPalette[["DarkRed2"]], 
                                                      "Positive - Low"=ColorPalette[["DarkGreen1"]], 
                                                      "Positive"=ColorPalette[["LightGreen1"]])) +
  theme_bw() +
  ylab("Telomere (log2 T/N)") 


#CNA panel
zip <- function(...) {
  mapply(list, ..., SIMPLIFY = FALSE)
}
alternating_labels <- c(unlist(zip(paste0("chr", seq(1,20,2)),"")),"","","chrX")


#Recoding A5_IDs to publication IDs 
  name_recode <- 
    a5_anno %>% 
    filter(`Patient ID`=="E143") %>% 
    dplyr::select(A5_ID, PublicationID) %>% 
    { set_names(x = .$PublicationID,nm = .$A5_ID) }
  name_recode_levels <- c(sort(name_recode[grepl("P",name_recode)]),
                          sort(name_recode[grepl("R",name_recode)]),
                          sort(name_recode[grepl("M",name_recode)]))
  
  cna_plot <- ggplot(a5_seg_keep %>% 
                                            mutate(Class=forcats:::fct_recode(.f=Class,  
                                                                              Loss="Loss + Subclonal CNLOH",
                                                                              Gain="WGD+Gain",
                                                                              Gain="Gain+LOH",
                                                                              #None="Minor Subclonal Loss", 
                                                                              None="Diploid/Haploid-X")) %>% 
                                            filter(A5_ID %in% names(name_recode)) %>% 
                                            mutate(A5_ID=dplyr::recode(gsub("T0","",A5_ID),!!!name_recode),
                                                   A5_ID=factor(A5_ID, levels = name_recode_levels)), 
                                          aes(x=A5_ID, xend=A5_ID, y=start_offset, yend=end_offset, color=Class)) + 
    geom_segment(size=2.5) +
    geom_hline(data = chr_offsets, mapping=aes(yintercept=offset), size=0.1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
          axis.ticks.y = element_blank(),
          panel.grid = element_blank()) +
    scale_y_reverse(breaks = chr_offsets$offset, labels = alternating_labels, expand=c(0,0)) +
    scale_color_manual(values = cna_palette)   + 
    ylab("Copy Number")

# Merge plots
  
  gg_tcr + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + 
    gg_vaf +
    cna_plot +
    plot_layout(heights = c(1,2), widths = c(3,1), design = "A#\nBC", nrow = 2, ncol = 2, guides = "collect")
