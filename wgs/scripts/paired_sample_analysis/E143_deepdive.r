

log2cpm  %>% 
  filter(symbol==g) %>%  
  left_join(a5_anno %>%  
              dplyr::select(A5_ID, `Patient ID`, TERT_ATRX_Mutation, 
                            TERT_ATRX_Mutation_Event, telhunter_log2_telcontentratio, 
                            Primary_Location_Simplified , differential_group_sampletype_strict, c_circle_result)) %>% 
  ggplot(aes(x=as.numeric(telhunter_log2_telcontentratio), 
             y=log2_cpm_z, 
             color=c_circle_result, 
             fill=TERT_ATRX_Mutation, 
             group=`Patient ID`)) + 
  geom_point(size=3, shape=21, stroke=2) + 
  geom_line( alpha=0.3, size=0.5) +
  #geom_text(aes(label=A5_ID)) + 
  scale_fill_manual(values = driver_cols) +
  scale_color_manual(values = c("No data"="darkgrey", 
                                "Negative"=ColorPalette[["DarkRed2"]], 
                                "Positive - Low"=ColorPalette[["DarkGreen1"]], 
                                "Positive"=ColorPalette[["LightGreen1"]])) +
  #scale_shape_manual(values = c(8,18,15,17,19,4)) +
  # scale_size_manual(values=c(1.3,2.2,2.2)) +
  theme_bw() +
  ylab(paste(g," Expression (z-scaled)")) +
  xlab("Tumour/Normal Telomere Content Ratio (log2)")


g="TERT"
pj <- position_jitter(width = 0.2, seed = 42)
log2cpm  %>% 
  filter(symbol==g) %>%  
  left_join(a5_anno %>%  
              dplyr::select(A5_ID, `Patient ID`, TERT_ATRX_Mutation, 
                            TERT_ATRX_Mutation_Event, telhunter_log2_telcontentratio, 
                            Primary_Location_Simplified, sample_purity, differential_group_sampletype_strict, c_circle_result)) %>% 
  filter(!(symbol == "TERT" & TERT_ATRX_Mutation_Event =="Missense")) %>% 
  ggplot(aes(x=TERT_ATRX_Mutation, 
             y=log2_cpm_z, 
             color=differential_group_sampletype_strict,
             shape=TERT_ATRX_Mutation_Event,
             group=`Patient ID`)) + 
  geom_point(position = pj) + 
  geom_path(position = pj,  alpha=0.1, size=0.5, color="black", linetype=2) +
  scale_color_manual(name="Sample Type", values = sampletype_strict_cols) +
  scale_shape_manual(name="Alteration Type", values = c(
    `Promotor Mutation`=4,
    `Structural Variant`=8, 
    `None`=15,
    `Stop/Frameshift`=17,
    `Splice`=19)) +
  theme_bw() +
  #geom_text(mapping=aes(label=A5_ID)) +
  ylab(paste(g," Expression (z-scaled)")) +
  xlab("TERT/ATRX Alteration")

E143_pileup <- pileup_summaries_annotated$E143 %>% 
  filter(variant_type == "SBS") %>% 
  pivot_longer(cols = c(-CHROM, -POS, -REF, -ALT, -variant_type), 
               names_to = c("A5_ID", "value_type"), 
               names_sep = "_",
               values_to="count") %>% 
  filter(value_type %in% c("altCount", "TotalCount")) %>% 
  distinct() %>% 
  pivot_wider(names_from = value_type, values_from = count) %>% 
  mutate(pileup_VAF=altCount/TotalCount) %>% 
  pivot_wider(id_cols = c(CHROM, POS, REF, ALT), names_from = A5_ID, values_from = pileup_VAF) %>% 
  filter(`E143-B01` == 0) %>% 
  dplyr::select(-`E143-B01`) %>% 
  pivot_longer(cols = c(`E143-T01`, `E143-T02`, `E143-T03`), names_to = "A5_ID",values_to = "pileup_VAF") %>% 
  mutate(A5_ID=gsub("T0", "", A5_ID))  

E143_pileup <- E143_pileup %>% inner_join(
  E143_pileup %>%  
  pivot_wider(names_from = A5_ID, values_from = pileup_VAF) %>%  
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
  mutate(end=ifelse(is.na(end), start, end)) %>% 
  inner_join(a5_anno %>% dplyr::select(A5_ID, PublicationID, sample_purity)) %>% 
  mutate(Tumour_AF=as.numeric(Tumour_AF),
         PublicationID=factor(PublicationID, levels=c("E143-P1","E143-M1","E143-M2")),
         color=ifelse(!is.na(PCGR_SYMBOL) & PCGR_SYMBOL=="TERT", "red", "black"),
         linetype=ifelse(!is.na(PCGR_SYMBOL) & PCGR_SYMBOL=="TERT", 1, 2),
         linewidth=ifelse(!is.na(PCGR_SYMBOL) & PCGR_SYMBOL=="TERT", 2, 1),
         alpha=ifelse(!is.na(PCGR_SYMBOL) & PCGR_SYMBOL=="TERT", 1, 0.1),
         varid=paste(seqnames, start, ALT)) %>% 
  arrange(PublicationID, Tumour_AF) %>% ungroup() %>% filter(pileup_VAF > 0)


source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/cna_event_labeller.r")
  
E143_segs <- list()
E143_segs[["E143-1"]] <- read.delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E143/umccrised/E143-1__E143-T01/purple/E143-1__E143-T01.purple.segment.tsv")
E143_segs[["E143-2"]] <- read.delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E143/umccrised/E143-2__E143-T02/purple/E143-2__E143-T02.purple.segment.tsv")
E143_segs[["E143-3"]] <- read.delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E143/umccrised/E143-3__E143-T03/purple/E143-3__E143-T03.purple.segment.tsv")

E143_segs <- map( E143_segs, .f =  \(x) { 
  x %>% filter(depthWindowCount > 20) %>% 
    mutate(cna_type= classify_cna_event(minorAlleleCopyNumber,majorAlleleCopyNumber, 
                                      "Male", chromosome, tumorCopyNumber, mean(tumorCopyNumber)))})

E143_segs_gr <- map(E143_segs, .f = \(x) GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = T))

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
  mcols(temp)[["cna_type"]][subjectHits(hits)] <- mcols(E143_segs_gr[[s]])[["cna_type"]][queryHits(hits)]
  
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

cna_palette <- c(`None`=ColorPalette[["LightGrey1"]],
                 Loss=ColorPalette[["DarkRed1"]],
                 `Subclonal Loss`="#fdadadff",
                 `Hom. Del.`=ColorPalette[["Yellow1"]], CNLOH="#f97344ff", 
                 Gain=ColorPalette[["DarkBlue3"]],
                 `Subclonal Gain`=ColorPalette[["LightBlue1"]],
                 WGD=ColorPalette[["Purple2"]],
                 Chromothripsis=ColorPalette[["LightGreen1"]],
                 Other=ColorPalette[["DarkGrey2"]],
                 `Gain+LOH`=ColorPalette[["LightBlue2"]],
                 `WGD+Gain`=ColorPalette[["DarkBlue1"]], 
                 `Loss + Subclonal CNLOH`=ColorPalette[["LightOrange2"]],
                 `Diploid/Haploid-X`=ColorPalette[["LightGrey1"]],
                 `Minor Subclonal Loss`= ColorPalette[["Salmon"]])

pj <- position_jitter(width = 0.2, seed = 42)
gg_vaf <- ggplot(data = E143_mut_seg_annotated, 
                 mapping=aes(x=PublicationID, 
                             y=pileup_VAF, 
                             group=varid, 
                             linetype=linetype,
                             linewidth=linewidth)) + 
  geom_point(mapping = aes(color=Privacy),position = pj, size=0.5) +
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
  # scale_color_identity() +
  # scale_color_manual(values=c(cna_palette, "TERT"="red")) +
  coord_cartesian(ylim = c(0,1)) + ylab("VAF") +
  guides(linewidth="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))


gg_dens <- ggplot(E143_mut_seg_annotated, aes(y=pileup_VAF, x=after_stat(count), color=`PublicationID`)) + 
  geom_density() + theme_bw() + coord_cartesian(ylim = c(0,1))


gg_vaf + (gg_dens + 
            theme(axis.text.y = element_blank(), 
                  axis.title.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.margin = margin(6,6,6,0))) + 
  plot_layout(guides="collect", widths = c(3,1))




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

library(patchwork)  
(gg_tcr + theme(axis.text.x = element_blank(), 
                axis.title.x = element_blank(), 
                plot.margin = margin(6,6,0,6))) / (
                  gg_vaf + theme(plot.margin = margin(0,6,6,6))) + 
  plot_layout(heights = c(1,2))

