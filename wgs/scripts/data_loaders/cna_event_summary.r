
cn_event_types <- c("Diploid/Haploid-X", "Loss",  "Subclonal Loss", "Minor Subclonal Loss",  
                    "Loss + Subclonal CNLOH", "CNLOH", "Gain", "Subclonal Gain", "Gain+LOH", "WGD", "WGD+Gain",
                    "Hom. Del.", "Other", "Chromothripsis")

A5_seg.files <- list.files("C:/ResearchData/RADIO/A5/Data/WGS/hg38/cn_segments", full.names = T)
names(A5_seg.files) <- gsub(".purple.segment.tsv", "", basename(A5_seg.files))
A5_seg <- lapply(A5_seg.files, read.delim)

A5_seg <- bind_rows(A5_seg, .id="A5_ID")
A5_seg <- A5_seg %>%  dplyr::filter(chromosome!="chrY")

A5_seg.keep <- A5_seg %>% group_by(A5_ID) %>% mutate(mean_tumorCopyNumber=mean(tumorCopyNumber)) %>%  ungroup()
A5_seg.keep <- A5_seg.keep %>%  filter(end-start > 5000, chromosome %in% c("chr1","chr3","chr6","chr7", "chr11"))
A5_seg.keep <- A5_seg.keep %>% mutate(A5_ID=gsub("T0","",A5_ID)) %>% left_join(a5_anno %>%  dplyr::select(A5_ID, Gender))
A5_seg.keep <- A5_seg.keep %>% filter(A5_ID %in% SampleOrder.A5_ID) %>% 
  mutate(A5_ID=factor(as.character(A5_ID), levels=SampleOrder.A5_ID)) #SampleOrder.A5_ID defined in master figure script


chr_offsets <- A5_seg.keep %>% 
  mutate(chromosome=factor(as.character(chromosome), levels=paste0("chr",c(1:22,"X")))) %>% 
  arrange(chromosome) %>% group_by(A5_ID, chromosome) %>% 
  summarise(max_end=max(end)) %>% 
  ungroup() %>% group_by(chromosome) %>% 
  summarise(chr_size=median(max_end)) %>% mutate(offset = cumsum(as.numeric(chr_size))-chr_size)


A5_seg.keep <- A5_seg.keep %>% left_join(chr_offsets) %>% mutate(start_offset=start+offset, end_offset=end+offset)


A5_seg.keep <- A5_seg.keep %>% mutate(Class=case_when(
  (majorAlleleCopyNumber > 1.85 & majorAlleleCopyNumber < 2.15) & (minorAlleleCopyNumber > 1.85 & minorAlleleCopyNumber < 2.15) & (mean_tumorCopyNumber > 3) ~ "WGD",
  (majorAlleleCopyNumber > 2.85) & (minorAlleleCopyNumber > 1.85 & minorAlleleCopyNumber < 2.15) & (mean_tumorCopyNumber > 3) ~ "WGD+Gain",
  (majorAlleleCopyNumber > 0.5 & majorAlleleCopyNumber < 1.05) & (minorAlleleCopyNumber > 0.95 & minorAlleleCopyNumber < 1.15) ~ "Diploid/Haploid-X",
  (majorAlleleCopyNumber > 0.5 & majorAlleleCopyNumber < 1.15) & (minorAlleleCopyNumber < 0.25) & (chromosome!="chrX" | (chromosome=="chrX" & Gender=="female")) ~ "Loss",
  (majorAlleleCopyNumber > 1.15 & majorAlleleCopyNumber < 1.85) & (minorAlleleCopyNumber < 0.25) & (chromosome!="chrX" | (chromosome=="chrX" & Gender=="female")) ~ "Loss + Subclonal CNLOH",
  (majorAlleleCopyNumber > 0.5 & majorAlleleCopyNumber < 1.15) & (minorAlleleCopyNumber < 0.25) & (chromosome=="chrX" & Gender=="male") ~ "Diploid/Haploid-X",
  (majorAlleleCopyNumber > 0.5 & majorAlleleCopyNumber < 1.15) & (minorAlleleCopyNumber < 0.85) & (chromosome!="chrX" | (chromosome=="chrX" & Gender=="female")) ~ "Subclonal Loss",
  (majorAlleleCopyNumber > 0.5 & majorAlleleCopyNumber < 1.15) & (minorAlleleCopyNumber < 0.95) & (chromosome!="chrX" | (chromosome=="chrX" & Gender=="female")) ~ "Minor Subclonal Loss",
  (majorAlleleCopyNumber > 1.85 & majorAlleleCopyNumber < 2.15) & (minorAlleleCopyNumber < 0.25) ~ "CNLOH",
  (majorAlleleCopyNumber > 1.85) & (minorAlleleCopyNumber > 0.85) ~ "Gain", #& minorAlleleCopyNumber < 1.15
  (majorAlleleCopyNumber > 1.05) & (minorAlleleCopyNumber > 0.85 ) ~ "Subclonal Gain", #& minorAlleleCopyNumber < 1.15
  (majorAlleleCopyNumber > 1.85) & (minorAlleleCopyNumber < 0.85) ~ "Gain+LOH",
  tumorCopyNumber < 0.5 ~ "Hom. Del.",
  TRUE ~ "Other"
)) %>% 
  mutate(Class=ifelse((Class=="Minor Subclonal Loss" & 
                         ((!is.na(lead(Class)) & lead(Class)=="Subclonal Loss")) | 
                         (!is.na(lag(Class))  & lag(Class)=="Subclonal Loss")), 
                      "Subclonal Loss", 
                      Class),
         Class=ifelse((Class=="Minor Subclonal Loss" & 
                         ((!is.na(lead(Class)) & lead(Class)=="Diploid/Haploid-X")) | 
                         (!is.na(lag(Class))  & lag(Class)=="Diploid/Haploid-X")), 
                      "Diploid/Haploid-X", 
                      Class)) %>% 
  mutate(Class=factor(Class,levels=cn_event_types))

A5_seg.keep <- A5_seg.keep %>% mutate(Class=case_when(
  ((A5_ID=="E123-1" | A5_ID=="E123-T01") & chromosome=="chr5") ~ "Chromothripsis",
  ((A5_ID=="E128-2" | A5_ID=="E128-T02") & (chromosome=="chr1" & start > 124*10^6)) ~ "Chromothripsis",
  ((A5_ID=="E138-1" | A5_ID=="E138-T01") & chromosome=="chr2") ~ "Chromothripsis",
  ((A5_ID=="E180-1" | A5_ID=="E180-T01") & (chromosome=="chr3" | chromosome=="chr11")) ~ "Chromothripsis",
  ((A5_ID=="E231-1" | A5_ID=="E231-T01") & (chromosome=="chr3" & start > 50*10^6)) ~ "Chromothripsis",
  ((A5_ID=="E198-1" | A5_ID=="E198-T01") & (chromosome=="chr17" & start < 10*10^6)) ~ "Chromothripsis",
  ((A5_ID=="E198-1" | A5_ID=="E198-T01") & (chromosome=="chr22" & ((start > 17*10^6) & (end < 50*10^6)))) ~ "Chromothripsis",
  TRUE ~ as.character(Class)
))

min_gap <-1000000
A5_seg.keep.merged <- A5_seg.keep %>% arrange(A5_ID, chromosome,start,end) %>% group_by(A5_ID, chromosome) %>% 
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
  group_by(A5_ID,chromosome,Gender,Class, new_seg_start, new_seg_end, offset) %>% 
  summarise(across(.cols = c(observedBAF,minorAlleleCopyNumber,minorAlleleCopyNumberDeviation,
                             observedTumorRatio, observedNormalRatio, majorAlleleCopyNumber,
                             tumorCopyNumber,fittedTumorCopyNumber,fittedBAF,tumorBAF),.fns = mean)) %>% 
  dplyr::rename(start=new_seg_start, end=new_seg_end) %>%  
  mutate(start_offset=start+offset, end_offset=end+offset) %>% 
  mutate(Class=factor(Class,levels=cn_event_types)) %>% 
  arrange(A5_ID, chromosome,start,end)

ggCNVSeg <- ggplot(A5_seg.keep, aes(x=A5_ID, xend=A5_ID, y=start_offset, yend=end_offset, color=Class)) + geom_segment(size=1.5) +
  geom_hline(data = chr_offsets, mapping=aes(yintercept=offset)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) +
  scale_y_reverse(breaks = chr_offsets$offset, labels = chr_offsets$chromosome, expand=c(0,0)) +
  scale_color_manual(values=c(`Diploid/Haploid-X`=ColorPalette[["LightGrey1"]], Loss=ColorPalette[["LightRed1"]], `Loss + Subclonal CNLOH`=ColorPalette[["LightRed1"]],
                              CNLOH=ColorPalette[["LightOrange1"]], Gain=ColorPalette[["DarkBlue3"]], `Subclonal Gain`=ColorPalette[["LightBlue1"]],
                              `Gain+LOH`=ColorPalette[["LightBlue2"]], WGD=ColorPalette[["LightBrown1"]], 
                              `WGD+Gain`=ColorPalette[["LightBrown2"]], 
                              `Subclonal Loss`="#fdadadff", "Minor Subclonal Loss"="#fdadadff", 
                              `Hom. Del.`=ColorPalette[["Yellow1"]], 
                              Chromothripsis=ColorPalette[["LightGreen1"]], Other=ColorPalette[["Purple3"]] )) + 
  ylab("Copy Number")




ggCNVSeg.merged <- ggplot(A5_seg.keep.merged %>% 
                            mutate(Class=forcats:::fct_recode(.f=Class,  
                                                              Loss="Loss + Subclonal CNLOH",
                                                              Gain="WGD+Gain",
                                                              Gain="Gain+LOH",
                                                              None="Minor Subclonal Loss", 
                                                              None="Diploid/Haploid-X")), 
                          aes(x=A5_ID, xend=A5_ID, y=start_offset, yend=end_offset, color=Class)) + 
  geom_segment(size=2.5) +
  geom_hline(data = chr_offsets, mapping=aes(yintercept=offset), size=0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) +
  scale_y_reverse(breaks = chr_offsets$offset, labels = chr_offsets$chromosome, expand=c(0,0)) +
  scale_color_manual(values=c(`None`=ColorPalette[["LightGrey1"]], Loss=ColorPalette[["DarkRed1"]], 
                              `Subclonal Loss`="#fdadadff", `Hom. Del.`=ColorPalette[["Yellow1"]], CNLOH="#f97344ff", 
                              Gain=ColorPalette[["DarkBlue3"]], `Subclonal Gain`=ColorPalette[["LightBlue1"]],
                              WGD=ColorPalette[["Purple2"]], Chromothripsis=ColorPalette[["LightGreen1"]], Other=ColorPalette[["DarkGrey2"]]
  )) + 
  ylab("Copy Number")
