

A5_seg.files <- list.files("C:/ResearchData/RADIO/A5/Data/WGS/hg38/cn_segments", full.names = T)
names(A5_seg.files) <- gsub(".purple.segment.tsv", "", basename(A5_seg.files))
A5_seg <- lapply(A5_seg.files, read.delim)

A5_seg <- bind_rows(A5_seg, .id="A5_ID")
A5_seg <- A5_seg %>%  dplyr::filter(chromosome!="chrY")

A5_seg.keep <- A5_seg %>%  filter(end-start > 5000, chromosome %in% c("chr1","chr3","chr6","chr7", "chr11"))
A5_seg.keep <- A5_seg.keep %>% mutate(A5_ID=gsub("T0","",A5_ID)) %>% left_join(A5_anno %>%  select(A5_ID, Gender))
A5_seg.keep <- A5_seg.keep %>% mutate(A5_ID=factor(as.character(A5_ID), levels=SampleOrder.A5_ID)) #SampleOrder.A5_ID defined in master figure script

chr_offsets <- A5_seg.keep %>% 
  mutate(chromosome=factor(as.character(chromosome), levels=paste0("chr",c(1:22,"X")))) %>% 
  arrange(chromosome) %>% group_by(A5_ID, chromosome) %>% 
  summarise(max_end=max(end)) %>% 
  ungroup() %>% group_by(chromosome) %>% 
  summarise(chr_size=median(max_end)) %>% mutate(offset = cumsum(as.numeric(chr_size))-chr_size)


A5_seg.keep <- A5_seg.keep %>% left_join(chr_offsets) %>% mutate(start_offset=start+offset, end_offset=end+offset)


A5_seg.keep <- A5_seg.keep %>% mutate(Class=case_when(
  (majorAlleleCopyNumber > 0.6 & majorAlleleCopyNumber < 1.4) & (minorAlleleCopyNumber > 0.6 & minorAlleleCopyNumber < 1.3) ~ "Diploid/Haploid-X",
  (majorAlleleCopyNumber > 0.6 & majorAlleleCopyNumber < 1.4) & (minorAlleleCopyNumber < 0.3) & (chromosome!="chrX" | (chromosome=="chrX" & Gender=="female")) ~ "Loss",
  (majorAlleleCopyNumber > 0.6 & majorAlleleCopyNumber < 1.4) & (minorAlleleCopyNumber < 0.3) & (chromosome=="chrX" & Gender=="male") ~ "Diploid/Haploid-X",
  (majorAlleleCopyNumber > 1.6 & majorAlleleCopyNumber < 2.4) & (minorAlleleCopyNumber < 0.3) ~ "CNLOH",
  (majorAlleleCopyNumber > 1.6 & majorAlleleCopyNumber < 2.4) & (minorAlleleCopyNumber > 1.7 & minorAlleleCopyNumber < 2.3) ~ "WGD",
  (majorAlleleCopyNumber > 1.6) & (minorAlleleCopyNumber > 0.6 & minorAlleleCopyNumber < 1.4) ~ "Gain",
  (majorAlleleCopyNumber > 1.6) & (minorAlleleCopyNumber < 0.6) ~ "Gain+LOH",
  (majorAlleleCopyNumber > 3.6) & (minorAlleleCopyNumber > 1.6 & minorAlleleCopyNumber < 2.4) ~ "WGD+Gain",
  tumorCopyNumber < 0.5 ~ "Hom. Del.",
  TRUE ~ "Other"
)) %>% mutate(Class=factor(Class,levels=c("Diploid/Haploid-X", "Loss", "CNLOH", "Gain", "Gain+LOH", "WGD", "WGD+Gain","Hom. Del.", "Other", "Chromothripsis")))

A5_seg.keep <- A5_seg.keep %>% mutate(Class=case_when(
  (A5_ID=="E123-1" & chromosome=="chr5") ~ "Chromothripsis",
  (A5_ID=="E128-2" & (chromosome=="chr1" & start > 124*10^6)) ~ "Chromothripsis",
  (A5_ID=="E138-1" & chromosome=="chr2") ~ "Chromothripsis",
  (A5_ID=="E180-1" & (chromosome=="chr3" | chromosome=="chr11")) ~ "Chromothripsis",
  (A5_ID=="E231-1" & (chromosome=="chr3" & start > 50*10^6)) ~ "Chromothripsis",
  (A5_ID=="E198-1" & (chromosome=="chr17" & start < 10*10^6)) ~ "Chromothripsis",
  (A5_ID=="E198-1" & (chromosome=="chr22" & ((start > 17*10^6) & (end < 50*10^6)))) ~ "Chromothripsis",
  TRUE ~ as.character(Class)
))

ggCNVSeg <- ggplot(A5_seg.keep, aes(x=A5_ID, xend=A5_ID, y=start_offset, yend=end_offset, color=Class)) + geom_segment(size=1.5) +
  geom_hline(data = chr_offsets, mapping=aes(yintercept=offset)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) +
  scale_y_reverse(breaks = chr_offsets$offset, labels = chr_offsets$chromosome, expand=c(0,0)) +
  scale_color_manual(values=c(`Diploid/Haploid-X`=ColorPalette[["LightGrey1"]], Loss=ColorPalette[["LightRed1"]], 
                              CNLOH=ColorPalette[["LightOrange1"]], Gain=ColorPalette[["DarkBlue3"]], 
                              `Gain+LOH`=ColorPalette[["LightBlue2"]], WGD=ColorPalette[["LightBrown1"]], 
                              `WGD+Gain`=ColorPalette[["LightBrown2"]],`Hom. Del.`=ColorPalette[["Yellow1"]], 
                              Chromothripsis=ColorPalette[["LightGreen1"]], Other=ColorPalette[["Purple3"]] )) + 
  ylab("Copy Number")
  

