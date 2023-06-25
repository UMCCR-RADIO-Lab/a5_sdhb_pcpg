library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tidyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(gridExtra)
mvals_save <- read.csv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/A5 methylation Avals.csv")
bVals2 <- read.csv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/A5 methylation Bvals.csv")

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

mvals_save.anno <- mvals_save %>% left_join(data.frame(annEPIC) %>% tibble::rownames_to_column("Probe") %>% select(Probe,chr,pos))
bvals_save.anno <- bVals2 %>% left_join(data.frame(annEPIC) %>% tibble::rownames_to_column("Probe") %>% select(Probe,chr,pos))
mvallist <- list()
bvallist <- list()
for (i in 2:96) {  
  s <- colnames(mvals_save.anno)[i]
  s <- gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0([0-9])_D$","\\1-\\2",s)
  s <- gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0[0-9]_D_([0-2])$","\\1-\\2",s)
  mvallist[[s]] <- mvals_save.anno[,c(1,97,98,i)]
  colnames(mvallist[[s]])[4] <- s
  
  s <- colnames(bvals_save.anno)[i]
  s <- gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0([0-9])_D$","\\1-\\2",s)
  s <- gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0[0-9]_D_([0-2])$","\\1-\\2",s)
  bvallist[[s]] <- data.frame(bvals_save.anno[,c(1,97,98,i)])
  colnames(bvallist[[s]])[4] <- s
}

purplecn.files <- list.files(path = "/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/WGS/Tool outputs/All UMCCRise outputs", 
                             recursive = T, pattern = ".purple.cnv",full.names = T)
purplecn.files <- purplecn.files[!grepl("work",purplecn.files)]
purplecn <- lapply(purplecn.files, read.delim)
snames <- basename(purplecn.files)
snames <- gsub("^(E[0-9]{3})__.+","\\1-1",snames)
snames <- gsub("^(E[0-9]{3})_([0-9]).+","\\1-\\2",snames)
snames <- gsub("^FL-2.+","E233-1",snames)
names(purplecn) <- snames

mbvals <- list()
for (s in names(mvallist))
{
  MGR <- GRanges(seqnames=mvallist[[s]]$chr,
                 ranges = IRanges(start=mvallist[[s]]$pos,
                                  end=mvallist[[s]]$pos),
                 Probe=mvallist[[s]]$Probe, 
                 Mval=mvallist[[s]][,s])
  
  BGR <- GRanges(seqnames=bvallist[[s]]$chr,
                 ranges = IRanges(start=bvallist[[s]]$pos,
                                  end=bvallist[[s]]$pos),
                 Probe=bvallist[[s]]$Probe, 
                 Bval=bvallist[[s]][,s])
  
  CNGR <- GRanges(seqnames=paste0("chr", purplecn[[s]]$X.chromosome),
                  ranges = IRanges(start=purplecn[[s]]$start,
                                   end=purplecn[[s]]$end),
                  CN=purplecn[[s]]$copyNumber)
  OVL <- findOverlaps(MGR,CNGR)
  mcols(MGR)[["CN"]] <- 0
  mcols(MGR)$CN[queryHits(OVL)] <- mcols(CNGR)$CN[subjectHits(OVL)]
  MGR <- data.frame(MGR) %>% rename(chr=seqnames) %>%  select(-width,-strand)
  MGR$Sample <- s
  mvallist[[s]] <- MGR
  
  OVL <- findOverlaps(BGR,CNGR)
  mcols(BGR)[["CN"]] <- 0
  mcols(BGR)$CN[queryHits(OVL)] <- mcols(CNGR)$CN[subjectHits(OVL)]
  BGR <- data.frame(BGR) %>% rename(chr=seqnames) %>%  select(-width,-strand)
  BGR$Sample <- s
  bvallist[[s]] <- BGR
  mbvals[[s]] <- MGR %>% full_join(BGR)
}
mbvals <- bind_rows(mbvals)

CheckSamples <- data.frame(Sample=c("E120-1","E121-1","E122-1","E122-1","E123-1","E126-1","E127-1","E135-1","E136-1","E143-1","E145-1","E145-1","E160-1","E160-1","E160-1","E171-1","E171-1","E188-1","E188-1"),
                           ROI=c("chr11p","chr11p","chr3q","chr11p","chr3q","chr3q","chr11p","chr11p","chr11p","chr14","chr3q","chr11p","chr3q","chr11p","chr14","chr11p","chr14","chr3q","chr11p"))
CheckSamples$Loss="Loss"
ggplot(mbvals %>% filter(Sample %in% CheckSamples$Sample) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
         filter(!is.na(ROI), 
                ,CN<5, CN>0.5)  %>%   
    inner_join(CheckSamples) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
         slice_sample(prop = 0.01) %>% 
         mutate(CN=factor(round(CN,0))), aes(x=Bval, color=Sample)) + 
  geom_density() + facet_wrap("chr") 
geom_jitter(size=0.05, alpha=0.1)


  
ggplot(mbvals %>% filter(Sample %in% CheckSamples$Sample, chr %in% c("chr3","chr11","chr14")) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
    filter(!is.na(ROI), 
           ,CN<5, CN>0.5)  %>%   
    left_join(CheckSamples, ) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
    #slice_sample(prop = 0.01) %>% 
    mutate(CN=factor(round(CN,0))), aes(x=Bval, color=Sample, linetype=Loss)) + 
  geom_density() + facet_grid(Loss~chr) 
geom_jitter(size=0.05, alpha=0.1)

ggplot(mbvals %>% filter(Sample %in% CheckSamples$Sample, chr %in% c("chr3","chr11","chr14")) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
    filter(!is.na(ROI), 
           ,CN<5, CN>0.5)  %>%   
    left_join(CheckSamples) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
    slice_sample(prop = 0.2) %>% 
    mutate(CN=factor(round(CN,0))), aes(x=Sample, y=Bval, color=Loss)) + 
  geom_jitter(size=0.05, alpha=0.1) + facet_wrap("chr")


mbvals %>% filter(Sample %in% CheckSamples$Sample, chr %in% c("chr3","chr11","chr14")) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
  filter(!is.na(ROI), 
         ,CN<5, CN>0.5)  %>%   
  left_join(CheckSamples) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
  pivot_wider(id_cols = c(Probe,chr,start), names_from=c(Sample), values_from=c(Mval,Loss)) %>% 
  slice_sample(prop = 0.2) 


plot.data <-mbvals %>% filter(Sample %in% CheckSamples$Sample, chr %in% c("chr3","chr11","chr14")) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
  filter(!is.na(ROI), 
         ,CN<5, CN>0.5)  %>%   
  left_join(CheckSamples) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
  mutate(ROI=ifelse((as.character(chr)=="chr14" & start >	100819176 &start < 101805942), "DLK-MEG3",ROI)) %>% 
  pivot_wider(id_cols = c(Probe,chr,start, ROI), names_from=c(Sample), values_from=c(Mval,Loss)) %>% 
  slice_sample(prop = 0.2) 
plot.list <- list()

for (s1 in c("E143-1", "E160-1", "E171-1"))
{
  for (s2 in c("E188-1", "E122-1", "E135-1"))
{
plot.list[[paste(s1,s2,sep="_")]] <- ggplot(plot.data,
  aes(x=!!sym(paste0("Mval_",s1)), y=!!sym(paste0("Mval_",s2)), color=gsub("No Loss/Loss","Loss/No Loss",paste(!!sym(paste0("Loss_",s1)),!!sym(paste0("Loss_",s2)),sep="/")))) + 
  geom_point(size=0.1, alpha=0.3) +
  facet_wrap("ROI", ncol=4) + labs(color="Combo") + guides(color = guide_legend(override.aes = list(size=5)))
  }
}

do.call("grid.arrange", plot.list)
