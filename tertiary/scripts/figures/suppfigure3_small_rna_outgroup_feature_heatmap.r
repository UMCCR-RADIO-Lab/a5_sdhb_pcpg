###########################################
###########################################
##                                       ##
##   __A5 small-RNA outliers heatmap__   ##
##                                       ##
##   Script for producing a heatmap      ##
##   highlighting the imprinted genes    ##
##   on chr14 in the small-RNA outlier   ##
##   group                               ##
##                                       ##
###########################################
###########################################


setwd("/g/data/pq08/projects/ppgl")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(googlesheets4)
library(patchwork)
library(parallel)
library(GenomicRanges)

source("/g/data/pq08/projects/ppgl/a5/scripts/geom_split_tile.r")

############################################
# Import data loaders and run data mergers #
############################################

#A5 CNA
source("./a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_cna_segs()

#Methylation
source("./a5/tertiary/scripts/data_mergers/combine_tcga_comete_a5_methylation_data.r")

#WTS
source("./a5/tertiary/scripts/data_mergers/combine_tcga_flynn_a5_wts_data.r")

#small-rna
source("./a5/tertiary/scripts/data_mergers/combine_tcga_comete_a5_smallrna_data.r")

#tcga cnv segments
source("./public_data/data_loaders/tcga_cnv_dataloader.r")
dataloader_tcga_ppgl_ascat()

#clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

##################
# Colour Palette #
##################

ColorPalette <- c(LightBlue1="#cadae8ff",LightBlue2="#abc4e2ff",LightBlue3="#8aa9d6ff",
                  DarkBlue1="#7884baff",DarkBlue2="#606d9fff",DarkBlue3="#4c5a88ff",
                  Purple1="#765b92ff",Purple2="#9370a5ff",Purple3="#b280a9ff",Purple4="#cc85b1ff",
                  LightRed1="#ee7474ff",LightRed2="#e2696aff",LightRed3="#de5353ff",
                  DarkRed1="#c25858ff",DarkRed2="#c04241ff",DarkOrange1="#c65a44ff",
                  DarkOrange2="#eb7b54ff",LightOrange1="#f18d73ff",LightOrange2="#f5a697ff",
                  LightBrown1="#f5b9b0ff",LightBrown2="#d6a097ff",DarkBrown1="#c17963ff",DarkBrown2="#96665aff",
                  DarkGreen1="#637b62ff",DarkGreen2="#6ea668ff",LightGreen1="#98ca8cff",LightGreen2="#e2eab5ff",
                  Yellow1="#fbf2adff",Yellow2="#fef8c6ff",Yellow3="#f4e764ff",Salmon="#f2c5a7ff",LightSalmon="#f2dec4ff",
                  LemonGrey1="#f7f0e2ff",LemonWhite1="#fefbe5ff",LemonWhite2="#f7f4dcff",LemonGrey2="#efecdcff",
                  LightGrey1="#e2e0d9ff",LightGrey2="#d6d6d4ff",DarkGrey1="#b4b4b4ff",DarkGrey2="#9a9999ff")

cna_palette <- c(`None`=ColorPalette[["LightGrey1"]],
                 Loss=ColorPalette[["DarkRed1"]],
                 `Subclonal Loss`="#fdadadff",
                 `Hom. Del.`=ColorPalette[["Yellow1"]], CNLOH="#f97344ff", 
                 Gain=ColorPalette[["DarkBlue3"]],
                 `Subclonal Gain`=ColorPalette[["LightBlue1"]],
                 WGD=ColorPalette[["Purple2"]],
                 Chromothripsis=ColorPalette[["LightGreen1"]],
                 Other=ColorPalette[["DarkGrey2"]],
                 `Gain+LOH`=ColorPalette[["LightBlue3"]],
                 `WGD+Gain`=ColorPalette[["DarkBlue1"]], 
                 `Loss + Subclonal CNLOH`=ColorPalette[["LightOrange2"]],
                 `Diploid/Haploid-X`=ColorPalette[["LightGrey1"]],
                 `Diploid`=ColorPalette[["LightGrey1"]],
                 `Minor Subclonal Loss`= ColorPalette[["Salmon"]])

###############
# Public Data #
###############

#Castro-Vega 2014, Supp Table 6, mir DE 
if(!file.exists("./public_data/annotation/castrovega_2014_supp6.xlsx")) {
  download.file(url = "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms7044/MediaObjects/41467_2015_BFncomms7044_MOESM371_ESM.xlsx",
                destfile ="./public_data/annotation/castrovega_2014_supp6.xlsx")
}

#Cytobands - HG38
if(!file.exists("/g/data/pq08/reference/GRCh38/cytobands_hg38.txt.gz")) {
  download.file(url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
                destfile ="/g/data/pq08/reference/GRCh38/cytobands_hg38.txt.gz")
}

cytoband <- read.delim("/g/data/pq08/reference/GRCh38/cytobands_hg38.txt.gz", header=F)
colnames(cytoband) <- c("chr","start","end","cytoband","stain")
cytoband <- cytoband %>% filter(chr %in% paste0("chr", c(1:22,"X","Y")))

cytoband_gr <- GenomicRanges::makeGRangesFromDataFrame(cytoband, 
                                                       seqnames.field ="chr", 
                                                       ignore.strand = T, 
                                                       keep.extra.columns = T)

####################
# Array Annotation #
####################

epic_array_annotation_hg38 <- as.data.frame(epic_array_annotation_hg38)
Ill450k_array_annotation_hg38 <- epic_array_annotation_hg38 %>% filter(Name %in% meth_tcga_comete_a5_450k_mval$Probe)

###################
# Reference data  #
###################

########
# miR  
########

castrovega_mi3vsmi4_7 <- readxl::read_xlsx("./public_data/annotation/castrovega_2014_supp6.xlsx", 
                                           sheet = 2, 
                                           range = "A3:E100",
                                           col_types = rep(x = c("text", "numeric"), times=c(1,4)))
castrovega_mi3vsmi4_7_use <- castrovega_mi3vsmi4_7 %>% 
  mutate(log_FC=log2(`Fold-change Mi4-7/Mi3`)) %>%
  filter(abs(log_FC)> 3)

mirbase_hsa_hg38 <- read.delim("/g/data/pq08/projects/ppgl/a5/small_rna/analysis/reference/Mirbase_hsa_hg38.gff3", 
                               comment.char = '#',
                               header=F)

colnames(mirbase_hsa_hg38) <- c("seqname", "source", "feature", "start", 
                                "end", "score", "strand", "frame", "attribute")

mirbase_hsa_hg38_mirs <- mirbase_hsa_hg38 %>%
  filter(feature != "miRNA_primary_transcript") %>%
  mutate(
    mir_name = gsub(".+Name=([A-Za-z0-9-]+);?.*", "\\1", attribute),
    mir_name = gsub("mir", "miR", mir_name)
  ) %>%
  dplyr::select(seqname, start, end, mir_name) %>%
  mutate(base_name = gsub("(hsa-(miR|let)-[0-9a-z]+)(-.+)?", "\\1", mir_name))

mirbase_hsa_hg38_loopstem <- mirbase_hsa_hg38 %>%
  filter(feature == "miRNA_primary_transcript") %>%
  mutate(
    mir_name = gsub(".+Name=([A-Za-z0-9-]+);?.*", "\\1", attribute),
    mir_name = gsub("mir", "miR", mir_name)
  ) %>%
  dplyr::select(seqname, start, end, mir_name) %>%
  mutate(base_name = gsub("(hsa-(miR|let)-[0-9a-z]+)(-.+)?", 
                          "\\1", 
                          mir_name)) %>%
  group_by(base_name, seqname) %>% 
  summarise(start = min(start), 
            end = max(end)) %>%
  dplyr::rename(loop_stem_start = start, 
                loop_stem_end = end)


mirs <- mirbase_hsa_hg38_mirs %>% 
  left_join(mirbase_hsa_hg38_loopstem, by=c("base_name","seqname")) %>% 
  inner_join(castrovega_mi3vsmi4_7_use %>% dplyr::select(id,log_FC), by=c("mir_name"="id"))

########
# Imprinted Genes 
########

GOI <- data.frame(Gene= c("DLK1", "DLK1", "MEG3", "MEG3", "BEGAIN", "RTL1", "RTL1", "DIO3", "DIO3"),
                  feature=c("promoter", "gene_body","promoter", "gene_body", "gene_body", "gene_body","promoter", "gene_body","promoter"),
                  chr=c("chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14","chr14"),
                  start=c(100726611,100726892, 100825867, 100826098,100537147, 100879753,100880668, 101561351,101561165),
                  end=c(100727237, 100738224, 100826215, 100861026, 100587417, 100903722, 100881000, 101563452, 101561510))

###########################
# small RNA UMAP outliers #
###########################

small_rna_outliers <- c(
  "TCGA-PR-A5PH-01", "TCGA-QR-A6GY-01", "TCGA-QR-A70H-01", "TCGA-QR-A70R-01", "TCGA-QT-A5XM-01", "TCGA-QT-A5XP-01", 
  "TCGA-QT-A69Q-01", "TCGA-RW-A67V-01", "TCGA-RW-A68G-01", "TCGA-S7-A7WR-01", "TCGA-S7-A7WT-01", "TCGA-SR-A6MS-01", 
  "TCGA-WB-A80L-01", "TCGA-WB-A80M-01", "TCGA-WB-A815-01", "TCGA-WB-A816-01", "TCGA-WB-A817-01", "TCGA-WB-A81A-01", 
  "TCGA-WB-A81K-01", "TCGA-WB-A81M-01", "E143-1", "E143-2", "E143-3", "E188-1")

###############################
# Process allele-specific CNV #
###############################

tcga_ppgl_cna_ascat_gr <- 
  GenomicRanges::makeGRangesFromDataFrame(tcga_ppgl_cna_ascat, 
                                          seqnames.field ="Chromosome", 
                                          ignore.strand = T,keep.extra.columns = T)

hits <- findOverlaps(subject = tcga_ppgl_cna_ascat_gr, query = cytoband_gr)
mcols(hits)[["cytoband"]] <- mcols(cytoband_gr)[["cytoband"]][queryHits(hits)]
mcols(hits)[["cytoband_start_pos"]] <- start(cytoband_gr)[queryHits(hits)]
mcols(hits)[["cytoband_end_pos"]] <- end(cytoband_gr)[queryHits(hits)]
mcols(hits)[["arm"]] <- stringr::str_extract(string = mcols(hits)[["cytoband"]], pattern = "p|q")
hits <- as.data.frame(hits) %>% 
  group_by(subjectHits) %>%
  arrange(cytoband_start_pos, cytoband_end_pos) %>% 
  filter(row_number()==1 | row_number()==n()) %>% 
  summarise(cytoband=paste(unique(cytoband), collapse = "-"), 
            arm=paste(unique(arm), collapse = "/"))

tcga_ppgl_cna_ascat[hits$subjectHits,"cytoband"] <- hits$cytoband
tcga_ppgl_cna_ascat[hits$subjectHits,"arm"] <- hits$arm


tcga_ppgl_cna_ascat_chr14 <- tcga_ppgl_cna_ascat %>% filter(Chromosome=="chr14")

#not true sizes just regions represented in segments
chrom_sizes <- tcga_ppgl_cna_ascat %>% 
  group_by(Chromosome) %>% 
  summarise(size=max(End)-min(Start)) %>%  
  {'names<-'(.$size, .$Chromosome)}

tcga_ppgl_cna_ascat_chr14.summarised <- tcga_ppgl_cna_ascat_chr14 %>%  
  mutate(cn_state=case_when(
    Major_Copy_Number == 1 & Minor_Copy_Number == 1 ~ "Diploid",
    Major_Copy_Number == 4 & Minor_Copy_Number == 2 ~ "WGD",
    Major_Copy_Number > 2 & Minor_Copy_Number == 0 ~ "Gain+LOH",
    Major_Copy_Number == 1 & Minor_Copy_Number == 0 ~ "Loss",
    Major_Copy_Number == 2 & Minor_Copy_Number == 0 ~ "CNLOH",
    Major_Copy_Number > 2 & Minor_Copy_Number > 0 ~ "Gain",
    Minor_Copy_Number == 0 ~ "Other-LOH",
    TRUE ~ "Other"),
    seg_size=End-Start) %>% 
  group_by(Sample, Chromosome, cn_state) %>% 
  summarise(seg_size=sum(seg_size)) %>% 
  mutate(prop=seg_size/chrom_sizes[Chromosome]) 

tcga_ppgl_cna_ascat_chr14.summarised <- 
  tcga_ppgl_cna_ascat_chr14.summarised %>% group_by(Sample) %>% 
  mutate(priority=case_when(
    cn_state=="Loss" & prop > 0.25 ~ 1,
    cn_state=="CNLOH" & prop > 0.25 ~ 1,
    cn_state=="Gain+LOH" & prop > 0.25 ~ 2,
    cn_state=="Gain" & prop > 0.25 ~ 3,
    cn_state=="Diploid" & prop > 0.25 ~ 4,
    prop > 0.5 ~ 5,
    prop > 0.25 ~ 6,
    TRUE ~ 7)) %>% 
  arrange(priority, desc(prop)) %>% 
  slice_head(n = 1) %>% 
  dplyr::select(Sample, cn_state) %>% 
  mutate(Sample=substr(x = Sample, start = 1,stop=15))

####
#A5 segment data
####

A5_seg_keep_gr <- 
  GenomicRanges::makeGRangesFromDataFrame(A5_seg_keep, 
                                          seqnames.field ="chromosome", 
                                          start.field =  "start", 
                                          end.field = "end", 
                                          ignore.strand = T,keep.extra.columns = T)

hits <- findOverlaps(subject = A5_seg_keep_gr, query = cytoband_gr)
mcols(hits)[["cytoband"]] <- mcols(cytoband_gr)[["cytoband"]][queryHits(hits)]
mcols(hits)[["cytoband_start_pos"]] <- start(cytoband_gr)[queryHits(hits)]
mcols(hits)[["cytoband_end_pos"]] <- end(cytoband_gr)[queryHits(hits)]
mcols(hits)[["arm"]] <- stringr::str_extract(string = mcols(hits)[["cytoband"]], pattern = "p|q")
hits <- as.data.frame(hits) %>% 
  group_by(subjectHits) %>%
  arrange(cytoband_start_pos, cytoband_end_pos) %>% 
  filter(row_number()==1 | row_number()==n()) %>% 
  summarise(cytoband=paste(unique(cytoband), collapse = "-"), 
            arm=paste(unique(arm), collapse = "/"))

A5_seg_keep[hits$subjectHits,"cytoband"] <- hits$cytoband
A5_seg_keep[hits$subjectHits,"arm"] <- hits$arm


A5_seg_keep_chr14 <- A5_seg_keep %>% filter(chromosome=="chr14")

#not true sizes just regions represented in segments
chrom_sizes <- A5_seg_keep %>% 
  group_by(chromosome) %>% 
  summarise(size=max(end)-min(start)) %>%  
  {'names<-'(.$size, .$chromosome)}

A5_seg_keep_chr14.summarised <- A5_seg_keep_chr14 %>%  
  mutate(cn_state=case_when(
      (majorAlleleCopyNumber >= 1.85 &
         majorAlleleCopyNumber <= 2.15) &
        (minorAlleleCopyNumber >= 1.85 &
           minorAlleleCopyNumber <= 2.15) &
        (mean_tumorCopyNumber >= 3) ~ "WGD",
      (majorAlleleCopyNumber >= 2.85) &
        (minorAlleleCopyNumber >= 1.85 &
           minorAlleleCopyNumber <= 2.15) &
        (mean_tumorCopyNumber >= 3) ~ "WGD+Gain",
      (majorAlleleCopyNumber >= 0.5 &
         majorAlleleCopyNumber <= 1.0) &
        (minorAlleleCopyNumber >= 0.97 &
           minorAlleleCopyNumber <= 1.15) ~ "Diploid/Haploid-X",
      (majorAlleleCopyNumber >= 0.5 &
         majorAlleleCopyNumber <= 1.15) &
        (minorAlleleCopyNumber <= 0.25) &
        (chromosome != "chrX" |
           (chromosome == "chrX" & Gender == "female")) ~ "Loss",
      (majorAlleleCopyNumber >= 1.15 &
         majorAlleleCopyNumber <= 1.85) &
        (minorAlleleCopyNumber <= 0.25) &
        (chromosome != "chrX" |
           (chromosome == "chrX" &
              Gender == "female")) ~ "Loss + Subclonal CNLOH",
      (majorAlleleCopyNumber >= 0.5 &
         majorAlleleCopyNumber <= 1.15) &
        (minorAlleleCopyNumber <= 0.25) &
        (chromosome == "chrX" &
           Gender == "male") ~ "Diploid/Haploid-X",
      (majorAlleleCopyNumber >= 0.5 &
         majorAlleleCopyNumber <= 1.15) &
        (minorAlleleCopyNumber <= 0.85) &
        (chromosome != "chrX" |
           (chromosome == "chrX" &
              Gender == "female")) ~ "Subclonal Loss",
      (majorAlleleCopyNumber >= 0.5 &
         majorAlleleCopyNumber <= 1.15) &
        (minorAlleleCopyNumber <= 0.97) &
        (chromosome != "chrX" |
           (chromosome == "chrX" &
              Gender == "female")) ~ "Minor Subclonal Loss",
      (majorAlleleCopyNumber >= 1.85 &
         majorAlleleCopyNumber <= 2.15) &
        (minorAlleleCopyNumber <= 0.25) ~ "CNLOH",
      (majorAlleleCopyNumber >= 1.85) &
        (minorAlleleCopyNumber >= 0.85) ~ "Gain",
      #& minorAlleleCopyNumber <= 1.15
      (majorAlleleCopyNumber >= 1.0) &
        (minorAlleleCopyNumber >= 0.85) ~ "Subclonal Gain",
      #& minorAlleleCopyNumber <= 1.15
      (majorAlleleCopyNumber >= 1.85) &
        (minorAlleleCopyNumber <= 0.85) ~ "Gain+LOH",
      tumorCopyNumber <= 0.5 ~ "Hom. Del.",
      TRUE ~ "Other"
    ),
    seg_size=end-start) %>% 
  group_by(A5_ID, chromosome, cn_state) %>% 
  summarise(seg_size=sum(seg_size)) %>% 
  mutate(prop=seg_size/chrom_sizes[chromosome]) 

A5_seg_keep_chr14.summarised <- 
  A5_seg_keep_chr14.summarised %>% 
  group_by(A5_ID) %>% 
  mutate(priority=case_when(
    cn_state=="Loss" & prop > 0.25 ~ 1,
    cn_state=="Subclonal Loss" & prop > 0.25 ~ 1,
    cn_state=="CNLOH" & prop > 0.25 ~ 1,
    cn_state=="Gain" & prop > 0.25 ~ 2,
    cn_state=="Diploid/Haploid-X" & prop > 0.25 ~ 3,
    prop > 0.5 ~ 4, 
    prop > 0.25 ~ 5, 
    TRUE ~ 6)) %>% 
  arrange(priority,desc(prop)) %>% 
  slice_head(n = 1) %>% 
  dplyr::select(A5_ID, cn_state) %>% 
  mutate(cn_state=dplyr::recode(cn_state,
                                "Diploid/Haploid-X"="Diploid",
                                "Minor Subclonal Loss"="Diploid"))

##############################
# Compute Methylation Values #
##############################

GOI$Probes <- vector(mode = "character", length=nrow(GOI))
for (i in 1:nrow(GOI))
{
  group_regex <- ifelse(GOI$feature[[i]]=="promoter", "TSS|UTR|1stExon","Body")
  
  GOI$Probes[[i]] <- paste(Ill450k_array_annotation_hg38 %>% filter(CHR_hg38==GOI[i,"chr"], 
                                                                    Start_hg38 > GOI[i,"start"],
                                                                    End_hg38 < GOI[i,"end"],
                                                                    grepl(group_regex, UCSC_RefGene_Group)) %>% pull(Name),
                           collapse=",")
}
GOI.long <- GOI %>% separate_rows(Probes, sep=',') %>% dplyr::rename(Probe=Probes)

mirs$Probes <- vector(mode = "list", length=nrow(mirs))
for (i in 1:nrow(mirs))
{
  mirs$Probes[[i]] <- paste(Ill450k_array_annotation_hg38 %>% filter(CHR_hg38==mirs[i,"seqname"], 
                                                                     Start_hg38 > mirs[i,"loop_stem_start"],
                                                                     End_hg38 < mirs[i,"loop_stem_end"]) %>% pull(Name),
                            collapse=",")
}
mirs <- mirs %>% filter(Probes != "", mir_name %in% rownames(smallrna_tcga_comete_a5_lcpm.batch_removed)) %>% slice_max(n=10, order_by = log_FC)
mirs.long <- mirs %>% separate_rows(Probes, sep=',') %>% dplyr::rename(Probe=Probes)

GOI_mir <- bind_rows(GOI.long %>% dplyr::select(Gene, feature, Probe),
                     mirs.long %>% dplyr::rename(Gene=mir_name) %>% dplyr::select(Gene, Probe) %>% mutate(feature="gene_body")) %>% 
  filter(Probe!="")

meth_tcga_comete_a5_450k_mval.batch_removed <- 
  data.frame(meth_tcga_comete_a5_450k_mval.batch_removed, check.names = F) %>% 
  tibble::rownames_to_column("Probe") 

GOI_mir_meth <- meth_tcga_comete_a5_450k_mval.batch_removed %>% 
  filter(Probe %in% GOI_mir$Probe ) %>% 
  pivot_longer(cols = -Probe, 
               names_to = "Sample", 
               values_to = "m_val") %>% 
  left_join(GOI_mir %>% dplyr::select(Probe, Gene, feature)) %>% 
  group_by(Gene,feature, Sample) %>% 
  summarise(mean_m_val=mean(m_val)) %>% 
  group_by(Gene) %>% 
  mutate(m_val_z=((mean_m_val-mean(mean_m_val))/sd(mean_m_val))) 

mir_expression <- data.frame(smallrna_tcga_comete_a5_lcpm.batch_removed, check.names = F) %>% 
  tibble::rownames_to_column("mir_name") %>% 
  filter(mir_name %in% mirs$mir_name) %>% 
  pivot_longer(cols = -mir_name, 
               names_to = "Sample", 
               values_to = "log2_cpm") %>% 
  group_by(mir_name) %>% 
  mutate(log2_cpm_z=((log2_cpm-mean(log2_cpm))/sd(log2_cpm))) 


ensgid_symbol_lookup <- data.frame(ensgid_symbol=rownames(a5_wts_counts)) %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), remove=F, extra = "merge", sep="_") %>% 
  mutate(ensgid=gsub("[.].+$", "", ensgid))

wts_tcga_comete_a5_lcpm.batch_removed <- data.frame(wts_tcga_flynn_a5_lcpm.batch_removed, check.names = F) %>% 
  tibble::rownames_to_column("ensgid") %>% inner_join(ensgid_symbol_lookup) %>% 
  dplyr::relocate(ensgid_symbol, symbol) %>% dplyr::select(-ensgid_symbol)


GOI_expression <- wts_tcga_comete_a5_lcpm.batch_removed %>% 
  filter(symbol %in% GOI$Gene) %>% 
  pivot_longer(cols = c(-symbol,-ensgid), 
               names_to = "Sample", 
               values_to = "log2_cpm") %>% 
  group_by(symbol) %>% 
  mutate(log2_cpm_z=((log2_cpm-mean(log2_cpm))/sd(log2_cpm))) 

# GOI_mir_meth %>%  mutate(feature=paste0(feature, "_methylation")) %>% 
#   dplyr::select(Gene, 
#                 Sample, 
#                 data_type=feature,   
#                 z_value=m_val_z) 

plot_data_expr_meth <- bind_rows(
  GOI_mir_meth %>% 
    dplyr::select(Gene, 
                  Sample, 
                  data_type=feature,   
                  value=m_val_z) %>%
    mutate(source="methylation"),
  GOI_expression %>% dplyr::select(Gene=symbol, 
                                   Sample, 
                                   value=log2_cpm_z) %>% 
    mutate(data_type="expression", source="mrna"),
  mir_expression %>% dplyr::select(Gene=mir_name, 
                                   Sample, 
                                   value=log2_cpm_z) %>% 
    mutate(data_type="expression", source="small_rna")
)

plot_data_expr_meth$Gene <- factor(plot_data_expr_meth$Gene, 
                                   levels = c(sort(unique(GOI$Gene)), sort(unique(mirs$mir_name))))
plot_data_expr_meth$Sample <- factor(plot_data_expr_meth$Sample)
plot_data_expr_meth$data_type <- factor(plot_data_expr_meth$data_type)
plot_data_expr_meth$source <- factor(plot_data_expr_meth$source)

plot_data_14q_cn <- bind_rows(
  tcga_ppgl_cna_ascat_chr14.summarised,
  A5_seg_keep_chr14.summarised %>% dplyr::rename(Sample=A5_ID))

complete_samples <- purrr::reduce(.x = list(colnames(smallrna_tcga_comete_a5_lcpm.batch_removed),
                        colnames(wts_tcga_comete_a5_lcpm.batch_removed),
                        plot_data_14q_cn$Sample,
                        colnames(meth_tcga_comete_a5_450k_mval.batch_removed)),.f = intersect)

plot_data_expr_meth <- plot_data_expr_meth %>% filter(Sample %in% complete_samples)
plot_data_14q_cn <- plot_data_14q_cn %>% filter(Sample %in% complete_samples)

set.seed(10)
use_samples <- wts_tcga_flynn_a5_anno %>% 
  filter(Sample %in% complete_samples, 
         new_naming != "C2C (Cortical admixture)", 
         !(Sample %in% small_rna_outliers)) %>% 
  filter(!(new_naming %in% c("A5 - NF1","A5 - VHL"))) %>% 
  mutate(new_naming=gsub("A5.+", "A5",new_naming)) %>% 
  split(x = ., f = .$new_naming) %>% 
  purrr::map(.f= ~sample(.x$Sample, 5, replace = T,)) %>% unlist()

use_samples <- unique(c(use_samples,small_rna_outliers))

plot_data_expr_meth.use  <- plot_data_expr_meth %>% filter(Sample %in% use_samples)
plot_data_14q_cn.use  <- plot_data_14q_cn %>% filter(Sample %in% use_samples)

# hc_input <- plot_data_expr_meth.use %>%  
#   pivot_wider(id_cols = c(Gene, data_type, source), 
#               names_from = Sample, values_from = value ) %>% 
#   unite(Gene, data_type, source, sep="_", col = "id")
# rownames(hc_input) <- hc_input$id
# hc_input <- hc_input %>%  dplyr::select(-id)
# hc_result.sample <- hclust(dist(hc_input))
# sample_order <- colnames(hc_input)[hc_result.sample$order]
# 


plot_data_annotation_bar <- wts_tcga_flynn_a5_anno %>% 
  dplyr::select(Sample, new_naming) %>%  
  filter(Sample %in% use_samples) %>% 
  mutate(smallrna_outlier=ifelse(Sample %in% small_rna_outliers, "Yes","No")) %>% 
  arrange(smallrna_outlier, new_naming) %>% 
  mutate(Sample=factor(Sample, levels=.$Sample)) %>% 
  pivot_longer(cols = -Sample, names_to="Feature", values_to = "Group")

gg_outlier <- ggplot(plot_data_annotation_bar %>%  
                       filter(Feature=="smallrna_outlier"), 
                     aes(x=Sample, y=Feature, fill=Group)) + 
  geom_tile() + 
  annotate("text", x=18, y=1, label="No") +
  annotate("text", x=48, y=1, label="Yes", color="white") +
  theme_bw() +
  scale_fill_manual(values=c("Yes"="Black",
                             "No"="white")) + 
  theme(panel.grid = element_blank(), axis.title.y = element_blank())

gg_subtype <- ggplot(plot_data_annotation_bar %>%  filter(Feature=="new_naming"), aes(x=Sample, y=Feature, fill=Group)) + geom_tile() + theme_bw() +
  scale_fill_manual(values=c(subtype_cols,
                             "small-RNA Outlier"="grey",
                             "A5 - Extraadrenal"="lightgoldenrod1", 
                             "A5 - Head_neck"="lightgoldenrod1",
                             "A5 - Adrenal"="lightgoldenrod1",
                             "A5 - Unspecified"="lightgoldenrod1")) + 
  theme(panel.grid = element_blank(), axis.title.y = element_blank())

plot_data_expr_meth.use$Sample <- factor(as.character(plot_data_expr_meth.use$Sample), levels = levels(plot_data_annotation_bar$Sample))


plot_data_expr_meth.use <- crossing(Gene=unique(plot_data_expr_meth.use$Gene),
         Sample=unique(plot_data_expr_meth.use$Sample),
         data_type=unique(plot_data_expr_meth.use$data_type)) %>% 
  full_join(plot_data_expr_meth.use)

gg_meth_exp <- ggplot(plot_data_expr_meth.use %>% 
                        filter(data_type!="promoter") %>% 
                        mutate(value=case_when(
                          value > 2 ~ 2,
                          value < -2 ~ -2,
                          TRUE ~ value))
                        , aes(x=Sample,y=Gene,fill=value,side=data_type)) + 
  geom_split_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white", high="red", midpoint = 0, name="Z-Score") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + ylab("Gene Expression/Methylation(Gene body)")

plot_data_14q_cn <- plot_data_14q_cn %>% 
  filter(Sample %in% use_samples) %>% 
  mutate(Sample=factor(as.character(Sample), 
                       levels = levels(plot_data_annotation_bar$Sample)))

gg_14qcn <- ggplot(plot_data_14q_cn, aes(x=Sample,y="14q CN", fill=cn_state)) + 
  geom_tile() + 
  scale_fill_manual(values=cna_palette) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.title.y = element_blank())

gg_promoter <- ggplot(plot_data_expr_meth.use %>% 
                        filter(data_type=="promoter", !is.na(value)) %>% 
                        mutate(value=case_when(
                          value > 2 ~ 2,
                          value < -2 ~ -2,
                          TRUE ~ value))
                      , aes(x=Sample,y=Gene,fill=value,side=data_type)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white", high="red", midpoint = 0, name="Z-Score") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + ylab("Promoter Methylation")

nox <-  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())

pdf(file =  "./a5/tertiary/results/plots/small_rna_outliers.pdf",width =12 ,height = 6)
gg_outlier + nox + theme(plot.margin = margin(0,6,0,6)) +
  gg_subtype + nox + theme(plot.margin = margin(0,6,0,6)) +
  gg_14qcn + nox + theme(plot.margin = margin(0,6,0,6)) +
  gg_promoter + nox + theme(plot.margin = margin(0,6,0,6)) +
  gg_meth_exp + theme(plot.margin = margin(0,6,0,6)) + 
  plot_layout(ncol = 1, heights=c(1,1,1,2,10), guides = "collect") 
  

x=0.91
y=0.15
height=0.025
width=0.02
grid::grid.rect(x=x,
                y=y,
                height=height,
                width = width, 
                draw = T, 
                gp = gpar(col = "black",
                          fill = "black",
                          lty = 1,
                          lwd = 1 * .pt,
                          linejoin = "mitre",
                          lineend = "butt"))
grid::grid.lines(x=c(x-(width/2),x+(width/2)),
                 y=c(y-(height/2),y+(height/2)), 
                 draw = T, 
                 gp = gpar(col = "white",
                           fill = "white",
                           lty = 1,
                           lwd = 1 * .pt,
                           linejoin = "mitre",
                           lineend = "butt"))

grid::grid.text(label="Exp/Meth", x=x+(width *2.1),y=y)
dev.off()
