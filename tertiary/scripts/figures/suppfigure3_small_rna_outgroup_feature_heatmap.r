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

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

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

a5_seg_keep_gr <- 
  GenomicRanges::makeGRangesFromDataFrame(a5_seg_keep, 
                                          seqnames.field ="chromosome", 
                                          start.field =  "start", 
                                          end.field = "end", 
                                          ignore.strand = T,keep.extra.columns = T)

hits <- findOverlaps(subject = a5_seg_keep_gr, query = cytoband_gr)
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

a5_seg_keep[hits$subjectHits,"cytoband"] <- hits$cytoband
a5_seg_keep[hits$subjectHits,"arm"] <- hits$arm


a5_seg_keep_chr14 <- a5_seg_keep %>% filter(chromosome=="chr14")

#not true sizes just regions represented in segments
chrom_sizes <- a5_seg_keep %>% 
  group_by(chromosome) %>% 
  summarise(size=max(end)-min(start)) %>%  
  {'names<-'(.$size, .$chromosome)}

a5_seg_keep_chr14.summarised <- a5_seg_keep_chr14 %>%  
  mutate(seg_size=end-start) %>% 
  group_by(A5_ID, chromosome, Class) %>% 
  summarise(seg_size=sum(seg_size)) %>% 
  mutate(prop=seg_size/chrom_sizes[chromosome]) 

a5_seg_keep_chr14.summarised <- 
  a5_seg_keep_chr14.summarised %>% 
  group_by(A5_ID) %>% 
  mutate(priority=case_when(
    Class=="Loss" & prop > 0.25 ~ 1,
    Class=="Subclonal Loss" & prop > 0.25 ~ 1,
    Class=="CNLOH" & prop > 0.25 ~ 1,
    Class=="Gain" & prop > 0.25 ~ 2,
    Class=="Diploid/Haploid-X" & prop > 0.25 ~ 3,
    prop > 0.5 ~ 4, 
    prop > 0.25 ~ 5, 
    TRUE ~ 6)) %>% 
  arrange(priority,desc(prop)) %>% 
  slice_head(n = 1) %>% 
  dplyr::select(A5_ID, Class) %>% 
  mutate(Class=dplyr::recode(Class,
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
  a5_seg_keep_chr14.summarised %>% dplyr::rename(Sample=A5_ID))

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
