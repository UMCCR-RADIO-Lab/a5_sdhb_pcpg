library(ggplot2)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MutationalPatterns)
library(tidyverse)

genie <- read.delim("/g/data/pq08/databases/genie_r15/data_mutations_extended.txt")
genie_meta <- read.delim("/g/data/pq08/databases/genie_r15/data_clinical_sample.txt", skip =4)
ppgl_muts <- genie %>%  inner_join(genie_meta %>%  filter(ONCOTREE_CODE %in% c("PGNG","PHC")), by=c("Tumor_Sample_Barcode"="SAMPLE_ID"))


dna_repair_genes <- read.delim("/g/data/pq08/reference/gene_lists/dna_damage_repair_pearl_etal_s3.txt")

mmr_genes <- dna_repair_genes %>%  
  filter(Pathway.2 %in% c("MutL homologs", "Mismatch and loop recognition factors")) %>% 
  pull(Gene.ID)

#############
# DFCI Case #
#############

dfci_case_muts <- ppgl_muts %>% 
  filter(Tumor_Sample_Barcode %in% c("GENIE-DFCI-001077-11300","GENIE-UHN-OCT385593-ARC1")) %>% 
  dplyr::select(Tumor_Sample_Barcode, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2) %>% 
  mutate(Chromosome = paste0("chr",Chromosome)) %>% 
  dplyr::rename("REF"=Reference_Allele, "ALT"=Tumor_Seq_Allele2) %>% 
  filter(ALT %in% c("A","C","T","G"), REF %in% c("A","C","T","G"))

dfci_case_muts.gr <- GenomicRanges::makeGRangesFromDataFrame(df = dfci_case_muts, seqnames.field = "Chromosome", 
                                                        start.field = "Start_Position", 
                                                        end.field = "End_Position", keep.extra.columns = T)  
genome(dfci_case_muts.gr) <- "hg19"

##############
# ctDNA case #
##############

cttso_muts <- read.delim("/g/data/pq08/projects/ppgl/a5/ctdna/PRJ230064_L2300271_TMB_Trace.tsv") 
cttso_muts <- cttso_muts %>%  
  dplyr::rename("REF"=RefCall, "ALT"=AltCall) %>% 
  filter(IncludedInTMBNumerator == "True", ALT %in% c("A","C","T","G"), REF %in% c("A","C","T","G"))
cttso_muts.gr <- GenomicRanges::makeGRangesFromDataFrame(df = cttso_muts, seqnames.field = "Chromosome", 
                                                         start.field = "position", 
                                                         end.field = "position", keep.extra.columns = T)
genome(cttso_muts.gr) <- "hg19"


###################
# mut matrix hg19 #
###################

snv_hg19 <- get_mut_type(list(
  DFCI=dfci_case_muts.gr[mcols(dfci_case_muts.gr)[["Tumor_Sample_Barcode"]] == "GENIE-DFCI-001077-11300",],
  UHN=dfci_case_muts.gr[mcols(dfci_case_muts.gr)[["Tumor_Sample_Barcode"]] == "GENIE-UHN-OCT385593-ARC1",],
  ctTSO=cttso_muts.gr), type = "snv")
mut_mat_hg19 <- mut_matrix(vcf_list = snv_hg19, ref_genome = BSgenome.Hsapiens.UCSC.hg19)

##########
# E167-1 #
##########

ref_genome <- ""

vcf_files <- list("E167-M1"="/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E167/umccrised/E167-2__E167-T02/small_variants/E167-2__E167-T02-somatic-PASS.vcf.gz",
                  "E167-M2"="/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E167/final/2021-10-28_E167/E167-1-ensemble-annotated.vcf.gz")

E167_grl <- read_vcfs_as_granges(vcf_files, names(vcf_files), genome=BSgenome.Hsapiens.UCSC.hg38, type = "all")

E167_snv <- get_mut_type(E167_grl, type = "snv")

mut_mat_hg38 <- mut_matrix(vcf_list = E167_snv, ref_genome = BSgenome.Hsapiens.UCSC.hg38)

#########
# merge #
#########

mut_mat <- cbind(mut_mat_hg19, mut_mat_hg38)

###############
# SNV 96 plot #
###############

plot_data <- mut_mat %>%
  as_tibble(rownames="context") %>%
  pivot_longer(cols = -context, names_to = "Sample", values_to = "count") %>%
  mutate(Alteration = stringr::str_extract(context,".>."),
         Alteration = factor(Alteration, levels=c("C>A","C>G","C>T","T>A","T>C","T>G")),
         context=gsub("\\[.>.\\]",".",context)) %>% 
  mutate(Sample = factor(Sample, levels = c("E167-M1","E167-M2","DFCI", "UHN","ctTSO")))

gg_snv96 <- ggplot(plot_data, aes(x=context, y=count,fill=Alteration)) + 
  geom_col() + 
  scale_fill_manual(values=c("#03bdee","#000000","#e52a25","#cdc9ca","#a3ce62","#edc5c5")) +
  facet_grid(Sample~Alteration, scale="free_y") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))


################
# GENIE counts #
################

genie_mmr <- ppgl_muts %>%  
  filter(Hugo_Symbol %in% mmr_genes) %>%  
  mutate(t_vaf = t_alt_count/ (t_ref_count + t_alt_count)) %>% 
  dplyr::select(Hugo_Symbol, Tumor_Sample_Barcode, SAMPLE_TYPE, ONCOTREE_CODE, Consequence, Polyphen_Prediction, t_vaf) %>% 
  filter(grepl("_damaging", Polyphen_Prediction)) %>%  #, t_vaf > 0.1
  distinct()
  

total_counts <- ppgl_muts %>% group_by(Tumor_Sample_Barcode) %>% dplyr::count() %>% left_join(genie_mmr %>%  group_by(Tumor_Sample_Barcode) %>%  summarise(Hugo_Symbol = paste(Hugo_Symbol, collapse = "+"))) %>% 
  arrange(n) %>%  ungroup() %>% mutate(Tumor_Sample_Barcode = factor(Tumor_Sample_Barcode, levels=Tumor_Sample_Barcode))

gg_geniecounts <- ggplot(total_counts %>% 
               mutate(Hugo_Symbol=replace_na(Hugo_Symbol,"None")), 
             aes(x=Tumor_Sample_Barcode, y=n, fill=Hugo_Symbol)) + 
  geom_col() + 
  scale_fill_manual(values=c(MSH6="red",MLH1="blue",MSH2="darkgreen", "None"="grey70")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


gg_geniecounts + theme(axis.text.x = element_blank(),
                       axis.title.x = element_blank(),
                       axis.ticks.x = element_blank()) +
  gg_snv96 + 
  plot_layout(nrow=2, heights=c(1,3))

