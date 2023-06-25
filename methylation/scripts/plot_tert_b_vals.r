library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

blank_theme <- theme_bw(base_size = 25)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

A5_bvals <- read_csv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/A5 methylation Bvals.csv")
A5_bvals[1:5,1:5]

# Read in the A5 clinical data
A5_clinical <- read_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Clinical data/A5 full clinical and genomic version 2 - A5 full clinical and genomic version 2.csv")

TERT_probe <- A5_bvals %>%
  filter(Probe == "cg11625005")%>%
  gather("A5 ID", "B value" ,-Probe)%>%
  mutate(`A5 ID` = gsub(".*\\.|_D", "", `A5 ID`))%>%
  mutate(`A5 ID` = gsub("_T0", "-", `A5 ID`))%>%
  mutate(`A5 ID` = gsub("-1_1", "-1", `A5 ID`))%>%
  mutate(`A5 ID` = gsub("-1_2", "-2", `A5 ID`))%>%
  left_join(A5_clinical)%>%
  mutate(T_A_other = replace(TERT_or_ATRX, !TERT_or_ATRX %in% c("TERT", "ATRX"), "Other"))%>%
  group_by(T_A_other)%>%
  mutate(med = median(`B value`))%>%
  arrange(-med)%>%
  ungroup()%>%
  mutate(T_A_other = factor(T_A_other, levels = unique(T_A_other)))

colours <- c("ATRX" ="#1F77B4FF","TERT" = "#FF7F0EFF", "Other" = "lightgrey")
               
pval <- "DE vs non-met primary\nFDR = 0.04\nlog2FC = 0.65"

# Plot the probe
ggplot(data = TERT_probe, aes(x = T_A_other, y = `B value`, fill = T_A_other))+
  geom_boxplot(outlier.alpha = 0)+
  geom_jitter(width = 0.3)+
  blank_theme+
  scale_fill_manual(values = colours)+
  labs(x = "TERT/ATRX status", y = "cg11625005 B value")+
  guides(fill = F)+
  annotate(geom="text", x=1, y=0.5, label=pval,
           color="black", size =5)+
  ggsave("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Probe_specific_plots/TERT_and_ATRX_cg11625005.pdf")


annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Get the methylation annotation
anno_df_meth <- data.frame(annEPIC)%>%
  rownames_to_column("Probe")%>%
  dplyr::select(Probe, chr, pos, SYMBOL = GencodeBasicV12_NAME, Regulatory_Feature_Group)%>%
  mutate(SYMBOL = gsub(";.*", "", SYMBOL))%>%
  filter(SYMBOL != "")

TERT_probes <-anno_df_meth %>%
  filter(SYMBOL == "TERT")

# Save a bed file of TERT probes
TERT_bed <- TERT_probes%>%
  mutate(end = pos+50)%>%
  select(chrom = chr, chromStart = pos, chromEnd = end, name = Probe)%>%
  write_tsv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/TERT_probes.bed",col_names = F)

# Read in RNA-Seq lcpms
lcpms <- read_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/RNA-Seq/Counts/RNA-Seq genewise log2 CPMs.csv")%>%
  dplyr::select(-ENSEMBL, -TXCHROM)%>%
  filter(!is.na(SYMBOL))%>%
  gather(Sample, lcpm, - SYMBOL)%>%
  mutate(Sample = gsub("\\.", "_", Sample))

bvals_long <- A5_bvals %>%
  filter(Probe %in% TERT_probes$Probe)%>%
  gather(Sample, B, -Probe)%>%
  mutate(Sample = gsub(".*\\.|_D|T0", "", Sample))

# Get an average B value per gene per sample
anno_df_meth_promo_TERT <- TERT_probes%>%
  left_join(bvals_long)%>%
  filter(!is.na(Sample))%>%
  # Get only promoter region
  filter(Probe == "cg11625005")%>%
  left_join(lcpms)%>%
  filter(!is.na(lcpm))%>%
  ungroup()%>%
  mutate(`A5 ID` = gsub("_", "-", Sample))%>%
  left_join(A5_clinical)

colours <- c("TERT" = "#FF7F0EFF","ATRX" =  "#1F77B4FF","Unknown_met_driver_met" = "dark red","Unknown_met_driver_primary" = "purple",
             "Short_follow_up_primary" = "grey","Non_met_primary" = "#2CA02CFF")

cor.test(anno_df_meth_promo_TERT$B, anno_df_meth_promo_TERT$lcpm)

ggplot(data = anno_df_meth_promo_TERT, aes(y = B, x = lcpm, colour = TERT_or_ATRX, label = `A5 ID`))+
  geom_point()+
  geom_text_repel(size =2)+
  blank_theme+
  labs(y = "cg11625005 B value", x = "Log2 CPM")+
  ggtitle("TERT")+
  coord_equal()+
  theme(aspect.ratio = 1)+
  scale_colour_manual(values = colours)+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Gene expression vs methylation/TERT cg11625005 vs gene expression.pdf", width = 16, height = 8.5)

ggplot(data = anno_df_meth_promo_TERT[anno_df_meth_promo_TERT$TERT_or_ATRX == "TERT",], aes(y = B, x = `sample_purity`, colour = TERT_or_ATRX, label = `A5 ID`))+
  geom_point()+
  geom_text_repel(size =2)+
  blank_theme+
  labs(y = "cg11625005 B value", x = "sample_purity")+
  ggtitle("TERT")+
  coord_equal()+
  theme(aspect.ratio = 1)+
  scale_colour_manual(values = colours)+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Gene expression vs methylation/TERT cg11625005 vs purity.pdf", width = 16, height = 8.5)
