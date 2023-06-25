library(tidyverse)
library(Glimma)
library(Limma)
library(edgeR)
library(ComplexHeatmap)

# Blang ggplot themeing 
blank_theme <- theme_bw()+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

#### TCGA and A5 pilot miRNA analysis ----

# Try with the counts mapped to mirbase transcripts
mirbase_ref_FC <- read.table("/Users/adpattison/Documents/Projects/Year_2018/A5study/small_RNA_seq/featureCounts/A5_subread_counts_miRbase_ref.tsv", header = T,sep = "\t")

counts_a5 <- mirbase_ref_FC[,7:ncol(mirbase_ref_FC)]

rownames(counts_a5) <- mirbase_ref_FC$Geneid

colnames(counts_a5) <- gsub("\\.sR.*|bams\\.|\\.subread.*", "", colnames(counts_a5))

# keep <- rowSums(counts>1)>=1
keep_a5 <- rowSums(cpm(counts_a5)>1)>=3

counts_a5 <- counts_a5[keep_a5,]

# Plot CPMs with batch effect still present
Heatmap(counts_a5,name = "Log2 CPM")

# TCGA counts mapped to mirbase transcripts
TCGA_ref_FC <- read.table("/Users/adpattison/Documents/Projects/Year_2018/A5study/small_RNA_seq/featureCounts/TGCA_subread_counts_miRbase_ref.tsv", header = T,sep = "\t")

# Both were aligned in the same way and have the same gene order so I can join them

all_counts <- cbind(mirbase_ref_FC, TCGA_ref_FC[,7:ncol(TCGA_ref_FC)])

annotation <-all_counts$Geneid
counts <- all_counts[,7:ncol(all_counts)]

rownames(counts) <- annotation
colnames(counts) <- gsub("\\.sR.*|bams\\.|\\.subread.*", "", colnames(counts))

# Comapre library sizes
groups <- c(rep("A5",6),rep("Brain",1), rep("TCGA",187))

lib_data <- data.frame(Sample = colnames(counts),Lib_sizes =  colSums(counts)/1E6,groups = groups)%>%
  arrange(-Lib_sizes)%>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

ggplot(data = lib_data, aes(x =Sample, y =  Lib_sizes, fill = groups))+
  geom_bar(stat = "identity")+
  blank_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Sample", y = "miRNA library size (M)")

# Keep rows with at least one value
# keep <- rowSums(counts>1)>=1
keep <- rowSums(cpm(counts)>1)>=3
# Filter counts and annotation
counts <- counts[keep,]
annotation <- annotation[keep]

# Keep only the samples that we have annotation for
rnaseq <- counts%>%
  DGEList(genes = annotation)

rnaseq <- calcNormFactors(rnaseq)

# Brain sample clusters away from the others which is good
# Clear plaform-specific bias
glMDSPlot(rnaseq,labels = colnames(rnaseq),groups = groups,  path = "/Users/adpattison/Documents/Projects/Year_2018/A5study/small_RNA_seq/featureCounts/")

# Get some idea of CPM from miRNA seq
cpm_miRNA <- cpm(rnaseq, log =T)

Heatmap(cpm_miRNA)

# Z score transform the log2 CPMs for plotting
scaled_cpm_miRNA <- cpm_miRNA  %>%
  t()%>%
  scale()%>%
  t()

Heatmap(scaled_cpm_miRNA)

# Plot TCGA mean miRNA expresion against A5
compare_df <- data.frame(miR = rownames(cpm_miRNA),A5_mean = rowMeans(cpm_miRNA[,groups == "A5"]),
                         TCGA_mean = rowMeans(cpm_miRNA[,groups == "TCGA"]))

correlation <- cor.test(compare_df$A5_mean, compare_df$TCGA_mean)
estimate <- paste("r =",round(correlation$estimate,2))

# Shows good correlation between TCGA and A5 pilot miRNAs
ggplot(data = compare_df,(aes(x = A5_mean, y = TCGA_mean)))+ geom_point()+
  annotate("text", label = estimate, x =15, y = 16)+
  blank_theme+
  labs(x = "A5 miRNA expression mean", y = "TCGA miRNA expression mean")

# Gather for later dplyr style use
gathered_CPM <- cpm_miRNA_A5 %>%
  data.frame()%>%
  mutate(Geneid = rownames(cpm_miRNA_A5))%>%
  gather(Sample, CPM, -Geneid)%>%
  arrange(-CPM)

