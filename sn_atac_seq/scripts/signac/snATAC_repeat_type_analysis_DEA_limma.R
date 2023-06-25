
rm(list=ls())
library(tidyverse)
library(patchwork)
library(limma)
library(edgeR)
library(ggsci)

#----
# Load Data and metadata into R
#----

rpt_counts <- read_csv("Data/repeat_analysis/counts/repeat_type_counts.csv")
rpt_counts <- rpt_counts %>%
  dplyr::rename(NAM018 = "Chromaffin_cell") %>%
  data.frame() %>% 
  column_to_rownames(var = "repeat_type") 

# store sample names in vector
sample.names <- colnames(rpt_counts)
# read in metadata
sampleinfo <- read_csv("Data/A5_singlenuclei_metadata.csv")
sampleinfo <- sampleinfo %>% 
  mutate(Sample.ID = gsub(pattern = "-", replacement = "_", Sample.ID)) # change underscores to hyphens for consistency

# reorder the info 
sample.order <- tibble(Sample.ID = sample.names)
sampleinfo <- sample.order %>% left_join(sampleinfo)

# make column describing primary and secondary driver mutation 
sampleinfo <- sampleinfo %>%
  mutate(
    group = paste(
      Secondary.driver,
      Primary.driver,
      sep = "_")) %>% 
  mutate(group = gsub("SDHA|SDHB", "SDHx",  group)) %>% 
  mutate(group = gsub ("Normal_", "", group)) %>% 
  mutate(TERT_mutation = recode(Secondary.driver, 
                                "TERT_subclonal" = "TERT",
                                "ATRX" = "WT"))
group <- sampleinfo$group

#----
# store counts and metadata in DGEList object
#----

# these library sizes were calculated using the sinto-filtered bams
# using the command 'samtools view -c'
# they represent the total number of alignments present in the neoplastic cells

rpt_dge <- DGEList(rpt_counts, group = group)

#----
# Data Normalisation and filtering 
#----

libdata <- rpt_dge$samples %>% 
  rownames_to_column(var = "sample")

ggplot(libdata) +
  aes(x = sample, y = lib.size,  fill = group) +
  geom_bar(stat = "identity") +
  ggtitle("library sizes") +
  theme_classic()

#filter out repeats with lowest counts
keep.exprs <- filterByExpr(rpt_dge, group=rpt_dge$samples$group)
rpt_dge <- rpt_dge[keep.exprs,, keep.lib.sizes=FALSE]
# TMM normalise gene expression 
rpt_dge <- calcNormFactors(rpt_dge, method = "TMM")

# get log-counts per million of raw counts 
lcpm <- cpm(rpt_dge, log=TRUE)
# Check distributions of samples using boxplots
boxplot(lcpm, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(lcpm),col="blue")
title("Boxplots of logCPMs (normalised)")

#----
# MDS plots for unsupervised clustering 
#----

# make MDS plot object
mds <- plotMDS(lcpm)

# make tbl which can be used for ggplot2 input
mds_ggplot <- tibble(
  x = mds$x,
  y =  mds$y,
  "Sample.ID" = colnames(mds$distance.matrix))%>%
  left_join(sampleinfo, by = "Sample.ID")

#organise distinct colour palette for each MDS
pal <- pal_d3(palette = "category20")(20)
pal2 <- pal[1:7]
pal3 <- pal[8:13]
pal4 <- pal[14:18]
pal5 <- pal_d3(palette = "category20c")(20)
pal6 <- pal5[1:3]
pal7 <- pal5[4:7]
pal8 <- pal5[9:10]

# Plot by group 
ggmds1 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= Sample.ID, colour = group))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal2) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

# Tumours also appear to cluster by secondary driver
ggmds3 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= Sample.ID, colour = Secondary.driver))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal4) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

ggmds4 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= Sample.ID, colour = Tumour.state))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal6 )+
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

ggmds5 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= Sample.ID, colour = Malignancy))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal7) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

ggmds6 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= Sample.ID, colour = Predicted.sex))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal8) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

(ggmds1 + ggmds3 + ggmds4 + ggmds5 + ggmds6) +
  plot_layout(ncol=3)

#----
# Differential expression analysis with the limma-voom pipeline
#----

# make group variable describing TERT mutation status and primary driver
design <- model.matrix(~0+group)
colnames(design)

contr.matrix <- makeContrasts(
  TERTvsATRX = groupTERT_SDHx - groupATRX_SDHx,
  levels = colnames(design))

#limma pipeline 
v <- voomWithQualityWeights(rpt_dge, design, plot=T)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
sum.fit <- decideTests(efit)

#check the number of DE genes 
summary(sum.fit)

# store the top diff. expressed genes 
# tt_TERTvsChromaffin <- topTable(efit, coef = "TERTvsChromaffin", p.value = 0.05, sort.by = "P", number =  282)
# head(tt_TERTvsChromaffin)
tt_TERTvsATRX <- topTable(efit, coef = "TERTvsATRX",sort.by = "P", number = 25)
head(tt_TERTvsATRX, 10)
tt_TERTvsATRX %>% rownames_to_column('gene') %>% as_tibble()

#----
# Plot specific region counts
#----

# convert data to longform for ggplot2 compatibility
data <- rownames_to_column(as.data.frame(lcpm), "gene")
data <- pivot_longer(
  data,
  cols = !c(gene),
  names_to = "Sample.ID",
  values_to = "Log2Cpm") %>% 
  left_join(sampleinfo, by = "Sample.ID")

ggplot(data)+
  geom_boxplot(aes(x =  Secondary.driver , y = Log2Cpm , fill = group)) + 
  geom_jitter(aes(x =  Secondary.driver , y = Log2Cpm))

# make a function that generates a barplot for gene expression across all samples 
# input will be the gene symbol
plot_cpm_bar <- function(genesymbol=NULL){
  gene.use <- genesymbol
  if (gene.use %in% data$gene) {
    genedata <- data %>% 
      dplyr::filter(gene == gene.use) %>%
      dplyr::select(Sample.ID, Log2Cpm, group) %>% 
      dplyr::group_by()
    plot <- ggplot(genedata) +
      aes(x = Sample.ID, y = Log2Cpm,  fill = group) +
      geom_bar(stat = "identity") +
      ggtitle(paste(gene.use, "repeat accessibility")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    return(plot)
  }
  else{
    print("error: Invalid gene symbol")
  }  
}

plot_cpm_bar(genesymbol = "Satellite/telo")

plot_cpm_box <- function(genesymbol=NULL){
  gene.use <- genesymbol
  if (gene.use %in% data$gene) {
    genedata <- data %>% 
      dplyr::filter(gene == gene.use) %>%
      dplyr::select(Secondary.driver, Log2Cpm, group) %>% 
      dplyr::group_by()
    plot <- ggplot(genedata) +
      aes(x = Secondary.driver, y = Log2Cpm,  fill = group) +
      geom_boxplot() +
      geom_jitter() + 
      ggtitle(paste(gene.use, "repeat accessibility")) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    return(plot)
  }
  else{
    print("error: Invalid gene symbol")
  }  
}

plot_cpm_box(genesymbol = "Satellite/telo")
