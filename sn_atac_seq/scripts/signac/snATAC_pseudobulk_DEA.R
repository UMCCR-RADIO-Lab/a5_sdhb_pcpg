
rm(list=ls())

library(tidyverse)
library(patchwork)
library(limma)
library(edgeR)
library(Homo.sapiens)
library(ChIPseeker)
library(qusage)
library(scales)
library(RColorBrewer)
library(ggsci)
library(ggrepel)

#----
# Load Data and metadata into R
#----

# read in the pseudo-bulked scRNA
counts <- read_csv("Data/snATAC-seq/processed/atac_pseudobulk.csv")

# make gene symbols rownames 
counts <- counts %>% column_to_rownames("X1")
# store sample names in vector
sample.names <- colnames(counts)
# read in metadata
sampleinfo <- read_csv("Data/A5_singlenuclei_metadata.csv")
sampleinfo <- sampleinfo %>% 
  mutate(Sample.ID = gsub(pattern = "-", replacement = "_", Sample.ID)) # change   hyphens to underscores for consistency

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

#----
# store counts and metadata in DGEList object
#----

x <- DGEList(counts)
#add sample metadata
group <- sampleinfo$group
sex <- sampleinfo$Predicted.sex
x$samples$group <- group 
x$samples$sex <- sex

#----
# organise peak annotations 
#----
# organise anotations 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)

#store peaks in a char vector
peaks <- rownames(x)
# convert peaks to granges format 
peaks.df <- data.frame(ranges = peaks)
peaks.df <- separate(
  data = peaks.df,
  col = "ranges",
  sep = paste0(":", "|", "-"),
  into = c("chr", "start", "end")
)
peaks.granges <- makeGRangesFromDataFrame(df = peaks.df)
#annotate peaks by the closest gene using chipseeker
# remove unmapped scaffolds
peaks.granges <- keepStandardChromosomes(peaks.granges, pruning.mode = 'coarse')
# make tagmatrix
tagMatrix <- getTagMatrix(peaks.granges, windows=promoter)
# make peak annotation
peakAnno <- annotatePeak(peaks.granges, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
# get annotation dataframe
peakAnno.df <- as.data.frame(peakAnno@anno)
# get DF of peak anno 
geneid <- unique(peakAnno.df$SYMBOL)
genes.df <- AnnotationDbi::select(Homo.sapiens, keys=geneid, columns=c("SYMBOL","GENENAME", "ENTREZID"), 
                                  keytype="SYMBOL")
# add gene info to the annotation DF
peakAnno.df <- peakAnno.df %>% left_join(genes.df, by = "SYMBOL")

# add peak string column 
peakAnno.df <- peakAnno.df %>% 
  unite(col = peaks,
        seqnames, start, end,
        sep = "-")
peakAnno.df$peaks

#remove all-but-first occurrence of any duplicated peaks
peakAnno.df <- peakAnno.df[!duplicated(peakAnno.df$peaks),]

peakAnno.df %>% filter(SYMBOL == "TERT")

#----
# Data Normalisation and filtering 
#----

libdata <- x$samples %>% 
  rownames_to_column(var = "sample")

ggplot(libdata) +
  aes(x = sample, y = lib.size,  fill = group) +
  geom_bar(stat = "identity") +
  ggtitle("library sizes") +
  theme_classic()

#filter out genes with low expression
keep.exprs <- filterByExpr(x, group=x$samples$group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
# TMM normalise gene expression 
x <- calcNormFactors(x, method = "TMM")

# get log-counts per million of raw counts 
lcpm <- cpm(x, log=TRUE)
# Check distributions of samples using boxplots
boxplot(lcpm, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(lcpm),col="blue")
title("Boxplots of logCPMs (normalised)")

#----
# Explore variation in the dataset with MDS
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

# Tumours appear to cluster by primary driver
ggmds2 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= Sample.ID, colour = Primary.driver))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal3) +
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

(ggmds1 + ggmds2 + ggmds3 + ggmds4 + ggmds5 + ggmds6) +
  plot_layout(ncol=3)

#----
# Differential expression analysis with the limma-voom pipeline
#----

# My contrast is TERT vs ATRX, and I will not include the subclonal TERT sample

group <- sampleinfo$Secondary.driver
sex <- sampleinfo$Predicted.sex
design <- model.matrix(~0+group)
colnames(design)

# make contrast matrix for tert vs atrx
contr.matrix <- makeContrasts(
  TERTvsATRX = groupTERT - groupATRX,
  levels = colnames(design))

#limma pipeline 
v <- voomWithQualityWeights(x, design, plot=T)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
sum.fit <- decideTests(efit)

#check the number of DE genes 
summary(sum.fit)

tt_TERTvsATRX <- topTable(efit, coef = "TERTvsATRX", sort.by = "P", number = Inf)
tt_TERTvsATRX <- rownames_to_column(tt_TERTvsATRX, "peaks")
# add closest gene info to the annotation DF
tt_TERTvsATRX <- tt_TERTvsATRX %>%
  left_join(peakAnno.df, by = "peaks") %>% 
  dplyr::select(peaks, logFC, adj.P.Val, SYMBOL, GENENAME.x)

tt_TERTvsATRX %>% head(30)

#write_csv(tt_TERTvsATRX, "results/limma_da_peaks/tt_TERTvsATRX.csv")

#----
# Plot the results 
#----







