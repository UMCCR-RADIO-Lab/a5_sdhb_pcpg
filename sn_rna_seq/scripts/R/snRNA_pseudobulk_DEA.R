# Differential expression and Gene-set enrichment Analysis for Pseudobulked snRNA-seq
rm(list=ls())

library(tidyverse)
library(patchwork)
library(limma)
library(edgeR)
library(Homo.sapiens)
library(qusage)
library(scales)
library(RColorBrewer)
library(ggsci)
library(ggrepel)

#----
# Load Data and metadata into R
#----

# read in the pseudo-bulked scRNA
#counts <- read_csv("/Users/blake/OneDrive - The University of Melbourne/Projects/Bioninformatics Project/PseudoBulk/pseudobulk_tightqc.csv")
counts <- read_csv("Data/snRNA-seq/processed/snRNA_pseudobulk.csv")
# remove E174 and E018 as they do no have SDHB- genotype 

# make gene symbols rownames 
counts <- counts %>% column_to_rownames("X1")

# reformat the colnames for the tumour samples
colnames(counts) <- if_else(grepl(pattern = "Tumour", x = colnames(counts)),
  sub(pattern = "\\.", replacement = "-", x = colnames(counts), colnames(counts)),
  colnames(counts))
colnames(counts) <- if_else(!grepl(pattern = "NAM0", x = colnames(counts)),
                            substr(colnames(counts), start = 1, stop = 6),
                            colnames(counts))

# store sample names in vector
sample.names <- colnames(counts)

# ----
# make a blacklist of adrenocortical genes
# ----

# make a matrix containing only adrenocortical counts
adrenocortical_counts <- as.matrix(counts[,colnames(counts) %in% c("NAM021_Adrenocortical.cells", "NAM025_Adrenocortical.cells")])

filterbyexpr()

# ---- 
# Organise the metadata 
# ---- 

# read in the clinical data 
sampleinfo <- read_csv("Data/A5_full_clinical_and_genomic_version_2.csv")
# filter to retain only the snRNA-seq samples
sampleinfo <- sampleinfo %>%
  dplyr::filter(`A5 ID` %in% sample.names) %>%
  mutate(sample = `A5 ID`)

# make new entry for the NAM samples 
sampleinfo <- sampleinfo %>%
  add_row(sample = "NAM021_Chromaffin.cells", Gender = "male", `Tumour metastasised` = "Normal", `Metastatic state of tumour sample` = "Normal") %>% 
  add_row(sample = "NAM021_Adrenocortical.cells", Gender = "male", `Tumour metastasised` = "Normal", `Metastatic state of tumour sample` = "Normal") #%>% 
  # add_row(sample = "NAM025_Chromaffin.cells", Gender = "male", `Tumour metastasised` = "Normal", `Metastatic state of tumour sample` = "Normal") %>% 
  # add_row(sample = "NAM025_Adrenocortical.cells", Gender = "male", `Tumour metastasised` = "Normal", `Metastatic state of tumour sample` = "Normal")

# add a column describing the 10x genomics chemistry that was used 
sampleinfo <- sampleinfo %>% 
  mutate(chemistry = if_else(sample %in% c("E140-1", "E143-1", "E171-1", "E225-1"), "3prime", "5prime"))

# make column describing primary and secondary driver mutation 
sampleinfo <- sampleinfo %>%
  mutate(metastasis_driver = case_when(
    sample %in% c("NAM021_Chromaffin.cells", "NAM025_Chromaffin.cells", "NPG103") ~ "Normal chromaffin",
    sample %in% c("NAM021_Adrenocortical.cells", "NAM025_Adrenocortical.cells") ~ "Normal adrenocortical",
    sample %in% c("E166-1", "E197-1", "E140-1", "E225-1") ~ "ATRX",
    sample %in% c("E146-1") ~ "TERT_subclonal",
    sample %in% c("E156-1","E123-1", "E143-1") ~ "TERT", 
    sample == "E171-1" ~ "Unknown_driver_metastatic"))

# reorder the DF
sample.order <- tibble(sample = sample.names)
sampleinfo <- sample.order %>%
  left_join(sampleinfo)

#----
# store counts and metadata in DGEList object
#----

# remove the adrenocotrical counts
counts <- counts[,!colnames(counts) %in% c("NAM021_Adrenocortical.cells", "NAM025_Adrenocortical.cells")]

x <- DGEList(counts)
# store the important stuff that i will later use in design matrix
group <- sampleinfo$metastasis_driver
sex <- sampleinfo$Gender
chemistry <- sampleinfo$chemistry
# add sample metadata to dgelist
x$samples$group <- group 
x$samples$sex <- sex

#----
# organise gene annotations 
#----

geneid <- rownames(x)
genes <- AnnotationDbi::select(Homo.sapiens, keys=geneid, columns=c("SYMBOL","GENENAME"), 
                               keytype="SYMBOL")
#remove all-but-first occurrence of any duplicated genes
genes <- genes[!duplicated(genes$SYMBOL),]
table(rownames(x$counts) == genes$SYMBOL)
#store genes in DGEList
x$genes <- genes

#----
# Data Normalisation and filtering 
#----

libdata <- x$samples %>% 
  rownames_to_column(var = "sample")

ggplot(libdata) +
  aes(x = sample, y = lib.size,  fill = group) +
  geom_bar(stat = "identity") +
  ggtitle("library sizes") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

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
# MDS plots for unsupervised clustering 
#----

# make MDS plot object
mds <- plotMDS(lcpm)

# make tbl which can be used for ggplot2 input
mds_ggplot <- tibble(
  x = mds$x,
  y =  mds$y,
  "sample" = colnames(mds$distance.matrix))%>%
  left_join(sampleinfo, by = "sample")

#organise distinct colour palette for each MDS
pal <- pal_d3(palette = "category20")(20)
pal2 <- pal[1:7]
pal3 <- pal[8:13]
pal4 <- pal[14:19]
pal5 <- pal_d3(palette = "category20c")(20)
pal6 <- pal5[1:3]
pal7 <- pal5[4:7]
pal8 <- pal5[9:10]

# Plot by group 
ggmds1 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= sample, colour = group))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal2) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

ggmds2 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= sample, colour = chemistry))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal3) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

# Tumours also appear to cluster by secondary driver
ggmds3 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= sample, colour = metastasis_driver))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal4) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

ggmds4 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= sample, colour = `Metastatic state of tumour sample`))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal6 )+
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

ggmds5 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= sample, colour = `Tumour metastasised`))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal7) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

ggmds6 <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= sample, colour = Gender))+
  geom_text() + theme_classic() + scale_colour_manual(values = pal8) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

(ggmds1 + ggmds2 + ggmds3 + ggmds4 + ggmds5 + ggmds6) +
  plot_layout(ncol=3)

#----
# Make a list of genes expressed above threshold in the adrenocortical cluster 
#----

# make list of genes likely to be in adrenocortical-derived ambient RNA
# These will be filtered out of the tumour v normal comparison









#----
# Differential expression analysis with the limma-voom pipeline
#----

chemistry <- sampleinfo$chemistry
# make group variable describing TERT mutation status and primary driver
# group <- gsub("ATRX", "WT", sampleinfo$group)
design <- model.matrix(~0+group+sex+chemistry)
colnames(design)

# make contrast matrix for each comparison I want to make:
# TERT mutant SDHB tumours  vs normal chromaffin cells
# TERT mutant SDHB tumours vs ATRX SDHB mutant tumour
# TERT mutant SDHB tumour vs wt SDHB

contr.matrix <- makeContrasts(
  TERTvsNormal = groupTERT - groupNormal,
  TERTvsATRX = groupTERT - groupATRX,
  ATRXvsNormal = groupATRX - groupNormal,
  levels = colnames(design))

#limma pipeline 
v <- voom(x, design, plot=T)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
sum.fit <- decideTests(efit)

#check the number of DE genes 
summary(sum.fit)

# store the top diff. expressed genes 
tt_TERTvsNormal <- topTable(efit, coef = "TERTvsNormal", p.value = 0.05, sort.by = "P", number = Inf)
tt_TERTvsATRX <- topTable(efit, coef = "TERTvsATRX", sort.by = "P", number = Inf)
tt_ATRXvsNormal <- topTable(efit, coef = "ATRXvsNormal", p.value = 0.05, sort.by = "P", number = Inf)


#----
# Plots after DE testing
#----

#venndiagram
par(mfrow=c(1,1))
vennDiagram(sum.fit)

# which genes are differentially expressed in both?
sum.fit[which(sum.fit[, "TERTvsChromaffin"] != 0 & sum.fit[, "TERTvsATRX"] != 0),]

# MDplot, highlight significant genes
par(mfrow=c(1,2))
plotMD(efit,coef="TERTvsChromaffin",
       status=sum.fit[,"TERTvsChromaffin"],
       values = c(-1, 1))
title("TERTvsChromaffin")
plotMD(efit,coef="TERTvsATRX",
       status=sum.fit[,"TERTvsATRX"],
       values = c(-1, 1))
title("TERTvsATRX")

#----
# Gene set enrichment analysis 
#----

#read in the MSigDB TFT gene sets 
tft_genes <- qusage::read.gmt("/Users/blake/OneDrive - The University of Melbourne/Projects/Honours Bioinformatics Project/bulk_rna/GSEA/c3.tft.v7.1.symbols.gmt")
indexed_for_camera <- ids2indices(gene.sets = tft_genes, identifiers = genes$SYMBOL)

TERTvsATRX_camera <- camera(
  y = v,
  index = indexed_for_camera,
  design = design,
  contrast = contr.matrix[,"TERTvsATRX"])
head(TERTvsATRX_camera)

par(mfrow=c(1,1))
barcodeplot(efit$t[,"TERTvsATRX"],
            index = indexed_for_camera[["ATM_TARGET_GENES"]],
            main = "TERTvsATRX: ATM_TARGET_GENES")

TERTvsChromaffin_camera <- camera(
  y = v,
  index = indexed_for_camera,
  design = design,
  contrast = contr.matrix[,"TERTvsChromaffin"])
head(TERTvsChromaffin_camera)

par(mfrow=c(1,1))
barcodeplot(efit$t[,"TERTvsChromaffin"],
            index = indexed_for_camera[["BAHD1_TARGET_GENES"]],
            main = "TERTvsChromaffin: BAHD1_TARGET_GENES")

#----
# Plot specific gene expression
#----

# convert data to longform for ggplot2 compatibility
data <- rownames_to_column(as.data.frame(lcpm), "gene")
data <- pivot_longer(
  data,
  cols = !c(gene),
  names_to = "sample",
  values_to = "Log2Cpm") %>% 
  left_join(sampleinfo, by = "sample")

# make a function that generates a barplot for gene expression across all samples 
# input will be the gene symbol
plot_cpm <- function(genesymbol=NULL){
  gene.use <- genesymbol
  if (gene.use %in% data$gene) {
    genedata <- data %>% 
      dplyr::filter(gene == gene.use) %>%
      dplyr::select(sample, Log2Cpm, group) %>% 
      dplyr::group_by()
    plot <- ggplot(genedata) +
      aes(x = sample, y = Log2Cpm,  fill = group) +
      geom_bar(stat = "identity") +
      ggtitle(paste(gene.use, "Pseudobulk gene-expression")) +
      theme_classic()
    return(plot)
  }
  else{
    print("Invalid gene symbol")
  }  
}

plot_cpm(genesymbol = "TERT")
plot_cpm(genesymbol = "TMEFF2")
plot_cpm(genesymbol = "SNAI2")
plot_cpm(genesymbol = "ATF3")
plot_cpm(genesymbol = "XIST")
plot_cpm(genesymbol = "TSIX")
plot_cpm(genesymbol = "AFF4")


#----
# Which DAR annotations intersect with the TERT-enriched genes?
#----

enriched_genes <- tt_TERTvsATRX_plot %>% filter(adj.P.Val < 0.05) %>% pull(SYMBOL)
tert_anno <- read_csv("/Users/blake/OneDrive - The University of Melbourne/Projects/Honours Bioinformatics Project/snATAC/Peaks/tert_annotated_peaks")
enriched_genes %in% tert_anno$SYMBOL

#----
# Make nice volcano plot 
#----

tt_TERTvsATRX_plot <- topTable(efit, coef = "TERTvsATRX", number = Inf, sort.by = "P")
# subtract genes differentially expressed in male and female samples 
genes.remove <- tt_FemalevsMale$SYMBOL
tt_TERTvsATRX_plot <- tt_TERTvsATRX_plot %>% 
  filter(!SYMBOL %in% genes.remove)
genes.highlight <- tt_TERTvsATRX_plot %>%
  filter(GENENAME != "<NA>") %>% 
  slice_min(n = 10, order_by = adj.P.Val) %>% pull(SYMBOL)
tt_TERTvsATRX_plot <- tt_TERTvsATRX_plot %>% 
  #flag the most significant genes to label
  mutate(sig = if_else(adj.P.Val < 0.05, "Adj. P Value < 0.05", "not significant")) %>% 
  mutate(highlight = if_else(SYMBOL %in% genes.highlight, TRUE, FALSE))


tert_atrx_volcano <- ggplot(tt_TERTvsATRX_plot, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black")) + theme_classic() +
  geom_text_repel(
    data = filter(
      tt_TERTvsATRX_plot, highlight == TRUE),
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.1, "lines"),
    segment.color = "grey",
    size = 2,
    aes(label=SYMBOL)) + ggtitle("")

tert_atrx_volcano

#colours <- sample(hue_pal()(12))
mdscolours <- c("#619CFF", "#B79F00", "#F8766D", "#00BA38", "#F564E3", "#00BFC4", "#FF64B0", "#C77CFF", "#DE8C00", "#00B4F0", "#00C08B", "#7CAE00")

secondarydriver_MDS <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= sample))+
  geom_point(aes(colour = Secondary.driver)) +
  geom_text_repel(size = 2)+ theme_classic() + scale_colour_manual(values = mdscolours[1:5]) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

tumorstate_MDS <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= sample))+
  geom_point(aes(colour = Tumour.state)) +
  geom_text_repel(size = 2)+ theme_classic() + scale_colour_manual(values = mdscolours[6:8] )+
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

malig_MDS <- ggplot(data = mds_ggplot, aes(x = x, y = y, label= sample))+
  geom_point(aes(colour = Malignancy)) +
  geom_text_repel(size = 2)+
  theme_classic() + scale_colour_manual(values = mdscolours[9:12]) +
  xlab("Leading LogFC Dim 1") + ylab("Leading LogFC Dim 2")

plots <- secondarydriver_MDS + tert_atrx_volcano +
  plot_layout(ncol = 2) & 
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size=8),
    legend.title = element_text(size=8),
    axis.text = element_text(size=8),
    legend.key.size = unit(8, "points"))

plots

ggsave(plot = plots, filename = "~/Desktop/MDS_Volcano_plots.pdf", width = 19, height = 7, units = "cm",  dpi = 300, useDingbats = FALSE)

degs_table <- tt_TERTvsATRX_plot %>% 
  filter(sig == "Adj. P Value < 0.05") %>% 
  dplyr::select(SYMBOL, logFC, adj.P.Val)

write_csv(x = degs_table, path = "~/Desktop/DEGs_pseudobulk.csv")

######

#----------------------------
#determine the top sex-linked genes 
sex <- sampleinfo$Predicted.sex
design <- model.matrix(~0+sex)
contr.matrix <- makeContrasts(
  FemalevsMale = sexF - sexM,
  levels = colnames(design))

#limma pipeline 
v <- voom(x, design, plot=F)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
sum.fit <- decideTests(efit)
summary(sum.fit)

tt_FemalevsMale <- topTable(efit, coef = "FemalevsMale", p.value = 0.05, sort.by = "P")

#--------------------------