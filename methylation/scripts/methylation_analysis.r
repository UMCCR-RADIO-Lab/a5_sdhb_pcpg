# Load packages required for analysis
library(tidyverse)
library(knitr)
library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(stringr)
library(minfiData)
library(maxprobes)
library(GSA)
library(egg)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(DMRcate)
library(Gviz)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(ggrepel)
library(umap)

# Set ggplot2 themes
blank_theme <- theme_bw(base_size = 18)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

# A5 Illumina infinum bisulfite sequencing analysis
# Following this guide modified for EPIC arrays:
# https://bioconductor.org/packages/devel/workflows/vignettes/methylationArrayAnalysis/inst/doc/methylationArrayAnalysis.html

# 'The manifest contains all of the annotation information for each of the CpG probes on the 450k array'

# Arrays used are infinium-methylationepic-v-1-0-b5
dataDirectory <- system.file("extdata", package = "methylationArrayAnalysis")
ins <- "/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Raw data/Methylation_arrays/ILMLEPIC-16614/"
inputs <- list.files(ins, recursive = TRUE)

# Take a loot at the annotation
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
head(annEPIC)

# Find the sample sheet
targets <- read.metharray.sheet(ins, pattern="SampleSheet.csv")%>%
  mutate(A5_ID = gsub("_D", "", Sample_Name))%>%
  mutate(A5_ID = gsub("_T0", "-", A5_ID))

# Read in the IDAT files
rgSet <- read.metharray.exp(targets=targets)
rgSet

# Give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID
rgSet

# P values that a probe was detected
detP <- detectionP(rgSet)
head(detP)

# Plot mean detection p-values 
# All samples seem to be way under 0.05 so looking good
pal <- rainbow(12)
p_val_df <- data.frame(ID= targets$ID, p_means = colMeans(detP),group =targets$Sample_Group, col=pal[factor(targets$Sample_Group)])

ggplot(data = p_val_df, aes(x = ID, y = p_means, fill = group))+
  geom_bar(stat= "identity")+
  scale_colour_manual(values = col)+
  blank_theme+
  labs(x = 'Sample', y = "Mean probe detection\n p value", fill = "Group")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/QC/Probe expression p values.pdf", width = 10, height = 7)

# Make a QC report
qcReport(rgSet, sampNames=targets$ID, sampGroups=targets$Sample_Group, 
         pdf="/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Outputs/QC/qcReport.pdf")

# Accoding to my guide preprocessQuantile is the normalisation I should 
# use as all the data is from a similar tissue type
#mSetSq <- preprocessQuantile(rgSet) 

# Save the object for easier loading later
#saveRDS(mSetSq, "/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Raw data/mSetSq_quantile.rds")

mSetSq <- readRDS("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Raw data/mSetSq_quantile.rds")

# create a MethylSet object from the raw data for plotting
#mSetRaw <- preprocessRaw(rgSet)
#saveRDS(mSetRaw, "/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Raw data/rgSet_plotting.rds")

mSetRaw <- readRDS("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Raw data/rgSet_plotting.rds")

# Look at the data before and after normalisation
# Looks like normalisation helped but there is still a bit of variation
pdf("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/QC/Quantile normalisation.pdf", width =10)
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE,
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE,
            cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
dev.off()

# ensure probes are in the same order in the mSetSq and detP objects
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 

# remove any probes that have failed in one or more samples
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

# Data to keep
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

# Remove probes with SNPs at CpG site
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt

# Get cross reactive probes from the maxprobes library
xReactiveProbes <- maxprobes::xreactive_probes(array_type = "EPIC")%>%
  unlist()
length(xReactiveProbes)

# Exclude cross reactive probes 
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes)
table(keep)

mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt

# Read in the A5 clinical data
A5_clinical <- read_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Clinical data/A5 full clinical and genomic version 2 - A5 full clinical and genomic version 2.csv")%>%
  dplyr::rename(A5_ID = `A5 ID`)%>%
  # Remove a sample with adrenocortical contamination
  filter(A5_ID != "E154-1")

# Remove that cortical admixture sample 
mSetSqFlt <- mSetSqFlt[,mSetSqFlt$Sample_Name != "E154_T01_D"] 

mvals_save <-  getM(mSetSqFlt)

mVals2 <- data.frame(mvals_save)%>%
  rownames_to_column("Probe")%>%
  write_csv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/A5 methylation Avals.csv")

set.seed(42)
# Plot a umap
umap_M <- mSetSqFlt%>%
  getM()%>%
  t()%>%
  umap()

to_plot_umap <- data.frame(umap_M$layout)%>%
  rownames_to_column("A5_ID")%>%
  mutate(A5_ID = gsub(".*\\.", "",A5_ID))%>%
  mutate(A5_ID = gsub("_", "-",A5_ID))%>%
  mutate(A5_ID = gsub("E158-T01-D-1", "E158-T01-D",A5_ID))%>%
  mutate(A5_ID = gsub("E158-T01-D-2", "E158-T02-D",A5_ID))%>%
  mutate(A5_ID = gsub("T0|-D", "", A5_ID))%>%
  left_join(A5_clinical)%>%
  write_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Dataset_integration/UMAP_coordinates/Methylation_probe_M_value_UMAP.csv")

ggplot(data = to_plot_umap, aes(x = X1, y = X2, colour = `is_head_and_neck`)) + 
  geom_point()+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  theme(aspect.ratio=1)

mds <- plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common", 
                col=pal[factor(targets$Sample_Group)])

mds_ggplot_1_2 <- tibble(x = mds$x,y =  mds$y, `A5 tumour ID` = colnames(mds$distance.matrix))%>%
  mutate(Group = gsub("\\..*", "",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub(".*\\.", "",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub("_", "-",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub("E158-T01-D-1", "E158-T01-D",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub("E158-T01-D-2", "E158-T02-D",`A5 tumour ID`))%>%
  left_join(A5_clinical)%>%
  dplyr::select(A5_ID,Group, x, y, predicted_sex, Is_HN)%>%
  mutate(Patient = gsub("-.*", "", A5_ID))

ggplot(data = mds_ggplot_1_2, aes(x = x, y = y,label= A5_ID))+
  geom_text()+
  blank_theme+
  aes(colour = predicted_sex)+
  labs(colour = "Gender", x= "MDS dim 1", y = "MDS dim 2")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/MDS plots/Dimension 1 and 2 gender.pdf")

# Run on dims 2+3 becuase dim 1 is clearly gender
mds <- plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
               col=pal[factor(targets$Sample_Group)],dim=c(2,3))

mds_ggplot <- tibble(x = mds$x,y =  mds$y, `A5 tumour ID` = colnames(mds$distance.matrix))%>%
  mutate(Group = gsub("\\..*", "",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub(".*\\.", "",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub("_", "-",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub("E158-T01-D-1", "E158-T01-D",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub("E158-T01-D-2", "E158-T02-D",`A5 tumour ID`))%>%
  left_join(A5_clinical)%>%
  dplyr::select(A5_ID,Group, x, y, predicted_sex, Metastatic, TERT_or_ATRX,redundant, miR_chr14_group, 
                ALT, Location_simple, Location_summarised,`TERT_or_ATRX_for_design`, Is_HN)%>%
  mutate(Patient = gsub("-.*", "", A5_ID))%>%
  mutate(Patient_plot = replace(Patient,!(duplicated(Patient)| duplicated(Patient,fromLast = T)),NA))

# Plot some possible and expected sources of variation
p <- ggplot(data = mds_ggplot, aes(x = x, y = y,label= A5_ID))+
  geom_text()+
  blank_theme+
  coord_equal()+
  aes(colour = predicted_sex)+
  labs(colour = "Gender", x= "MDS dim 2", y = "MDS dim 3")+
  theme(aspect.ratio = 1)

# From this I will include location and gender as batch variables
Batch <- p + 
  aes(colour = Group)+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/MDS plots/Dimension 2 and 3 group.pdf", useDingbats = F, width =10, height = 8)
Location <- p + 
  aes(colour = Location_simple)+
  labs(colour = "Location")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/MDS plots/Dimension 2 and 3 location.pdf", useDingbats = F, width =10, height = 8)
Patient <- p + 
  aes(colour = Patient_plot)+
  labs(colour = "Patient")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/MDS plots/Dimension 2 and 3 patient.pdf", useDingbats = F, width =10, height = 8)
p + aes(colour = Location_summarised)
# TERT and ATRX form their own clustering as well
TERT_ATRX <- p + 
  aes(colour = TERT_or_ATRX)+
  labs(colour = "TERT or ATRX")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/MDS plots/Dimension 2 and 3 TERT or ATRX.pdf", useDingbats = F, width =10, height = 8)
p + aes(colour = TERT_or_ATRX == "TERT")
p + aes(colour = TERT_or_ATRX == "ATRX")
p + aes(colour = ALT)
Met <- p + 
  aes(colour = Metastatic)+
  labs(colour = "Metastatic")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/MDS plots/Dimension 2 and 3 metastatic.pdf", useDingbats = F, width =10, height = 8)
p +aes(colour = miR_chr14_group) 

plot_list <- list(Batch, Patient, Location, Met, TERT_ATRX)
pdf("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/MDS plots/Combined_mds.pdf",width = 30, height = 20)
ggarrange(plots = plot_list,newpage = F,ncol = 3)
dev.off()

# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
head(mVals[,1:5])

# Plot a UMAP as well
umap <- umap(t(mVals))

umap_ggplot <- data.frame(umap$layout)%>%
  rownames_to_column("A5 tumour ID")%>%
  mutate(Group = gsub("\\..*", "",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub(".*\\.", "",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub("_", "-",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub("E158-T01-D-1", "E158-T01-D",`A5 tumour ID`))%>%
  mutate(`A5 tumour ID` = gsub("E158-T01-D-2", "E158-T02-D",`A5 tumour ID`))%>%
  left_join(A5_clinical)

ggplot(data = umap_ggplot, aes(x = X1, y = X2, colour = Location_simple)) + 
  geom_point()+
#  geom_text_repel(color ="black",nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  theme(aspect.ratio=1)

ggplot(data = umap_ggplot, aes(x = X1, y = X2, colour = TERT_or_ATRX)) + 
  geom_point()+
  #  geom_text_repel(color ="black",nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  theme(aspect.ratio=1)

ggplot(data = umap_ggplot, aes(x = X1, y = X2, colour = predicted_sex)) + 
  geom_point()+
  #  geom_text_repel(color ="black",nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  theme(aspect.ratio=1)

# B vals for visualisation
bVals <- getBeta(mSetSqFlt)
head(bVals[,1:5])

hist(bVals)

# Save the beta values table 
bVals2 <- data.frame(bVals)%>%
  rownames_to_column("Probe")%>%
  write_csv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/A5 methylation Bvals.csv")


pdf("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/QC/Beta and M-values filtered.pdf", width =10)
par(mfrow=c(1,2))
densityPlot(bVals, sampGroups=targets$Sample_Group, main="Beta values", 
            legend=FALSE, xlab="Beta values")
densityPlot(mVals, sampGroups=targets$Sample_Group, main="M-values", 
            legend=FALSE, xlab="M values") 
dev.off()

# Design matrix with sex and location as confounders
# Set up a design matrix
design <- model.matrix(~0 + TERT_or_ATRX_for_design + predicted_sex + Is_HN, data = mds_ggplot)
colnames(design) <- gsub("TERT_or_ATRX_for_design| ","", colnames(design))

# Find the correlation between samples from the same patient
# Warnings are apparently ok for only a few genes:
# https://support.bioconductor.org/p/6618/
#corfit <- duplicateCorrelation(mVals, design, block = mds_ggplot$Patient)
# Save the RDS becuase this step takes forever
#saveRDS(corfit, "/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/corfit.RDS")
corfit <- readRDS("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/corfit.RDS")


# fit the linear model with blocking where the same sample is from multiple patients
fit <- lmFit(mVals, design, block = mds_ggplot$Patient, correlation =
               corfit$consensus)

contMatrix <- makeContrasts(ATRX_all = (ATRX_MET_Metastatic + ATRX_MET_Primary + ATRX_MET_Recurrent)/3 - Non_met_primary_MET_Primary,
                             TERT_all = (TERT_MET_Primary + TERT_MET_Metastatic)/2 - Non_met_primary_MET_Primary,
                             ATRX_met = ATRX_MET_Metastatic - Non_met_primary_MET_Primary,
                             TERT_met = TERT_MET_Metastatic - Non_met_primary_MET_Primary,
                             Unknown_driver_met = Unknown_met_driver_met_MET_Metastatic - Non_met_primary_MET_Primary,
                             Unknown_driver_all = (Unknown_met_driver_met_MET_Metastatic + Unknown_met_driver_primary_MET_Primary)/2 - Non_met_primary_MET_Primary,
                             levels= design)

fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summa.fit <- decideTests(fit2)

# Seems like there is a clear pattern in the TERT+ ATRX samples but not the others
summary(summa.fit)

# Get the results for the ATRX comparison
annEPICSub <- annEPIC[match(rownames(mVals),annEPIC$Name),
                      c(1:4,12:19,24:ncol(annEPIC))]

collections <- list.files("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Reference_data/MSIGdb/", full.names = T)

contrasts_to_plot <- c("ATRX_all","TERT_all", "ATRX_met", "TERT_met")

# Save each toptable and run GSEA analysis on each contrast
for(i in 1:length(contrasts_to_plot)){
  
  contrast <- contrasts_to_plot[i]
  output <- paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Outputs/DM toptables/", contrast,"_toptable.csv")
  toptable <- topTable(fit2, num=Inf, coef=contrast, genelist=annEPICSub, p.value = 0.05)%>%
    write_csv(output)
  
  # Do some gene set enrichment analysis
  
  # Get the table of results for the first contrast (ATRX - non met primary)
  DMPs <- topTable(fit2, num=Inf, coef=contrast, genelist=annEPICSub)
  
  # Get the significant CpG sites at less than 5% FDR
  sigCpGs <- DMPs$Name[DMPs$adj.P.Val<0.05]
  # Get all the CpG sites used in the analysis to form the background
  all <- DMPs$Name
  # Look for bias in number of CPGs to a gene being called as signifcant
  # and run GO analysis
  # This is an overrepresentation analysis
  
  # Run for each gene set in the msigDB
  
  # Do region by region DE now
  myAnnotation <- cpg.annotate(object = mVals, datatype = "array", what = "M", 
                               analysis.type = "differential", design = design, 
                               contrasts = TRUE, cont.matrix = contMatrix, 
                               coef = contrast, arraytype = "EPIC")
  
  # Run a function to calculate differentially methylated regions
  DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
  results.ranges <- extractRanges(DMRs)
  
  results_df <- data.frame(results.ranges)%>%
    write_csv(paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Outputs/DMR tables/", contrast,"_dmrcate_table.csv"))
  pal <- brewer.pal(8,"Dark2")
  # set up the grouping variables and colours
  groups <- pal[1:length(unique(mds_ggplot$TERT_or_ATRX))]
  names(groups) <- levels(factor(mds_ggplot$TERT_or_ATRX))
  cols <- groups[as.character(factor(mds_ggplot$TERT_or_ATRX))]
  
  bVals_plot <- bVals
  colnames(bVals_plot) <- gsub(".*\\.|_D", "", colnames(bVals_plot))
  
  # Plot dmrcate plots for the top 10 regions
  for(number in 1:30){
    
    save <- paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/DMR plots/", contrast, "_", number, ".pdf")
    pdf(save, height = 30, width =20)
    par(mfrow=c(1,1))
    
    skip_to_next <- FALSE
    
    # Note that print(b) fails since b doesn't exist
    
    tryCatch(DMR.plot(ranges = results.ranges, dmr = number, CpGs = bVals_plot, phen.col = cols, 
                          what = "Beta", arraytype = "EPIC", genome = "hg19"), 
                 error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next) { 
      dev.off()
      next 
      }

    
    dev.off()
    
}
 #GO term analysis of differentially methylated probes
 for(collection in collections){

   gene_set <- GSA.read.gmt(collection)

   gene_set_formatted <- gene_set$genesets

   names(gene_set_formatted) <- gene_set$geneset.names

   collection_name <- gsub(".entrez.gmt","", basename(collection))

   # Perform the gene set test
   print(contrast)
   print(collection)
   gsa <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=gene_set_formatted, array.type ="EPIC", plot.bias=F)

   top_gene_set <- topGSA(gsa, number=Inf)%>%
     rownames_to_column("Gene set")%>%
     #filter(FDR < 0.05)%>%
     mutate(collection = collection_name)%>%
     write_csv(paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Outputs/gometh_GO/", contrast,"_",collection_name , "_signf_go_terms.csv"))
 }
}

# plot the top 4 most significantly differentially methylated CpGs 
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals, cpg=cpg, pheno=mds_ggplot$TERT_or_ATRX, ylab = "Beta values")
})

# Compare methylation and RNA-Seq counts
# Get the beta values (% methylation) to plot against the CPMs
bvals_long <- bVals %>%
  data.frame()%>%
  rownames_to_column("Probe")%>%
  gather(Sample, B, -Probe)%>%
  mutate(Sample = gsub(".*\\.|_D|T0", "", Sample))

bvals_long[1:5,]

# Correlate the b values and the RNA-Seq log2 CPMs
lcpms <- read_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/RNA-Seq/Counts/RNA-Seq genewise log2 CPMs.csv")%>%
  dplyr::select(-ENSEMBL, -TXCHROM)%>%
  filter(!is.na(SYMBOL))%>%
  gather(Sample, lcpm, - SYMBOL)%>%
  mutate(Sample = gsub("\\.", "_", Sample))

# Get the methylation annotation
anno_df_meth <- data.frame(annEPIC)%>%
  rownames_to_column("Probe")%>%
  dplyr::select(Probe, chr, pos, SYMBOL = GencodeBasicV12_NAME, Regulatory_Feature_Group)%>%
  mutate(SYMBOL = gsub(";.*", "", SYMBOL))%>%
  filter(SYMBOL != "")

anno_df_meth

# Get an average B value per gene per sample
anno_df_meth_promo_B <- anno_df_meth %>%
  left_join(bvals_long)%>%
  filter(!is.na(Sample))%>%
  # Average M vals for a given gene
  group_by(Sample, SYMBOL)%>%
  summarise(Mean_B = mean(B, na.rm = T))%>%
  left_join(lcpms)%>%
  filter(!is.na(lcpm))%>%
  ungroup()

# Try correlating some key genes
correlate_genes <-function(anno_df_meth_promo_B, gene){
  
  gene_filt <- anno_df_meth_promo_B%>%
    filter(SYMBOL == gene)%>%
    mutate(A5_ID = gsub("_", "-", Sample))%>%
    left_join(A5_clinical)
  
  colours <- c("TERT" = "red","ATRX" =  "orange","Unknown_met_driver_met" = "dark red","Unknown_met_driver_primary" = "purple",
               "Short_follow_up_primary" = "grey","Non_met_primary" = "dark green")
  
  plt <- ggplot(data = gene_filt, aes(x = Mean_B, y = lcpm, colour = TERT_or_ATRX, label = A5_ID))+
    geom_point()+
    geom_text_repel(size =2)+
    blank_theme+
    labs(x = "Mean methylation B value", y = "Log2 CPM")+
    ggtitle(gene)+
    coord_equal()+
    theme(aspect.ratio = 1)+
    scale_colour_manual(values = colours)+
    ggsave(paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Gene expression vs methylation/",gene, ".pdf"), width = 15, height = 8, useDingbats = F)
  
  print(plt)
  
}

genes <- c("TERT", "CDKN2A")

for(gene in genes){
  correlate_genes(anno_df_meth_promo_B, gene)
}

TERT_probes <-anno_df_meth %>%
  filter(SYMBOL == "TERT")

# TERT specific analysis
# Get an average B value per gene per sample
anno_df_meth_promo_TERT <- TERT_probes%>%
  left_join(bvals_long)%>%
  filter(!is.na(Sample))%>%
  # Get only promoter region
  filter(Probe == "cg02545192")%>%
  left_join(lcpms)%>%
  filter(!is.na(lcpm))%>%
  ungroup()%>%
  mutate(A5_ID = gsub("_", "-", Sample))%>%
  left_join(A5_clinical)

colours <- c("TERT" = "red","ATRX" =  "orange","Unknown_met_driver_met" = "dark red","Unknown_met_driver_primary" = "purple",
             "Short_follow_up_primary" = "grey","Non_met_primary" = "green")

ggplot(data = anno_df_meth_promo_TERT, aes(y = B, x = lcpm, colour = TERT_or_ATRX, label = A5_ID))+
  geom_point()+
  geom_text_repel(size =2)+
  blank_theme+
  labs(y = "cg02545192 B value", x = "Log2 CPM")+
  ggtitle("TERT")+
  coord_equal()+
  theme(aspect.ratio = 1)+
  scale_colour_manual(values = colours)+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Gene expression vs methylation/TERT COMETE probe recapitulation.pdf")

  mvals_long <- mVals %>%
  data.frame()%>%
  rownames_to_column("Probe")%>%
  gather(Sample, M, -Probe)%>%
  mutate(Sample = gsub(".*\\.|_D|T0", "", Sample))

mvals_long[1:5,]

# Get just the promoter regions and average per gene per sample
# Then get Z score
anno_df_meth_promo <- anno_df_meth %>%
  filter(Regulatory_Feature_Group == "Promoter_Associated")%>%
  left_join(mvals_long)%>%
  filter(!is.na(Sample))%>%
  # Average M vals for a given gene
  group_by(Sample, SYMBOL)%>%
  summarise(Mean_M = mean(M, na.rm = T))%>%
  ungroup()%>%
  group_by(SYMBOL)%>%
  # Get a Z score of the mean M values
  mutate(z_score_promoter = (Mean_M - mean(Mean_M)) / sd(Mean_M))%>%
  filter(abs(z_score_promoter) >3)

anno_df_meth_promo[1:5,]

robz <- function(row){
  # Compute median and mean absolute deviation for row
  m <- median(row,na.rm = T)
  s <- mad(row,na.rm = T)
  
  # If the MAD is 0, set it to a very small number
  if(s == 0){
    s <- 1E-100
  }
  robzscore <- (row - m) / (s)
  return(robzscore)
}

# Get just the promoter regions and average per gene per sample
# Then get robust Z score
anno_df_meth_promo_robz <- anno_df_meth %>%
  filter(Regulatory_Feature_Group == "Promoter_Associated")%>%
  left_join(mvals_long)%>%
  filter(!is.na(Sample))%>%
  # Average M vals for a given gene
  group_by(Sample, SYMBOL)%>%
  summarise(Mean_M = mean(M, na.rm = T))%>%
  ungroup()%>%
  group_by(SYMBOL)%>%
  # Get a Z score of the mean M values
  mutate(rob_z_score_promoter = robz(Mean_M))%>%
  filter(abs(rob_z_score_promoter) >3)

# Look for things that are z score outliers by gene expression (>3) and have a change in methylation (>3)
z_scores_RNA_long <- read_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/RNA-Seq/Counts/RNA-Seq genewise log2 CPM Z-scores.csv")%>%
  dplyr::select(-ENSEMBL, -TXCHROM)%>%
  gather(Sample, Z_score, -SYMBOL)%>%
  filter(abs(Z_score) >3)%>%
  mutate(Sample = gsub("\\.", "_", Sample))%>%
  left_join(anno_df_meth_promo)%>%
  # Remove NA values
  filter(!is.na(Mean_M))%>%
  # Remove a sample with adrenocortical contamination
  filter(Sample != "E154_1")%>%
  dplyr::select(A5_ID = Sample,Gene =SYMBOL, `RNA-Seq Z score` = Z_score, `Mean array promoter methylation Z score` = z_score_promoter)%>%
  write_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Outputs/Z score analysis/RNA-Seq and Methylation Z score outliers.csv")

# Look for things that are robust z score outliers (>3) and have a change in methylation
rob_z_scores_RNA_long <- read_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/RNA-Seq/Counts/RNA-Seq genewise log2 CPM robust Z-scores.csv")%>%
  dplyr::select(-ENSEMBL, -TXCHROM)%>%
  gather(Sample, Robust_Z_score, -SYMBOL)%>%
  filter(abs(Robust_Z_score) >3)%>%
  mutate(Sample = gsub("\\.", "_", Sample))%>%
  left_join(anno_df_meth_promo_robz)%>%
  # Remove NA values
  filter(!is.na(Mean_M))%>%
  # Remove a sample with adrenocortical contamination
  filter(Sample != "E154_1")

# Most things that are expressed seem to be mostly unmethylated 
ggplot(data = z_scores_RNA_long, aes(x = `RNA-Seq Z score`, y = `Mean array promoter methylation Z score`))+
  geom_point()+
  labs(x = "RNA-Seq Z score", y = "Mean array promoter\nmethylation Z score")+
  blank_theme+
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_hline(yintercept = -3, linetype = "dashed")+
  geom_vline(xintercept = 3, linetype = "dashed")+
  geom_vline(xintercept = -3, linetype = "dashed")+
  coord_equal()+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Z score comparison/RNA-Seq vs methylation z scores.pdf")

# Most things that are expressed seem to be mostly unmethylated 
ggplot(data = rob_z_scores_RNA_long, aes(x = Robust_Z_score, y = rob_z_score_promoter))+
  geom_point()+
  labs(x = "Robust Z score RNA", y = "Mean Z robust score methylation")+
  blank_theme

# Get the contribution to outliers
contribution <- table(z_scores_RNA_long$Sample)%>%
  data.frame()

# Set up reusable tracks
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
#Ensure that the methylation data is ordered by chromosome and base position.
annEPICOrd <- annEPIC[order(annEPIC$chr,annEPIC$pos),]
bValsOrd <- bVals[match(annEPICOrd$Name,rownames(bVals)),]
# Create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(annEPICOrd$chr),
                   ranges=IRanges(start=annEPICOrd$pos, end=annEPICOrd$pos),
                   strand=Rle(rep("*",nrow(annEPICOrd))),
                   betas=bValsOrd)
islandHMM <- read.csv(paste0(dataDirectory,
                             "/model-based-cpg-islands-hg19-chr17.txt"),
                      sep="\t", stringsAsFactors=FALSE, header=FALSE)
islandData <- GRanges(seqnames=Rle(islandHMM[,1]), 
                      ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))
dnase <- read.csv(paste0(dataDirectory,"/wgEncodeRegDnaseClusteredV3chr17.bed"),
                  sep="\t",stringsAsFactors=FALSE,header=FALSE)
dnaseData <- GRanges(seqnames=dnase[,1],
                     ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                     strand=Rle(rep("*",nrow(dnase))),
                     data=dnase[,5])

gene_info <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)%>%
  data.frame()

plot_region <- function(gene, sample){
  
  # Get gene aliases
  gene <- select(org.Hs.eg.db,keys = gene, columns  = c("SYMBOL","ENTREZID"), keytype = "SYMBOL")
  
  location <- filter(gene_info, gene_id == gene$ENTREZID[1])
  
  chromosome <- location$seqnames
  
  if(location$strand == "-"){
    start <- location$end-500
    end <- location$end+5000
  }
  
  else{
    start <- location$start-5000
    end <- location$start+500
    
  }
  
  gen <- "hg19"
  # add 25% extra space to plot
  minbase <- start - (0.25*(end-start))
  maxbase <- end + (0.25*(end-start))
  
  df <- data.frame(chr=chromosome, start=minbase, end=maxbase,strand="*")
  ranges <- makeGRangesFromDataFrame(df) 
  
  iTrack <- IdeogramTrack(genome = gen, chromosome = chromosome)
  
  rTrack <- UcscTrack(genome=gen, chromosome=chromosome, track="NCBI RefSeq", 
                      from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                      rstarts="exonStarts", rends="exonEnds", gene="name", 
                      symbol="name2", transcript="name", strand="strand", 
                      fill="darkblue",stacking="squish", name="RefSeq", 
                      showId=TRUE, geneSymbol=TRUE)
  
  # extract data on CpGs in DMR
  cpgData_subset <- subsetByOverlaps(cpgData, ranges)
  
  groups <- left_join(targets, sample)
  
  # Set up the sample to plot
  groups <- factor(groups$TERT_mutant_or_expressed,levels = unique(groups$TERT_mutant_or_expressed))
  
  # Methylation data track
  methTrack <- DataTrack(range=cpgData_subset, genome = gen, groups = groups,
                         chromosome=chromosome, ylim=c(-0.05,1.05), col=pal,
                         type=c("a","p"), name="DNA Meth.\n(beta value)",
                         background.panel="white", legend=TRUE, cex.title=0.8,
                         cex.axis=0.8, cex.legend=0.8)
  
  # CpG island track
  islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                                 chromosome=chromosome,fill="darkgreen")
  
  # DNaseI hypersensitive site data track
  dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                          type="gradient", chromosome=chromosome)
  
  # DMR position data track
  dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                              chromosome=chromosome,fill="darkred")
  
  tracks <- list(iTrack, gTrack, methTrack,rTrack)
  sizes <- c(2,2,4,3) # set up the relative sizes of the tracks
  
  # Plot the tracks
  #pdf("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Gviz plots/Gviz_test.pdf")
  plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
             add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
  #dev.off()
}

# Try and plot some of the genes
sample <- A5_clinical%>%
  dplyr::select(A5_ID, TERT_or_ATRX,TERT_expression)%>%
  mutate(TERT_mutant_or_expressed = "No")%>%
  mutate(TERT_mutant_or_expressed = replace(TERT_mutant_or_expressed, TERT_or_ATRX == "TERT", "TERT mutant"))%>%
  mutate(TERT_mutant_or_expressed = paste0("Mutated: ", TERT_mutant_or_expressed,", ", TERT_expression))%>%
  dplyr::select(A5_ID, TERT_mutant_or_expressed)

gene <- "TERT"

plot_region(gene = gene, sample = sample)

# Plot RNA gene expression vs methylation

# Read in the RNA-Seq counts and join on B vals
RNAseq_CPM_long <- read_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/RNA-Seq/Counts/RNA-Seq genewise log2 CPMs.csv")%>%
  dplyr::select(-ENSEMBL, -TXCHROM)%>%
  gather(Sample, lcpm, -SYMBOL)%>%
  mutate(Sample = gsub("\\.", "_", Sample))%>%
  left_join(anno_df_meth_promo)

test <- filter(RNAseq_CPM_long, Sample == "E140_1")
# Most things that are expressed seem to be mostly unmethylated
ggplot(data = test, aes(x = lcpm, y = Mean_M))+
  geom_point()
cor.test(test$lcpm, test$Mean_B)