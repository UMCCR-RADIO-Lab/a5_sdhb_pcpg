# Generate pseudobulk atac samples 
rm(list=ls())
library(Seurat)
library(Signac)
library(tidyverse)


set.seed(1234)
setwd("/data/cephfs/punim0648/Projects/Blake/A5_R_project/")

atac <-  readRDS("Data/snATAC-seq/processed/snATAC_AGG_01_processed.RDS")
DefaultAssay(atac) <- "MACS2"

#----
# Create Pseudobulk dataset
#----

# Make data frame of annotation for each cell type, filter out unwanted barcodes
ann_df <- data.frame(Clusters = atac$cell_type, Sample = atac$sampleID, Tumour = atac$predicted.id, stringsAsFactors = F)%>%
  rownames_to_column("Barcode")%>%
  # calculate the percentage of each sample within each celltype cluster 
  group_by(Clusters)%>%
  mutate(Count_per_cluster = n())%>%
  ungroup()%>%
  group_by(Clusters, Sample)%>%
  mutate(Count_sample_per_cluster = n())%>%
  mutate(Percent_sample = Count_sample_per_cluster/Count_per_cluster*100)%>%
  # We are only interested in the tumour and normal chromaffin clusters for now
  filter(Tumour %in% c("Tumour", "Chromaffin_cell"))%>%
  ungroup()

# store sample names
samples <- unique(ann_df$Sample)

# List to add each cluster df
sample_list <- list()
# Loop to generate a pseudo bulk column for each celltype cluster
for (i in 1:length(samples)){
  # The name of the cluster
  sample_name <- samples[[i]]
  # Filter the annotation to keep only clusters where the sample of interest is the most frequent sample of the cluster
  anno_sample <- filter(ann_df, Sample == sample_name)%>%
    mutate(Percent_sample_max = max(Percent_sample))%>%
    dplyr::filter(Percent_sample == Percent_sample_max)
  # Subset the seurat object to get the cluster
  #cluster <- subset(x = atac, subset = orig.ident == sample_name)
  cluster <- atac[,anno_sample$Barcode]
  # Get the counts from the cluster
  counts <- as.matrix(cluster@assays$MACS2@counts)
  # Sanity check
  sane <- sum(colnames(counts) == anno_sample$Barcode) == nrow(anno_sample)
  print(paste("Sanity check passed",sane))
  counts_cluster <- rowSums(counts)  
  df <- data.frame(counts_cluster)
  colnames(df) <- sample_name
  # Add the counts to a list
  sample_list [[i]] <-  df
}

PB <- bind_cols(sample_list)
head(PB)

# write.csv(PB, "Data/snATAC-seq/processed/atac_pseudobulk.csv")

