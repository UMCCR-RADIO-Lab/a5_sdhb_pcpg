# Comparison of fetal adrenal and normal adrenal cell type annotation
# perform label transfer to annotate normal adrenal medulla and normal paraganglia cells
# using a previously published reference dataset

rm(list=ls())

library(Seurat)
library(tidyverse)
library(ggsci)

# samples F2, F7, F106 AND F107 are human fetal adrenal samples 

# read in the gene-cell matrix
F2.table <- read.delim("Data/published_datasets/Dong_2020/GSM4088785_F2_gene_cell_exprs_table.xls.gz", sep = "\t", )
# remove second instance of the non-unique gene Symbols
F2.table <- F2.table[!duplicated(F2.table$Symbol),]
# remove the ensemble gene IDs 
F2.data <- F2.table[,2:ncol(F2.table)]
# convert to matrix format
F2.data <- F2.data %>%
  remove_rownames %>%
  column_to_rownames(var = "Symbol")
F2.data <- as.matrix(F2.data)
# create seurat object
F2_rna <- CreateSeuratObject(counts = F2.data, project = "F2")

# read in the gene-cell matrix
F7.table <- read.delim("Data/published_datasets/Dong_2020/GSM4088786_F7_gene_cell_exprs_table.xls.gz", sep = "\t", )
# remove second instance of the non-unique gene Symbols
F7.table <- F7.table[!duplicated(F7.table$Symbol),]
# remove the ensemble gene IDs 
F7.data <- F7.table[,2:ncol(F7.table)]
# convert to matrix format
F7.data <- F7.data %>%
  remove_rownames %>%
  column_to_rownames(var = "Symbol")
F7.data <- as.matrix(F7.data)
# create seurat object
F7_rna <- CreateSeuratObject(counts = F7.data, project = "F7")

# read in the gene-cell matrix
F106.table <- read.delim("Data/published_datasets/Dong_2020/GSM4088787_F106_gene_cell_exprs_table.xls.gz", sep = "\t", )
# remove second instance of the non-unique gene Symbols
F106.table <- F106.table[!duplicated(F106.table$Symbol),]
# remove the ensemble gene IDs 
F106.data <- F106.table[,2:ncol(F106.table)]
# convert to matrix format
F106.data <- F106.data %>%
  remove_rownames %>%
  column_to_rownames(var = "Symbol")
F106.data <- as.matrix(F106.data)
# create seurat object
F106_rna <- CreateSeuratObject(counts = F106.data, project = "F106")

# read in the gene-cell matrix
F107.table <- read.delim("Data/published_datasets/Dong_2020/GSM4088788_F107_gene_cell_exprs_table.xls.gz", sep = "\t", )
# remove second instance of the non-unique gene Symbols
F107.table <- F107.table[!duplicated(F107.table$Symbol),]
# remove the ensemble gene IDs 
F107.data <- F107.table[,2:ncol(F107.table)]
# convert to matrix format
F107.data <- F107.data %>%
  remove_rownames %>%
  column_to_rownames(var = "Symbol")
F107.data <- as.matrix(F107.data)
# create seurat object
F107_rna <- CreateSeuratObject(counts = F107.data, project = "F107")

# merge all of the samples
fetal_rna <- merge(
  F2_rna, y = c(F7_rna, F106_rna, F107_rna),
  add.cell.ids = c("F2","F7","F106","F107"),
  project = "fetal_adrenal_gland")

# remove the (poor-quality), non-annotated cells that are not in the annotation csv
adrenal_annotation <- read.csv("Data/published_datasets/Dong_2020/GSE137804_Adrenal_gland_annotation.csv")
# reformat the cell barcodes to match mine
adrenal_annotation <- adrenal_annotation %>%
  mutate(barcode = substr(cell_id,1,18)) %>% 
  unite(barcode, sample, barcode)
cells.use <- adrenal_annotation$barcode

fetal_rna <- subset(fetal_rna, cells = cells.use)
# the number of cells for each sample matches the number that passed QC for each sample
table(fetal_rna@meta.data$orig.ident)

# add the cell type metadata
cell.type <- adrenal_annotation$annotation
names(cell.type) <- adrenal_annotation$barcode
fetal_rna <- AddMetaData(fetal_rna, cell.type, col.name = "cell.type")
#saveRDS(fetal_rna, "Data/published_datasets/Dong_2020/fetal_adrenal_scRNA.RDS")

# read in the processed scRNA-seq data 
a5_rna <- readRDS("Data/snRNA-seq/processed/snRNA_cs_regressed_annotated_doublets_removed.RDS")
# subset the normal adrenal cells
nam_rna <- subset(a5_rna, subset = orig.ident %in% c("NAM021","NAM025", "NPG103"))
rm(a5_rna) # free up memory 
# saveRDS(nam_rna, "Data/snRNA-seq/processed/snRNA_nam_cs_regressed_annotated.RDS")

fetal_rna <- readRDS("Data/published_datasets/Dong_2020/fetal_adrenal_scRNA.RDS")

# Run the seurat pipeline on the fetal dataset 
fetal_rna <- NormalizeData(fetal_rna)
fetal_rna <- FindVariableFeatures(fetal_rna, selection.method = "vst", nfeatures = 2000)
fetal_rna <- ScaleData(fetal_rna)
fetal_rna <- RunPCA(fetal_rna)
fetal_rna <- RunUMAP(fetal_rna, dims = 1:20)
DimPlot(fetal_rna, group.by = "cell.type")

# Run the seurat pipeline on the nam dataset 
nam_rna <- NormalizeData(nam_rna)
nam_rna <- FindVariableFeatures(nam_rna, selection.method = "vst", nfeatures = 2000)
nam_rna <- ScaleData(nam_rna)
nam_rna <- RunPCA(nam_rna)
nam_rna <- RunUMAP(nam_rna, dims = 1:20)
DimPlot(nam_rna, group.by = "cell_type")

# Use label transfer to classify the NAM data 
anchors <- FindTransferAnchors(
  reference = fetal_rna,
  query = nam_rna)

predictions <- TransferData(
  anchorset = anchors,
  refdata = fetal_rna$cell.type)

nam_rna <- AddMetaData(nam_rna, metadata = predictions)

colour.pal <- pal_d3()(10)
colour.pal2 <- pal_igv()(10)

DimPlot(nam_rna, group.by = "predicted.id", cols = colour.pal) +
  DimPlot(nam_rna, group.by = "cell_type", cols = colour.pal2) +
  FeaturePlot(nam_rna, features = "prediction.score.max")

table(nam_rna$cell_type)

hist(nam_rna$prediction.score.max)

table(fetal_rna$cell.type)

DimPlot(nam_rna, group.by = "orig.ident", cols = colour.pal) 

x <- DotPlot(nam_rna, features = c("CARTPT", "PHOX2A", "PHOX2B", "NPY","SOX10", "CDH19"), group.by = 'predicted.id') + ggtitle("NAM")

VlnPlot(fetal_rna, features = c("CARTPT", "PHOX2A", "PHOX2B", "SOX10", "CDH19"), group.by = "cell.type")

FeaturePlot(nam_rna, features = c("CARTPT", "PHOX2B", "PNMT", "TH"))
FeaturePlot(fetal_rna, features = c("CARTPT", "PHOX2B", "PNMT", "TH"))


y <- DotPlot(fetal_rna, features = c("CARTPT", "PHOX2A", "PHOX2B", "NPY","SOX10", "CDH19"), group.by = 'cell.type') + ggtitle("fetal")
x + y
