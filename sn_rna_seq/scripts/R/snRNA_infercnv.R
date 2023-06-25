# inferCNV for A5
rm(list=ls())

library(infercnv)
library(Seurat)
library(tidyverse)
library(future)

plan("multicore", workers = 8) # change to 'multisession' if running in Rstudio
options(future.globals.maxSize = 200 * 1000 * 1024^2) # 200gb max RAM

# ----
# infercnv inputs
# ----

# read in seurat object
rna <- readRDS("Data/snRNA-seq/processed/snRNA_cs_regressed_annotated.RDS")

# genes with start and end genomic coordinates
# made this list using the "gtf_to_position_file.py" script supplied with infercnv
gene.positions.file <- "Data/snRNA-seq/infercnv/hg38_gene_positions.tsv"
gene.positions <- read.delim(gene.positions.file, sep = "\t")

length(setdiff(rownames(rna), gene.positions[,1])) # only 10 genes are different, will just remove them 

# subset to only have the genes that are in the sn-RNA-seq
gene.positions <- gene.positions[gene.positions[,1] %in% rownames(rna), ]
table(gene.positions[,1] %in% rownames(rna))

# ----
# run infercnv 
# ----

# run infercnv on each sample, using the normal cells from all samples as a reference
all.samples <- unique(rna$A5.ID)
normal.cell.types <- unique(rna$cell_type)[!unique(rna$cell_type) %in% c("Tumor", "SCLCs")]

for (i in seq_along(all.samples)) {
  sample <- all.samples[i]
  print(sample)
  # subset the seurat object to contain sample tumor and sustentactular cells from a sample
  # plus all of the normal cells
  rna_sample <- subset(rna, subset = (A5.ID == sample &
                                        cell_type %in% c("Tumor", "SCLCs")) |
                         (cell_type %in% normal.cell.types))
  mat <- as.matrix(rna_sample@assays$RNA@counts)

  mat <- mat[gene.positions[,1],]
  
  print(mat[1:5,1:5])
  print(dim(mat))
  print(table(rownames(mat) == gene.positions[,1]))
  annot <- rna_sample@meta.data[,"cell_type", drop=F]
  rm(rna_sample)
  gc()
  # get the raw counts matrix, match the order to the gene.positions file 
  # mat = as.matrix(x@assays$RNA@counts[match(gene.positions[,1], rownames(rna_sample@assays$RNA@counts)),])
  infercnv_obj = CreateInfercnvObject(raw_counts_matrix=mat,
                                      annotations_file=annot,
                                      delim="\t",
                                      gene_order_file=gene.positions.file,
                                      ref_group_names=normal.cell.types)
  
  outdir = file.path("results/infercnv", sample)
  if(!dir.exists(outdir)){dir.create(outdir)}
  infercnv_obj = infercnv::run(infercnv_obj,
                               cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                               out_dir=outdir,
                               cluster_by_groups=FALSE, 
                               denoise=TRUE,
                               HMM=TRUE,
                               analysis_mode = "subclusters",
                               tumor_subcluster_partition_method = "qnorm",
                               tumor_subcluster_pval = 0.01,
                               num_threads=8,
                               no_plot=TRUE,
                               no_prelim_plot=TRUE
                               )
  saveRDS(infercnv_obj, file=file.path(outdir, "output.RDS"))
  print(file.path(outdir, "output.RDS"))
}





# ERROR:

# Error in UseMethod("as.phylo") : 
#   no applicable method for 'as.phylo' applied to an object of class "NULL"
# In addition: Warning message:
#   In max(nchar(obs_annotations_names)) :
#   no non-missing arguments to max; returning -Inf