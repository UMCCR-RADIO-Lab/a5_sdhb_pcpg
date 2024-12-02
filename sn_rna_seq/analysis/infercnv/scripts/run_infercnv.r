# inferCNV for A5

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
n_threads = as.numeric(args[2])

library(infercnv)
library(Seurat)
library(tidyverse)
library(future)

setwd("/g/data/pq08/projects/ppgl/")

plan("multicore", workers = n_threads) # change to 'multisession' if running in Rstudio
options(future.globals.maxSize = 200 * 1000 * 1024^2) # 200gb max RAM

# ----
# infercnv inputs
# ----

# read in seurat object
source("./a5/sn_rna_seq/scripts/data_loaders/a5_snrna_dataloader.r")
data_loader_a5_snrna()
a5_snrna <- snrna_annotate_cell_types(a5_snrna)

# genes with start and end genomic coordinates
# made this list using the "gtf_to_position_file.py" script supplied with infercnv
gene.positions.file <- "./a5/sn_rna_seq/analysis/infercnv/gencode.v36.gene_positions.tsv"
gene.positions <- read.delim(gene.positions.file, sep = "\t", header = F, col.names = c("gene_name", "chr", "pos_start", "pos_end"))

length(setdiff(rownames(a5_snrna), gene.positions$gene_name)) # 710 genes are different, will just remove them 

# subset to only have the genes that are in the sn-RNA-seq
gene.positions <- gene.positions[gene.positions$gene_name %in% rownames(a5_snrna), ]
#gene.positions <- gene.positions %>% filter(chr %in% paste0("chr", c(1:22,"X","Y")))
#write_delim(x = gene.positions, file="./a5/sn_rna_seq/analysis/infercnv/gencode.v36.gene_positions.tsv", delim = "\t", col_names = F)
table(gene.positions$gene_name %in% rownames(a5_snrna))

# ----
# run infercnv 
# ----

# run infercnv on each sample, using the normal cells from all samples as a reference
normal.cell.types <- unique(a5_snrna$cell_type)[!unique(a5_snrna$cell_type) %in% c("Tumor", "SCLCs")]
  
# subset the seurat object to contain sample tumor and sustentactular cells from a sample
# plus all of the normal cells
rna_sample <- subset(a5_snrna, subset = (A5_ID == sample & cell_type %in% c("Tumor", "SCLCs")) |
                        (cell_type %in% normal.cell.types))
mat <- as.matrix(rna_sample@assays$RNA@counts)

mat <- mat[gene.positions$gene_name,]

#print(mat[1:5,1:5])
#print(dim(mat))
#print(table(rownames(mat) == gene.positions$gene_name))
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

outdir = file.path("./a5/sn_rna_seq/analysis/infercnv/output", sample)
if(!dir.exists(outdir)){dir.create(outdir,recursive = T)}
infercnv_obj = infercnv::run(infercnv_obj,
                              cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                              out_dir=outdir,
                              cluster_by_groups=FALSE, 
                              denoise=TRUE,
                              HMM=TRUE,
                              analysis_mode = "subclusters",
                              tumor_subcluster_partition_method = "qnorm",
                              tumor_subcluster_pval = 0.01,
                              num_threads=n_threads,
                              no_plot=FALSE,
                              no_prelim_plot=FALSE
                              )
saveRDS(infercnv_obj, file=file.path(outdir, "output.RDS"))
print(file.path(outdir, "output.RDS"))