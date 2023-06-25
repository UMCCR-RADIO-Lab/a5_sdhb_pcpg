rm(list=ls())
library(Seurat)
library(Signac)
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)

#----
# Make a cells tsv for sinto filterbarcodes
#----

# I first converted the repeatmask annotation file to bed using awk

# make a cells file for sinto containing the bardoces that i want to keep 
atac <- readRDS("Data/snATAC-seq/temp/snATAC_AGG_01.RDS")
# make a vector with all the cell types that I want to keep - neoplastic and chromaffin cells 

cells.keep <- c(
  "E123_1",
  "E146_1", 
  "E156_1", 
  "E166_1",
  "E188_1",
  "E197_1",
  "E200_1",
  "E201_1",
  "Chromaffin_cell")

barcodes.table.all <- atac@meta.data %>%
  rownames_to_column(var = "barcode") %>% 
  dplyr::select(
    "barcode",
    "sequencing_library",
    "sampleID",
    "cell_type") %>% 
  dplyr::filter(
    cell_type %in% cells.keep
  )

# make barcode tables for sinto
for (i in 1:length(cells.keep)) {
  cell_type.use <- cells.keep[i]
  # if chromaffin cells, remove cells not in NAM018 sample
  if (cell_type.use == "Chromaffin_cell") {
    cells.table <- barcodes.table.all %>%
      dplyr::filter(
        cell_type == cell_type.use &
          sampleID == "NAM018") %>% 
      dplyr::select(barcode, cell_type) %>%
      mutate(barcode = sub(".$", "1", barcode))
  }
  # for other samples, filter out cells not in the tumour sample (cell type is the same as sample name )  
  else {
    cells.table <- barcodes.table.all %>%
      dplyr::filter(
        sampleID == cell_type.use &
          sampleID == cell_type.use) %>%
      dplyr::select(barcode, cell_type) %>%
      mutate(barcode = sub(".$", "1", barcode))
  }
  # remove column names and save as a TSV
  print(paste0(cell_type.use, ".table"))
  print(head(cells.table))
  assign(paste0(cell_type.use, ".table"), cells.table)
  file.path.use <- paste0("Data/repeat_analysis/cell_tables/", cell_type.use, "_cells.tsv")
  #write.table(x = cells.table, file = file.path.use, sep = '\t', row.names = F, col.names = F, quote = F)
}

#---- 
# Format the repeat annotation table 
#----

repeat_annotation.tbl <- read_tsv("Data/repeat_analysis/repeatmasker/awk.hg19.fa.out.bed", col_names = F)

# first I use this conversion table from GATK to rename the hg19 names to GRCH37 names
hg19_GRCh37_chr_names <- read_tsv("Data/repeat_analysis/repeatmasker/hg19_GRCh37_comparison_table.tsv")
hg19_GRCh37_chr_names <- hg19_GRCh37_chr_names %>%
  filter(!is.na(HG19_Contig) & !is.na(GRCh37_Contig))
hg19_to_GRCh37 <- setNames(hg19_GRCh37_chr_names$GRCh37_Contig, hg19_GRCh37_chr_names$HG19_Contig)
repeat_annotation.tbl <- repeat_annotation.tbl %>%
  mutate(X1 = recode(X1, !!!hg19_to_GRCh37)) %>%  # use the hg19_toGRCh37 named vector to rename the chromosome contigs
  filter(X1 != "chrM") # the mitochondrial chromosome is different in hg19 compared to GRCh37 so I'll remove repeat annotations located on chrM 

repeat_annotation.tbl$X1 %>% unique()
write_tsv(repeat_annotation.tbl, col_names = F, file = "Data/repeat_analysis/repeatmasker/awk.GRCh37.hg19.fa.out.bed")

#---- 
# Mosdepth-quantified repeats analysis 
#----

# read in the repeat coverage files (quantified using mosdepth)

E123_1.repeats.gr <- import.bed("Data/repeat_analysis/mosdepth/E123_1.regions.bed.gz")
E123_1.repeats.gr <- sort(E123_1.repeats.gr)
E146_1.repeats.gr <- import.bed("Data/repeat_analysis/mosdepth/E146_1.regions.bed.gz")
E146_1.repeats.gr <- sort(E146_1.repeats.gr)
E156_1.repeats.gr <- import.bed("Data/repeat_analysis/mosdepth/E156_1.regions.bed.gz")
E156_1.repeats.gr <- sort(E156_1.repeats.gr)
E166_1.repeats.gr <- import.bed("Data/repeat_analysis/mosdepth/E166_1.regions.bed.gz")
E166_1.repeats.gr <- sort(E166_1.repeats.gr)
E188_1.repeats.gr <- import.bed("Data/repeat_analysis/mosdepth/E188_1.regions.bed.gz")
E188_1.repeats.gr <- sort(E188_1.repeats.gr)
E197_1.repeats.gr <- import.bed("Data/repeat_analysis/mosdepth/E197_1.regions.bed.gz")
E197_1.repeats.gr <- sort(E197_1.repeats.gr)
E200_1.repeats.gr <- import.bed("Data/repeat_analysis/mosdepth/E200_1.regions.bed.gz")
E200_1.repeats.gr <- sort(E200_1.repeats.gr)
E201_1.repeats.gr <- import.bed("Data/repeat_analysis/mosdepth/E201_1.regions.bed.gz")
E201_1.repeats.gr <- sort(E201_1.repeats.gr)
Chromaffin_cell.repeats.gr <- import.bed("Data/repeat_analysis/mosdepth/Chromaffin_cell.regions.bed.gz")
Chromaffin_cell.repeats.gr <- sort(Chromaffin_cell.repeats.gr)

#  combine into a list
repeats.gr.list <- list(
  E123_1 = E123_1.repeats.gr,
  E146_1 = E146_1.repeats.gr,
  E156_1 = E156_1.repeats.gr,
  E166_1 = E166_1.repeats.gr,
  E188_1 = E188_1.repeats.gr,
  E197_1 = E197_1.repeats.gr,
  E200_1 = E200_1.repeats.gr,
  E201_1 = E201_1.repeats.gr,
  Chromaffin_cell = Chromaffin_cell.repeats.gr)

# get the repeat type from the annotation bed file 
extra_cols <-c(repeat_type = "character") 
repeat_annotation.gr <- import.bed(
  "Data/repeat_analysis/repeatmasker/awk.hg19.fa.out.bed",
  extraCols = extra_cols)
# the unmapped chromosomes are different - might need to realign to a reference with the same names for unmapped contigs
# these were removed by mosdepth 
# will need to rename these, and rerun mosdepth but for now i'll just remove the extras
repeat_annotation.gr <- keepStandardChromosomes(repeat_annotation.gr, pruning.mode = "coarse")
# make the chromosome levels the same order
seqlevels(repeat_annotation.gr) <- seqlevels(E123_1.repeats.gr)
repeat_annotation.gr <- sort(repeat_annotation.gr, ignore.strand = T)

# check that the names and ranges match 
for (i in 1:length(repeats.gr.list)) {
  gr <- repeats.gr.list[[i]]
  matching_names <- table(repeat_annotation.gr$name == gr$name)
  print(paste("names match:", matching_names))
}

counts_list <- list()
# get the name for each repeat
counts_list[[1]] <-  as_tibble(repeat_annotation.gr$name)
names(counts_list[[1]]) <- "repeat_name"
#get the type for each repeat
counts_list[[2]] <- as_tibble(repeat_annotation.gr$repeat_type)
names(counts_list[[2]]) <- "repeat_type"
# get the range of the repeat 
counts_list[[3]] <- as_tibble(GRangesToString(repeat_annotation.gr))
names(counts_list[[3]]) <- "repeat_range"

# extract the counts for each region and put them in a dataframe  
for (i in 1:length(repeats.gr.list)) {
  sample.name <- names(repeats.gr.list)[i]
  gr <- repeats.gr.list[[i]]
  repeat_counts <- gr$score
  repeat_counts <- as_tibble(repeat_counts)
  counts_list[[3+i]] <- repeat_counts
  names(counts_list[[3+i]]) <- c(paste(sample.name))
}
repeat_counts <- bind_cols(counts_list)

# sum the counts for each repeat and type
counts_per_repeat <- repeat_counts %>% 
  group_by(repeat_name) %>% 
  summarise(across(c(E123_1, E146_1, E156_1, E166_1, E188_1, E197_1, E200_1, E201_1, Chromaffin_cell), sum))
counts_per_repeat_type <- repeat_counts %>% 
  group_by(repeat_type) %>% 
  summarise(across(c(E123_1, E146_1, E156_1, E166_1, E188_1, E197_1, E200_1, E201_1, Chromaffin_cell), sum))

# write_csv(counts_per_repeat, "Data/repeat_analysis/counts/repeat_counts.csv")
# write_csv(counts_per_repeat_type, "Data/repeat_analysis/counts/repeat_type_counts.csv")
# write_csv(repeat_counts, "Data/repeat_analysis/counts/repeat_type_regions.csv")

x <- letters[c(1,2,2,1)]
y <- c(a = "Apple", b = "Banana")
dplyr::recode(x, !!!y)
