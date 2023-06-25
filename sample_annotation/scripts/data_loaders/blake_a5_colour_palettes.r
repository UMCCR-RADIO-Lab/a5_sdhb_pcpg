# colour palette scripts for the A5 single cell 

library(tidyverse)
library(ggsci)
library(scales)
library(RColorBrewer)
library(ggthemes)

# ----
# 1. sample colour palette:
# ----

all_samples <- c("E019-1",
                 "E225-1",
                 "E197-1",
                 "E171-1",
                 "E166-1",
                 "E201-1",
                 "E188-1",
                 "E200-1",
                 "E123-1",
                 "E156-1",
                 "E140-1",
                 "E146-1",
                 "E143-1",
                 "NAM018",
                 "E240-1",
                 "E243-1",
                 "NPG103",
                 "P018-PGL1",
                 "P018-PGL3")

#show_col(c(calc_pal()(12), "#A65628" , "#F781BF", "#999999"))

sample_colours <- setNames(c(calc_pal()(12), "#A65628" , "#F781BF", "#999999", "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF"), all_samples)

# ----
# 2. Cell type colour palette
# ----

all_cell_types <- c("Tumor",
                    "Endothelial cells" ,
                    "Chromaffin cells" ,
                    "Myeloid cells",
                    "Lymphocytes",
                    "Adrenocortical cells",
                    "Fibroblasts",
                    "SCLCs")

cell_type_colours <- setNames(brewer.pal(name = "Set1", n = 8), all_cell_types)

# # This script sets up consistent colour scheme for the pheo groups/cell types
# all_genotypes = c("ATRX",
#                   "TERT",
#                   "Non-met primary",
#                   "other/unknown"#,
#                   #"Normal"
# )
# # all.cell.types = c("Tumour", "Chromaffin cells", "Adrenocortical cells", "Endothelial cells", "Fibroblasts", "Sustentacular cells", "Macrophages/Monocytes", "Lymphocytes")
# all_cell_types = c("Tumour",
#                    "Chromaffin cells",
#                    "Sustentacular cells",
#                    "Endothelial cells",
#                    "Fibroblasts",
#                    "Myeloid cells",
#                    "Lymphocytes",
#                    "Adrenocortical cells")
# 
# # colours from pal_d3
# genotype_cols <- setNames(c(pal_d3()(3), "lightgrey"),
#                           nm=all_genotypes)
# gender_cols <- setNames(pal_d3("category20")(10)[9:10], nm = c("male", "female"))
# ALT_cols <- setNames(pal_d3()(10)[5:6], c("ALT", "No ALT"))

# old code from snRNA_analysis ----

# ----
# setup nice colour palettes 
# ----

# TODO: make a colour palette for the cell types and genotypes
# genotypes_and_celltypes <- c()

# palettes for genotypes and cell types 
pal_genotypes_cell_types <- c(pal_d3("category20")(20), "lightgrey") 

# store all the samples 
rna_samples <- c("E019",
                 "E123-1",
                 "E140-1",
                 "E143-1",
                 "E146-1",
                 "E156-1",
                 "E166-1",
                 "E171-1",
                 "E197-1",
                 "E225-1",
                 "NAM021",
                 "NAM025",
                 "NPG103")

pal_samples <- setNames(c(calc_pal()(12), "#7F7F7FFF"), rna_samples)


