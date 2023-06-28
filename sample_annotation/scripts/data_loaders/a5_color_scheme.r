library(RColorBrewer)
library(ggplot2)

source("/g/data/pq08/projects/ppgl/public_data/annotation/zethoven_sn_rna_cluster_naming.r")

if (!exists("differential_group_levels"))
{
  source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_add_differential_groups.r")
}

ColorPalette <- c(LightBlue1="#cadae8ff",LightBlue2="#abc4e2ff",LightBlue3="#8aa9d6ff",
                  DarkBlue1="#7884baff",DarkBlue2="#606d9fff",DarkBlue3="#4c5a88ff",
                  Purple1="#765b92ff",Purple2="#9370a5ff",Purple3="#b280a9ff",Purple4="#cc85b1ff",
                  LightRed1="#ee7474ff",LightRed2="#e2696aff",LightRed3="#de5353ff",
                  DarkRed1="#c25858ff",DarkRed2="#c04241ff",DarkOrange1="#c65a44ff",
                  DarkOrange2="#eb7b54ff",LightOrange1="#f18d73ff",LightOrange2="#f5a697ff",
                  LightBrown1="#f5b9b0ff",LightBrown2="#d6a097ff",DarkBrown1="#c17963ff",DarkBrown2="#96665aff",
                  DarkGreen1="#637b62ff",DarkGreen2="#6ea668ff",LightGreen1="#98ca8cff",LightGreen2="#e2eab5ff",
                  Yellow1="#fbf2adff",Yellow2="#fef8c6ff",Yellow3="#f4e764ff",Salmon="#f2c5a7ff",LightSalmon="#f2dec4ff",
                  LemonGrey1="#f7f0e2ff",LemonWhite1="#fefbe5ff",LemonWhite2="#f7f4dcff",LemonGrey2="#efecdcff",
                  LightGrey1="#e2e0d9ff",LightGrey2="#d6d6d4ff",DarkGrey1="#b4b4b4ff",DarkGrey2="#9a9999ff")

ColorPalette.show <- function() {ggplot(data.frame(Color=factor(names(ColorPalette), levels=names(ColorPalette)),Hex=ColorPalette),
                                        mapping = aes(x=Color,y=1,fill=Hex)) +
    geom_tile() +
    scale_fill_identity() +
    theme_void() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
# This script sets up consistent colour scheme for the pheo groups/cell types

all.genotypes = c("HRAS", "NF1","EPAS1","SDHD","MAX","FH","RET","Normal","SDHB","MAML3","TMEM127","H3F3A","VHL","Unknown", "SDHA","SDHC","EGLN1","SETD2","KIF1B","MAML","NGFR","BRAF","CSDE1","IDH1")
all.cell.types = c("Tumour","Chromaffin cells",  "Adrenocortical cells", "Endothelial cells", "Fibroblasts", "Sustentacular cells", "Myeloid cells", "T/NK cells", "B cells", "Mast cells")

# Stolen from https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
c25 = c("dodgerblue2","#E31A1C", # red
        "green4",
        "#6A3D9A", # purple
        "#FF7F00", # orange
        "gold1",
        "skyblue2","#FB9A99", # lt pink
        "palegreen2",
        "#CAB2D6", # lt purple
        "#FDBF6F", # lt orange
        "gray70", "khaki2",
        "maroon","orchid1","deeppink1","blue1","steelblue4",
        "darkturquoise","green1","yellow4","yellow3",
        "darkorange4","brown",
        "black")


genotype_cols = setNames(c25[-25], all.genotypes)
genotype_cols["Unknown"] = alpha('grey50', 0.25) # fade out missing genotypes
genotype_cols = genotype_cols[c(setdiff(names(genotype_cols), "Unknown"), "Unknown")]
genotype_cols["VHL"] = "darkorchid1" # remove bright yellow for white backgrounds
genotype_cols["H3F3A"] = "grey10"


cell_cols = setNames(brewer.pal(n = length(all.cell.types), name = "Set3")[c(9,10, 1:8)], all.cell.types) # reorder to make Tumour grey
cell_cols[c("Endothelial cells")] = "yellow3" # replace bright yellow for white backgrounds
cell_cols[make.names(names(cell_cols))] = unname(cell_cols) # include fixed R column names

subtype_cols <- c("yellow3", "#A65628" , "#FF7F00", "#984EA3", "#E41A1C", 
                  "#377EB8","#4DAF4A",  "light blue", "#FB9A99", "black")
subtype_cols <- setNames(rep(subtype_cols,5), unlist(zethoven_cluster_naming))
subtype_cols <- subtype_cols[!duplicated(names(subtype_cols))]
subtype_cols[["C2B2 (MAML3)"]] <- subtype_cols[["C2B2 (MAML)"]]


differential_group_colors=c("green","darkgreen","burlywood1", "burlywood4",  "red", "coral","coral", "coral3",
                            "lightblue", "darkblue", "grey", "grey", "lightgreen"," lightgreen", "green", "purple","palevioletred2")
names(differential_group_colors) <- differential_group_levels

# ggplot scales

genotype_color_scale = scale_color_manual(values=genotype_cols, name="Genotype")
genotype_fill_scale = scale_fill_manual(values=genotype_cols, name="Genotype")

cell_color_scale = scale_color_manual(values=cell_cols, name="Cell Type")
cell_fill_scale2 = scale_fill_manual(values=cell_cols, name="Cell Type")

subtype_color_scale = scale_color_manual(values=subtype_cols, name="")
subtype_fill_scale = scale_color_manual(values=subtype_cols, name="")

continuous_colours = c("blue", "white", "yellow")
cont_color_scale = scale_color_gradientn(colours = continuous_colours)
cont_fill_scale2 = scale_fill_gradientn(colours = continuous_colours)

specimen_type_cols = setNames(c(ColorPalette[["DarkGreen2"]], ColorPalette[["LightOrange1"]],
                              ColorPalette[["LightRed2"]], ColorPalette[["LightGrey2"]]),
                              c("Non-malignant Primary", "Malignant primary", 
                                "Metastasis", "Local Recurrance"))
specimen_type_cols[["MetastaticPrimary"]] <- specimen_type_cols[["Malignant primary"]]
specimen_type_cols[["NonMetastaticPrimary"]] <- specimen_type_cols[["Non-malignant Primary"]]
specimen_type_cols[["Other"]] <- "#999999"

location_cols = c("Adrenal"=ColorPalette[["LightGreen1"]], 
                  "Extraadrenal"=ColorPalette[["DarkBlue1"]], 
                  "Head_neck"=ColorPalette[["LightOrange1"]], 
                  "Extraadrenal_aortic"=ColorPalette[["DarkRed1"]], 
                  "Metastasis"=ColorPalette[["Purple3"]], 
                  "Unspecified"=ColorPalette[["DarkGrey1"]])
location_cols[["Head and neck"]] <- location_cols[["Head_neck"]] 
location_cols[["Extraadrenal_bladder"]] = "#99cccc"

driver_cols <- c(ATRX=ColorPalette[["DarkBlue1"]], 
                TERT=ColorPalette[["LightOrange1"]], 
                WT=ColorPalette[["LightGreen1"]])

all_cell_types <- c("Tumor",
                    "Endothelial cells" ,
                    "Chromaffin cells" ,
                    "Myeloid cells",
                    "Lymphocytes",
                    "Adrenocortical cells",
                    "Fibroblasts",
                    "SCLCs")

cell_type_colours <- setNames(brewer.pal(name = "Set1", n = 8), all_cell_types)
