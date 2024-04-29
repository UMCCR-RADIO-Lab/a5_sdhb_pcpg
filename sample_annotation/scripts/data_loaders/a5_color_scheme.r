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


differential_group_colors=
  c("NMP_WT_VAFOK_PurityLow" = "darkolivegreen",
    "NMP_WT_VAFOK_PurityOK" = "darkgreen",
    "SFU_ATRX_VAFOK_PurityOK" = "burlywood1",
    "SFU_TERT_VAFOK_PurityLow" = "burlywood4",
    "SFU_WT_VAFOK_PurityOK" = "yellow",
    "MetatastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityLow" = "lightblue",
    "MetatastaticPrimaryUnconfirmed_ATRX_VAFOK_PurityOK" = "lightblue",
    "MetatastaticPrimaryUnconfirmed_TERT_VAFLow_PurityLow" = "coral",
    "MetatastaticPrimaryUnconfirmed_TERT_VAFLow_PurityOK" = "coral3",
    "MetatastaticPrimaryUnconfirmed_TERT_VAFOK_PurityOK" = "coral3",
    "MetatastaticPrimaryUnconfirmed_WT_VAFOK_PurityLow" = "purple",
    "MetatastaticPrimaryUnconfirmed_WT_VAFOK_PurityOK" = "purple",
    "MetatastaticPrimaryConfirmed_ATRX_VAFOK_PurityOK" = "blue",
    "MetatastaticPrimaryConfirmed_TERT_VAFLow_PurityOK" = "chocolate3",
    "MetatastaticPrimaryConfirmed_TERT_VAFOK_PurityOK" = "chocolate3",
    "MetatastaticPrimaryConfirmed_WT_VAFOK_PurityOK" = "red1",
    "Metastasis_ATRX_VAFOK_PurityOK" = "darkblue",
    "Metastasis_TERT_VAFOK_PurityLow" = "orange",
    "Metastasis_TERT_VAFOK_PurityOK" = "orange",
    "Metastasis_WT_VAFOK_PurityOK" = "red3",
    "Normal_WT_VAFOK_PurityOK" = "palevioletred2")

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

specimen_type_cols = c("Non-malignant Primary"=ColorPalette[["DarkGreen2"]], 
                       "Malignant primary"="orange",
                       `SFU` = ColorPalette[["Yellow3"]],
                       "Metastasis"= ColorPalette[["LightRed2"]], 
                       "Local Recurrance"= ColorPalette[["LightGrey2"]])
                              
specimen_type_cols[["MetastaticPrimary"]] <- specimen_type_cols[["Malignant primary"]]
specimen_type_cols[["NonMetastaticPrimary"]] <- specimen_type_cols[["Non-malignant Primary"]]
specimen_type_cols[["Other"]] <- "#999999"


sampletype_strict_cols <- c(`Non-metastatic primary` = ColorPalette[["DarkGreen2"]],
  `Non-metastatic local recurrence` = ColorPalette[["LightBlue1"]],
  `Primary (short follow up)` = ColorPalette[["Yellow3"]],
  `Primary (metastasis reported)` = ColorPalette[["DarkBrown2"]],
  `Local recurrence (metastasis reported)` = ColorPalette[["LightBlue3"]],
  `Metastatic primary` = ColorPalette[["LightOrange2"]],
  `Metastatic local recurrence` = ColorPalette[["DarkBlue3"]],
  `Metastasis` = ColorPalette[["LightRed2"]],
  `SFU` = ColorPalette[["DarkBlue3"]])


location_cols = c("Adrenal"=ColorPalette[["LightGreen1"]], 
                  "Extraadrenal"=ColorPalette[["DarkBlue1"]], 
                  "Head_neck"=ColorPalette[["LightOrange1"]], 
                  "Thoracic_non_chromaffin"=ColorPalette[["DarkRed1"]], 
                  "Metastasis"=ColorPalette[["Purple3"]], 
                  "Ambiguous"=ColorPalette[["DarkGrey1"]])
location_cols[["Head and neck"]] <- location_cols[["Head_neck"]] 
location_cols[["Extraadrenal_bladder"]] = "#99cccc"
location_cols[["Extraadrenal (abdominal/thoracic)"]] <-  location_cols[["Extraadrenal"]]
location_cols[["Abdominal_Thoracic"]] <-  location_cols[["Extraadrenal"]]
location_cols[["Thoracic (non-chromaffin)"]] <-  location_cols[["Thoracic_non_chromaffin"]]
location_cols[["Extraadrenal (bladder)"]] <-  location_cols[["Extraadrenal_bladder"]]

driver_cols <- c(ATRX=ColorPalette[["DarkBlue1"]], 
                TERT=ColorPalette[["LightOrange1"]], 
                WT=ColorPalette[["LightGreen1"]])

cell_of_origin_cols <- c(Chromaffin="#95665a",
                         Non_chromaffin="#f18d73ff",
                         Ambiguous="#b4b4b4ff")

all_cell_types <- c("Tumor",
                    "Endothelial cells" ,
                    "Chromaffin cells" ,
                    "Myeloid cells",
                    "Lymphocytes",
                    "Adrenocortical cells",
                    "Fibroblasts",
                    "SCLCs")

cell_type_colours <- setNames(brewer.pal(name = "Set1", n = 8), all_cell_types)


genomic_alteration_cols <- c(Missense=ColorPalette[["DarkBlue1"]], 
                             `Stop gained`=ColorPalette[["DarkRed1"]],
                             `Stop lost`=ColorPalette[["DarkOrange2"]], 
                             Frameshift=ColorPalette[["LightBrown2"]], 
                             `Promoter mutation`=ColorPalette[["Purple3"]],  
                             `Structural variant`=ColorPalette[["LightGreen1"]], 
                             `Splice acceptor`=ColorPalette[["DarkGrey2"]], 
                             `Splice donor`=ColorPalette[["DarkGrey1"]], 
                             `Splice region`=ColorPalette[["LightGrey1"]], 
                             `Homozyg. del.`=ColorPalette[["DarkGreen1"]],
                             `In-frame deletion`=ColorPalette[["DarkBrown1"]])

cna_palette <- c(`None`=ColorPalette[["LightGrey1"]],
                 Loss=ColorPalette[["DarkRed1"]],
                 `Subclonal Loss`="#fdadadff",
                 `Hom. Del.`=ColorPalette[["Yellow1"]], CNLOH="#f97344ff", 
                 Gain=ColorPalette[["DarkBlue3"]],
                 `Subclonal Gain`=ColorPalette[["LightBlue1"]],
                 WGD=ColorPalette[["Purple2"]],
                 Chromothripsis=ColorPalette[["LightGreen1"]],
                 Other=ColorPalette[["DarkGrey2"]],
                 `Gain+LOH`=ColorPalette[["LightBlue2"]],
                 `WGD+Gain`=ColorPalette[["DarkBlue1"]], 
                 `Loss + Subclonal CNLOH`=ColorPalette[["LightOrange2"]],
                 `Diploid/Haploid-X`=ColorPalette[["LightGrey1"]],
                 `Minor Subclonal Loss`= ColorPalette[["Salmon"]])


go_biotype_cols <- c("Cell adhesion" = c25[[1]],
                     "Cell cycle/proliferation" = c25[[2]],
                     "Cell migration" = c25[[3]],
                     "Cellular developmental process" = c25[[4]],
                     "Cellular response to DNA damage stimulus" = c25[[5]],
                     "Chromatin/chromosome organization" = c25[[6]],
                     "Nervous system development/neurogenesis" = c25[[7]],
                     "None" = "gray30",
                     "Other" = "gray70",
                     "Programmed cell death" = c25[[10]],
                     "Pseudogene" = c25[[11]],
                     "Regulation of DNA-templated transcription" = c25[[23]],
                     "lncRNA" = c25[[18]],
                     "Angiogenesis" = c25[[14]],
                     "DNA replication"  = c25[[15]])

bio_types <-  c(
   "protein_coding","processed_pseudogene",
   "transcribed_unprocessed_pseudogene", "snRNA","misc_RNA",
   "unitary_pseudogene","unprocessed_pseudogene", "miRNA",
   "Not available", "TEC", "pseudogene","snoRNA", "transcribed_processed_pseudogene",
   "IG_V_pseudogene","transcribed_unitary_pseudogene","rRNA_pseudogene", 
   "IG_D_gene","lncRNA", "TR_V_pseudogene","TR_V_gene")

bio_type_cols <- setNames(c25[1:length(bio_types)], bio_types)
rm(bio_types)


