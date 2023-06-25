setwd("/g/data/pq08/projects/ppgl")
renv::activate("/g/data/pq08/projects/ppgl/a5")

library(tidyverse)
library(limma)
library(edgeR)
library(umap)
library(ggrepel)
library(ComplexHeatmap)



################
# Data Loaders #
################

source("./a5/scripts/a5_clinical_annotation_dataloader.R")
source("./public_data/data_loaders/comete_smallrna_data_loader.r")
source("./public_data/data_loaders/tcga_smallrna_data_loader.r")
source("./a5/small_rna/scripts/data_loaders/a5_smallrna_seq_dataloader.r")

######################
# Colours and Themes #
######################

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

##################
# Load data sets #
##################


blank_theme <- theme_bw(base_size = 15)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

# Read in the featurecounts reads
TCGA_small_RNA_Seq <- read.table("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/counts/TCGA_subread_counts_miRbase_grch37_ref.tsv", header = T,sep = "\t")
colnames(TCGA_small_RNA_Seq) <- gsub("TCGA_bams\\.|\\.subread_results.bam", "", colnames(TCGA_small_RNA_Seq))
colnames(TCGA_small_RNA_Seq) <- gsub('.{8}$', '', colnames(TCGA_small_RNA_Seq))

A5_raw <- read.table("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/counts/A5_subread_counts_miRbase_grch37_ref.tsv", header = T,sep = "\t")
A5_raw <- A5_raw[,!grepl("Human.Brain.Total.RNA.subread_results.bam",colnames(A5_raw))]
colnames(A5_raw) <- gsub("bams\\.|\\.subread_results.bam", "", colnames(A5_raw))
colnames(A5_raw) <- gsub("\\.T0", "-", colnames(A5_raw))

COMETE_raw <-  read.table("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/counts/COMETE_subread_counts_miRbase_grch37_ref.tsv", header = T,sep = "\t")
colnames(COMETE_raw) <- gsub("bams\\.|\\.subread_results.bam|COMETE_", "", colnames(COMETE_raw))

TCGA_small_RNA_Seq[1:5,1:10]
A5_raw[1:5,1:10]
COMETE_raw[1:5,1:10]

# Sanity check labels are the same
sum(TCGA_small_RNA_Seq$Geneid == A5_raw$Geneid) == nrow(A5_raw)
sum(COMETE_raw$Geneid == A5_raw$Geneid) == nrow(A5_raw)

# Get the comete sample download info
comete_download <- read_tsv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Reference data/COMETE_small_RNA_refernce_E-MTAB-2833.sdrf.txt")%>%
  select(Sample = `Comment[ENA_RUN]`, Alias = `Comment[SUBMITTED_FILE_NAME]`)%>%
  mutate(Alias = gsub("_fastq.txt.gz", "", Alias))

COMETE_raw <- COMETE_raw[,7:ncol(COMETE_raw)]

rename_df <- data.frame(Sample = colnames(COMETE_raw), check.names = F, stringsAsFactors = F)%>%
  left_join(comete_download)

colnames(COMETE_raw) <- rename_df$Alias

# Combine the counts with the same refernce annotation
all_counts <- bind_cols(TCGA_small_RNA_Seq, A5_raw[,7:ncol(A5_raw)], COMETE_raw)

# Get the batch and clinical data

# A5_clinical <- read_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Clinical data/A5 full clinical and genomic version 2 - A5 full clinical and genomic version 2.csv")%>%
#   rename(A5_ID = `A5 ID`)

gs4_auth()
A5_clinical <- read_sheet("1hnXdXI29KvvuLxsaTBID1-EbE7mTSk6bG05HcSFfhgo", col_types = "c")
A5_clinical <- A5_clinical %>% dplyr::rename(`A5_ID`=`A5 ID`)

# Read in the compendium annotation that magnus generated
subtypes <- read_tsv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/RNA-Seq/Reference data/Pheo-atlas/bulk_metadata.tsv")%>%
  # Keep TCGA and COMETE samples
  filter(Dataset %in% c("TCGA","E-MTAB-733"))%>%
  dplyr::select(Barcode = "Sample.raw", Cluster, Dataset, Alias, Sex,Castro_Vega_miRNA_Classification, Malignancy)%>%
  mutate(Barcode = gsub("\\.", "-", Barcode))%>%
  # Cut down the TCGA barcodes
  mutate(Barcode = replace(Barcode, Dataset == "TCGA", gsub('.{8}$', '', Barcode[Dataset == "TCGA"])))%>%
  mutate(Barcode = replace(Barcode, Dataset == "E-MTAB-733", Alias [Dataset == "E-MTAB-733"]))%>%
  select(-Alias)%>%
  mutate(Dataset = replace(Dataset, Dataset == "E-MTAB-733", "COMETE"))

old_naming <- c("SDHx","SDHx (H&N)","VHL", "PH-NOS","Kinase", "IDC", "MAML3", "Normal", "Cortical admixture")
new_naming <- c("C1Ai (SDHx)","C1Aii (SDHx-HN)","C1Bi (VHL)", "C1Bii (PH-NOS)","C2A (Kinase)", 
                "C2Bi (IDC)", "C2Bii (MAML)","Normal", "C2C (Cortical admixture)")
names <- data.frame(Cluster = old_naming, new_naming)
# Set the subtype colours
subtpye_cols <- c("yellow3", "#A65628" , "#FF7F00", "#984EA3", "#E41A1C","#377EB8","#4DAF4A", "#FB9A99", "light blue")
names(subtpye_cols) <- c("C1Ai (SDHx)","C1Aii (SDHx-HN)","C1Bi (VHL)", "C1Bii (PH-NOS)","C2A (Kinase)", 
                         "C2Bi (IDC)", "C2Bii (MAML)", "Normal", "C2C (Cortical admixture)")

A5_HN <- A5_clinical%>%
  filter(`is_head_and_neck` == "H_N")
A5_outgroup <- A5_clinical%>%
  filter(chr_14_miRNA_outgroup != "Main_group")
A5_female <- A5_clinical%>%
  filter(Gender == "female")

mat <- all_counts[,7:ncol(all_counts)]

colnames(mat) <- gsub("\\.", "-", colnames(mat) )

rownames(mat) <- all_counts$Geneid

# Remove samples that were specific to comete in the future UMAP and seem
# to be some sort of technical effect
samples_to_remove <- readRDS("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Intermediate data/Comete batch effect samples.RDS")

mat <- mat[,!colnames(mat) %in% samples_to_remove]

annotation <- data.frame(Barcode = colnames(mat), check.names = F)%>%
  left_join(subtypes)%>%
  mutate(Cluster = replace(Cluster, Barcode == "TCGA-QR-A6GZ-01A-11R", "VHL"))%>%
  # Set the IDs SDHx for the A5
  # Add in the A5 samples that are known to be the wrong genotype
  #mutate(Cluster = replace(Cluster, Barcode == "E124-1", "Kinase"))%>%
  #mutate(Cluster = replace(Cluster, Barcode == "E145-1", "VHL"))%>%
  #mutate(Cluster = replace(Cluster, Barcode == "E154-1", "Cortical admixture"))%>%
  mutate(Dataset = replace(Dataset, Barcode %in% A5_clinical$A5_ID, "A5"))%>%
  mutate(Dataset = replace(Dataset, grepl("E1", Barcode) & is.na(Dataset), "A5"))%>%
  mutate(Dataset = replace(Dataset, grepl("TCGA", Barcode), "TCGA"))%>%
  mutate(Cluster = replace(Cluster, Dataset == "A5", "SDHx"))%>%
  # Updated with the known A5 head and necks
  mutate(Cluster = replace(Cluster, Barcode %in% A5_HN$A5_ID, "SDHx (H&N)"))%>%
  mutate(Sex = replace(Sex, Dataset == "A5", "Male"))%>%
  mutate(Sex = replace(Sex, Barcode %in% A5_female, "Female"))%>%
  mutate(Sex = replace(Sex, Sex == "F", "Female"))%>%
  mutate(Sex = replace(Sex, Sex == "M", "Male"))%>%
  # Add in samples where we didn't have gender for some reason that I found on the GDC data portal
  mutate(Sex = replace(Sex, Barcode %in%  c("TCGA-P7-A5NX-01A-11R", "TCGA-P8-A5KD-11A-11R"), "Female"))%>%
  mutate(Sex = replace(Sex, Barcode %in%  c("TCGA-P8-A5KC-11A-11R","TCGA-SQ-A6I4-11A-11R", "TCGA-QR-A6GZ-01A-11R"), "Male"))%>%
  left_join(names)%>%
  mutate(Label = Barcode)%>%
  mutate(Label = replace(Label,! Barcode %in% c("E124-1", "E154-1", "E145-1", "E183-1", "E230-1"), NA))%>%
  mutate(Size = 4)%>%
  # Change font size for known outliers
  mutate(Label = replace(Label,Barcode %in% A5_outgroup$A5_ID, "Chr14 outgroup"))%>%
  filter(!is.na(Dataset))

rnaseq <- DGEList(mat[,annotation$Barcode])

# Filter lowly expressed genes
keep.exprs <- filterByExpr(rnaseq, group = annotation$new_naming)
rnaseq <- rnaseq[keep.exprs,, keep.lib.sizes=FALSE]

# Apply TMM normalisation to the DGElist
rnaseq <- calcNormFactors(rnaseq)
# Convert the TMM counts to log2 CPMs
lcpm <- cpm(rnaseq, log = T)

# Sanity check for ordering
sum(annotation$Barcode == rownames(rnaseq$samples)) == length(rownames(rnaseq$samples))

design <- model.matrix(~0 + new_naming, data = annotation)
colnames(design) <- gsub("new_naming| ","", colnames(design))
rownames(design) <- rownames(rnaseq$samples)

# Remove the library specific batch effects
batch_removed <- removeBatchEffect(lcpm, batch=annotation$Dataset, batch2 = annotation$Sex, design = design)

# Plot an MDS of the data
mds <- plotMDS(batch_removed)

# Join MDS info to annotation
mds_ggplot <- tibble(x = mds$x,y =  mds$y, Barcode = colnames(mds$distance.matrix))%>%
  left_join(annotation)

ggplot(data = mds_ggplot, aes(x = x, y = y, colour = new_naming,shape= Dataset, label = Label)) + 
  geom_point()+
  geom_text_repel(color ="black", size = annotation$Size,nudge_y = -0.2)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "MDS dim 1", y = "MDS dim 2", colour = "Subtype")+
  scale_color_manual(values = subtpye_cols)+
  theme(aspect.ratio=1)+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/Integrated dataset plots/A5-TCGA MDS.pdf", useDingbats = F, width =10, height = 8)

set.seed(42)

# Plot a UMAP as well
umap_config <- umap.defaults
umap_config$n_neighbors=10
umap_config$spread=3

umap <- umap(t(batch_removed), config = umap_config)


saveRDS(umap, "/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Intermediate data/Small RNA-Seq UMAP.pdf")
umap <- readRDS("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Intermediate data/Small RNA-Seq UMAP.pdf")

to_plot_umap <- data.frame(umap$layout)%>%
  rownames_to_column("Barcode")%>%
  left_join(annotation)

write_delim(to_plot_umap, 
            file =  "/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/Integrated dataset plots/umap_coord_smallrna_com+tcga+a5_allsamples_nn10_seed42_spread3.txt", 
            delim = "\t")

ggplot(data = to_plot_umap, aes(x = X1, y = X2, colour = new_naming,shape= Dataset, label = Label)) + 
  geom_point()+
  geom_text_repel(color ="black", size = annotation$Size,nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  scale_color_manual(values = subtpye_cols)+
  theme(aspect.ratio=1)+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/Integrated dataset plots/A5-TCGA-COMETE UMAP.pdf", useDingbats = F, width =10, height = 8)

# Label the outgroups
outgroups <- to_plot_umap %>%
  mutate(miRNA_cluster = "Main")%>%
  mutate(miRNA_cluster = replace(miRNA_cluster, X2 < -3, "Cluster 1"))%>%
  mutate(miRNA_cluster = replace(miRNA_cluster, X1 > 10, "Cluster 2"))%>%
  mutate(Label = Castro_Vega_miRNA_Classification)%>%
  mutate(Label = replace(Label, miRNA_cluster == "Main", NA))%>%
  mutate(Alias = Barcode)%>%
  left_join(comete_download)

ggplot(data = outgroups, aes(x = X1, y = X2, colour = new_naming, shape= Dataset, label = Label)) + 
  geom_point()+
  geom_text_repel(color ="black", size = annotation$Size,nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  scale_color_manual(values = subtpye_cols)+
  theme(aspect.ratio=1)+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/Integrated dataset plots/A5-TCGA COMETE UMAP COMETE miR labels.pdf", useDingbats = F, width =10, height = 8)
