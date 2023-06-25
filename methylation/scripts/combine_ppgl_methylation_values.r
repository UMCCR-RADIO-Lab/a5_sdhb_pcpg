# library(TCGAbiolinks)
# library(DT)
# library(SummarizedExperiment)
# library(Homo.sapiens)
# library(umap)
# library(ggrepel)
# library(GEOquery)
# library(googlesheets4)

setwd("/g/data/pq08/projects/ppgl")
renv::activate("/g/data/pq08/projects/ppgl/a5")

library(limma)
library(edgeR)
library(umap)
library(tidyverse)
library(ggplot2)


#######################
# Import data loaders #
#######################

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
source("./public_data/data_loaders/comete_methylation_dataloader.r")
source("./public_data/data_loaders/tcga_methylation_dataloader.r")
source("./a5/methylation/scripts/data_loaders/a5_methylation_dataloader.r")

######################
# Colours and Themes #
######################

blank_theme <- theme_bw(base_size = 15)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

##################
# Load data sets #
##################

data_loader_comete_methylation(remove_27k_availaible_on_450k = T, only_zethoven_annotated = T)
data_loader_tcga_methylation()
data_loader_a5_methylation_array(quickload = T)
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.org", use_cache = T)

a5_methylation_beta <- 
  data.frame(getBeta(a5_methylation_filtered), check.names = F) %>% 
  tibble::rownames_to_column("Probe")
  
a5_methylation_mval <-  
  data.frame(getM(a5_methylation_filtered), check.names = F) %>% 
  tibble::rownames_to_column("Probe")

#####################
# Combine data sets #
#####################

mval_all_27k <- comete_methylation_27k_mval %>% 
  inner_join(comete_methylation_450k_mval) %>% 
  inner_join(tcga_methylation_450k_mval) %>% 
  inner_join(a5_methylation_mval)

beta_all_27k <- comete_methylation_27k_beta %>% 
  inner_join(comete_methylation_450k_beta) %>% 
  inner_join(tcga_methylation_450k_beta) %>% 
  inner_join(a5_methylation_beta)

mval_all_450k <- comete_methylation_450k_mval %>% 
  inner_join(tcga_methylation_450k_mval) %>% 
  inner_join(a5_methylation_mval)

beta_all_450k <- comete_methylation_450k_beta %>% 
  inner_join(tcga_methylation_450k_beta) %>% 
  inner_join(a5_methylation_beta)



# old_naming <- c("SDHx","SDHx (H&N)","VHL", "PH-NOS","Kinase", "IDC", "MAML3", "Normal", "Cortical admixture")
# new_naming <- c("C1Ai (SDHx)","C1Aii (SDHx-HN)","C1Bi (VHL)", "C1Bii (PH-NOS)","C2A (Kinase)", 
#                 "C2Bi (IDC)", "C2Bii (MAML)","Normal", "C2C (Cortical admixture)")
# names <- data.frame(Cluster = old_naming, new_naming)
# # Set the subtype colours
# subtpye_cols <- c("yellow3", "#A65628" , "#FF7F00", "#984EA3", "#E41A1C","#377EB8","#4DAF4A", "#FB9A99", "light blue")
# names(subtpye_cols) <- c("C1Ai (SDHx)","C1Aii (SDHx-HN)","C1Bi (VHL)", "C1Bii (PH-NOS)","C2A (Kinase)", 
#                          "C2Bi (IDC)", "C2Bii (MAML)", "Normal", "C2C (Cortical admixture)")

a5_anno.compatible <- a5_anno %>% dplyr::select(A5_ID, Gender, tumour_metastasised) %>% 
  dplyr::rename(Sample=A5_ID, Sex=Gender) %>% 
  filter(Sample %in% colnames(a5_methylation_beta)) %>% 
  mutate(Dataset="A5", 
         Genotype="SDHB",
         Malignancy=ifelse(tumour_metastasised=="Yes", "Malignant", "Benign")) %>% 
  dplyr::select(-tumour_metastasised)
  

annotation_27 <- bind_rows(comete_methylation_27k_anno %>% mutate(across(everything(),as.character)),
          comete_methylation_450k_anno %>% mutate(across(everything(),as.character)),
          tcga_methylation_450k_anno %>% mutate(across(everything(),as.character))) %>% 
  dplyr::select(Sample, Cluster, Dataset, Sex, Castro_Vega_miRNA_Classification, Malignancy, Genotype) %>%
  bind_rows(a5_anno.compatible)

#Remove 27k only samples to make annotation for 450K
annotation_450 <- annotation_27 %>% filter(Sample %in% colnames(beta_all_450k))  

#Reorder annotation to match data
annotation_27 <- annotation_27[match(colnames(beta_all_27k)[-1], annotation_27$Sample),]
annotation_450 <- annotation_450[match(colnames(beta_all_450k)[-1], annotation_450$Sample),]

# Sanity check for ordering
if(!(all(colnames(beta_all_27k)[-1]==annotation_27$Sample) & 
     all(colnames(beta_all_450k)[-1]==annotation_450$Sample))) {
  stop("Annotation ordering sanity check failed")
}

#Convert to matrix objects
mval_all_27k.matrix <- as.matrix(mval_all_27k[,-1])
mval_all_450k.matrix <- as.matrix(mval_all_450k[,-1])
rownames(mval_all_27k.matrix) <- mval_all_27k$Probe
rownames(mval_all_450k.matrix) <- mval_all_450k$Probe

#Create design matrices
design_27 <- model.matrix(~0 + Genotype, data = annotation_27)
colnames(design_27) <- gsub("Genotype| ","", colnames(design_27))
rownames(design_27) <- colnames(mval_all_27k.matrix)
design_450 <- model.matrix(~0 + Genotype, data = annotation_450)
colnames(design_450) <- gsub("Genotype| ","", colnames(design_450))
rownames(design_450) <- colnames(mval_all_450k.matrix)


# weights=ifelse(is.na(annotation_27$Genotype) | annotation_27$Genotype=="Unknown", 0, 1)

# Remove the library specific batch effects
mval_all_27k.batch_removed <- removeBatchEffect(mval_all_27k.matrix, batch=annotation_27$Dataset,batch2 =annotation_27$Sex, design = design_27)
mval_all_450k.batch_removed <- removeBatchEffect(mval_all_450k.matrix, batch=annotation_450$Dataset,batch2 =annotation_450$Sex, design = design_450)

# Plot an MDS of the data
mds <- plotMDS(mval_all_27k.batch_removed, plot = F)

# Join MDS info to annotation
mds_ggplot <- tibble(x = mds$x,y =  mds$y, Sample = colnames(mds$distance.matrix))%>%
  left_join(annotation_27)

ggplot(data = mds_ggplot, aes(x = x, y = y, colour = Genotype,shape= Dataset)) + 
  geom_point()+
  #geom_text_repel(color ="black", size = annotation$Size,nudge_y = -0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "MDS dim 1", y = "MDS dim 2", colour = "Subtype")+
  scale_color_manual(values = genotype_cols)+
  theme(aspect.ratio=1)
  #+ggsave("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Integrated dataset plots/Methylation MDS.pdf", useDingbats = F, width =10, height = 8)

set.seed(42)

umap_config <- umap.defaults
umap_config$n_neighbors=10
umap_config$spread=3

# Plot a UMAP as well
umap <- umap(t(mval_all_27k.batch_removed), config = umap_config)

to_plot_umap <- data.frame(umap$layout)%>%
  rownames_to_column("Sample")%>%
  left_join(annotation_27)

# write_delim(to_plot_umap, 
#             file =  "/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Integrated dataset plots/umap_coord_meth_com+tcga+a5_allsamples_nn10_seed42_spread3.txt", 
#             delim = "\t")

ggplot(data = to_plot_umap, aes(x = X1, y = X2, colour = Genotype,shape= Dataset)) + 
  geom_point()+
#  geom_text_repel(color ="black", size = annotation$Size,nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  scale_color_manual(values = genotype_cols)+
  theme(aspect.ratio=1)
#+  ggsave("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Integrated dataset plots/Methylation UMAP.pdf", useDingbats = F, width =10, height = 8)

ggplot(data = to_plot_umap, aes(x = X1, y = X2, colour = Genotype =="FH",shape= Dataset, label = Label)) + 
  geom_point()+
  geom_text_repel(color ="black",nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Genotype is FH")+
  theme(aspect.ratio=1)+
  ggsave("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Integrated dataset plots/FH plot.pdf", useDingbats = F, width =10, height = 8)

# Replot the UMAP excluding SDHx and normal 
no_SDHX <- annotation %>%
  filter(!new_naming %in% c("C1Ai (SDHx)", "C1Aii (SDHx-HN)", "Normal", "C2C (Cortical admixture)"))

umap <- umap(t(batch_removed[,no_SDHX$Barcode]), config = umap_config)

to_plot_umap <- data.frame(umap$layout)%>%
  rownames_to_column("Barcode")%>%
  left_join(annotation)

ggplot(data = to_plot_umap, aes(x = X1, y = X2, colour = new_naming,shape= Dataset, label = Label)) + 
  geom_point()+
  geom_text_repel(color ="black",nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  scale_color_manual(values = subtpye_cols)+
  theme(aspect.ratio=1)+
  ggsave("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Integrated dataset plots/Methylation UMAP no SDHx.pdf", useDingbats = F, width =10, height = 8)

write_delim(to_plot_umap, 
            file =  "/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Integrated dataset plots/umap_coord_meth_com+tcga+a5_nosdhx_nn10_seed42_spread3.txt", 
            delim = "\t")


SDHX_only <- annotation %>%
  filter(new_naming %in% c("C1Ai (SDHx)", "C1Aii (SDHx-HN)"))

umap <- umap(t(batch_removed[,SDHX_only$Barcode]), config = umap_config)

to_plot_umap <- data.frame(umap$layout)%>%
  rownames_to_column("Barcode")%>%
  left_join(annotation)

ggplot(data = to_plot_umap, aes(x = X1, y = X2, colour = new_naming,shape= Dataset, label = Label)) + 
  geom_point()+
  geom_text_repel(color ="black",nudge_y = 0.5)+
  blank_theme+
  guides(label= F)+
  labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = "Subtype")+
  scale_color_manual(values = subtpye_cols)+
  theme(aspect.ratio=1)+
  ggsave("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Integrated dataset plots/Methylation UMAP no SDHx.pdf", useDingbats = F, width =10, height = 8)

write_delim(to_plot_umap, 
            file =  "/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Integrated dataset plots/umap_coord_meth_com+tcga+a5_sdhxonly_nn10_seed42_spread3.txt", 
            delim = "\t")
