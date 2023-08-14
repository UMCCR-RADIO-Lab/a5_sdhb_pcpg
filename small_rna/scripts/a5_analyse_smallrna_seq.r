library(tidyverse)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(scales)
library(ggrepel)
library(RColorBrewer)
library(egg)
library(circlize)
library(umap)

setwd("/g/data/pq08/projects/ppgl")

#######################
# Import data loaders #
#######################

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
source("./a5/small_rna/scripts/data_loaders/a5_smallrna_seq_dataloader.r")

######################
# Colours and Themes #
######################

source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#############
# Load data #
#############

data_loader_a5_smallrna(genome_version="hg38")

################
# Prepare data #
################

# Convert the TMM counts to log2 CPMs
a5_smallrna_lcpm_list <- map(.x = a5_smallrna_dge_list, \(x) cpm(x, log = T))

#############
# UMAP data #
#############

set.seed(10)
umap_config <- umap.defaults
umap_config$n_neighbors=10
#umap_config$spread=3
umap_list <- lapply(a5_wts_lcpm_list, function (lcpm, umc) { umap(t(lcpm), config = umc) }, umc=umap_config) 

annotate_umap <- function (umap_result, annotation)
{
  umap_layout <- data.frame(umap_result$layout)
  return(
    tibble(UMAP1 = umap_layout$X1,
           UMAP2 = umap_layout$X2,
           A5_ID = rownames(umap_layout),
           Dim = "1+2") %>% 
      left_join(annotation) %>%
      # Select key columns
      dplyr::select(A5_ID,
                    UMAP1,
                    UMAP2,
                    Gender,
                    differential_group,
                    # chr_14_miRNA_outgroup, 
                    `is_head_and_neck`,
                    TERT_ATRX_Mutation,
                    `sample_purity`,
                    #Run,
                    Dim,
                    `Patient ID`,
                    `Year of birth`,
                    Primary_Location_Simplified) %>% 
      tidyr::separate(differential_group, 
                      into=c("SampleType"), 
                      extra="drop", 
                      sep="_", 
                      remove=F) %>% 
      mutate(Primary_Location_Simplified=gsub("_abdominal|_thoracic|_bladder|_[Ll]eft|_[Rr]ight",
                                              "",
                                              Primary_Location_Simplified)) %>% 
      mutate(`sample_purity`=as.numeric(`sample_purity`))
  )
}


umap_annotated_list <- lapply(umap_list,annotate_umap, annotation=a5_anno.wts)

colour_by <- c("differential_group", "Primary_Location_Simplified", "SampleType", "sample_purity",
               "TERT_ATRX_Mutation", "is_head_and_neck", "Gender")

umap_plot_list <- lapply(umap_annotated_list, function (umap_tbl) {
  umap_ggplot <- ggplot(data = umap_tbl, aes(x = UMAP1, y = UMAP2))
  
  plot_list <- lapply(colour_by, function (colour_variable, plot_prototype) { 
    plot_prototype + 
      geom_point(aes(colour = !!sym(colour_variable))) +
      ggtitle(colour_variable)
  }, 
  plot_prototype=umap_ggplot )
  
  plot_list[[length(plot_list)+1]] <-  umap_ggplot + geom_text(aes(label=A5_ID))
  
  return(plot_list)
} )

#######
# MDS #
#######

# Look at various mds plot dims
mds <- plotMDS(lcpm)
mds_dims_3_4 <- plotMDS(lcpm, dim.plot=c(3,4))
mds_dims_5_6 <- plotMDS(lcpm, dim.plot=c(5,6))
mds_dims_7_8 <- plotMDS(lcpm, dim.plot=c(7,8))

# Join MDS info to annotation
mds_ggplot <- tibble(x = mds$x, y = mds$y, `A5 ID` = colnames(mds$distance.matrix), Dim = "1+2")
mds_ggplot_3_4 <- tibble(x = mds_dims_3_4$x,y =  mds_dims_3_4$y, `A5 ID` = colnames(mds$distance.matrix), Dim = "3+4")
mds_ggplot_5_6 <- tibble(x = mds_dims_5_6$x,y =  mds_dims_5_6$y, `A5 ID` = colnames(mds$distance.matrix), Dim = "5+6")
mds_ggplot_7_8 <- tibble(x = mds_dims_7_8$x,y =  mds_dims_7_8$y, `A5 ID` = colnames(mds$distance.matrix), Dim = "7+8")

mds_ggplot <- bind_rows(mds_ggplot, mds_ggplot_3_4, mds_ggplot_5_6, mds_ggplot_7_8)%>%
  left_join(annotation)
set.seed(20)
colours <- c(brewer.pal(8, "Dark2"), brewer.pal(4, "Set2"))

# Plot some possible and expected sources of variation
p <- ggplot(data = mds_ggplot, aes(x = x, y = y,label= `A5 ID`))+
  facet_wrap(~Dim, scales = "free", ncol = 2)+
  geom_text()+
  blank_theme+
  scale_colour_manual(values = colours)

# Seq batch
batch <- p + aes(colour = batch)+
  labs(colour = "Sequencing batch")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MDS plots/Batch samples.pdf")
# From this I will include location and gender as batch variables
Patient <- p + aes(colour = Patient_plot)+
  labs(colour = "Patient")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MDS plots/Patient samples.pdf")
Location <- p + aes(colour = `is_head_and_neck`)+
  labs(colour = "Sample location")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MDS plots/Sample location.pdf")
Gender <- p + aes(colour = Gender)+
  labs(colour = "Gender")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MDS plots/Gender.pdf")
# TERT and ATRX form their own clustering as well
TERT_ATRX <- p + aes(colour = TERT_or_ATRX)+
  labs(colour = "TERT or ATRX\nmutation status")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MDS plots/TERT or ATRX.pdf")
Met <- p + aes(colour = is_primary_or_met)+
  labs(colour = "Tumour state")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MDS plots/Met status.pdf")
chr14 <- p + aes(colour = chr_14_miRNA_outgroup)+
  labs(colour = "Chr14 outgroup")+
  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MDS plots/Chr14 status.pdf")

plot_list <- list(batch, Patient, Location, Gender, Met, TERT_ATRX, chr14)
pdf("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MDS plots/Combined_mds.pdf",width = 30, height = 20)
ggarrange(plots = plot_list,newpage = F,ncol = 3)
dev.off()

# Get the MDS plot back to dims 1 and 2 for further analysis 
mds_ggplot <- filter(mds_ggplot, Dim == "1+2")

# Limma pipeline with blocking

# Check ordering
annotation$`A5 ID` == rownames(rnaseq$samples)

# Set up a design matrix
design <- model.matrix(~0 + TERT_or_ATRX_for_design + `is_head_and_neck` + Gender, data = annotation)
colnames(design) <- gsub("TERT_or_ATRX_for_design| |`","", colnames(design))
rownames(design) <- rownames(rnaseq$samples)

# Estimate correlation with double voom for tumours from the same patient
# Solves the issue of how to deal with that
# https://support.bioconductor.org/p/59700/

# Voom and estimate correlation
v <- voom(rnaseq,design,plot = TRUE)
corfit <- duplicateCorrelation(v, design, block = annotation$`Patient ID`)
# Voom again with correlations known
v <- voom(rnaseq, design, block = annotation$`Patient ID`, correlation =
                              corfit$consensus)
# Fit a lineat mode to the design matrix
fit <- lmFit(v, design, block = annotation$`Patient ID`, correlation =
               corfit$consensus)

# Make a kind of complicated contrast matrix
cont.matrix <- makeContrasts(ATRX_all = (ATRX_MET_Shortfollowup + ATRX_MET_Yes)/2 - Non_met_primary_MET_No,
                             TERT_all = (TERT_MET_Shortfollowup + TERT_MET_No + TERT_MET_Yes)/3 - Non_met_primary_MET_No,
                             ATRX_met = ATRX_MET_Yes - Non_met_primary_MET_No,
                             TERT_met = TERT_MET_Yes - Non_met_primary_MET_No,
                             Tert_or_ATRX = (ATRX_MET_Shortfollowup + ATRX_MET_Yes + TERT_MET_Shortfollowup + TERT_MET_No + TERT_MET_Yes)/5 - Non_met_primary_MET_No,
                             Unknown_driver_met = Unknown_met_driver_met_MET_Yes - Non_met_primary_MET_No,
                             Unknown_driver_all = (Unknown_met_driver_met_MET_Yes + Unknown_met_driver_primary_MET_Yes)/2 - Non_met_primary_MET_No,
                             levels= design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont)

# Seems like there is a clear pattern in the TERT+ ATRX samples but not the others
summary(summa.fit)

for(contrast in colnames(summa.fit)){
  toptable <- topTable(fit.cont, coef= contrast, sort.by="p", number = Inf)%>%
    rownames_to_column("SYMBOL")%>%
    filter(adj.P.Val <= 0.1)%>%
    write_csv(paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Outputs/Toptables/", contrast, ".csv"))
}


contrasts <- c("ATRX_met", "TERT_met", "Tert_or_ATRX", "ATRX_all", "TERT_all")

for(contrast in contrasts){
  
  toptable <- topTable(fit.cont,coef=contrast,sort.by="p",number = Inf)
  
  plot_input <- toptable %>%
    rownames_to_column("SYMBOL")%>%
    mutate(label = paste0("italic('",SYMBOL, "')"))%>%
    mutate(colour = ifelse(adj.P.Val <= 0.05, "FDR <= 0.05", "FDR > 0.05"))%>%
    mutate(FC = ifelse(sign(logFC) > 0, "Up", "Down"))%>%
    mutate(colour = paste(colour, FC))%>%
    mutate(colour = replace(colour, colour %in% c("FDR > 0.05 Up", "FDR > 0.05 Down"),"FDR > 0.05"))%>%
    mutate(colour = factor(colour, levels = c("FDR > 0.05","FDR <= 0.05 Down","FDR <= 0.05 Up")))%>%
    arrange(adj.P.Val)%>%
    mutate(label = replace(label,adj.P.Val > 0.05, NA))%>%
    mutate(label = factor(label, levels = unique(label)))%>%
    arrange(colour)
  
  # Make a volcano and MA plot for each contrast
  volcano <- ggplot(data = plot_input, aes (x  = logFC, y = -log10(adj.P.Val), colour = colour))+
    geom_point( size = 0.5, alpha = 0.8)+
    geom_text_repel(parse = T, size =5,colour ="black", aes(x  = logFC, y = -log10(adj.P.Val),label = label))+
    scale_colour_manual(values = c("grey","blue", "red"))+
    labs(colour = "Significance", x = expression('Log'[2]*'FC'), y = expression('-Log'[10]*'FDR'))+
    geom_vline(xintercept = 0,alpha = 0.5,linetype="dashed")+
    #scale_y_continuous(expand = c(0, 0), limits = c(0, 8))+
    blank_theme+
    #ggtitle(contrast)+
    guides(colour = F)+
    ggsave(paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/Volcano plots/", contrast, "_volcano.pdf"),useDingbats = F)+
    ggsave(paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/Volcano plots/", contrast, "_volcano.png"))
  
  # MA plot
  MA <- ggplot(data = plot_input, aes (x  = AveExpr, y = logFC, colour = colour))+
    geom_point(size = 0.5, alpha = 0.8)+
    geom_text_repel(parse = T, size =5, colour ="black", aes(x  = AveExpr, y = logFC,label = label))+
    scale_colour_manual(values = c("grey","blue", "red"))+
    labs(x = expression('Mean expression (log'[2]*'CPM)'), y = expression('Log'[2]*'FC'), colour = "Significance")+
    geom_hline(yintercept = 0)+
    blank_theme+
    guides(colour = F)+
    ggsave(paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MA plots/", contrast, "_MA.pdf"),useDingbats= F)+
    ggsave(paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/MA plots/", contrast, "_MA.png"))
}


# Do an outlier analysis to find genes that are way out in ceretain samples

# Z score transform the log2 CPMs for the cohort
z_scores <- lcpm%>%
  t()%>%
  scale()%>%
  t()

hist(z_scores)

# Find outliers
outliers <- z_scores %>%
  data.frame()%>%
  rownames_to_column("miR")%>%
  gather(Sample, z_score, -miR)%>%
  filter(abs(z_score) > 4)

outlier_mirs <- z_scores[unique(outliers$miR),]

# Set the colours for the top anno
colours <- c(TERT = "red",ATRX =  "orange", "Non_met_primary" = "green","Unknown_met_driver_met" = "dark red", "Short_follow_up_primary" = "grey","Unknown_met_driver_primary" = "purple")
Location_cols <- c("H_N" = "purple", "PC_or_PGL" = "pink")

mds_ggplot$`A5 ID` == colnames(outlier_mirs)

# Set the top annotation
top_anno <- HeatmapAnnotation(`Sample group` = mds_ggplot$TERT_or_ATRX,
                              Location = mds_ggplot$`is_head_and_neck`,
                              col= list(`Sample group` = colours, 
                                        Location = Location_cols),
                              show_annotation_name = T)

pdf("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/Heatmaps/miRNA outliers.pdf", width =20, height = 10)
Heatmap(outlier_mirs, top_annotation = top_anno, column_split = mds_ggplot$TERT_or_ATRX,
        name = "Z score", column_title_rot = 90)
dev.off()

# Robust z score scaled lcpm, scaled by median and MAD stolen from data camp
# and applied using apply
# Lines up with the RobScale function from DescTools
robz <- function(row){
  # Compute median and mean absolute deviation for row
  m <- median(row,na.rm = T)
  s <- mad(row,na.rm = T)
  
  # If the MAD is 0, set it to a very small number
  if(s == 0){
    s <- 1E-100
  }
  robzscore <- (row - m) / (s)
  return(robzscore)
}

# Compute robust z-score for each observation
robust_z_score <- apply(lcpm,1,robz)%>%
  t()

# Find outliers by robust z score
outliers <- robust_z_score %>%
  data.frame()%>%
  rownames_to_column("miR")%>%
  gather(Sample, rob_z_score, -miR)%>%
  arrange(rob_z_score)

# Top 10 up and down by mad
mirs_to_plot <- c(unique(outliers$miR)[1:10], tail(unique(outliers$miR), 10))

# top 20 MAD outliers
outlier_mirs <- robust_z_score[mirs_to_plot,]

hist(robust_z_score)

col_fun = colorRamp2(c(-1E6, 1E6), c("green","red"))

# Robust Z scores have a crazy range
pdf("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Small RNA-Seq/Plots/Heatmaps/miRNA robust outliers.pdf", width =20, height = 10)
Heatmap(outlier_mirs, top_annotation = top_anno, column_split = annotation$TERT_or_ATRX,
        name = "Z score", column_title_rot = 90, col = col_fun)
dev.off()

# # miRNA database format:https://psicquic.github.io/MITAB27Format.html
# # Database from analysing the literature: https://rnajournal.cshlp.org/content/24/8/1005.full#sec-19
# miRNA_database <- read_tsv("~/Documents/Projects/Year_2018/A5study/Small_RNA_Seq/Databases/miRNA_go_terms.txt",col_names = F)%>%
#   select(miRNA = X5,Gene = X6,mir_species = X10,gene_species=X11, paper = X8)%>%
#   filter(mir_species == "taxid:9606(Homo sapiens)")%>%
#   filter(Gene != "-")%>%
#   mutate(Gene = gsub("uniprotkb:|\\(gene name\\)", "",Gene))%>%
#   mutate(miRNA = gsub('.*\\(human\\) |".*', "",miRNA))
# 
# down_mir <- filter(tt_tert, logFC < 0)
# 
# # Look at ATRX mir targets
# TERT_mir_down_genes <- miRNA_database%>%
#   filter(miRNA %in% down_mir$miR)%>%
#   filter(!duplicated(Gene))%>%
#   write_csv("~/Documents/Projects/Year_2018/A5study/Small_RNA_Seq/miRNA_targets_analysis/TERT_mir_down_genes.csv")
# 
# 
# 
# # Look at ATRX mir targets
# ATRX_mir_up_genes <- miRNA_database%>%
#   filter(miRNA %in% up_mir$miR)%>%
#   filter(!duplicated(Gene))%>%
#   write_csv("~/Documents/Projects/Year_2018/A5study/Small_RNA_Seq/miRNA_targets_analysis/ATRX_mir_up_genes.csv")
# 
# ATRX_mir_down_genes <- miRNA_database%>%
#   filter(miRNA %in% down_mir$miR)%>%
#   filter(!duplicated(Gene))%>%
#   write_csv("~/Documents/Projects/Year_2018/A5study/Small_RNA_Seq/miRNA_targets_analysis/ATRX_mir_down_genes.csv")
# 
# # Plot the CPMS
# CPM <- cpm(counts.keep,log = T)
# 
# # Plot a TERT heatmap
# pdf("~/Documents/Projects/Year_2018/A5study/Small_RNA_Seq/Plots/Results/TERT_genes_heatmap.pdf")
# mirs <- CPM[tt_tert$miR,]
# top_anno = HeatmapAnnotation(df = data.frame(TERT = TERT), col = list(TERT= c("Yes" = "dark red", "No" = "dark green", "Unknown_or_redundant" = "grey")))
# Heatmap(mirs-rowMeans(mirs),top_annotation = top_anno,name = "Log2 CPM-\nmean(log2 CPM)")
# dev.off()
# 
# # Do a boxplot of TERT expression 
# to_boxplot <- CPM[c("hsa-miR-181c-5p", "hsa-miR-181d-5p"),]%>%
#   data.frame(check.names = F)%>%
#   rownames_to_column("miR")%>%
#   gather(`A5 ID`, CPM, -miR)%>%
#   left_join(A5_clinical)%>%
#   mutate(redundant = replace(redundant, is.na(redundant), FALSE))%>%
#   mutate(TERT_mutant= replace(TERT_mutant, is.na(TERT_mutant)| redundant ==T, "Unknown or redundant"))%>%
#   mutate(ATRX_mutant= replace(ATRX_mutant, is.na(ATRX_mutant)| redundant ==T, "Unknown or redundant"))%>%
#   filter(TERT_mutant != "Unknown or redundant")%>%
#   select(miR, CPM, TERT_mutant,ATRX_mutant, `A5 ID`, Metastatic)
# 
# ggplot(data = to_boxplot, aes(x = TERT_mutant,y= CPM, fill = TERT_mutant))+geom_boxplot(outlier.alpha = 0)+
#   geom_jitter(size = 0.5, alpha = 0.5)+
#   facet_wrap(~miR)+
#   blank_theme+
#   labs(x = "TERT mutant", y = "Log2 CPM")+
#   guides(fill = F)+
#   scale_fill_manual(values = c("grey", "red"))+
#   ggsave("~/Documents/Projects/Year_2018/A5study/Small_RNA_Seq/Plots/Results/TERT_mirs_CPM_boxplot.pdf")
# 
# ggplot(data = to_boxplot, aes(x = ATRX_mutant,y= CPM, fill = ATRX_mutant))+geom_boxplot(outlier.alpha = 0)+
#   geom_jitter(size = 0.5, alpha = 0.5)+
#   facet_wrap(~miR)+
#   theme_bw(base_size = 25)
# 
# # Plot a weird sample heatmap
# mirs <- CPM[tt_weird$miR,]
# top_anno = HeatmapAnnotation(df = data.frame(TERT = TERT), col = list(Metastatic= c("Yes" = "dark red", "No" = "dark green")))
# Heatmap(mirs-rowMeans(mirs),top_annotation = top_anno,cluster_columns = F)
# 
# # Known pheo mirs
# pheo_mirs <- c("hsa-miR-21-3p", "hsa-miR-183-5p","hsa-miR-182-5p", "hsa-miR-96-5p", "hsa-miR-551b-3p","hsa-miR-202-5p")
# 
# # Plot a weird sample heatmap
# mirs <- CPM[rownames(CPM)%in% pheo_mirs,]
# top_anno = HeatmapAnnotation(df = data.frame(TERT = TERT), col = list(Metastatic= c("Yes" = "dark red", "No" = "dark green")))
# Heatmap(mirs-rowMeans(mirs),top_annotation = top_anno)
# 
# 
# # Plot the library sizes of each batch
# batch_plot <- data.frame(counts.keep)%>%
#   rownames_to_column(var = "gene")%>%
#   gather(`A5 ID`, count, -gene)%>%
#   mutate(`A5 ID` = gsub("\\.", "-", `A5 ID`))%>%
#   left_join(annotation)%>%
#   group_by(`A5 ID`,batch)%>%
#   summarise(total = sum(count))
# 
# ggplot(data = batch_plot, aes(x = batch, y = log2(total)))+ geom_boxplot()+
#   geom_jitter()+
#   blank_theme+
#   ggsave("~/Documents/Projects/Year_2018/A5study/Small_RNA_Seq/Plots/QC/Batch_assigned_total.pdf")
# 
# 
# mirs <- CPM[c("hsa-miR-181c-5p", "hsa-miR-181d-5p"),]
# 
# # Shorten TERT 
# TERT<- replace(TERT,TERT == "Unknown_or_redundant", "UK/R")
# 
# colnames(mirs) <- rep("", ncol(mirs))
# 
# top_anno = HeatmapAnnotation(df = data.frame(TERT = TERT), col = list(TERT= c("Yes" = "dark red", "No" = "dark green", "UK/R" = "grey")))
# pdf("~/Documents/Projects/Year_2018/A5study/Small_RNA_Seq/Plots/Results/TERT_genes_heatmap.pdf",width = 7,height = 3)
# Heatmap(mirs-rowMeans(mirs),top_annotation = top_anno,name = "Log2 CPM -\nrow means", 
#         row_title_gp = gpar(fontsize = 20))
# dev.off()
# 
# 
# # Known pheo mirs
# pheo_mirs <- c("hsa-miR-21-3p", "hsa-miR-183-5p","hsa-miR-182-5p", "hsa-miR-96-5p", "hsa-miR-551b-3p","hsa-miR-202-5p")
# 
# raw_all<- A5_raw[A5_raw$Geneid %in% pheo_mirs,]
# 
# anno_met <- select(miRNA_clinical, `A5 ID`, Metastatic, `Tumour state`)%>%
#   mutate(Metastatic = factor(Metastatic, levels = c("No", "Short follow up","Yes")))%>%
#   mutate(`Tumour state` = factor(`Tumour state`, levels = c("Primary", "Metastatic","Recurrent","Unknown")))%>%
#   arrange(`Tumour state`,Metastatic)
# 
# CPM_pheo <- CPM[rownames(CPM) %in% pheo_mirs,]
# CPM_pheo <- CPM_pheo - rowMeans(CPM_pheo)
# CPM_pheo <- CPM_pheo[,anno_met$`A5 ID`]
# colnames(CPM_pheo)== anno_met$`A5 ID`
# 
# top_anno = HeatmapAnnotation(df = data.frame(Metastatic = anno_met$Metastatic, State = anno_met$`Tumour state`), col = list(Metastatic= c("Yes" = "dark red", "Short follow up" = "grey", "No" = "dark green"),
#                                                                                                                             State = c("Primary" = "green", "Metastatic" = "red", "Recurrent" = "Yellow", "Unknown" = "grey")))
# Heatmap(CPM_pheo,cluster_columns = F,top_annotation = top_anno,name = "Difference from row mean (log2 CPM)")
# 
