setwd("/g/data/pq08/projects/ppgl/a5/")
#renv::activate("./")

#library(knitr)



#library(stringr)
#library(minfiData)

#library(GSA)
#library(egg)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#library(IlluminaHumanMethylationEPICmanifest)
#library(Gviz)
#library(org.Hs.eg.db)
#library(ggrepel)


library(tidyverse)
library(RColorBrewer)
library(limma)
library(umap)
library(ggplot2)
library(DMRcate)
library(missMethyl)
library(patchwork)


###########
# Imports #
###########

source("./methylation/scripts/go_meth_offline.r")

###############
# DataLoaders #
###############

#load clinical annotation
source("./sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#load array data
source("./methylation/scripts/data_loaders/a5_methylation_dataloader.r")
data_loader_a5_methylation_array(quickload = T, output_qc = F, normalisation="functional")

#######################
# Colours and Themes #
#######################

source("./sample_annotation/scripts/data_loaders/a5_color_scheme.r")

blank_theme <- theme_bw(base_size = 18)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

###################
# Make annotation #
###################

#compute age at surgery
a5_anno <- a5_anno %>% mutate(age_at_resection=
                                floor(lubridate::interval(lubridate::my(paste("01", a5_anno$`Year of birth`,sep="-")),
                                                          lubridate::mdy(a5_anno$`Date of resection (DD/MM/YYYY)`)) 
                                      / lubridate::years(1)))

# reorder clinical data in same order as arrays
a5_anno.meth <- a5_anno %>% 
  mutate(A5_ID=factor(A5_ID, 
                      levels=colnames(a5_methylation_filtered))) %>% 
  arrange(A5_ID) %>% 
  filter(!is.na(A5_ID))

#Update annotation based on imaging review
a5_anno.meth <- a5_anno.meth %>% 
  mutate(Primary_Location_Simplified=replace(
    Primary_Location_Simplified, A5_ID=="E185-1", "Extraadrenal_thoracic"))

a5_anno.meth <- a5_anno.meth %>% 
  mutate(Primary_Location_Base=gsub("_abdominal|_thoracic|_bladder|_cardiac|_[Ll]eft|_[Rr]ight",
                                    "",
                                    Primary_Location_Simplified))


samples.exclude <- a5_anno.meth %>% filter(Exclude=="Y") %>% pull(A5_ID)
samples.hn <- a5_anno.meth %>% filter(Primary_Location_Base == "Head_neck") %>% pull(A5_ID)

a5_anno.meth.noex <- a5_anno.meth %>% filter(!(A5_ID %in% samples.exclude))

a5_anno.meth.nohn_noex <-  a5_anno.meth %>% filter(!(A5_ID %in% c(samples.hn, samples.exclude)))

#####################################
# Make contrast and design matrices #
#####################################

source("./sample_annotation/scripts/data_loaders/a5_contrast_groups.r")

make_hn_vs_abdominothoracic_contrasts(sample_anno = a5_anno.meth.noex,
                                      exclude_samples = c("E185-1", "E167-2", "E135-1", "E166-1",  "E166-2", "E188-1"))

make_genotype_sampletype_contrasts(sample_anno = a5_anno.meth.nohn_noex, 
                                   exclude_samples = c("E185-1", "E167-2", "E135-1", "E166-1",  "E166-2", "E188-1"))


####################
# Compute M/B-values #
####################

# calculate M-values for statistical analysis
m_vals <- getM(a5_methylation_filtered)
b_vals <- getBeta(a5_methylation_filtered)

################
# MDS QC plots #
################


pdf(file = "methylation/results/plots/mds/a5_methylation_mds_dim1_2_3.pdf", onefile = T, width=10, height = 10)
for (dim_pair in list(c(1,2), c(2,3)))
{
  mds <- plotMDS(m_vals, top=10000, gene.selection="common",
                 col=pal[factor(targets$Sample_Group)],dim=dim_pair, plot=F)
  
  plot.data <- tibble(x = mds$x, 
                      y =  mds$y, 
                      A5_ID = colnames(mds$distance.matrix)) %>%
    left_join(a5_anno) %>%
    left_join(a5_methylation_targets %>% dplyr::select(A5_ID=Sample_Name, Batch=Sample_Group)) %>% 
    dplyr::select(A5_ID,Batch, x, y, Gender, Batch, Primary_Location_Simplified, 
                  TERT_ATRX_Mutation,is_primary_or_met, age_at_resection,
                  chr_14_miRNA_outgroup) %>%
    mutate(Patient = gsub("-.*", "", A5_ID)) 
  
  for (covar in c("Gender", "Batch", "Primary_Location_Simplified", "TERT_ATRX_Mutation",  "is_primary_or_met", "chr_14_miRNA_outgroup"))
  {
    print(ggplot(data = plot.data, aes(x = x, y = y,label= A5_ID)) +
            geom_text() +
            blank_theme +
            theme(aspect.ratio=1) +
            aes(colour = !!sym(covar)) +
            labs(colour = covar, x= "MDS dim 1", y = "MDS dim 2") + ggtitle(paste(covar, "- Dims", dim_pair[[1]],"+", dim_pair[[2]])))
  }
}
dev.off()

#################
# UMAP QC plots #
#################



seed_val <- 10
set.seed(seed_val)
for (nn in c(5,10,15,20))
{
  umap_config <- umap.defaults
  umap_config$n_neighbors=nn
  umap_config$metric="euclidean"
  umap <- umap(t(m_vals),config = umap_config)
  
  to_plot_umap <- data.frame(umap$layout) %>%
    rownames_to_column("A5_ID") %>%
    left_join(a5_anno) %>% 
    left_join(a5_methylation_targets %>% dplyr::select(A5_ID=Sample_Name, Batch=Sample_Group)) %>% 
    dplyr::select(
      A5_ID, Batch, X1, X2,
      Gender, Batch, Primary_Location_Simplified, TERT_ATRX_Mutation, 
      is_primary_or_met, chr_14_miRNA_outgroup, Catecholamine_profile
    ) %>%
    mutate(Patient = gsub("-.*", "", A5_ID))
  
  
  write.table(to_plot_umap,file=paste0("./methylation/results/plots/umap/umap_mvals_allsamples_nn",nn,"seed", seed_val,".tsv"), sep="\t", row.names = F)
  
  pdf(paste0("./methylation/results/plots/umap/umap_mvals_allsamples_nn",nn,"seed", seed_val,".pdf"), width = 15, height = 15, onefile = T)
  for (covar in c("Gender", "Batch", "Primary_Location_Simplified", "TERT_ATRX_Mutation",  "is_primary_or_met", "chr_14_miRNA_outgroup", "Catecholamine_profile"))
  {
    print(   ggplot(data = to_plot_umap, aes(x = X1, y = X2, label=A5_ID, colour = !!sym(covar))) + 
               #geom_point() +
               geom_text()+
               #  geom_text_repel(color ="black",nudge_y = 0.5)+
               blank_theme+
               guides(label= F)+
               labs(shape ="", x= "UMAP 1", y = "UMAP 2", colour = covar)+
               theme(aspect.ratio=1) + ggtitle(paste0(covar, " - ", "nNeighbours=",umap_config$n_neighbors)))
  }
  dev.off()
}


########################
# B/M Density QC plots #
########################

#QC now generated by data-loader script

# pdf("./methylation/qc/Beta_M_values_filtered.pdf", width =10)
# par(mfrow=c(1,2))
# densityPlot(getBeta(a5_methylation_filtered), sampGroups=a5_methylation_targets$Sample_Group, main="Beta values",
#             legend=FALSE, xlab="Beta values")
# densityPlot(m_vals, sampGroups=a5_methylation_targets$Sample_Group, main="M-values",
#             legend=FALSE, xlab="M values")
# dev.off()

#########################################
### Differential Methylation Analysis ###
#########################################

fit_objects <- list()
design_matrices <- list()
contrast_matrices <- list()

##############
# Head and Neck vs Adrenal/Abdominal/Thoracic
##############

count_contrast_members(contrast_matrix_hn, design_matrix_hn)

# Find the correlation between samples from the same patient
# Warnings are apparently ok for only a few genes:
# https://support.bioconductor.org/p/6618/
file_hashes <- purrr::map(list(m_vals, design_matrix_hn, a5_anno.meth.noex), digest::digest, algo = "md5")
suffix <- paste(purrr::map(file_hashes, stringr::str_sub, start=25, end=32), collapse = "_")
checkpoint_file <- paste0("./methylation/quickload_checkpoints/corfit_mvals_all_samples_",suffix,".rds")
if(file.exists(checkpoint_file)) {
  corfit <- readRDS(checkpoint_file)
} else {
  corfit <- duplicateCorrelation(m_vals, design_matrix_hn, block = a5_anno.meth.noex$`Patient ID`)
  # Save the RDS as this step has a long processing time
  saveRDS(corfit, checkpoint_file)
}

# fit the linear model with blocking where the same sample is from multiple patients
hn_vs_abdo_fit <- lmFit(m_vals, design_matrix_hn, block = a5_anno.meth.noex$`Patient ID`, correlation = corfit$consensus)

hn_vs_abdo_fit2 <- contrasts.fit(hn_vs_abdo_fit, contrast_matrix_hn)
hn_vs_abdo_fit2 <- eBayes(hn_vs_abdo_fit2)

hn_vs_abdo_testresult <- decideTests(hn_vs_abdo_fit2)

summary(hn_vs_abdo_testresult)

fit_objects[["hn_vs_abdo"]] <- hn_vs_abdo_fit2
contrast_matrices[["hn_vs_abdo"]] <- contrast_matrix_hn
design_matrices[["hn_vs_abdo"]] <- design_matrix_hn

#####
# BACON inflation adjustment
#####
# library(bacon)
# current_coef <- "Non_chromaffin_vs_Chromaffin"
# z_score_t_geno <- zscoreT(x = geno_fit2$t, 
#                           df=max(geno_fit2$df.total))
# bc_geno <- bacon(teststatistics = z_score_t_geno, niter = as.integer(5000),nburnin = as.integer(1000))
# estimates(bc_geno)
# #plot(bc_geno, type="qq")
# bacon_corrected_pval <- pval(bc_geno)
# 
# test <- limma::topTable(geno_fit2, coef = current_coef, number = nrow(geno_fit2)) %>%  data.frame(check.names = F) %>% 
#   tibble::rownames_to_column("probe_id") %>% 
#   inner_join(bacon_corrected_pval %>%  data.frame(check.names = F) %>% 
#                tibble::rownames_to_column("probe_id") %>% select(probe_id, bacon_p=!!sym(current_coef))) %>% 
#   mutate(bacon_p_adj=p.adjust(bacon_p))


##############
# Genotype/Clinical-course comparisons
##############

m_val.nohn <- m_vals[,a5_anno.meth.nohn_noex$A5_ID]

count_contrast_members(contrast_matrix_genosampletype, design_matrix_genosampletype)

# Find the correlation between samples from the same patient
# Warnings are apparently ok for only a few genes:
# https://support.bioconductor.org/p/6618/
# m_val.nohn="4ff191b0eff9be2836737c5057ac8217"
# design_matrix_genosampletype="c24c546980b968a9d8d455ce77f05f39"
# a5_anno.meth.nohn_noex="2abc07df3f550ce0eb4cc7fba32ee5a6"
file_hashes <- purrr::map(list(m_val.nohn, design_matrix_genosampletype, a5_anno.meth.nohn_noex), digest::digest, algo = "md5")
suffix=paste(purrr::map(file_hashes, stringr::str_sub, start=25, end=32), collapse = "_")
checkpoint_file=paste0("./methylation/quickload_checkpoints/corfit_mvals_nohn_",suffix,".rds")
if(file.exists(checkpoint_file)) {
  corfit <- readRDS(checkpoint_file)
} else {
  corfit <- duplicateCorrelation(m_val.nohn, design_matrix_genosampletype, block = a5_anno.meth.nohn_noex$`Patient ID`)
  # Save the RDS as this step has a long processing time
  saveRDS(corfit, checkpoint_file)
}


# fit the linear model with blocking where the same sample is from multiple patients
geno_fit <- lmFit(m_val.nohn, design_matrix_genosampletype, block = a5_anno.meth.nohn_noex$`Patient ID`, correlation = corfit$consensus)


geno_fit2 <- contrasts.fit(geno_fit, contrast_matrix_genosampletype)
geno_fit2 <- eBayes(geno_fit2)

geno_testresult <- decideTests(geno_fit2)

# Seems like there is a clear pattern in the TERT+ ATRX samples but not the others
summary(geno_testresult)

fit_objects[["genotype_sampletype"]] <- geno_fit2
contrast_matrices[["genotype_sampletype"]] <- contrast_matrix_genosampletype
design_matrices[["genotype_sampletype"]] <- design_matrix_genosampletype

#####
# BACON inflation adjustment
#####
# library(bacon)
# current_coef <- "ATRX_AllvsNonATRX"
# z_score_t_geno <- zscoreT(x = geno_fit2$t, 
#                      df=max(geno_fit2$df.total))
# bc_geno <- bacon(teststatistics = z_score_t_geno, niter = as.integer(5000),nburnin = as.integer(1000))
# estimates(bc_geno)
# #plot(bc_geno, type="qq")
# bacon_corrected_pval <- pval(bc_geno)
# 
# test <- limma::topTable(geno_fit2, coef = current_coef, number = nrow(geno_fit2)) %>%  data.frame(check.names = F) %>% 
#   tibble::rownames_to_column("probe_id") %>% 
#   inner_join(bacon_corrected_pval %>%  data.frame(check.names = F) %>% 
#                tibble::rownames_to_column("probe_id") %>% select(probe_id, bacon_p=!!sym(current_coef))) %>% 
#   mutate(bacon_p_adj=p.adjust(bacon_p))

########
# DMRs #
########

label_me <- c("E134-1","E144-1","E186-1","E143-3")

contrasts_to_use <- NULL #c("Metastasis_AllvsNonMetPri_WT", "ATRX_AllvsTERT_All")

dmr_lists <- list()
plots_summary_data <- list()
plots_detailed <- list()
for (comparison in c("hn","genosampletype")) {
  contrast_matrix <- get(paste0("contrast_matrix_", comparison))
  design_matrix <- get(paste0("design_matrix_", comparison))
  available_contrasts <- dimnames(contrast_matrix)$Contrasts
  
  if(!is.null(contrasts_to_use))
  {
    available_contrasts <- intersect(available_contrasts, contrasts_to_use)
  }
  
  top_lists[[comparison]] <- list()
  
  for(contrast in available_contrasts){
    
    contrast_membership <- contrastdesign_to_memberlist(contrast_name = contrast, 
                                                        contrast_matrix = contrast_matrix,
                                                        design_matrix = design_matrix)
    
    
    myAnnotation <- cpg.annotate(object = m_vals[,rownames(design_matrix)], datatype = "array", what = "M",
                                 analysis.type = "differential", design = design_matrix,
                                 contrasts = TRUE, cont.matrix = contrast_matrix,
                                 coef = contrast, arraytype = "EPIC")
    
    # Run a function to calculate differentially methylated regions
    DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
    results.ranges <- extractRanges(DMRs,genome = "hg19")
    
    dmr_lists[[comparison]][[contrast]] <- results.ranges
    
    # results.ranges <- dmr_lists[[fit]][[contrast]]
    
    write_delim(data.frame(results.ranges), 
                paste0("./methylation/results/tables/dmr/",
                       fit,"_", 
                       contrast,"_dmrcate_table.tsv"))
    
    # Plot dmrcate plots for the top 10 regions
    for (i in 1:30){
      
      overlapping.genes <- mcols(results.ranges)[i,"overlapping.genes"]
      
      if(is.na(overlapping.genes)) { next }
      
      region_probes <- epic_array_annotation_hg38 %>% GenomicRanges::as.data.frame() %>% 
        filter(chr==as.character(seqnames(results.ranges[i,])),
               pos >= start(results.ranges[i,]),
               pos <= end(results.ranges[i,])) %>% 
        filter(Name %in% rownames(m_vals))
      
      probe_data_b <- b_vals[region_probes$Name,] %>% data.frame(check.names = F) %>% 
        tibble::rownames_to_column("probe_id") %>% 
        pivot_longer(cols=-probe_id, names_to = "A5_ID", values_to = "b_val")
      
      probe_data_m <- m_vals[region_probes$Name,] %>% data.frame(check.names = F) %>% 
        tibble::rownames_to_column("probe_id") %>% 
        pivot_longer(cols=-probe_id, names_to = "A5_ID", values_to = "m_val")
      
      probe_data <-  probe_data_b %>% inner_join(probe_data_m)
      
      plot.data <- probe_data %>%  
        inner_join(region_probes %>% select(probe_id=Name, chr, pos)) %>%  
        inner_join(a5_anno.meth %>%  select(A5_ID, Primary_Location_Base, differential_group)) %>% 
        inner_join(contrast_membership) 
      
      
      plot.data <- plot.data %>% 
        bind_rows(plot.data %>% group_by(A5_ID, 
                                         Primary_Location_Base, differential_group, 
                                         group) %>% 
                    summarise(probe_id="mean_all_probes", 
                              pos=max(pos)+1, 
                              b_sd=sd(b_val, na.rm = T),
                              m_sd=sd(m_val, na.rm = T),
                              b_val=mean(b_val, na.rm = T), 
                              m_val=mean(m_val, na.rm = T)))
      
      
      plot.data <- plot.data %>% arrange(pos) %>% mutate(probe_id = factor(probe_id, levels=unique(.$probe_id)))
      
      plots_summary_data[[comparison]][[contrast]][[i]] <- plot.data %>%  filter(probe_id=="mean_all_probes") %>% group_by(group) %>% 
        mutate(m_z=(m_val-mean(m_val, na.rm = T))/sd(m_val), 
               label=ifelse(abs(m_z) > 2 | n() <= 5 | A5_ID %in% label_me , A5_ID, NA),
               overlapping_genes=overlapping.genes) 
      
      
      plots_detailed[[comparison]][[contrast]][[i]] <-
        ggplot(plot.data, aes(x=group, y=b_val)) + 
        geom_boxplot() +
        geom_jitter(mapping=aes(color=differential_group, shape=Primary_Location_Base), width = 0.2) + 
        scale_color_manual(values = differential_group_colors) + facet_wrap("probe_id") + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1)) + 
        ggtitle(paste(contrast,overlapping.genes))
      
      # jitter_width=0.4
      # plots_summary[[comparison]][[contrast]][[i]] <- 
      #   ggplot(plot.data.summary, aes(x=group, y=b_val)) + 
      #   geom_boxplot(outlier.alpha = 0) +
      #   geom_segment(mapping=aes(xend=group, y=b_val+b_sd,yend=b_val-b_sd), 
      #                position = position_jitter(seed = 10, width = jitter_width), alpha=0.05) +
      #   scale_color_manual(values = differential_group_colors) + 
      #   geom_point(mapping=aes(color=differential_group, shape=Primary_Location_Base),
      #              position = position_jitter(seed = 10, width = jitter_width)) +
      #   ggrepel::geom_text_repel(aes(label=label), position = position_jitter(seed = 10, width = jitter_width)) +
      #   theme_bw() + 
      #   theme(axis.text.x = element_blank()) + 
      #   ggtitle(overlapping.genes, subtitle = contrast)
      
      
    }
    
    plot_file <- paste0("./methylation/results/plots/dmr/",comparison,"_",contrast,".pdf")
    
    pdf(file = plot_file, width = 20, height = 15, onefile = T)
    
    summary_data <- bind_rows(plots_summary_data[[comparison]][[contrast]])
    
    print(ggplot(summary_data, aes(x=group, y=b_val)) + 
            geom_boxplot() +
            geom_jitter(mapping=aes(color=differential_group, shape=Primary_Location_Base), width = 0.2) + 
            scale_color_manual(values = differential_group_colors) + 
            facet_wrap("overlapping_genes", ncol=5) +
            theme_bw() + 
            theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1)) + 
            ggtitle(paste(contrast, " - mean probe betas")))
    
    purrr::walk(plots_detailed[[comparison]][[contrast]], print)
    
    dev.off()
  }
}


########
# GSEA #
########

probe_gene_annotation <- epic_array_annotation_hg38[match(rownames(m_val.nohn),epic_array_annotation_hg38$Name),
                                                    c(1:4,12:19,24:ncol(epic_array_annotation_hg38))]

collections <- list.files("/g/data/pq08/reference/msig_db/")

# Get the table of results for the first contrast
DMPs <- limma::topTable(fit_objects[[fit]], num=Inf, coef=contrast, genelist=probe_gene_annotation)

# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs$Name[DMPs$adj.P.Val < 0.05]

# Get all the CpG sites used in the analysis to form the background
all <- DMPs$Name

#GO term analysis of differentially methylated probes
for(collection in collections){
  
  gene_set <- GSA.read.gmt(collection)
  
  gene_set_formatted <- gene_set$genesets
  
  names(gene_set_formatted) <- gene_set$geneset.names
  
  collection_name <- gsub(".entrez.gmt","", basename(collection))
  
  # Perform the gene set test
  print(contrast)
  print(collection)
  gsa <- gsameth(sig.cpg=sigCpGs, all.cpg=all, collection=gene_set_formatted, array.type ="EPIC", plot.bias=F)
  
  top_gene_set <- topGSA(gsa, number=Inf)%>%
    rownames_to_column("Gene set")%>%
    #filter(FDR < 0.05)%>%
    mutate(collection = collection_name)%>%
    write_csv(paste0("results/tables/gometh/", contrast,"_",collection_name , "_signf_go_terms.csv"))
}
}
}
# plot the top 4 most significantly differentially methylated CpGs 
pdf("results/plots/dmr/top4_dmps.pdf", width = 10, height = 10)
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(bVals[,colnames(bVals) %in% colnames(m_vals.noHN)], cpg=cpg, pheno=design_anno$TERT_or_ATRX, ylab = "Beta values")
})
dev.off()

# Compare methylation and RNA-Seq counts

# Get the methylation annotation
anno_df_meth <- data.frame(annEPIC)%>%
  rownames_to_column("Probe")%>%
  dplyr::select(Probe, chr, pos, SYMBOL = GencodeBasicV12_NAME, Regulatory_Feature_Group)%>%
  mutate(SYMBOL = gsub(";.*", "", SYMBOL))%>%
  filter(SYMBOL != "")

# Get the beta values (% methylation) to plot against the CPMs
bvals_anno <- bVals %>% rownames_to_column(var = "Probe") %>% inner_join(anno_df_meth) %>% 
  rename_with(.fn = ~gsub(".*\\.|_D|T0", "", .x)) %>%  rename_with(.fn = ~gsub("(E[0-9]{3})_([0-9])", "\\1-\\2", .x))

bvals_gene_mean <- bvals_anno %>% 
  group_by(SYMBOL) %>% summarise(across(.cols = matches("^E[0-9]{3}-[0-9]"), .fns = mean, na.rm = T)) 

bvals_genepromoter_mean <- bvals_anno %>% filter(Regulatory_Feature_Group=="Promoter_Associated") %>% 
  group_by(SYMBOL) %>% summarise(across(.cols = matches("^E[0-9]{3}-[0-9]"), .fns = mean, na.rm = T))

bvals_gene_mean_long <- bvals_gene_mean %>%
  pivot_longer(cols = -SYMBOL, names_to = "Sample", values_to = "B") 

bvals_genepromoter_mean_long <- bvals_genepromoter_mean %>%
  pivot_longer(cols = -SYMBOL, names_to = "Sample", values_to = "B") 


# Correlate the b values and the RNA-Seq log2 CPMs
lcpms <- read_csv("../WTS/results/differential_gene_expression/bulk_RNA_seq/bulk_counts_hg38_tmm_normalized_log2cpm_All.csv") %>%
  dplyr::rename(EnsGeneID_Symbol=X1) %>% 
  filter(!is.na(EnsGeneID_Symbol))%>%
  pivot_longer(cols = -EnsGeneID_Symbol, names_to = "Sample", values_to = "log2_cpm") %>% 
  separate(EnsGeneID_Symbol, into = c("EnsGeneID", "SYMBOL"), sep="_", extra = "merge")

mean_bval_lcpm <- bvals_gene_mean_long %>%  dplyr::rename(gene_mean_b=B) %>% 
  left_join(
    bvals_genepromoter_mean_long %>%  dplyr::rename(genepromoter_mean_b=B)) %>% 
  inner_join(lcpms %>% select(-EnsGeneID)) %>% dplyr::rename(A5_ID=Sample)


# Try correlating some key genes
correlate_genes <-function(mean_bval_lcpm, gene){
  
  gene_filt <- mean_bval_lcpm %>%
    filter(SYMBOL == gene)%>%
    left_join(A5_clinical)
  
  colours <- c("TERT" = "red","ATRX" =  "orange","Unknown_met_driver_met" = "dark red","Unknown_met_driver_primary" = "purple",
               "Short_follow_up_primary" = "grey","Non_met_primary" = "dark green")
  
  plot.list <- list()
  
  plotbase  <- ggplot(data = gene_filt, aes(y = log2_cpm, colour = TERT_or_ATRX, label = A5_ID))+
    geom_point()+
    geom_text_repel(size =2)+
    blank_theme+
    labs(x = "Mean gene methylation B value", y = "Log2 CPM")+
    ggtitle(gene)+
    coord_equal()+
    theme(aspect.ratio = 1)+
    scale_colour_manual(values = colours) 
  
  
  plot.list[[1]] <- plotbase + aes(x=gene_mean_b)
  plot.list[[2]] <- plotbase + aes(x=genepromoter_mean_b)
  
  ggsave(filename = paste0("results/plots/expression_vs_methylation/",gene, ".pdf"), plot = plot.list[[1]], width = 25, height = 15, useDingbats = F)
  ggsave(filename = paste0("results/plots/expression_vs_methylation/",gene, "_promoter_annotated_only.pdf"), plot = plot.list[[2]], width = 25, height = 15, useDingbats = F)
  
  return(plot.list)
  
}

genes <- c("TERT", "CDKN2A")

plot.list <- list()
for(gene in genes){
  plot.list[[gene]]  <- correlate_genes(mean_bval_lcpm, gene)
}

# TERT specific analysis
# Get an average B value per gene per sample
anno_df_meth_promo_TERT <- bvals_anno %>% 
  filter(Probe == "cg02545192") %>%
  pivot_longer(cols = c(-Probe, -SYMBOL, -Regulatory_Feature_Group, -chr,-pos), names_to = "Sample", values_to = "B") %>% 
  left_join(lcpms) %>%
  filter(!is.na(log2_cpm)) %>%
  ungroup() %>%
  dplyr::rename(A5_ID=Sample) %>%
  left_join(A5_clinical)

colours <- c("TERT" = "red","ATRX" =  "orange","Unknown_met_driver_met" = "dark red","Unknown_met_driver_primary" = "purple",
             "Short_follow_up_primary" = "grey","Non_met_primary" = "green")

ggplot(data = anno_df_meth_promo_TERT, aes(y = B, x = log2_cpm, colour = TERT_or_ATRX, label = A5_ID))+
  geom_point()+
  geom_text_repel(size =2)+
  blank_theme+
  labs(y = "cg02545192 B value", x = "Log2 CPM")+
  ggtitle("TERT")+
  coord_equal()+
  theme(aspect.ratio = 1)+
  scale_colour_manual(values = colours)+
  ggsave("results/plots/expression_vs_methylation/TERT_cg02545192_beta_vs_cpm.pdf")


mvals_anno <- m_vals %>% rownames_to_column(var = "Probe") %>% inner_join(anno_df_meth) %>% 
  rename_with(.fn = ~gsub(".*\\.|_D|T0", "", .x)) %>%  rename_with(.fn = ~gsub("(E[0-9]{3})_([0-9])", "\\1-\\2", .x))

mvals_gene_mean <- mvals_anno %>% 
  group_by(SYMBOL) %>% summarise(across(.cols = matches("^E[0-9]{3}-[0-9]"), .fns = mean, na.rm = T)) 

mvals_genepromoter_mean <- mvals_anno %>% filter(Regulatory_Feature_Group=="Promoter_Associated") %>% 
  group_by(SYMBOL) %>% summarise(across(.cols = matches("^E[0-9]{3}-[0-9]"), .fns = mean, na.rm = T))

mvals_gene_mean_long <- mvals_gene_mean %>%
  pivot_longer(cols = -SYMBOL, names_to = "Sample", values_to = "M") 

mvals_genepromoter_mean_long <- mvals_genepromoter_mean %>%
  pivot_longer(cols = -SYMBOL, names_to = "Sample", values_to = "M") 

robz <- function(row){
  # Compute median and mean absolute deviation for row
  
  if(all(is.na(row))) { return (row) }
  
  m <- median(row,na.rm = T)
  s <- mad(row,na.rm = T)
  
  # If the MAD is 0, set it to a very small number
  if(s == 0){
    s <- 1E-100
  }
  robzscore <- (row - m) / (s)
  return(robzscore)
}

# Get just the promoter regions and average per gene per sample
# Then get Z score
mean_mval_lcpm_z <- mvals_gene_mean_long %>%  dplyr::rename(gene_mean_m=M) %>% 
  left_join(
    mvals_genepromoter_mean_long %>%  dplyr::rename(genepromoter_mean_m=M)) %>% 
  inner_join(lcpms %>% select(-EnsGeneID)) %>% dplyr::rename(A5_ID=Sample) %>% 
  group_by(SYMBOL) %>%
  mutate(gene_mval_z = (gene_mean_m - mean(gene_mean_m)) / sd(gene_mean_m),
         genepromoter_mval_z = (genepromoter_mean_m - mean(genepromoter_mean_m)) / sd(genepromoter_mean_m),
         log2_cpm_z = (log2_cpm - mean(log2_cpm)) / sd(log2_cpm),
         gene_mval_robustz = robz(gene_mean_m),
         genepromoter_mval_robustz = robz(genepromoter_mean_m),
         log2_cpm_robustz = robz(log2_cpm))


mean_mval_lcpm_z %>% 
  filter(abs(log2_cpm_z) >2 & abs(gene_mval_z) > 2) %>%  
  write_delim("results/tables/logcpm_plus_methylation_Z_gt2.csv", delim = "\t")

# Most things that are expressed seem to be mostly unmethylated 
ggplot(data = mean_mval_lcpm_z  %>% mutate(Outlier=ifelse(abs(log2_cpm_z) >3 & abs(genepromoter_mval_z) > 3, "Outlier", "No")), 
       aes(x = log2_cpm_z, y = genepromoter_mval_z, color=Outlier))+
  geom_point(size=0.3)+
  labs(x = "RNA-Seq Z score", y = "Mean array promoter\nmethylation Z score")+
  blank_theme+
  geom_hline(yintercept = 3, linetype = "dashed")+
  geom_hline(yintercept = -3, linetype = "dashed")+
  geom_vline(xintercept = 3, linetype = "dashed")+
  geom_vline(xintercept = -3, linetype = "dashed")+
  scale_color_manual(values=list(Outlier="black", No="grey")) +
  coord_equal()+
  guides(color=F) +
  ggsave("results/plots/expression_vs_methylation/logcpm_vs_methylation_Zscore.png", device = png())


# Get the contribution to outliers
contribution <- mean_mval_lcpm_z %>% filter(abs(log2_cpm_z) >2 & abs(gene_mval_z) > 2) %>% 
  group_by(A5_ID) %>% 
  dplyr::count() %>% 
  arrange(desc(n))



# Set up reusable tracks
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
#Ensure that the methylation data is ordered by chromosome and base position.
annEPICOrd <- annEPIC[order(annEPIC$chr,annEPIC$pos),]
bValsOrd <- bVals[match(annEPICOrd$Name,rownames(bVals)),]
# Create genomic ranges object from methylation data
cpgData <- GRanges(seqnames=Rle(annEPICOrd$chr),
                   ranges=IRanges(start=annEPICOrd$pos, end=annEPICOrd$pos),
                   strand=Rle(rep("*",nrow(annEPICOrd))),
                   betas=bValsOrd)
islandHMM <- read.csv(paste0(dataDirectory,
                             "/model-based-cpg-islands-hg19-chr17.txt"),
                      sep="\t", stringsAsFactors=FALSE, header=FALSE)
islandData <- GRanges(seqnames=Rle(islandHMM[,1]), 
                      ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                      strand=Rle(strand(rep("*",nrow(islandHMM)))))
dnase <- read.csv(paste0(dataDirectory,"/wgEncodeRegDnaseClusteredV3chr17.bed"),
                  sep="\t",stringsAsFactors=FALSE,header=FALSE)
dnaseData <- GRanges(seqnames=dnase[,1],
                     ranges=IRanges(start=dnase[,2], end=dnase[,3]),
                     strand=Rle(rep("*",nrow(dnase))),
                     data=dnase[,5])

gene_info <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)%>%
  data.frame()

plot_region <- function(gene, sample){
  
  # Get gene aliases
  gene <- select(org.Hs.eg.db,keys = gene, columns  = c("SYMBOL","ENTREZID"), keytype = "SYMBOL")
  
  location <- filter(gene_info, gene_id == gene$ENTREZID[1])
  
  chromosome <- location$seqnames
  
  if(location$strand == "-"){
    start <- location$end-500
    end <- location$end+5000
  }
  
  else{
    start <- location$start-5000
    end <- location$start+500
    
  }
  
  gen <- "hg19"
  # add 25% extra space to plot
  minbase <- start - (0.25*(end-start))
  maxbase <- end + (0.25*(end-start))
  
  df <- data.frame(chr=chromosome, start=minbase, end=maxbase,strand="*")
  ranges <- makeGRangesFromDataFrame(df) 
  
  iTrack <- IdeogramTrack(genome = gen, chromosome = chromosome)
  
  rTrack <- UcscTrack(genome=gen, chromosome=chromosome, track="NCBI RefSeq", 
                      from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                      rstarts="exonStarts", rends="exonEnds", gene="name", 
                      symbol="name2", transcript="name", strand="strand", 
                      fill="darkblue",stacking="squish", name="RefSeq", 
                      showId=TRUE, geneSymbol=TRUE)
  
  # extract data on CpGs in DMR
  cpgData_subset <- subsetByOverlaps(cpgData, ranges)
  
  groups <- left_join(targets, sample)
  
  # Set up the sample to plot
  groups <- factor(groups$TERT_mutant_or_expressed,levels = unique(groups$TERT_mutant_or_expressed))
  
  # Methylation data track
  methTrack <- DataTrack(range=cpgData_subset, genome = gen, groups = groups,
                         chromosome=chromosome, ylim=c(-0.05,1.05), col=pal,
                         type=c("a","p"), name="DNA Meth.\n(beta value)",
                         background.panel="white", legend=TRUE, cex.title=0.8,
                         cex.axis=0.8, cex.legend=0.8)
  
  # CpG island track
  islandTrack <- AnnotationTrack(range=islandData, genome=gen, name="CpG Is.", 
                                 chromosome=chromosome,fill="darkgreen")
  
  # DNaseI hypersensitive site data track
  dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                          type="gradient", chromosome=chromosome)
  
  # DMR position data track
  dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                              chromosome=chromosome,fill="darkred")
  
  tracks <- list(iTrack, gTrack, methTrack,rTrack)
  sizes <- c(2,2,4,3) # set up the relative sizes of the tracks
  
  # Plot the tracks
  #pdf("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Plots/Gviz plots/Gviz_test.pdf")
  plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
             add35=TRUE, grid=TRUE, lty.grid=3, sizes = sizes, length(tracks))
  #dev.off()
}

# Try and plot some of the genes
sample <- A5_clinical%>%
  dplyr::select(A5_ID, TERT_or_ATRX,TERT_expression)%>%
  mutate(TERT_mutant_or_expressed = "No")%>%
  mutate(TERT_mutant_or_expressed = replace(TERT_mutant_or_expressed, TERT_or_ATRX == "TERT", "TERT mutant"))%>%
  mutate(TERT_mutant_or_expressed = paste0("Mutated: ", TERT_mutant_or_expressed,", ", TERT_expression))%>%
  dplyr::select(A5_ID, TERT_mutant_or_expressed)

gene <- "TERT"

plot_region(gene = gene, sample = sample)

# Plot RNA gene expression vs methylation

####
# Effect of Copy-number on methylation?
####

library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tidyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(gridExtra)
mvals_save <- read.csv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/A5 methylation Avals.csv")
bVals2 <- read.csv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/A5 methylation Bvals.csv")

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

mvals_save.anno <- mvals_save %>% left_join(data.frame(annEPIC) %>% tibble::rownames_to_column("Probe") %>% select(Probe,chr,pos))
bvals_save.anno <- bVals2 %>% left_join(data.frame(annEPIC) %>% tibble::rownames_to_column("Probe") %>% select(Probe,chr,pos))
mvallist <- list()
bvallist <- list()
for (i in 2:96) {  
  s <- colnames(mvals_save.anno)[i]
  s <- gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0([0-9])_D$","\\1-\\2",s)
  s <- gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0[0-9]_D_([0-2])$","\\1-\\2",s)
  mvallist[[s]] <- mvals_save.anno[,c(1,97,98,i)]
  colnames(mvallist[[s]])[4] <- s
  
  s <- colnames(bvals_save.anno)[i]
  s <- gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0([0-9])_D$","\\1-\\2",s)
  s <- gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0[0-9]_D_([0-2])$","\\1-\\2",s)
  bvallist[[s]] <- data.frame(bvals_save.anno[,c(1,97,98,i)])
  colnames(bvallist[[s]])[4] <- s
}

purplecn.files <- list.files(path = "/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/WGS/Tool outputs/All UMCCRise outputs", 
                             recursive = T, pattern = ".purple.cnv",full.names = T)
purplecn.files <- purplecn.files[!grepl("work",purplecn.files)]
purplecn <- lapply(purplecn.files, read.delim)
snames <- basename(purplecn.files)
snames <- gsub("^(E[0-9]{3})__.+","\\1-1",snames)
snames <- gsub("^(E[0-9]{3})_([0-9]).+","\\1-\\2",snames)
snames <- gsub("^FL-2.+","E233-1",snames)
names(purplecn) <- snames

mbvals <- list()
for (s in names(mvallist))
{
  MGR <- GRanges(seqnames=mvallist[[s]]$chr,
                 ranges = IRanges(start=mvallist[[s]]$pos,
                                  end=mvallist[[s]]$pos),
                 Probe=mvallist[[s]]$Probe, 
                 Mval=mvallist[[s]][,s])
  
  BGR <- GRanges(seqnames=bvallist[[s]]$chr,
                 ranges = IRanges(start=bvallist[[s]]$pos,
                                  end=bvallist[[s]]$pos),
                 Probe=bvallist[[s]]$Probe, 
                 Bval=bvallist[[s]][,s])
  
  CNGR <- GRanges(seqnames=paste0("chr", purplecn[[s]]$X.chromosome),
                  ranges = IRanges(start=purplecn[[s]]$start,
                                   end=purplecn[[s]]$end),
                  CN=purplecn[[s]]$copyNumber)
  OVL <- findOverlaps(MGR,CNGR)
  mcols(MGR)[["CN"]] <- 0
  mcols(MGR)$CN[queryHits(OVL)] <- mcols(CNGR)$CN[subjectHits(OVL)]
  MGR <- data.frame(MGR) %>% rename(chr=seqnames) %>%  select(-width,-strand)
  MGR$Sample <- s
  mvallist[[s]] <- MGR
  
  OVL <- findOverlaps(BGR,CNGR)
  mcols(BGR)[["CN"]] <- 0
  mcols(BGR)$CN[queryHits(OVL)] <- mcols(CNGR)$CN[subjectHits(OVL)]
  BGR <- data.frame(BGR) %>% rename(chr=seqnames) %>%  select(-width,-strand)
  BGR$Sample <- s
  bvallist[[s]] <- BGR
  mbvals[[s]] <- MGR %>% full_join(BGR)
}
mbvals <- bind_rows(mbvals)

CheckSamples <- data.frame(Sample=c("E120-1","E121-1","E122-1","E122-1","E123-1","E126-1","E127-1","E135-1","E136-1","E143-1","E145-1","E145-1","E160-1","E160-1","E160-1","E171-1","E171-1","E188-1","E188-1"),
                           ROI=c("chr11p","chr11p","chr3q","chr11p","chr3q","chr3q","chr11p","chr11p","chr11p","chr14","chr3q","chr11p","chr3q","chr11p","chr14","chr11p","chr14","chr3q","chr11p"))
CheckSamples$Loss="Loss"
ggplot(mbvals %>% filter(Sample %in% CheckSamples$Sample) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
    filter(!is.na(ROI), 
           ,CN<5, CN>0.5)  %>%   
    inner_join(CheckSamples) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
    slice_sample(prop = 0.01) %>% 
    mutate(CN=factor(round(CN,0))), aes(x=Bval, color=Sample)) + 
  geom_density() + facet_wrap("chr") 
geom_jitter(size=0.05, alpha=0.1)



ggplot(mbvals %>% filter(Sample %in% CheckSamples$Sample, chr %in% c("chr3","chr11","chr14")) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
    filter(!is.na(ROI), 
           ,CN<5, CN>0.5)  %>%   
    left_join(CheckSamples, ) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
    #slice_sample(prop = 0.01) %>% 
    mutate(CN=factor(round(CN,0))), aes(x=Bval, color=Sample, linetype=Loss)) + 
  geom_density() + facet_grid(Loss~chr) 
geom_jitter(size=0.05, alpha=0.1)

ggplot(mbvals %>% filter(Sample %in% CheckSamples$Sample, chr %in% c("chr3","chr11","chr14")) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
    filter(!is.na(ROI), 
           ,CN<5, CN>0.5)  %>%   
    left_join(CheckSamples) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
    slice_sample(prop = 0.2) %>% 
    mutate(CN=factor(round(CN,0))), aes(x=Sample, y=Bval, color=Loss)) + 
  geom_jitter(size=0.05, alpha=0.1) + facet_wrap("chr")


mbvals %>% filter(Sample %in% CheckSamples$Sample, chr %in% c("chr3","chr11","chr14")) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
  filter(!is.na(ROI), 
         ,CN<5, CN>0.5)  %>%   
  left_join(CheckSamples) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
  pivot_wider(id_cols = c(Probe,chr,start), names_from=c(Sample), values_from=c(Mval,Loss)) %>% 
  slice_sample(prop = 0.2) 


plot.data <-mbvals %>% filter(Sample %in% CheckSamples$Sample, chr %in% c("chr3","chr11","chr14")) %>%  mutate(ROI=case_when(
  (as.character(chr)=="chr11" &	start < 53700000) ~ "chr11p",
  (as.character(chr)=="chr3" & start >	91000000) ~ "chr3q",
  as.character(chr)=="chr14" ~ "chr14")) %>% 
  filter(!is.na(ROI), 
         ,CN<5, CN>0.5)  %>%   
  left_join(CheckSamples) %>% mutate(Loss=replace_na(Loss,"No Loss")) %>% 
  mutate(ROI=ifelse((as.character(chr)=="chr14" & start >	100819176 &start < 101805942), "DLK-MEG3",ROI)) %>% 
  pivot_wider(id_cols = c(Probe,chr,start, ROI), names_from=c(Sample), values_from=c(Mval,Loss)) %>% 
  slice_sample(prop = 0.2) 
plot.list <- list()

for (s1 in c("E143-1", "E160-1", "E171-1"))
{
  for (s2 in c("E188-1", "E122-1", "E135-1"))
  {
    plot.list[[paste(s1,s2,sep="_")]] <- ggplot(plot.data,
                                                aes(x=!!sym(paste0("Mval_",s1)), y=!!sym(paste0("Mval_",s2)), color=gsub("No Loss/Loss","Loss/No Loss",paste(!!sym(paste0("Loss_",s1)),!!sym(paste0("Loss_",s2)),sep="/")))) + 
      geom_point(size=0.1, alpha=0.3) +
      facet_wrap("ROI", ncol=4) + labs(color="Combo") + guides(color = guide_legend(override.aes = list(size=5)))
  }
}

do.call("grid.arrange", plot.list)

bvals_save.anno %>% filter(Probe %in% c("cg12434587", "cg12981137")) %>%  #chr=="chr5", pos > 1294819, pos < 1295869
  pivot_longer(cols=c(-Probe,-chr,-pos), 
               names_to="A5_ID", values_to="bval") %>% 
  mutate(A5_ID=gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0([0-9])_D$","\\1-\\2",A5_ID),
         A5_ID=gsub("Group[0-9][0-9]?.(E[0-9]{3})_T0[0-9]_D_([0-2])$","\\1-\\2",A5_ID)) %>% 
  left_join(A5_clinical %>% select(A5_ID,`Assumed driver of metastais`)) %>% rename(Annotation=`Assumed driver of metastais`) %>% 
  mutate(Annotation=ifelse(grepl("E167",A5_ID),"E167",Annotation)) %>% 
  ggplot(aes(x=Probe,y=bval, color=Annotation)) + 
  geom_jitter(size=0.8, width=0.2) + 
  #geom_segment(x=1295321, xend=1295753, y=0.5, yend=0.5, size=2, color="orange") +
  #geom_vline(xintercept = 1295180, linetype=2) +
  #geom_vline(xintercept = c(1295228,1295250), linetype=3) +
  #geom_text(aes(label=A5_ID)) + 
  scale_color_manual(values = c("blue", "red", "green", "grey")) 

1295321-1295753
pos > 1294819, pos < 1295869


m_vals_z <- scale(t(m_vals), center = T, scale = T)
outlier_p <- m_vals_z > 2 
outlier_n <- m_vals_z < -2 
outliers <- bind_rows(data.frame(A5_ID=names(colSums(outlier_p)), n_outliers=colSums(outlier_p), direction="positive"),
                      data.frame(A5_ID=names(colSums(outlier_n)), n_outliers=colSums(outlier_n), direction="negative")) %>% 
  left_join(a5_anno.meth)  %>% group_by(A5_ID) %>% mutate(Total=sum(n_outliers)) %>% arrange(Total) %>% 
  mutate(A5_ID = factor(A5_ID,levels=unique(.$A5_ID))) %>% 
  mutate(percent=(Total/nrow(m_vals_z))*100) 

gg_pos_neg <- ggplot(outliers, 
                     aes(x=A5_ID,y=n_outliers, fill=Primary_Location_Base)) + 
  geom_col() + 
  facet_wrap("direction", nrow = 2, scales = "free_y") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=7)) 

gg_total <- ggplot(outliers, aes(x=A5_ID,y=Total, fill=Primary_Location_Base)) + 
  geom_col() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=7))

gg_percent <- ggplot(outliers, aes(x=A5_ID,y=percent, fill=Primary_Location_Base)) + 
  geom_col() + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=7))

nox <- theme(axis.text.x = element_blank(), axis.title.x = element_blank())
gg_pos_neg + nox + 
  gg_total + nox +
  gg_percent + plot_layout(ncol=1, guides = "collect")

contrasts_to_use <- NULL
gsea_result <- list()
kegg_offline=T
for (comparison in c("hn","genosampletype")) {
  message("Processing comparison: ", comparison)
  
  gsea_result[[comparison]] <- list()
  
  contrast_matrix <- get(paste0("contrast_matrix_", comparison))
  design_matrix <- get(paste0("design_matrix_", comparison))
  available_contrasts <- dimnames(contrast_matrix)$Contrasts
  
  if(!is.null(contrasts_to_use))
  {
    available_contrasts <- intersect(available_contrasts, contrasts_to_use)
  }
  
  for(contrast in available_contrasts){
    
    message("Processing contrast: ", contrast)
    
    gsea_result[[comparison]][[contrast]] <- list()
    
    contrast_membership <- contrastdesign_to_memberlist(contrast_name = contrast, 
                                                        contrast_matrix = contrast_matrix,
                                                        design_matrix = design_matrix)
    
    participants <- contrast_membership[contrast_membership$group != "non_participant",]
    
    m_vals_active <- m_vals[, participants$A5_ID]
    b_vals_active <- b_vals[, participants$A5_ID]
    
    if(!all(participants$A5_ID == colnames(m_vals_active)) | 
       !all(participants$A5_ID == colnames(b_vals_active))) {
      stop("Methylation data column order does not match expected sample order")
    }
    
    m_vals_active_nc <- bind_rows(data.frame(m_vals_active, 
                                             check.names = F),
                                  data.frame(a5_methylation_neg_cont_probes[, participants$A5_ID], 
                                             check.names = F))
    
    is_neg_ctl_probe <- rownames(m_vals_active_nc) %in% rownames(a5_methylation_neg_cont_probes)
    message("Using ", sum(is_neg_ctl_probe), " primary control probes")
    
    
    rfit_control1 <- RUVfit(Y = m_vals_active_nc, X = participants$group, ctl = is_neg_ctl_probe) # Stage 1 analysis
    rfit_control1_adj <- RUVadj(Y = m_vals_active_nc, fit = rfit_control1)
    
    suffix <- dimnames(rfit_control1_adj$R)[[2]][[1]]
    pBH_control2_cutoff <- 0.5
    top_control1_adj <- topRUV(rfit_control1_adj, number=Inf, p.BH = 1)
    is_low_p_probe <- rownames(m_vals_active) %in% rownames(top_control1_adj[top_control1_adj[[paste0("p.BH_",suffix)]] > pBH_control2_cutoff,])
    message("Identified ", table(is_low_p_probe)[["TRUE"]], " secondary control probes with and p-adj above ", pBH_control2_cutoff)
    
    rfit_control2 <- RUVfit(Y = m_vals_active, X = participants$group, ctl = is_low_p_probe) # Stage 2 analysis
    rfit_control2_adj <- RUVadj(Y = m_vals_active, fit = rfit_control2)
    # Look at table of top results
    top_control2_adj <- topRUV(rfit_control2_adj,  number=Inf, p.BH = 1)
    
    gsea_result[[comparison]][[contrast]][["top_control2_adj"]] <- top_control2_adj
    
    grp_a <- unique(participants$group)[[1]]
    grp_b <- unique(participants$group)[[2]]
    
    gsea_bhp_threshold <- 0.01
    min_beta_delta <- 0.25
    gsea_max_cpgs <- 5000
    gsea_min_cpgs <- 1000
    
    b_vals_active <- b_vals_active[match(rownames(top_control2_adj),rownames(b_vals_active)),]
    mean_beta_a <- rowMeans(b_vals_active[,participants$group==grp_a])
    mean_beta_b <- rowMeans(b_vals_active[,participants$group==grp_b])
    beta_delta <- mean_beta_a - mean_beta_b
    sigDM <- top_control2_adj[[paste0("p.BH_",suffix)]] < gsea_bhp_threshold & abs(beta_delta) > min_beta_delta
    message("Identifed ", sum(sigDM), " significant CpGs for GSEA using p-adj threshold of ", 
            gsea_bhp_threshold, " and minimum delta beta of ", min_beta_delta) 
    selection_mode <- "threshold"
    
    sigCpGs <- names(sigDM)[sigDM]
    
    if(length(sigCpGs) > gsea_max_cpgs)
    {
      message("Number of significant CpGs exceeds threshold (", gsea_max_cpgs,"). Using first ", 
              gsea_max_cpgs," from top table with p-adj below ", gsea_bhp_threshold," ranked by delta-beta")
      
      #CpGs below FDR threshold
      cpgs_pthreshold <- rownames(top_control2_adj)[top_control2_adj[[paste0("p.BH_",suffix)]] < gsea_bhp_threshold]
      
      #CpGs above delta threshold, rank by delta
      cpgs_deltathreshold <- sort(beta_delta[abs(beta_delta) > min_beta_delta], decreasing = T)
      
      #Delta ranked CpGs filtered for those below FDR threshold 
      sig <- cpgs_deltathreshold[names(cpgs_deltathreshold) %in% cpgs_pthreshold]
      
      #limit to max CpGs
      sigCpGs <- names(sig[1:min(length(sigDM),gsea_max_cpgs)])
      
      selection_mode <- paste0("toptable_top",gsea_max_cpgs,"_padj",gsea_bhp_threshold, "_ranked_delta")
    }
    
    if(length(sigCpGs) < gsea_min_cpgs)
    {
      message("Number of significant CpGs is below threshold (", 
              gsea_min_cpgs,"). Using all CpGs below p-adj threshold of ", 
              gsea_bhp_threshold)
      
      
      sigCpGs <- rownames(topRUV(rfit_control2_adj, number=gsea_max_cpgs, p.BH = gsea_bhp_threshold))
      
      message("Found ", length(sigCpGs), " below p-adj threshold of ", gsea_bhp_threshold)
      
      selection_mode <- paste0("toptable_padj",gsea_bhp_threshold)
    }
    
    #Record CpGs and selection method
    gsea_result[[comparison]][[contrast]][["sigCpGs"]] <- sigCpGs
    gsea_result[[comparison]][[contrast]][["selection_mode"]] <- selection_mode
    
    check <- getMappedEntrezIDs(sig.cpg = sigCpGs)
    message("Significant CpGs map to ", length(check$sig.eg), " genes")
    
    
    gsea_result[[comparison]][[contrast]][["gsa"]][["body"]] <- list()
    gsea_result[[comparison]][[contrast]][["gsa"]][["promoter"]] <- list()
    
    for (current_collection in c("GO", "KEGG")) #, 
    {
      
      gsea_result[[comparison]][[contrast]][["gsa"]][["promoter"]][[current_collection]] <-
        gometh(
          sig.cpg = sigCpGs,
          all.cpg = rownames(top_control2_adj),
          collection = current_collection,
          plot.bias = TRUE,
          genomic.features = c("TSS200",
                               "TSS1500",
                               "1stExon"),
          offline = kegg_offline,
          offline_cache_dir = "/g/data/pq08/projects/ppgl/a5/methylation/offline_cache"
          
        )
      
      gsea_result[[comparison]][[contrast]][["gsa"]][["body"]][[current_collection]] <-
        gometh(
          sig.cpg = sigCpGs,
          all.cpg = rownames(top_control2_adj),
          collection = current_collection,
          plot.bias = TRUE,
          genomic.features =  c("Body"),
          offline = kegg_offline,
          offline_cache_dir = "/g/data/pq08/projects/ppgl/a5/methylation/offline_cache"
        )
    }
  }
}
saveRDS(gsea_result, "./methylation/quickload_checkpoints/gsea_result.rds")

abdo_thoracic_include <- design_matrix_hn[,"Adrenal"] | design_matrix_hn[,"Extraadrenal"]
head_neck_include <- as.logical(design_matrix_hn[,"Head_neck"])

group <- ifelse(abdo_thoracic_include, "abdo_thoracic",  ifelse(head_neck_include, "head_neck", "exclude"))

group <- group[abdo_thoracic_include | head_neck_include]

m_vals_hnat <- m_vals[, abdo_thoracic_include | head_neck_include]

m_vals_hnat_nc <- bind_rows(data.frame(m_vals_hnat, 
                                       check.names = F),
                            data.frame(a5_methylation_neg_cont_probes[, abdo_thoracic_include | head_neck_include], 
                                       check.names = F))
# create vector marking negative controls in data matrix
is_neg_ctl_probe <- rownames(m_vals_hnat_nc) %in% rownames(a5_methylation_neg_cont_probes)



rfit1 <- RUVfit(Y = m_vals_hnat_nc, X = group, ctl = is_neg_ctl_probe) # Stage 1 analysis
rfit2 <- RUVadj(Y = m_vals_hnat_nc, fit = rfit1)

top1 <- topRUV(rfit2, num=Inf, p.BH = 1)
ctl2 <- rownames(m_vals_hnat) %in% rownames(top1[top1$p.BH_X1.head_neck > 0.5,])
table(ctl2)

rfit3 <- RUVfit(Y = m_vals_hnat, X = group, ctl = ctl2) # Stage 2 analysis
rfit4 <- RUVadj(Y = m_vals_hnat, fit = rfit3)
# Look at table of top results
topRUV(rfit4)

Madj <- getAdj(Y = m_vals_hnat, fit = rfit3) # get adjusted values


par(mfrow=c(1,2))
plotMDS(m_vals_hnat, labels=names(group), col=as.integer(factor(group)),
        main="Unadjusted", gene.selection = "common", dim.plot = c(1,2))
legend("right",legend=c("Cancer","Normal"),pch=16,cex=1,col=1:2)
plotMDS(Madj, labels=names(group), col=as.integer(factor(group)),
        main="Adjusted: RUV-inverse", gene.selection = "common", dim.plot = c(1,2))



probe_data_b <- b_vals[rownames(topRUV(rfit4)),] %>% data.frame(check.names = F) %>% 
  tibble::rownames_to_column("probe_id") %>% 
  pivot_longer(cols=-probe_id, names_to = "A5_ID", values_to = "b_val")

probe_data_m <- m_vals[rownames(topRUV(rfit4)),] %>% data.frame(check.names = F) %>% 
  tibble::rownames_to_column("probe_id") %>% 
  pivot_longer(cols=-probe_id, names_to = "A5_ID", values_to = "m_val")

probe_data <-  probe_data_b %>% inner_join(probe_data_m) 



plot.data <- probe_data %>%  
  inner_join(a5_anno.meth %>%  select(A5_ID, Primary_Location_Base, differential_group)) %>% 
  inner_join(data.frame(A5_ID=names(group), group=group)) %>% mutate(group=replace_na(group, "other"))

ggplot(plot.data, aes(x=group, y=b_val)) + 
  geom_boxplot() +
  geom_jitter(mapping=aes(color=differential_group, shape=Primary_Location_Base), width = 0.2) + 
  scale_color_manual(values = differential_group_colors) + facet_wrap("probe_id") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust =0.5, hjust=1))

jitter_width=0.4
plots_summary[[fit]][[contrast]][[i]] <- 
  ggplot(plot.data.summary, aes(x=participation, y=b_val)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_segment(mapping=aes(xend=participation, y=b_val+b_sd,yend=b_val-b_sd), 
               position = position_jitter(seed = 10, width = jitter_width), alpha=0.05) +
  scale_color_manual(values = differential_group_colors) + 
  geom_point(mapping=aes(color=differential_group, shape=Primary_Location_Base),
             position = position_jitter(seed = 10, width = jitter_width)) +
  ggrepel::geom_text_repel(aes(label=label), position = position_jitter(seed = 10, width = jitter_width)) +
  theme_bw() + 
  theme(axis.text.x = element_blank()) + 
  ggtitle(overlapping.genes, subtitle = contrast)
