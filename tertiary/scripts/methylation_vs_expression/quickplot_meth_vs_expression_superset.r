setwd("/g/data/pq08/projects/ppgl/a5")

library(tidyverse)
library(ggplot2)
library(patchwork)

###########
# Imports #
###########

#Modified version of gometh from missmethyl to use an offline cache for KEGG pathways
source("./methylation/scripts/go_meth_offline.r")

#Plotting helper functions
source("./methylation/scripts/helpers/plot_methylation.r")

###############
# DataLoaders #
###############

#load clinical annotation
source("./sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#load methylation array data
#Methylation
source("./tertiary/scripts/data_mergers/combine_tcga_comete_a5_methylation_data.r")

#load WTS data
source("./tertiary/scripts/data_mergers/combine_tcga_flynn_a5_wts_data.r")

##############
# Sample Set #
##############

sample_anno <- wts_tcga_flynn_a5_anno %>% 
  filter(Sample %in% intersect(wts_tcga_flynn_a5_anno$Sample, meth_tcga_comete_a5_450k_anno$Sample)) %>%
  mutate(subtype=dplyr::recode(new_naming,
                               "A5 - Extraadrenal"="C1A1 (SDHx)",
                               "A5 - NF1"="C2A (Kinase)",
                               "A5 - Adrenal"="C1A1 (SDHx)",
                               "A5 - Head_neck"="C1A2 (SDHx-HN)",
                               "A5 - VHL"="C1B1 (VHL)",
                               "A5 - Extraadrenal_aortic"="C1A2 (SDHx-HN)",
                               "A5 - Unspecified"="C1A1 (SDHx)",
                               "C2B2 (MAML)"="C2B2 (MAML3)"),
         Malignancy=replace_na(Malignancy, "Unknown")) 

######################
# Compute M/B-values #
######################

# calculate M-values for statistical analysis
m_vals <- meth_tcga_comete_a5_450k_mval.batch_removed[, intersect(colnames(meth_tcga_comete_a5_450k_mval.batch_removed), sample_anno$Sample)]
b_vals <- lumi::m2beta(meth_tcga_comete_a5_450k_mval.batch_removed[, intersect(colnames(meth_tcga_comete_a5_450k_mval.batch_removed), sample_anno$Sample)])

#################
# Relabel genes #
#################

ensgid_symbol_lookup <- data.frame(ensgid_symbol=rownames(a5_wts_counts)) %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), remove=F, extra = "merge", sep="_") %>% 
  mutate(ensgid=gsub("[.].+$", "", ensgid))

rownames(wts_tcga_flynn_a5_lcpm.batch_removed) <- ensgid_symbol_lookup$symbol[match(rownames(wts_tcga_flynn_a5_lcpm.batch_removed), ensgid_symbol_lookup$ensgid)]


############
# Plotting #
############

plot_methylation_vs_expr(#gene_symbol = "MGMT",
                         probes= c("cg12434587","cg12981137","cg09993459", "cg04473030","cg13405636"),
                         samples_to_label = c("E169-1", "E169-2"),
                         label_outliers = F,
                         sample_annotation = (sample_anno %>% dplyr::select(Sample_ID=Sample, subtype, Sex, Malignancy, Genotype, PCPG, ATRX_mutation)), 
                         b_vals = b_vals, 
                         summarise_probes_by_region = F,
                         m_vals = m_vals,
                         log2_cpm =wts_tcga_flynn_a5_lcpm.batch_removed, 
                         plot_mode = "beta", 
                         array_type = "450K", 
                         colour_column = "subtype", 
                         grouping_column = "subtype",
                         shape_column = "Sex", 
                         colour_scale = subtype_cols 
)

plot_methylation(#gene_symbol = "MGMT", 
                 region = "chr10:131265043-131265625",
                 show_mean_summarised_only = F,
                 label_outliers = F,
                 sample_annotation = (sample_anno %>% dplyr::select(Sample_ID=Sample, subtype, Malignancy, Genotype, PCPG, ATRX_mutation)), 
                 b_vals = b_vals, 
                 m_vals = m_vals, 
                 plot_mode = "m", 
                 array_type = "450K", 
                 colour_column = "subtype", 
                 grouping_column = "subtype",
                 colour_scale = subtype_cols,
                 shape_column = "Malignancy")  + coord_cartesian(ylim=c(-8,2))
