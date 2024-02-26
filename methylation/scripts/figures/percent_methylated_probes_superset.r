######################################
# This script generates a box-jitter #
# plot of the proportion of probes   #
# above a beta-value cutoff across   #
# all samples in the methylation     #
# superset                           #
#                                    #
# Author: Aidan Flynn                #
# Date: 11/08/2023                   #
######################################

setwd("/g/data/pq08/projects/ppgl/")

##################
# Color Definitions
##################

#Color pallets and themes
source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")


############################################
# Import data loaders and run data mergers #
############################################

#clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#Methylation
source("./a5/tertiary/scripts/data_mergers/combine_tcga_comete_a5_methylation_data.r")

#Convert m_vals to beta
meth_tcga_comete_a5_450k_bval.batch_removed <- lumi::m2beta(meth_tcga_comete_a5_450k_mval.batch_removed)

####################################
# Mark and count methylated probes #
####################################

#Mark probes above cutoff
high_meth_cutoff <-  0.7
high_meth_probes <- meth_tcga_comete_a5_450k_bval.batch_removed > high_meth_cutoff

#Count probes above cutoff
n_high_meth <- colSums(high_meth_probes)

############
# Plotting #
############

plot_data <- tibble(Sample=names(n_high_meth), 
                    n_probes=n_high_meth, 
                    pcnt_probes=(n_high_meth/nrow(meth_tcga_comete_a5_450k_bval.batch_removed))*100) %>% 
  inner_join(meth_tcga_comete_a5_450k_anno) %>% 
  mutate(subtype=dplyr::recode(new_naming,
                               "A5 - Extraadrenal"="C1A1 (SDHx)",
                               "A5 - NF1"="C2A (Kinase)",
                               "A5 - Adrenal"="C1A1 (SDHx)",
                               "A5 - Head_neck"="C1A2 (SDHx-HN)",
                               "A5 - VHL"="C1B1 (VHL)",
                               "A5 - Extraadrenal_cardiac"="C1A1 (SDHx)",
                               "A5 - Extraadrenal_mediastinal"="C1A2 (SDHx-HN)",
                               "A5 - Unspecified"="C1A1 (SDHx)",
                               "C2B2 (MAML)"="C2B2 (MAML3)")) %>% 
  filter(!(subtype %in% c("Normal","C2C (Cortical admixture)"))) %>% 
  mutate(label = ifelse(Sample %in% c("E229-1","E229-2"), Sample, NA))

ggplot(plot_data, 
       aes(x=subtype, 
           y=pcnt_probes)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(mapping=aes(color=Dataset),
              width=0.2) + 
  theme_bw() +
  geom_text(mapping = aes(label=label)) +
  theme(axis.text.x = element_text(angle=45, 
                                   hjust=1, 
                                   vjust=1), 
        axis.title.x = element_blank()) +
  scale_color_manual(values=c(TCGA=ColorPalette[["LightGreen1"]],
                              `E-MTAB-733`=ColorPalette[["DarkBlue2"]],
                              A5=ColorPalette[["DarkRed2"]])) +
  ylab(paste("Probes with beta > ", high_meth_cutoff, "(%)")) 



t.test(pcnt_probes~subtype , plot_data %>%  filter(subtype %in% c("C1A1 (SDHx)","C1A2 (SDHx-HN)"))) 
