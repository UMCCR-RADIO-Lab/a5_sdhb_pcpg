#############################################
# Script to generate a plot of 
# polyclonality scores from purple          #
# Author: Aidan Flynn                       #
# Date: 13/11/2024                          #
# Languages: R                              #
#############################################

library(ggplot2)

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#############
# Load data #
#############

#####
# Annotation
#####

data_loader_a5_clinical_anno("aidan.flynn@umccr-radio-lab.page", use_cache = T)

#####
# PURPLE
#####
purple_output_files <- list.files(path = "/g/data/pq08/projects/ppgl/a5/wgs/analysis/purple/", pattern = "*.purple.purity.tsv", recursive = T, full.names = T)
names(purple_output_files) <- gsub(".purple.purity.tsv","",basename(purple_output_files))
purple_output <- purrr::map(purple_output_files,read.delim) %>% bind_rows(.id="A5_ID")

plot_data <- purple_output %>% 
  mutate(A5_ID=gsub("-T0","-", A5_ID)) %>%  
  left_join(a5_anno %>%  dplyr::select(A5_ID, PublicationID, differential_group_sampletype_strict, cell_of_origin, TERT_ATRX_Mutation)) %>% 
  filter(!is.na(PublicationID)) %>% 
  filter(cell_of_origin != "Ambiguous")

ggplot(plot_data %>% 
         mutate(differential_group_sampletype_strict=dplyr::recode(
           differential_group_sampletype_strict,
           "Metastatic local recurrence"="Metastatic primary",
           "Local recurrence (metastasis reported)"="Primary (metastasis reported)",)), 
       aes(x=differential_group_sampletype_strict, y=polyclonalProportion)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(mapping = aes(color=TERT_ATRX_Mutation), width = 0.2, height=0) + 
  scale_color_manual(values = driver_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90,  vjust=0.5, hjust = 1)) +
  facet_grid(~cell_of_origin, space="free_x", scales = "free_x") 


##t-test

t_test_data <- plot_data %>% filter(cell_of_origin == "Chromaffin") 

t_result_NMP_vs_MP <- t.test(polyclonalProportion~differential_group_sampletype_strict,  
                              t_test_data %>% filter(differential_group_sampletype_strict %in% c("Non-metastatic primary","Metastatic primary")), alternative = c("greater"))  


t_result_NMP_vs_Met <- t.test(polyclonalProportion~differential_group_sampletype_strict,  
                                t_test_data %>% filter(differential_group_sampletype_strict %in% c("Non-metastatic primary","Metastasis")), alternative = c("greater"))  

