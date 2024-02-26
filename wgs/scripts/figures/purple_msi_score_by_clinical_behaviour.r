#############################################
# Script to generate a plot of MSI scores   # 
# from purple or MSI Sensor                 #
# Author: Aidan Flynn                       #
# Date: 18/12/2023                          #
# Languages: R                              #
#############################################

library(ggplot2)

source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")

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

#####
# MSI sensor
#####

msi_output_files <- list.files(path = "/g/data/pq08/projects/ppgl/a5/wgs/analysis/msi_sensor_pro", pattern = "*.msi.txt$", recursive = T, full.names = T)
names(msi_output_files) <- gsub(".msi.txt","",basename(msi_output_files))
msi_output <- purrr::map(msi_output_files,read.delim) %>% bind_rows(.id="A5_ID")
colnames(msi_output)[4] <- "pcnt_sites"

#####
# Merge data and annotation
#####

msi <- msi_output %>% left_join(purple_output %>% dplyr::select(A5_ID, msIndelsPerMb))

msi <- msi %>%  
  mutate(A5_ID = gsub("-T0","-", A5_ID)) %>% 
  left_join(a5_anno) %>% 
  filter(Exclude == "N") %>% 
  mutate(`Treatment exposure`=dplyr::case_match(Resection_post_dna_damaging_treatment,
                                                "No" ~ "No",
                                                "Uncertain" ~ "Uncertain",
                                                "Yes (I131-MIBG - 7.9 GBq)" ~ "I131-MIBG",
                                                "Yes (CVD - 2 Cycles)" ~ "CVD <10 cycles",
                                                "Yes (CVD - 14 Cycles)" ~ "CVD >10 cycles",
                                                "Yes (CVD - 26 Cycles)" ~ "CVD >20 cycles",                        
                                                "Possibly (CVD)" ~ "Uncertain",
                                                "Yes (CVD - 22 cycles;I131-MIBG - 5.3 Gbq)" ~ "CVD >20 cycles + I131-MIBG",  
                                                "Yes (I131-MIBG - 5.3 Gbq)"  ~ "I131-MIBG",                      
                                                "Yes (CVD - 23 cycles)" ~ "CVD >20 cycles", 
                                                "Yes (I131-MIBG - 33 GBq;Carboplatin - 3 cycles)"  ~ "I131-MIBG",
                                                "Yes (CVD - Unknown;I131-MIBG - Unknown)" ~ "CVD Unknown Cycles + I131-MIBG")) %>% 
  mutate(Resection_post_dna_damaging_treatment=case_when(
    grepl("Yes",Resection_post_dna_damaging_treatment) ~ "Yes",
    Resection_post_dna_damaging_treatment == "No" ~ "No",
    .default = "Uncertain"))


############
# Plotting #
############

# p1 <- ggplot(msi, aes(x=Resection_post_dna_damaging_treatment, y=pcnt_sites, color=treatment_type)) + 
#   geom_jitter(width = 0.15, height=0, alpha=0.5) + 
#   theme(axis.text.x = element_text(angle=90,  vjust=0.5, hjust = 1)) + 
#   scale_color_manual(values =  c25[c(2:5,8,17,18)])
# 
# 
# p2 <- ggplot(msi, aes(x=TERT_ATRX_Mutation, y=pcnt_sites, color=treatment_type)) + 
#   geom_jitter(width = 0.15, height=0, alpha=0.5) + 
#   theme(axis.text.x = element_text(angle=90,  vjust=0.5, hjust = 1)) + 
#   scale_color_manual(values =  c25[c(2:5,8,17,18)])
# 
# p3 <- ggplot(msi, aes(x=differential_group_sampletype_strict, y=pcnt_sites, color=treatment_type)) + 
#   geom_jitter(width = 0.15, height=0, alpha=0.5) + 
#   theme(axis.text.x = element_text(angle=90,  vjust=0.5, hjust = 1)) + 
#   scale_color_manual(values = c25[c(2:5,8,17,18)])
# 
# p1 + p2 + p3 + plot_layout(guides="collect")

ggplot(msi %>% 
         mutate(differential_group_sampletype_strict=dplyr::recode(
           differential_group_sampletype_strict,
           "Metastatic local recurrence"="Metastatic primary",
           "Local recurrence (metastasis reported)"="Primary (metastasis reported)",)), 
       aes(x=differential_group_sampletype_strict, y=msIndelsPerMb, color=`Treatment exposure`)) + 
  geom_jitter(width = 0.2, height=0) + 
  scale_color_manual(values = c25[c(2:5,8,16,25,21)]) +
  theme_bw() +
  geom_hline(yintercept = 4, linetype=2, color="red") +
  theme(axis.text.x = element_text(angle=90,  vjust=0.5, hjust = 1)) 

