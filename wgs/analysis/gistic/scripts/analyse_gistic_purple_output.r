library(tidyverse)
library(ComplexHeatmap)
library(ggplot2)

# Read in GISTIC results at the confidence level that I set (0.90)
# Explanation of options: brlen if the proportion of a chromosome required to be altered to be considered a 'broad peak'
# -genegistic 1 means to calculate the significance of deletions at a gene level instead of a marker level. 
# Conf was increased to 90 as opposed to default 75 (more stringent)
# gcm: Method for reducing marker-level copy number data to the gene-level copy number data in the gene tables. 
# Markers contained in the gene are used when available, otherwise the flanking marker or markers are used. 
# The extreme method chooses whichever of min or max is furthest from diploid. 
# -broad 1 means do an additional broad level analysis
# -rx 0 means do not remove sex chromosome data
# GISTIC user guide: 
# ftp://ftp.broadinstitute.org/pub/genepattern/modules_public_server_doc/GISTIC2.pdf
# See this blog post for more info 
# http://crazyhottommy.blogspot.com/2017/11/run-gistic2-with-sequenza-segmentation.html

# Raw copy numbr values used were purple: 'Fitted absolute copy number of segment adjusted for purity and ploidy'

# Gistic command
#$GIS_folder/gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 1 -smallmem 0 -broad 1 -brlen 0.5 \
# -conf 0.90 -armpeel 1 -savegene 1 -gcm extreme -rx 0

setwd("/g/data/pq08/projects/ppgl/")

#######################
# Clinical Annotation #
#######################

source(paste0(basedir,"/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r"))
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

###############
# GISTIC Data #
###############

##########
# GISTIC Narrow Peaks 
##########

gistic_output_dir="./a5/wgs/analysis/gistic/output"

gistic_score_files <- list.files(gistic_output_dir, pattern="scores.gistic", recursive = T, full.names = T)
names(gistic_score_files) <- basename(gsub("/scores.gistic","",gistic_score_files))

gistic_scores <- map(.x=gistic_score_files, .f=read_tsv) %>% bind_rows(.id = "Group")


ggplot(gistic_scores %>% 
         group_by(Chromosome) %>% 
         mutate(ChrStart=min(Start), ChrEnd=max(End)) %>% 
         filter(Start > ChrStart + 10^6, End < ChrEnd-10^6, Type=="Del"), #, Group %in% c("malignant_all", "non_malignant_primary_nohn")
       aes(x=Start, xend=Start, y=0,yend=`-log10(q-value)`, color=Group)) + 
  geom_segment( alpha=0.5) + 
  geom_hline(yintercept = -log10(0.05), linetype=2, color="orange") +
  ylab("-log10(q-value)`") +
  facet_grid(Group~Chromosome, scales = "free_x") + theme_bw() + 
  theme(panel.grid = element_blank(), panel.spacing = unit(0, "lines"), axis.text.x = element_text(angle=90)) + 
  scale_color_manual(values=c("red","black", "darkgreen", "darkblue")) + coord_cartesian(ylim=c(0,8))

gistic_scores %>% filter(Group %in% c("malignant_tert_atrx", "non_malignant_primary_nohn"), Type=="Amp") %>% 
  pivot_wider(id_cols = c(Chromosome, Start, End), names_from = Group, values_from = `-log10(q-value)`)

##########
# GISTIC Broad Peaks 
##########

broad_peaks <- read.delim(file.path(gistic_output_dir,"broad_values_by_arm.txt"))

broad_peaks <- broad_peaks %>% 
  pivot_longer(cols = -`Chromosome Arm`, names_to="A5_ID", values_to = "value") %>% 
  mutate(A5_ID = gsub("-T0", "-", A5_ID))

broad_peaks <- broad_peaks %>% 
  mutate(value_catagory = 
           case_when(
             value > 0.25 ~ "Gain",
             value < -0.25 ~ "Loss",
             TRUE ~ "Neutral"
           ))

broad_peaks <- broad_peaks %>% inner_join(a5_anno)

######################
# ANDREW LEGACY CODE #
######################

# Ggplot2 blank theme
blank_theme <- theme_bw(base_size = 15)+
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

# Have a look at the broad changes in the context of TERT, ATRX ane metastasis

# Read in the A5 clinical data
gs4_auth("aidan.flynn@umccr-radio-lab.page")
a5_anno <- read_sheet("1hnXdXI29KvvuLxsaTBID1-EbE7mTSk6bG05HcSFfhgo", col_types = "c")
a5_anno <- a5_anno %>% dplyr::rename(`A5_ID`=`A5_ID`)
a5_anno <- a5_anno %>% filter(!is.na(`A5_ID`))

a5_anno <-  a5_anno %>%  mutate(TERT_and_ATRX = ifelse(TERT_or_ATRX %in% c("TERT", "ATRX"), "T or A", "Non T or A"))%>%
  mutate(Non_met_primary = ifelse(TERT_or_ATRX %in% c("Non_met_primary", "Short_follow_up_primary","Unknown_driver_relapse"), "Non met primary", "TERT/ATRX/Met"))

# This script is for exploring and analysing GISTIC broad and focal results
# All lesions file (this file is the focal lesions GISTIC has called)
all_lesions <- read_tsv(paste(gistic_output_dir,"all_samples", "all_lesions.conf_90.txt", sep="/"))
colnames(all_lesions) <- gsub("-T0","-",colnames(all_lesions))

# Gistic CNV values for every gene and every sample
all_data_by_genes <- read_tsv(paste(gistic_output_dir,"all_samples", "all_data_by_genes.txt", sep="/"))

# Read in the broadly significant regions by chromosome arm
broad_signif <- read_tsv(paste(gistic_output_dir,"all_samples", "broad_significance_results.txt", sep="/"))

# Broad chromosome changes by chromosome arm for the full cohort
broad_vals_by_arm_all <- read_tsv(paste(gistic_output_dir,"all_samples", "broad_values_by_arm.txt", sep="/")) %>%
  gather(key, value,-`Chromosome Arm`)%>%
  filter(value != 0)%>%
  dplyr::select(key,value,`Chromosome Arm`)%>%
  mutate(`Gainloss` = ifelse(value >0 ,"Gain","Loss")) %>%
  unite(`Arm change`,c(`Chromosome Arm`,Gainloss), sep = "-")%>%
  dplyr::rename(A5_ID = key,`Arm change` = `Arm change`, `GISTIC copy number change` = value)%>%
  group_by(A5_ID) %>%
  summarise(`GISTIC arm changes`= paste(`Arm change`, collapse=", "),
            `GISTIC copy number changes` = paste(`GISTIC copy number change`, collapse=", "))
colnames(broad_vals_by_arm_all)

# Broad chromosome changes by chromosome arm for the full cohort in a format that is easy to analyse
include_samples <- a5_anno %>% filter(Exclude=="N") %>% pull(A5_ID)
broad_vals_by_arm_all <- read_tsv(paste(gistic_output_dir,"all_samples", "broad_values_by_arm.txt", sep="/"))
colnames(broad_vals_by_arm_all) <- gsub("-T0","-",colnames(broad_vals_by_arm_all))
broad_vals_by_arm_all <- broad_vals_by_arm_all[,colnames(broad_vals_by_arm_all) %in% c(include_samples, "Chromosome Arm")]

chrmat <- as.matrix(broad_vals_by_arm_all[,2:ncol(broad_vals_by_arm_all)])
rownames(chrmat) <- broad_vals_by_arm_all$`Chromosome Arm`

# Get rid of chr X as it is annotated as 'missing' in men
chrmat <- chrmat[!grepl("X", rownames(chrmat)),]

colours <- c("TERT" = "red","ATRX" =  "orange","Unknown_met_driver_met" = "darkred","Unknown_met_driver_primary" = "purple",
             "Short_follow_up_primary" = "limegreen","Non_met_primary" = "darkgreen", "Unknown_driver_relapse"="yellow")

met_colours = c("No" = "lime green", "Short follow up" = "light green", "Yes" = "red","No Data"="grey")

anno <- a5_anno %>% filter(A5_ID %in% colnames(chrmat)) %>% 
  mutate(TERT_or_ATRX = factor(TERT_or_ATRX, 
                               levels = c("TERT", "ATRX", "Unknown_met_driver_met", 
                                          "Unknown_met_driver_primary","Short_follow_up_primary", "Non_met_primary", "Unknown_driver_relapse")))

# Top annotation
gen_anno1 <- data.frame(T_AT_stat = anno$TERT_or_ATRX, Met_status = anno$`Tumour metastasised`)
top_ano <- HeatmapAnnotation(df = gen_anno1, col = list(T_AT_stat = colours, Met_status = met_colours))

pdf("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/WGS/Plots/Copy_number/All samples copy number.pdf", width = 16, height = 8)
Heatmap(chrmat, top_annotation = top_ano,cluster_rows = F, column_split = anno$TERT_or_ATRX, column_title_rot = 90, cluster_columns = F,
        cluster_column_slices = T, name = "PURPLE copy\nnumber change")
dev.off()

# Show that chr 11 is associated with TERT
all_q <- broad_vals_by_arm_all %>%
  gather(`A5_ID`, arm_change, -`Chromosome Arm`)%>%
  left_join(a5_anno)%>%
  # Only keep arms that changed more than a few times
  filter(!grepl("X",`Chromosome Arm`))

hist(all_q$arm_change, breaks = 50)

pq_list <- list()
# Loop t tests for chromosome arms
for(i in 1:length(unique(all_q$`Chromosome Arm`))){
  
  pq <- unique(all_q$`Chromosome Arm`)[i]
  arm <- filter(all_q, `Chromosome Arm` == pq)%>%
    dplyr::select(Non_met_primary, arm_change)%>%
    mutate(Gainloss = ifelse(arm_change > 0, "Gain", "Loss"))%>%
    mutate(Gainloss = replace(Gainloss, arm_change  == 0, "No change"))%>%
    mutate(Gainloss = factor(Gainloss, levels = c("Gain", "Loss", "No change")))%>%
    dplyr::select(-arm_change)%>%
    table()%>%
    t()
  
  chisq <- chisq.test(arm)
  
  p <- chisq$p.value
  
  p_df <- data.frame(arm =pq, p = p)
  
  pq_list[[i]] <- p_df
  
}

all_pq <- bind_rows(pq_list)%>%
  mutate(FDR = p.adjust(p))%>%
  arrange(-FDR)%>%
  filter(p < 0.05)
  #%>% write_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/WGS/Outputs/Copy_number_statistical_tests/Non met vs met chromosome changes chi square.csv")

signif <- all_q %>%
  filter(`Chromosome Arm` %in% all_pq$arm)%>%
  dplyr::select(`A5_ID`, Non_met_primary, arm_change, TERT_or_ATRX, `Chromosome Arm`, `Metastatic state of tumour sample`)

ggplot(data = signif, aes(x = Non_met_primary, y = arm_change))+
  geom_violin(aes(fill = Non_met_primary))+
  geom_boxplot(outlier.alpha = 0, width =0.3, aes(fill = Non_met_primary))+
  facet_wrap(~`Chromosome Arm`, ncol = 1)+
  geom_jitter(aes(colour = TERT_or_ATRX), width =0.3)+
  #blank_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(x = "Tumour state", y = "Chromosome arm change\n(purtiy/plody adjusted)", colour = "TERT or ATRX status")+
  coord_flip()+
  scale_colour_manual(values =  colours)+
  scale_fill_manual(values = c("light green", "coral1"))
  #+  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/WGS/Plots/Copy_number/All sample arm level copy change vs met state.pdf", width =10, height =10)

# Have a look at if copy number changes in paired tumours
dups <- unique(all_q$`A5_ID`)[!grepl("-1",unique(all_q$`A5_ID`))]
dups <- gsub("-.*", "", dups)
# Remove E167 

paired_arm <- all_q %>%
  filter(`Patient ID` %in% dups)%>%
  filter(arm_change != 0)%>%
  mutate(`Metastatic state of tumour sample` = factor(`Metastatic state of tumour sample`, levels = c("Primary", "Recurrent", "Metastatic")))%>%
  arrange(`Metastatic state of tumour sample`)%>%
  mutate(ID_state = paste(`A5_ID`, `Metastatic state of tumour sample`, "n=", round(as.numeric(`Sample ploidy`),1)))%>%
  mutate(`Chromosome Arm` = factor(`Chromosome Arm`, levels = broad_vals_by_arm_all$`Chromosome Arm`))%>%
  mutate(ID_state = factor(ID_state, levels = unique(ID_state)))

#saveRDS(paired_arm, "/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/WGS/Intermediate_files/copy_number/Clonal_arm_changes.rds")

ggplot(paired_arm, aes(x = ID_state, y = `Chromosome Arm`, fill = arm_change))+
  facet_wrap(~`Patient ID`, scales = "free", nrow = 1)+
  geom_tile()+
  scale_fill_gradient2(low = "blue", high = "red")+
  #blank_theme+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
#+  ggsave("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/WGS/Plots/Copy_number/Clonal copy number.pdf", width = 20, height = 8)

# Have a look at the focal events called by gistic
lesions_order <- all_lesions%>%
  mutate(`q values` = as.numeric(`q values`))%>%
  arrange(`q values`)%>%
  filter(`q values` < 0.05)

# Make plots of the focal regions that gistic has called
for(i in 1:nrow(lesions_order)){
  sample_cols <- lesions_order[i,10:ncol(lesions_order)]
  info_cols <- lesions_order[i,1:9]
  
  sample_gathered <- sample_cols %>%
    gather(`A5_ID`, CNV)%>%
    left_join(a5_anno)%>%
    filter(!is.na(CNV))
  
  save <- paste0("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/WGS/Plots/Copy_number/Gistic_focal_", info_cols$Descriptor,".pdf")
  
  ggplot(data = sample_gathered, aes(x = Non_met_primary, y = CNV))+
    geom_violin(aes(fill = Non_met_primary))+
    geom_boxplot(outlier.alpha = 0, width =0.3, aes(fill = Non_met_primary))+
    geom_jitter(aes(colour = TERT_or_ATRX), width =0.3)+
    #blank_theme+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    labs(x = "Tumour state", y = "Chromosome arm change\n(purtiy/plody adjusted)", colour = "TERT or ATRX status")+
    coord_flip()+
    scale_colour_manual(values =  colours)+
    scale_fill_manual(values = c("light green", "coral1"))+
    ggtitle(paste0("Gistic lesion\n", info_cols$Descriptor, " ", gsub("\\(.*", "",info_cols$`Wide Peak Limits`)))
    #+ ggsave(save)
}

# Function to parse gistic outputs
fix_gisitc_genes <- function (gistic_file){
  
  # Look at the gene changes called by GISTIC
  deleted_or_amp_genes <- read_tsv(gistic_file, col_names = F)
  
  del_trans <- data.frame(t(deleted_or_amp_genes),stringsAsFactors = F)
  
  united <- unite(del_trans, new, 5:ncol(del_trans), sep=",")
  
  colnames(united) <- c(as.character(united[1,1:4]),"genes in wide peak")
  
  united <- united [-1,]%>%
    filter(cytoband != "NA")
  
  united$`genes in wide peak` <- gsub(",NA.*","",united$`genes in wide peak`)
  
  united$`genes in wide peak` <- gsub("\\[|\\]","",united$`genes in wide peak`)
  
  united <- united %>%
    mutate(`q value` = as.numeric(`q value`))%>%
    arrange(`q value`)
  
  return(united)
}

# Get a table of the focally deleted and amplified genes from the CNV profiles of the cohort
deleted_genes <- fix_gisitc_genes(paste(gistic_output_dir,"all_samples", "del_genes.conf_90.txt", sep="/")) %>%
  mutate(State = "Deleted")
amplified_genes <- fix_gisitc_genes (paste(gistic_output_dir,"all_samples", "amp_genes.conf_90.txt", sep="/"))%>%
  mutate(State = "Amplified")

amp_del <- bind_rows(deleted_genes, amplified_genes)
#%>% write_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/WGS/Outputs/Copy_number_tables/Purple unique sample amplification and deletion genes.csv")



