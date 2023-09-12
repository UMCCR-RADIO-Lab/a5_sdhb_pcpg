# Modify purple output so that it can be used by GISTIC
library(tidyverse)
library(googlesheets4)

####################
# Data Directories #
####################

base_dir="/g/data/pq08/projects/ppgl/a5"
purple_dir=file.path(base_dir,"wgs/analysis/purple")
gistic_dir=file.path(base_dir,"wgs/analysis/gistic")

################
# Data Loaders #
################

source(paste0(base_dir, "/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")) 
if (!exists("a5_anno"))
{
  data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = TRUE)
}

#####################
# Read PURPLE input #
#####################

purple_outputs <- list.files(purple_dir, recursive = T,
                             full.names = T, pattern = "*.purple.segment.tsv")
names(purple_outputs) <- gsub(".purple.segment.tsv", "", basename(purple_outputs))

combined <- map(.x=purple_outputs, .f=read_tsv) %>% bind_rows(.id = "A5_ID")


#######################
# Format GISTIC input #
#######################

# Gisitc input format 
# Sample	Chrom	Start	Stop	#Mark	Seg.CN
#   S1	   1	   61735	77473	10	-0.1234
# Seg.CN = log2()-1 of copy number

# The purple copy number is:
# 'Fitted absolute copy number of segment adjusted for purity and ploidy'

# Make a .seg file
gistic_format <- combined %>%
  # Set number of marks to BafCount. That would be equivalent No. of spots on a copy number array
  dplyr::select(A5_ID, Chrom = chromosome, Start= start, Stop = end, `#Mark` = bafCount, Seg.CN = tumorCopyNumber)%>%
  # Remove some below zero copy numbers
  filter(Seg.CN >=0)%>%
  # Center around 0 for GISTIC
  mutate(Seg.CN = log2(as.numeric(Seg.CN)+0.01)-1) %>% 
  mutate(A5_ID = gsub("-T0","-", A5_ID))

# Histogram of values. Should be centered around 0
hist(gistic_format$Seg.CN,breaks = 1000)

###########
# Write all sample GISTIC input 
###########

# Save a seg file to view in IGV
write_tsv(gistic_format, file.path(gistic_dir, "input", "a5_purple_segs_gistic_format.seg"))


###########
# Write sub-group sample GISTIC input 
###########

bio_groups <- list()

bio_groups[["all_hnpgl"]] <- a5_anno %>% 
  filter(differential_group_anatomy == "Head_Neck") %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)

bio_groups[["non_metastatic_atpgl"]] <- a5_anno %>% 
  filter(`differential_group_sampletype_strict` == "Non-metastatic primary",
         differential_group_anatomy == "Abdominal_Thoracic") %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)


bio_groups[["metastases_atpgl"]] <- a5_anno %>% 
  filter(`differential_group_sampletype_strict` %in% c("Metastatic primary", "Metastatic local recurrence", "Metastasis"),
         differential_group_anatomy == "Abdominal_Thoracic") %>% 
  group_by(`Patient ID`) %>% slice_max(n=1, order_by = differential_group_sampletype_strict) %>% pull(A5_ID)


bio_groups[["tert_atpgl"]] <- a5_anno %>% 
  filter(TERT_ATRX_Mutation == "TERT",
         differential_group_anatomy == "Abdominal_Thoracic") %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)

bio_groups[["atrx_atpgl"]] <- a5_anno %>% 
  filter(TERT_ATRX_Mutation == "ATRX",
         differential_group_anatomy == "Abdominal_Thoracic") %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)


# Make a file for gistic from the PURPLE outputs. 
for (g in names(bio_groups))
{
  member_data <- gistic_format %>%  filter(A5_ID %in% bio_groups[[g]])
  
  write_tsv(member_data,file.path(gistic_dir, "input",paste0("a5_purple_segs_gistic_format_", g,".seg")))  
}

