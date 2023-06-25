# Modify purple output so that it can be used by GISTIC
library(tidyverse)
library(googlesheets4)

base_dir="/media/gadi/projects/A5"
#base_dir="/g/data/pq08/projects/A5"

purple_outputs <- list.files(paste(base_dir,"WGS/analysis/purple/",sep="/"), recursive = T,
                             full.names = T, pattern = "*.purple.segment.tsv")
names(purple_outputs) <- gsub(".purple.segment.tsv", "", basename(purple_outputs))

combined <- map(.x=purple_outputs, .f=read_tsv) %>% bind_rows(.id = "A5_ID")

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
  mutate(Seg.CN = log2(as.numeric(Seg.CN)+0.01)-1)

# Histogram of values. Should be centered around 0
hist(gistic_format$Seg.CN,breaks = 1000)

# Save a seg file to view in IGV
write_tsv(gistic_format,paste(base_dir,"WGS/analysis/gistic/input/A5_purple_segs_gistic_format.seg",sep="/"))

# Read in the A5 clinical data
gs4_auth("aidan.flynn@umccr-radio-lab.page")
a5_anno <- read_sheet("1hnXdXI29KvvuLxsaTBID1-EbE7mTSk6bG05HcSFfhgo", col_types = "c")
a5_anno <- a5_anno %>% dplyr::rename(`A5_ID`=`A5 ID`)
a5_anno <- a5_anno %>% filter(!is.na(`A5_ID`))

#Select non-malignant primaries
exclude_samples <- c(a5_anno %>% filter(is.na(PublicationID)) %>% pull(A5_ID), #quality excluded samples
                     "E124-1", "E145-1") #not SDHB samples

keep_samples.non_malignant_primary_nohn <- a5_anno %>% 
  filter(`Metastatic state of tumour sample` == "Primary", 
         `Tumour metastasised` == "No", 
         `Tumour is a head and neck` != "H_N",
         !(A5_ID %in% exclude_samples)) %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)

keep_samples.malignant_all_nohn <- a5_anno %>% 
  filter(`Metastatic state of tumour sample`=="Metastatic", 
         `Tumour metastasised`=="Yes", 
         `Tumour is a head and neck` != "H_N",
         !(A5_ID %in% exclude_samples)) %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)

keep_samples.malignant_tert_atrx_nohn <- a5_anno %>% 
  filter(`Metastatic state of tumour sample`=="Metastatic", 
         `Tumour metastasised`=="Yes", 
         `Tumour is a head and neck` != "H_N",
         `Assumed driver of metastasis`!="Unknown",
         !(A5_ID %in% exclude_samples)) %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)

# Make a file for gistic from the PURPLE outputs. 
for (g in c("keep_samples.non_malignant_primary_nohn", "keep_samples.malignant_all_nohn", "keep_samples.malignant_tert_atrx_nohn"))
{
  group_name <- gsub("keep_samples.","",g)
  group_members <- get(g)
  
  member_data <- gistic_format %>%  filter(A5_ID %in% group_members)
  
  write_tsv(gistic_format,paste0(base_dir, "/WGS/analysis/gistic/input/","A5_purple_segs_gistic_format_", group_name,".seg"))  
}

