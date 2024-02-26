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
  data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)
}

##################
# Load cytobands #
##################

#These cytoband match those in the hg38.UCSC.add_miR.160920.refgene.mat reference file used by GISTIC

# Cytobands - HG38
# if(!file.exists("/g/data/pq08/reference/GRCh38/cytobands_hg38.txt.gz")) {
#   download.file(url = "https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz",
#                 destfile ="/g/data/pq08/reference/GRCh38/cytobands_hg38.txt.gz")
# }
# 
# cytoband <- read.delim("/g/data/pq08/reference/GRCh38/cytobands_hg38.txt.gz", header=F)
# colnames(cytoband) <- c("chr","start","end","cytoband","stain")
# cytoband <- cytoband %>% filter(chr %in% paste0("chr", c(1:22,"X","Y")))
# 
# chr_mid_point <- cytoband %>%
#   filter(chr %in% paste0("chr", c(1:22,"X","Y"))) %>%
#   mutate(arm=stringr::str_extract(cytoband, "^[pq]")) %>%
#   filter(arm == "p") %>%
#   group_by(chr) %>% summarise(midpoint = max(end))

#####################
# Read PURPLE input #
#####################

purple_segments <- list.files(purple_dir, recursive = T,
                              full.names = T, pattern = "*.purple.cnv.somatic.tsv")
names(purple_segments) <- gsub(".purple.cnv.somatic.tsv", "", basename(purple_segments))

seg_combined <- map(.x=purple_segments, .f=read_tsv) %>% bind_rows(.id = "A5_ID")

purple_purity <- list.files(purple_dir, recursive = T,
                            full.names = T, pattern = "*.purple.purity.tsv")
names(purple_purity) <- gsub(".purple.purity.tsv", "", basename(purple_purity))

purity_combined <- map(.x=purple_purity, .f=read_tsv) %>% bind_rows(.id = "A5_ID")

seg_combined <- seg_combined %>% left_join(purity_combined %>% dplyr::select(A5_ID, ploidy))

# seg_combined_with_limits <- seg_combined %>%  inner_join(chr_mid_point, by=c("chromosome"="chr"))
# 
# #Trim segments so they don't cross the midpoint of chromosomes
# for (i in 1:nrow(seg_combined_with_limits)) {
#   if (seg_combined_with_limits[i,]$start < seg_combined_with_limits[i,]$midpoint &
#       seg_combined_with_limits[i,]$end > seg_combined_with_limits[i,]$midpoint) {
#     if (seg_combined_with_limits[i,]$midpoint - seg_combined_with_limits[i,]$start <
#         seg_combined_with_limits[i,]$end - seg_combined_with_limits[i,]$midpoint) {
#       seg_combined_with_limits[i,]$start = seg_combined_with_limits[i,]$midpoint
#     } else {
#       seg_combined_with_limits[i,]$end = seg_combined_with_limits[i,]$midpoint
#     }
#   }
# }

########
# Recenter WGD samples around diploid
########

seg_combined_wgd_norm <- seg_combined %>% 
  rowwise() %>% 
  mutate(copyNumber=case_when(
    ploidy > 3 ~ min(max(1,copyNumber-2),3),
    #ploidy > 3 & copyNumber < 0.8 ~ max(0,copyNumber-2),
    #ploidy > 3 & copyNumber > 1.2 ~ max(1,copyNumber-2),
    .default = copyNumber
  ))

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
gistic_format <- seg_combined_wgd_norm %>%
  # Set number of marks to BafCount. That would be equivalent No. of spots on a copy number array
  dplyr::select(A5_ID, Chrom = chromosome, Start= start, Stop = end, `#Mark` = depthWindowCount, Seg.CN = copyNumber)%>%
  # Remove some below zero copy numbers or no markers
  filter(Seg.CN >=0, `#Mark` > 0)%>%
  # Center around 0 for GISTIC
  mutate(Seg.CN = log2(as.numeric(Seg.CN)+0.01)-1) %>% 
  mutate(A5_ID = gsub("-T0","-", A5_ID)) 

# Histogram of values. Should be centered around 0
hist(gistic_format$Seg.CN,breaks = 1000)

###########
# Filter out excluded samples
###########

gistic_format <- gistic_format %>% filter(A5_ID %in% a5_anno$A5_ID)

###########
# Write all sample GISTIC input 
###########

# Save a seg file to view in IGV
write_tsv(gistic_format, file.path(gistic_dir, "input", "a5_purple_segs_gistic_format_all.seg"))


###########
# Write sub-group sample GISTIC input 
###########

bio_groups <- list()

bio_groups[["parasympathetic"]] <- a5_anno %>% 
  filter(differential_group_anatomy %in% c("Head_Neck","Mediastinum")) %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)

bio_groups[["sympathetic"]] <- a5_anno %>% 
  filter(differential_group_anatomy %in% c("Abdominal_Thoracic")) %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)


bio_groups[["non_metastatic_sympathetic"]] <- a5_anno %>% 
  filter(`differential_group_sampletype_strict` == "Non-metastatic primary",
         differential_group_anatomy == "Abdominal_Thoracic") %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)


bio_groups[["metastatic_sympathetic"]] <- a5_anno %>% 
  filter(`differential_group_sampletype_strict` %in% c("Primary (metastasis present)","Metastatic primary", "Local recurrence (metastasis present)", "Metastatic local recurrence", "Metastasis"),
         differential_group_anatomy == "Abdominal_Thoracic") %>% 
  group_by(`Patient ID`) %>% 
  slice_max(n=1, order_by = differential_group_sampletype_strict, with_ties = F) %>% 
  pull(A5_ID)


bio_groups[["tert"]] <- a5_anno %>% 
  filter(TERT_ATRX_Mutation == "TERT",
         differential_group_anatomy == "Abdominal_Thoracic") %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)

bio_groups[["atrx"]] <- a5_anno %>% 
  filter(TERT_ATRX_Mutation == "ATRX",
         differential_group_anatomy == "Abdominal_Thoracic") %>% 
  group_by(`Patient ID`) %>% slice_head(n=1) %>% pull(A5_ID)


# Make a file for gistic from the PURPLE outputs. 
for (g in names(bio_groups))
{
  member_data <- gistic_format %>%  filter(A5_ID %in% bio_groups[[g]])
  
  write_tsv(member_data,file.path(gistic_dir, "input",paste0("a5_purple_segs_gistic_format_", g,".seg")))  
}

