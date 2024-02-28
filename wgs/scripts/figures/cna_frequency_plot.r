library(ggplot2)
#library(patchwork)
library(plyranges)
library(tidyr)

################
# Data Loaders #
################

if(!exists("a5_anno"))
{
source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)
}

if(!exists("a5_seg_keep"))
{
source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_cna_segs()
}

if(!exists("ColorPalette"))
{
  source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")
}

####################
# Data preparation #
####################

# cna_simplify_coding <- c(
#   "Loss"="Loss",
#   "Loss"="Loss + Subclonal CNLOH",
#   "Subclonal Loss"="Subclonal Loss",
#   "Subclonal Loss"="Minor Subclonal Loss",
#   "CNLOH"="CNLOH",
#   "Diploid/Haploid-X"="Diploid/Haploid-X",
#   "Gain"="Gain+LOH",
#   "Gain"="Gain",
#   "Gain"="WGD+Gain",              
#   "Subclonal Gain"="Subclonal Gain",
#   "Other"="Hom. Del.",
#   "Diploid/Haploid-X"="Other",
#   "WGD"="WGD",
#   "Chromothripsis"="Chromothripsis")

cna_simplify_coding <- c(
  "Loss"~"Loss",
  "Loss + Subclonal CNLOH"~"Loss",
  "Subclonal Loss"~"Subclonal Loss",
  "Minor Subclonal Loss"~"Diploid/Haploid-X",
  "CNLOH"~"CNLOH",
  "Diploid/Haploid-X"~"Diploid/Haploid-X",
  "Gain+LOH"~"Gain",
  "Gain"~"Gain",
  "WGD+Gain"~"Gain",              
  "Subclonal Gain"~"Subclonal Gain",
  "Hom. Del."~"Other",
  "Other"~"Diploid/Haploid-X",
  "WGD"~"WGD",
  "Chromothripsis"~"Chromothripsis")

########
# Recenter WGD samples around diploid
########

a5_seg_keep_wgd_norm <- a5_seg_keep %>% 
  rowwise() %>% 
  mutate(
    copyNumber=ifelse(mean_copyNumber > 2.5, copyNumber - 2, copyNumber),
    minorAlleleCopyNumber=ifelse(mean_copyNumber > 2.5, max(0,minorAlleleCopyNumber - 1), minorAlleleCopyNumber),
    majorAlleleCopyNumber=ifelse(mean_copyNumber > 2.5 & majorAlleleCopyNumber > 1.5, max(0,majorAlleleCopyNumber - 1), majorAlleleCopyNumber),
    mean_copyNumber=ifelse(mean_copyNumber > 2.5, mean_copyNumber - 2, mean_copyNumber)
    ) %>% 
  ungroup() %>% 
  mutate(
      Class = classify_cna_event(minorAlleleCopyNumber,majorAlleleCopyNumber, Gender,chromosome, copyNumber, mean_copyNumber) #function sourced from cna_event_labeller.R
    ) %>%
  mutate(Class = factor(Class, levels = cn_event_types)) %>% 
  mutate(Class = case_match(Class, 
                            !!!cna_simplify_coding,
                            .default = Class)) 

for (cr in 1:nrow(chromothripsis_regions))
{
  a5_seg_keep_wgd_norm <- a5_seg_keep_wgd_norm %>% mutate(Class=ifelse(
    (A5_ID==chromothripsis_regions$A5_ID[[cr]] | A5_ID==gsub("T0", "", chromothripsis_regions$A5_ID[[cr]])) & 
      (chromosome==chromothripsis_regions$chromosome[[cr]] & 
         ((start >= chromothripsis_regions$start[[cr]]) & (end <= chromothripsis_regions$end[[cr]]))), "Chromothripsis", as.character(Class)))
}

######
# Make windows 
######

chrom_gr <- GRanges(seqnames = chr_offsets$chromosome,  ranges = IRanges(start = 1, end = chr_offsets$chr_size))
chrom_windows <- tile(chrom_gr,width = 1000000)
names(chrom_windows) <- purrr::map_chr(.x=chrom_windows, .f=\(x) { as.character(seqnames(x))[1] })

######
# Prepare CNA 
######

seg_gr <- a5_seg_keep_wgd_norm %>%  
  group_by(chromosome) %>%  
  group_split() %>% 
  purrr::set_names(purrr::map_chr(., ~.x$chromosome[1])) %>% 
  purrr::map(.f = \(segs) {
    return(GenomicRanges::makeGRangesFromDataFrame(segs %>% dplyr::select(chromosome, start,end, A5_ID, Class), 
                                        start.field = "start", 
                                        end.field = "end", 
                                        seqnames.field = "chromosome",
                                        keep.extra.columns = T)) 

  })
  
#sort order
seg_gr <- seg_gr[paste0("chr", c(1:22,"X"))]

#order sanity check
if (!all(names(seg_gr) == names(chrom_windows))) { stop("Name sanity check failed") }

#######################
# Join windows to CNA #
#######################


seg_window_overlap <- purrr::map2(.x = chrom_windows, 
            .y = seg_gr, 
            .f = \(windows, segs) {
              
              return(windows %>% plyranges::join_overlap_inner(segs))
              
            }) 

#####################
# Summarise windows #
#####################

seg_window_overlap_df <- 
  purrr::map(.x=seg_window_overlap, .f=GenomicRanges::as.data.frame) 


count_members <- function(chrom_segs, grouping_feature, permitted_anatomical_groups, permitted_cell_of_origin) {
  segs_annotated <- 
    chrom_segs %>%
    mutate(Class=factor(as.character(Class, levels=cn_event_types)),
           Class=forcats::fct_drop(Class)) %>% 
    dplyr::select(seqnames,  start, end, A5_ID, Class) %>% 
    inner_join(a5_anno %>%  
                 dplyr::select(A5_ID, differential_group_anatomy, cell_of_origin, !!sym(grouping_feature))) 
  
  if(!is.null(permitted_anatomical_groups))  {
    segs_annotated <-   segs_annotated %>% filter(differential_group_anatomy %in% permitted_anatomical_groups) }
  
  if(!is.null(permitted_cell_of_origin))  {
    segs_annotated <-   segs_annotated %>% filter(cell_of_origin %in% permitted_cell_of_origin) }
    
  segs_annotated <- segs_annotated %>% 
    distinct() %>% 
    group_by(seqnames, start, end, A5_ID) %>% 
    slice_max(order_by = Class) %>% #select one classification for the segment
    group_by(seqnames, start, end, !!sym(grouping_feature)) %>% 
    mutate(total_samples_window=n()) 
    
    total_samples_per_class <- segs_annotated %>% 
      group_by(!!sym(grouping_feature)) %>% 
      summarise(total_samples_class=max(total_samples_window))
    
    count_table <- segs_annotated %>% group_by(seqnames, start, end, !!sym(grouping_feature), total_samples_window, Class) %>% 
      dplyr::count() %>% 
      left_join(total_samples_per_class) 
    
    return(count_table)
}

############
# Plotting #
############


######
# Anatomy
######

seg_counts_anatomy <- seg_window_overlap_df %>%  
  purrr:::map(.f = ~count_members(.x, 
                                  grouping_feature="cell_of_origin",
                                  permitted_anatomical_groups=NULL,
                                  permitted_cell_of_origin=NULL)
                                  )  

seg_counts_anatomy <- bind_rows(seg_counts_anatomy)

seg_counts_anatomy$Class <- factor(as.character(seg_counts_anatomy$Class), levels=intersect(cn_event_types, unique(seg_counts_anatomy$Class)))

#Require 75% of samples in a class to have a classification for a window
plot_data <- seg_counts_anatomy %>% 
  filter(total_samples_window/total_samples_class > 0.75, total_samples_class > 10) %>% 
  mutate(proportion=n/total_samples_window) %>% 
  inner_join(chr_offsets ,by=c("seqnames"="chromosome")) %>% 
  mutate(start_offset=start+offset)

x_breaks <- chr_offsets %>% mutate(offset_midway=offset+((lead(offset, default = (155689000+2863281000))-offset)/2)) %>% pull(offset_midway)
names(x_breaks) <- chr_offsets$chromosome

gg_cna_freq_anatomy <- ggplot(plot_data, 
       aes(x=start_offset, y=proportion, fill= Class,color= Class)) + geom_col() + 
  facet_wrap("cell_of_origin", ncol=1) +
  geom_vline(data=chr_offsets, aes(xintercept=offset)) +
  scale_color_manual(values=cna_palette) +
  scale_fill_manual(values=cna_palette) +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1, vjust = 1)) +
  scale_x_continuous(breaks = x_breaks, labels = names(x_breaks), expand = c(0,0))

plot_data_flip <- plot_data %>% mutate(proportion=ifelse(Class %in% c("Loss", "CNLOH"), proportion * -1, proportion)) #, "Subclonal Loss"

gg_cna_freq_anatomy_flip <- ggplot(plot_data_flip  %>%  filter(Class %in% c("Gain","Loss", "CNLOH")), #, "Subclonal Loss"
                                    aes(x=start_offset, y=proportion, fill= Class,color= Class)) + geom_col() + 
  facet_wrap("cell_of_origin", ncol=1) +
  geom_vline(data=chr_offsets, aes(xintercept=offset)) +
  scale_color_manual(values=cna_palette) +
  scale_fill_manual(values=cna_palette) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle=45, hjust=1, vjust = 1),
        strip.text = element_text(margin = margin(0,0,0,0))) +
  scale_x_continuous(breaks = x_breaks, labels = names(x_breaks), expand = c(0,0)) +
  geom_hline(yintercept = 0)

######
# Clinical Behaviour
######


seg_counts_sampletype <- seg_window_overlap_df %>%  
  purrr:::map(.f = ~count_members(.x, 
                                  grouping_feature="differential_group_sampletype_strict",
                                  permitted_anatomical_groups=NULL,
                                  permitted_cell_of_origin="Chromaffin"))  

seg_counts_sampletype <- bind_rows(seg_counts_sampletype)

seg_counts_sampletype$Class <- factor(as.character(seg_counts_sampletype$Class), levels=intersect(cn_event_types, unique(seg_counts_sampletype$Class)))

#Require 75% of samples in a class to have a classification for a window
plot_data <- seg_counts_sampletype %>% 
  filter(total_samples_window/total_samples_class > 0.75, total_samples_class > 10) %>% 
  mutate(proportion=n/total_samples_window) %>% 
  inner_join(chr_offsets ,by=c("seqnames"="chromosome")) %>% 
  mutate(start_offset=start+offset) %>% 
  mutate()

x_breaks <- chr_offsets %>% mutate(offset_midway=offset+((lead(offset, default = (155689000+2863281000))-offset)/2)) %>% pull(offset_midway)
names(x_breaks) <- chr_offsets$chromosome

gg_cna_freq_clinical <- ggplot(plot_data %>%  filter(differential_group_sampletype_strict %in% c("Non-metastatic primary", "Metastasis")), 
       aes(x=start_offset, y=proportion, fill= Class,color= Class)) + geom_col() + 
  facet_wrap("differential_group_sampletype_strict", ncol=1) +
  geom_vline(data=chr_offsets, aes(xintercept=offset)) +
  scale_color_manual(values=cna_palette) +
  scale_fill_manual(values=cna_palette) +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1, vjust = 1)) +
  scale_x_continuous(breaks = x_breaks, labels = names(x_breaks), expand = c(0,0))

plot_data_flip <- plot_data %>% mutate(proportion=ifelse(Class %in% c("Loss", "CNLOH", "Subclonal Loss"), proportion * -1, proportion))

gg_cna_freq_clinical_flip <- ggplot(plot_data_flip %>%  
                                      filter(differential_group_sampletype_strict %in% c("Non-metastatic primary", "Metastasis")) %>%  
                                      filter(Class %in% c("Gain","Loss", "CNLOH")), 
       aes(x=start_offset, y=proportion, fill= Class,color= Class)) + geom_col() + 
  facet_wrap("differential_group_sampletype_strict", ncol=1) +
  geom_vline(data=chr_offsets, aes(xintercept=offset)) +
  scale_color_manual(values=cna_palette) +
  scale_fill_manual(values=cna_palette) +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1, vjust = 1), ) +
  scale_x_continuous(breaks = x_breaks, labels = names(x_breaks), expand = c(0,0)) +
  geom_hline(yintercept = 0)


