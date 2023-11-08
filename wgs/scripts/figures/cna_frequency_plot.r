library(ggplot2)
#library(patchwork)
library(plyranges)
library(tidyr)

################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_cna_segs()


source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

################
# Make windows #
################

chrom_gr <- GRanges(seqnames = chr_offsets$chromosome,  ranges = IRanges(start = 1, end = chr_offsets$chr_size))
chrom_windows <- tile(chrom_gr,width = 1000000)
names(chrom_windows) <- purrr::map_chr(.x=chrom_windows, .f=\(x) { as.character(seqnames(x))[1] })

###############
# Prepare CNA #
###############

seg_gr <- a5_seg_keep %>%  
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


cna_simplify_coding <- c(
  "Loss"="Loss",
  "Loss"="Loss + Subclonal CNLOH",
  "Subclonal Loss"="Subclonal Loss",
  "Subclonal Loss"="Minor Subclonal Loss",
  "CNLOH"="CNLOH",
  "Diploid/Haploid-X"="Diploid/Haploid-X",
  "Gain"="Gain+LOH",
  "Gain"="Gain",
  "Gain"="WGD+Gain",              
  "Gain"="Subclonal Gain",
  "Other"="Hom. Del.",
  "Diploid/Haploid-X"="Other",
  "WGD"="WGD",
  "Chromothripsis"="Chromothripsis")

count_members <- function(chrom_segs, grouping_feature, permitted_anatomical_groups) { 
  segs_annotated <- 
    chrom_segs %>%
    mutate(Class=forcats::fct_recode(Class, !!!cna_simplify_coding),
           Class=factor(as.character(Class, levels=cn_event_types)),
           Class=forcats::fct_drop(Class)) %>% 
    dplyr::select(seqnames,  start, end, A5_ID, Class) %>% 
    inner_join(a5_anno %>%  
                 dplyr::select(A5_ID, differential_group_anatomy, !!sym(grouping_feature)) %>% 
                 mutate(differential_group_anatomy=ifelse(A5_ID=="E185-1", "Head_Neck",differential_group_anatomy))) %>% 
    filter(differential_group_anatomy %in% permitted_anatomical_groups) %>% 
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


#################
# Antatomy
#################

seg_counts_anatomy <- seg_window_overlap_df %>%  
  purrr:::map(.f = ~count_members(.x, 
                                  grouping_feature="differential_group_anatomy", 
                                  permitted_anatomical_groups=c("Head_Neck", "Abdominal_Thoracic")))  

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
ggplot(plot_data, 
       aes(x=start_offset, y=proportion, fill= Class,color= Class)) + geom_col() + 
  facet_wrap("differential_group_anatomy", ncol=1) +
  geom_vline(data=chr_offsets, aes(xintercept=offset)) +
  scale_color_manual(values=cna_palette) +
  scale_fill_manual(values=cna_palette) +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1, vjust = 1)) +
  scale_x_continuous(breaks = x_breaks, labels = names(x_breaks), expand = c(0,0))


#################
# Clinical Behaviour
#################

seg_counts_sampletype <- seg_window_overlap_df %>%  
  purrr:::map(.f = ~count_members(.x, 
                                  grouping_feature="differential_group_sampletype_strict",
                                  permitted_anatomical_groups=c("Abdominal_Thoracic")))  

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
ggplot(plot_data %>%  filter(differential_group_sampletype_strict %in% c("Non-metastatic primary", "Metastasis")), 
       aes(x=start_offset, y=proportion, fill= Class,color= Class)) + geom_col() + 
  facet_wrap("differential_group_sampletype_strict", ncol=1) +
  geom_vline(data=chr_offsets, aes(xintercept=offset)) +
  scale_color_manual(values=cna_palette) +
  scale_fill_manual(values=cna_palette) +
  theme_bw() + 
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle=45, hjust=1, vjust = 1)) +
  scale_x_continuous(breaks = x_breaks, labels = names(x_breaks), expand = c(0,0))




