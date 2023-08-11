#########################################
# This script makes a variant blacklist #
# based on summarised counts of mpileup #
# data from normal samples              #
#                                       #
# Author: Aidan Flynn                   #
# Date: 26/05/2023                      #
#########################################

library(tidyverse)

########################################
# Read in count summaries from pileups #
########################################

pileup_summary_files <- list.files("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/count_summaries/paired_samples/", 
                                   pattern = "_pileup_summary.txt", 
                                   full.names = T)
names(pileup_summary_files) <- gsub("_pileup_summary.txt","",basename(pileup_summary_files))
pileup_summaries <- map(.x = pileup_summary_files,.f = read_delim)


##########################
# Read in position files #
##########################

position_files <- list.files("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/positions/paired_samples", 
                                   pattern = ".txt", 
                                   full.names = T)
names(position_files) <- gsub(".txt","",basename(position_files))
positions <- map(.x = position_files,.f = \(x) read_delim(file = x, 
                                                          col_names = c("CHROM","POS","POS_END","REF","ALT","variant_id")))



################################
# Join pilup to position files #
################################

if (!all(names(pileup_summaries) == names(positions))) {stop("Sanity check failed")}
    
pileup_summaries_annotated <- map2(.x = pileup_summaries,
     .y = positions,
     .f = function(pileup_summary, position_file)
          {
            pileup_summary %>% 
            inner_join(position_file %>% 
                         dplyr::select(-POS_END), 
                       by = c("CHROM", "POS", "REF", "ALT")) %>% 
            dplyr::relocate(variant_id)
          })



################
# Summarise INDELs VAFs #
################

pileup_summaries_annotated$E143 %>% group_by(variant_id) %>% 
  summarise(across(.cols = dplyr::matches("_altCount|_refCount|_TotalCount"),
                   .fns = \(x) { ceiling(mean(x)) }))



above_rs_threshold_indel <- above_rs_threshold_indel %>% 
  group_by(variant_id) %>% 
  summarise(across(.cols=everything(), .fns = \(x) all(x==T)))

################
# Compute VAFs #
################

pileup_summaries_annotated$E143 %>%  pivot_longer(cols = c(-variant_id, -CHROM, -POS, REF, ALT, ))

vafs <- map(.x = pileup_summaries,
           .f = function(pileup_summary, position_file)
           {
             
             pileup_summary %>% 
               inner_join(position_file %>% 
                            dplyr::select(-POS_END), 
                          by = c("CHROM", "POS", "REF", "ALT")) %>% 
               dplyr::relocate(variant_id)
           })