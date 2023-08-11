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

pileup_summary_files <- list.files("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/paired_samples/count_summaries", 
                                   pattern = ".txt", 
                                   full.names = T)
names(pileup_summary_files) <- gsub(".txt","",basename(pileup_summary_files))
pileup_summaries <- map(.x = pileup_summary_files,.f = read_delim, show_col_types = FALSE)


##########################
# Read in position files #
##########################

position_files <- list.files("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/paired_samples/positions/", 
                                   pattern = ".txt", 
                                   full.names = T)
names(position_files) <- gsub(".txt","",basename(position_files))
positions <- map(.x = position_files,.f = \(x) read_delim(file = x, 
                                                          col_names = c("CHROM","POS","POS_END","REF","ALT","variant_type"),
                                                          show_col_types = FALSE))



#################################
# Join pileup to position files #
#################################

if (!all(names(pileup_summaries) == names(positions))) {stop("Sanity check failed")}
    
pileup_summaries_annotated <- map2(.x = pileup_summaries,
     .y = positions,
     .f = function(pileup_summary, position_file)
          {
            pileup_summary %>% 
            inner_join(position_file %>% 
                         dplyr::select(-POS_END), 
                       by = c("CHROM", "POS", "REF", "ALT"))
          })



################
# Compute VAFs #
################

pileup_to_vaf <- function (pileup) 
{
  pileup %>% 
  pivot_longer(cols = c(-CHROM, -POS, -REF, -ALT, -variant_type), 
               names_to = c("A5_ID", "value_type"), 
               names_sep = "_",
               values_to="count") %>% 
  filter(value_type %in% c("altCount", "TotalCount")) %>% 
  distinct() %>% 
  pivot_wider(names_from = value_type, values_from = count) %>% 
  mutate(tumour_pileup_VAF=altCount/TotalCount) %>% 
  pivot_wider(id_cols = c(CHROM, POS, REF, ALT), names_from = A5_ID, values_from = tumour_pileup_VAF) %>% 
    dplyr::rename_with(.fn=\(x) ifelse(grepl("-[BT]0[0-9]", x), paste0(x,"_pileup_vaf"), x))
}

pileup_vafs <- map(.x = pileup_summaries_annotated,
           .f = pileup_to_vaf)
