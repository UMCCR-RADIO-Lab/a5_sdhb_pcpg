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

pileup_summary_files <- list.files("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/count_summaries/", 
                                   pattern = "_pileup_summary.txt", 
                                   full.names = T)

pileup_summaries <- map(.x = pileup_summary_files,.f = read_delim)

pileup_summary <- bind_rows(pileup_summaries)
rm(pileup_summaries)

#Keep only ALT count columns
pileup_summary_alt_only <- pileup_summary %>%  dplyr::select(CHROM, POS, REF, ALT, ends_with("_altCount"))

##########################
# Read in position files #
##########################

position_files <- list.files("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/positions", 
                                   pattern = "chr.+_part.+.txt", 
                                   full.names = T)

positions <- map(.x = position_files,.f = \(x) read_delim(file = x, 
                                                          col_names = c("CHROM","POS","POS_END","REF","ALT","variant_id")))

positions <- bind_rows(positions)


################################
# Join pilup to position files #
################################

pileup_summary_alt_only <- pileup_summary_alt_only %>% 
  inner_join(positions %>% 
               dplyr::select(-POS_END), 
             by = c("CHROM", "POS", "REF", "ALT")) %>% 
  dplyr::relocate(variant_id)


#########################################
# Assess attrition at different cutoffs #
#########################################

#Threshold ranges
sample_support <- 2:10
read_support <- 3:10

##############
# Iterated over thresholds
##############

observed_normals <- list()

proportion_blacklisted <- tibble(n_observed_reads=vector(mode="numeric"),
                                 n_observed_normals=vector(mode="numeric"),
                                 prop_blacklisted=vector(mode="numeric"))

key_counts <- pileup_summary_alt_only %>% group_by(variant_id) %>% dplyr::count()
indels <- key_counts %>% filter(n>1) %>% pull(variant_id)
snvs <- key_counts %>% filter(n==1) %>% pull(variant_id)

for(rs in read_support){
  above_rs_threshold <- pileup_summary_alt_only[, -c(1,2,3,4,5)] >= rs
  
  above_rs_threshold <- bind_cols(variant_id=pileup_summary_alt_only$variant_id, above_rs_threshold)
  
  above_rs_threshold_snv <- above_rs_threshold %>% filter(variant_id %in% snvs)
  
  above_rs_threshold_indel <- above_rs_threshold %>% filter(variant_id %in% indels)
  
  above_rs_threshold_indel <- above_rs_threshold_indel %>% 
    group_by(variant_id) %>% 
    summarise(across(.cols=everything(), .fns = \(x) all(x==T)))
  
  above_rs_threshold <- bind_rows(above_rs_threshold_indel, above_rs_threshold_snv)
  
  observed_normals[[rs]] <-  bind_cols(above_rs_threshold[,1, drop=F], 
                                       n_above_threshold=rowSums(above_rs_threshold[,-1]))
  
  total_mutations <- nrow(observed_normals[[rs]])
  
  for (ss in sample_support)
  {
    
    n_filtered <- observed_normals[[rs]] %>% 
      mutate(threshold_category=cut(x = n_above_threshold, 
                                    breaks=c(-1,ss,100), 
                                    labels=c("below_threshold", "above_threshold"))) %>% 
      group_by(threshold_category) %>% dplyr::count() %>% 
      filter(threshold_category=="above_threshold") %>% 
      pull(n)
    prop_filtered <- n_filtered/total_mutations
    
    proportion_blacklisted <- add_case(proportion_blacklisted, 
                                       n_observed_reads=rs,
                                       n_observed_normals=ss,
                                       prop_blacklisted=prop_filtered)
  }
  
}
proportion_blacklisted <- proportion_blacklisted %>%  mutate(n_observed_normals=factor(n_observed_normals,levels=sample_support))


##############
# Plot blacklist occupancy 
##############

ggplot(proportion_blacklisted, aes(x=n_observed_reads, y=prop_blacklisted*100, color=n_observed_normals, group=n_observed_normals)) +
  geom_point() +
  geom_line() +
  xlab("Read support threshold") +
  ylab("% blacklisted") +
  scale_color_discrete(name="nNormals Threshold") + theme_bw()

###################
# Write blacklist #
###################

#Select thresholds
read_support_cutoff <- 3
normal_support_cutoff <- 3

blacklist <- observed_normals[[read_support_cutoff]] %>%  filter(n_above_threshold >= normal_support_cutoff)

blacklist <- blacklist %>%  separate(variant_id, into=c("seqnames","start","end", "ALT"))

write.table(blacklist, 
            file = "/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/blacklists/blacklist_readsupport_gteq3_samplesupport_gteq3.tsv",
            sep="\t",
            quote = F, 
            row.names = F)
