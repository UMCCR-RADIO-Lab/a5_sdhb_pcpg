#############################
#                           #
# This scripts generates a  #
# variant blacklist based   #
# on variant sites that are # 
# recurrently rejected in   #
# multiple samples          #
#                           #
# Author: Aidan Flynn       #
# Date 16/10/2024           #
#                           #
#############################

library(tidyverse)

################
# Data loading #
################

######
# Read in position file with all PASS call variant sites
######

pos_files <- list.files("/g/data/pq08/projects/ppgl/a5/wgs/analysis/blacklist_from_vcf/input/positions", 
                        full.names = T, 
                        pattern = "*.txt")

names(pos_files) <- gsub(".txt", "", basename(pos_files))
var_pos <- map(pos_files, ~read_tsv(.x,
                                    col_names = c("Chr","Pos_Start", "Pos_End", "Ref", "Alt","Type") ))

######
# Read in observation files (see make_extract_calls_qsub.sh for creation code)
######

var_obs_files <- list.files("/g/data/pq08/projects/ppgl/a5/wgs/analysis/blacklist_from_vcf/output/variant_observations", 
                            full.names = T, 
                            pattern = "*.tsv")
names(var_obs_files) <- gsub(".tsv", "", basename(var_obs_files))

if(!all(names(var_pos) == names(var_obs_files))) stop("Failed order check")

var_obs <- map(var_obs_files, 
               ~read_tsv(.x,
                         col_names = c("Chr","Pos","Alt","Filter","Allele_Depth", 
                                       "Total_Depth","Caller", "Sample")))

#######
# Process mutect, vardict, and strelka2 call info to uniform format
#######

obs <- map(names(var_pos), .f = function(var_set_name) {
  
  var_set_pos <- distinct(var_pos[[var_set_name]])
  var_set_obs <- var_obs[[var_set_name]]
  
  strelka_obs <- var_set_obs %>% filter(Caller=="strelka2") %>% 
    inner_join(var_set_pos %>%  dplyr::select(-Pos_End) %>% dplyr::rename(Pos=Pos_Start), by=c("Chr", "Pos", "Alt"))
    
  strelka_obs <- strelka_obs %>% separate_wider_delim(cols = Allele_Depth, 
                                                      names = c("A_count", "T_count", "C_count", "G_count", 
                                                                "Indel_count", "Alt_count"), 
                                                      delim=":", 
                                                      too_few = c("error"),
                                                      too_many = c("error"),
                                                      cols_remove = T) %>% 
    mutate(across(ends_with("_count"),
                  .fns = \(x) {gsub(",.+", "", x)})) %>% 
    mutate(Alt_Depth = case_when(
      Type != "SBS" ~ Indel_count,
      Alt == "A" ~ A_count,
      Alt == "C" ~ C_count,
      Alt == "T" ~ T_count,
      Alt == "G" ~ G_count
    ),
    Ref_Depth = case_when(
      Type != "SBS" ~ Alt_count,
      Ref == "A" ~ A_count,
      Ref == "C" ~ C_count,
      Ref == "T" ~ T_count,
      Ref == "G" ~ G_count
    )) %>% 
    dplyr::select(Chr, Pos, Alt, Filter, Ref_Depth, Alt_Depth, Total_Depth, Caller, Sample)
  
  varmut_obs <- var_set_obs %>% filter(Caller != "strelka2")

  varmut_obs <- varmut_obs %>% 
    separate_wider_delim(cols = Allele_Depth, names = c("Ref_Depth","Alt_Depth"), delim =",",
                         too_few = "align_start",
                         too_many = "merge") %>% 
    separate_longer_delim(cols = c(Alt, Alt_Depth), delim = ",")
  
  return_set <- var_set_pos %>%  dplyr::select(-Pos_End) %>% dplyr::rename(Pos=Pos_Start) %>% filter(!is.na(Type)) %>% 
  left_join(bind_rows(varmut_obs, strelka_obs), by=c("Chr", "Pos", "Alt")) %>% 
    dplyr::select(Chr, Pos, Ref, Alt, Filter, Ref_Depth, Alt_Depth, Total_Depth, Caller, Sample, Type)
  
  
  return(return_set)
  
})

obs <- bind_rows(obs)

obs <- distinct(obs) %>%  filter(!is.na(Caller)) #Remove a small number of variants that are altered by UMCCRISE

obs <- obs %>%  mutate(Filter_Binary = ifelse(Filter=="PASS", "PASS", "REJECT"))

#######
# Count PASS/REJECT observations for each variant
#######

obs_counts <- obs %>%  group_by(Chr, Pos, Alt, Caller, Filter_Binary) %>% dplyr::count() 

obs_counts_wide <- obs_counts %>% 
  pivot_wider(id_cols = c(Chr, Pos, Alt), 
              names_from = c(Caller, Filter_Binary), values_from = n) %>% 
  mutate(across(matches("PASS|REJECT"), ~replace_na(.x, 0)))

#########
# Build blacklist based on thresholds
#########

max_rejected_samples_multi = 3
max_rejected_samples_single = 5

# blacklist <- obs_counts_wide %>%  
#   mutate(PASS_REJECT_blacklist = ifelse(
#     (mutect2_REJECT > max_rejected_samples_multi & strelka2_REJECT> max_rejected_samples_multi) |
#     (mutect2_REJECT > max_rejected_samples_multi & vardict_REJECT > max_rejected_samples_multi) |
#     (strelka2_REJECT > max_rejected_samples_multi & vardict_REJECT > max_rejected_samples_multi) |
#       mutect2_REJECT > max_rejected_samples_single |
#     strelka2_REJECT > max_rejected_samples_single | 
#       vardict_REJECT > max_rejected_samples_single,
#     T, F)
#   )


blacklist <- obs_counts_wide %>% 
  mutate(pass_reject_ratio = (mutect2_REJECT + strelka2_REJECT + vardict_REJECT)/(mutect2_PASS + strelka2_PASS + vardict_PASS),
         pass_reject_blacklist = pass_reject_ratio > 1)

#########
# Write blacklist
#########

readr::write_tsv(blacklist,
                 file="/g/data/pq08/projects/ppgl/a5/wgs/analysis/blacklist_from_vcf/output/blacklist_reject_pr_ratio_gt1.tsv")
