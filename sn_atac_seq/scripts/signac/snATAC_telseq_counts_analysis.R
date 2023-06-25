# analysis of telomeric read-counts (computed using telseq)

library(tidyverse)
library(patchwork)

telseq_table <- read_tsv("Data/repeat_analysis/telseq/output/telseq_output.tsv")
metadata <- read_csv("Data/A5_singlenuclei_metadata.csv") %>%
  mutate(External.ID = gsub(x = External.ID, pattern = "_", replacement = "-"))
tel_lengths_wgs <- read_tsv("Data/repeat_analysis/wgs_tel_lengths_table.tsv") %>%
  dplyr::rename(Sample.ID = A5_ID)

# make a column for the number of telomeric reads (reads with 7+ telomeric repeats?)
# might change this later 
telseq_table <- telseq_table %>%
  mutate(telomeric_read_count = (TEL7 + TEL8 + TEL9 + TEL10 + TEL11 + TEL12 + TEL13 + TEL14 + TEL15 + TEL16)) %>% 
  mutate(non_duplicate_count = (Total - Duplicates))

telo_counts <- telseq_table %>%
  dplyr::select(Sample, telomeric_read_count, non_duplicate_count) 

telo_counts <- telo_counts %>%
  group_by(Sample) %>% 
  summarise(telomeric_read_count = sum(telomeric_read_count),
            non_duplicate_count = sum(non_duplicate_count)) %>% 
  rename(External.ID = Sample) %>% 
  mutate(External.ID = recode(
    External.ID, 
    E156_1_combined = "E156-1",
    E166_1_combined = "E166-1", 
    S126 = "S126-1")
  ) %>%
  left_join(metadata) %>% 
  left_join(tel_lengths_wgs)

# Calculate TPM (sort of) using tumour telomere length (WGS-derived) instead of gene length 
telo_counts <- telo_counts %>% 
  mutate(tel_RPK = telomeric_read_count/Tumour_tel_length_kb) %>% 
  mutate(per_million_scale_factor = non_duplicate_count/1000000) %>% # normally you would sum the RPK values to get the total counts for scaling
  mutate(tel_TPM = tel_RPK/per_million_scale_factor) # Technically not TPM but close enough for now 
  

# normalise to telomere length of the tumour 
# normalise to total non-duplicate # of reads?

telomeric_counts <- ggplot(telo_counts) +
  geom_bar(stat = "identity", aes(x = Sample.ID, y = tel_TPM, fill = Secondary.driver)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

non_duplicate_counts <- ggplot(telo_counts) +
  geom_bar(stat = "identity", aes(x = Sample.ID, y = non_duplicate_count, fill = Secondary.driver)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

telomeric_counts +
  non_duplicate_counts + plot_layout(ncol = 1)
