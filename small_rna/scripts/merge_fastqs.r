library(tidyverse)

samples <- list.files("/Users/adpattison/igv/mounts/Spartan/A5_study/Complete_dataset/small_RNA_seq/fastqs_all/")

# Apparently it's fine to just concatenate gzipped files
samples_df <- data.frame(samples)%>%
  mutate(Tumour = gsub("-sR.*", "", samples))%>%
  mutate(to_write = paste0(Tumour, "_merged_R1_001.fastq.gz"))%>%
  mutate(command = paste0("cat ", samples, " >> ", to_write))%>%
  select(command)%>%
  write_tsv("/Users/adpattison/igv/mounts/Spartan/A5_study/Complete_dataset/small_RNA_seq/fastqs_all/merge_fastqs.sh")
