library(tidyverse)
library(TCGAbiolinks)


# Find the LUSC primary RNA-Seq
query_rseq <- GDCquery(project = c("TCGA-PCPG"),
                       data.category = "Sequencing Reads",  
                       experimental.strategy = "miRNA-Seq")
# Get results and keep only primary samples
results_rseq <- getResults(query_rseq)%>%
  filter(!duplicated(substr(cases,1,12)))

# Get the manifest of the full panel of normals for PON creation
miRNA_maifest <- getManifest(query_rseq,save = "~/igv/mounts/Spartan_2/download-manifests/PCPG_mir_seq_manifest.tsv")
