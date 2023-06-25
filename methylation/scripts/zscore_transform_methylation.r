library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Read in the A5 clinical data
A5_clinical <- read_csv("/data/cephfs/punim0648/Pattison_projects/A5/A5-paper/Clinical data/A5 full clinical and genomic version 2 - A5 full clinical and genomic version 2.csv")


# Read in the A5 alpha values
A5_avals <- read_csv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Intermediate data/A5 methylation Avals.csv")
A5_avals[1:5,1:5]

annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
annEPIC[1:5,]

genes_tab <- annEPIC%>%
  data.frame()%>%
  # Select out the annotation columns I want
  dplyr::select(Probe = Name, Regulatory_Feature_Group, GencodeBasicV12_NAME)%>%
  # Remove extra gene names
  mutate(Gene = gsub(";.*" ,"", GencodeBasicV12_NAME))

unique(genes_tab$Regulatory_Feature_Group)

join_genes_tab <- genes_tab%>%
  select(Probe, Gene)

A5_avals_gene <- A5_avals %>%
  left_join(join_genes_tab)%>%
  select(-Probe)%>%
  filter(Gene != "")%>%
  gather(`A5 ID`, A_value, -Gene)%>%
  mutate(`A5 ID` = gsub(".*\\.|_D", "", `A5 ID`))%>%
  mutate(`A5 ID` = gsub("_T0", "-", `A5 ID`))%>%
  mutate(`A5 ID` = gsub("-1_1", "-1", `A5 ID`))%>%
  mutate(`A5 ID` = gsub("-1_2", "-2", `A5 ID`))%>%
  group_by(Gene,`A5 ID`)%>%
  summarise(Mean_A = mean(A_value))%>%
  pivot_wider(names_from = `A5 ID`, values_from = Mean_A)

write_csv(A5_avals_gene, "/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Outputs/Methylation_tables/Mean methylation Avals gene level all probes.csv")

meth_matrix <- A5_avals_gene%>%
  ungroup()%>%
  select(-Gene)%>%
  as.matrix()

rownames(meth_matrix) <- A5_avals_gene$Gene
meth_z_score <- meth_matrix%>%
  t()%>%
  scale()%>%
  t()%>%
  data.frame(check.names = F)%>%
  rownames_to_column("Gene")%>%
  write_csv("/data/gpfs/projects/punim0648/Pattison_projects/A5/A5-paper/Illumina-EPIC-Methylation/Outputs/Methylation_tables/Mean methylation Avals gene level Z score all probes.csv")

# Check a gene I know should be methylated in TERT or ATRX
TWIST <-meth_z_score%>%
  filter(Gene == "TWIST1")%>%
  gather(`A5 ID`, Z, -Gene)%>%
  left_join(A5_clinical)%>%
  select(`A5 ID`, Z, Gene, TERT_or_ATRX)




