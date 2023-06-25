rm(list=ls())

library(edgeR)
library(tidyverse)

#----
# Load Data and metadata into R
#----

rpt_counts <- read_csv("Data/repeat_analysis/counts/repeat_type_counts.csv")
rpt_counts <- rpt_counts %>%
  dplyr::rename(NAM018 = "Chromaffin_cell") %>%
  data.frame() %>% 
  column_to_rownames(var = "repeat_type") 

# store sample names in vector
sample.names <- colnames(rpt_counts)
# read in metadata
sampleinfo <- read_csv("Data/A5_singlenuclei_metadata.csv")
sampleinfo <- sampleinfo %>% 
  mutate(Sample.ID = gsub(pattern = "-", replacement = "_", Sample.ID)) # change underscores to hyphens for consistency

# reorder the info 
sample.order <- tibble(Sample.ID = sample.names)
sampleinfo <- sample.order %>% left_join(sampleinfo)

# make column describing primary and secondary driver mutation 
sampleinfo <- sampleinfo %>%
  mutate(
    group = paste(
      Secondary.driver,
      Primary.driver,
      sep = "_")) %>% 
  mutate(group = gsub("SDHA|SDHB", "SDHx",  group)) %>% 
  mutate(group = gsub ("Normal_", "", group)) %>% 
  mutate(TERT_mutation = recode(Secondary.driver, 
                                "TERT_subclonal" = "TERT",
                                "ATRX" = "WT"))
group <- sampleinfo$group

#----
# EdgeR pipeline 
#----

# these library sizes were calculated using the sinto-filtered bams
# using the command 'samtools view -c'
# they represent the total number of alignments present in the neoplastic cells

rpt_dge <- DGEList(counts = rpt_counts, group = group)
keep <- filterByExpr(rpt_dge)
rpt_dge <- rpt_dge[keep,,keep.lib.sizes=FALSE]
rpt_dge <- calcNormFactors(rpt_dge)

design <- model.matrix(~0+group)

TERTvsATRX <- makeContrasts(groupTERT_SDHx-groupATRX_SDHx, levels=design)
rpt_dge <- estimateDisp(rpt_dge, design)

#----
# quasi-likelihood F-tests
#----

fit <- glmQLFit(rpt_dge, design)
qlf <- glmQLFTest(fit, contrast = TERTvsATRX)
toptags <- topTags(qlf, adjust.method = "BH")
toptags
