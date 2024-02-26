##################################################
# Run Mutational Patterns on A5 somatic variants #
# for de novo signature analysis                 #
# Author: Aidan Flynn                            #
# Date: 20/10/2023                               #
# Languages: R                                   #
##################################################


library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(purrr)
library(ggplot2)

setwd("/g/data/pq08/projects/ppgl/")

#load clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#Load cosmic signatures
signatures = get_known_signatures(genome = "GRCh38", muttype="snv")

#Load black listed variants
blacklist="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/blacklist/blacklists/blacklist_readsupport_gteq3_samplesupport_gteq3.tsv"
blacklisted_variants <- read.delim(blacklist)
blacklisted_variants_gr <- makeGRangesFromDataFrame(blacklisted_variants, seqnames.field = "CHROM", start.field = "POS", end.field = "POS", keep.extra.columns = T)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

vcf_files <- list.files("./a5/wgs/symlinks/somatic_vcf",pattern = "-somatic-PASS.vcf.gz$", full.names = TRUE,recursive = T)

sample_names <- gsub("__.*","", basename(vcf_files))

all_grl <- read_vcfs_as_granges(vcf_files = vcf_files, sample_names = sample_names, genome=ref_genome, type = "all")


#Remove blacklisted SNVs
all_grl_noblacklist <- map(.x = all_grl, .f = function(vcf, blacklist_gr=blacklisted_variants_gr) {
  hits <- findOverlaps(subject=vcf, query=blacklist_gr)
  mcols(vcf)[["FILTER"]][subjectHits(hits)] <- "blacklist"
  return(vcf[mcols(vcf)[["FILTER"]]=="PASS",])
})

##SNVs
grl <- get_mut_type(all_grl_noblacklist, type = "snv")


mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
mut_mat <- mut_mat + 0.0001


#Sample sets
all_samples <- a5_anno %>% filter(Exclude == "N") %>% pull(A5_ID)
no_e167_1 <- a5_anno %>% filter(Exclude == "N", A5_ID != "E167-1") %>% pull(A5_ID)
treated <- a5_anno %>% filter(Exclude == "N", grepl("Yes",Resection_post_dna_damaging_treatment)) %>% pull(A5_ID)
primary <- a5_anno %>% filter(Exclude == "N", differential_group_sampletype_strict=="Non-metastatic primary") %>% pull(A5_ID)
met <- a5_anno %>% filter(Exclude == "N", differential_group_sampletype_strict=="Metastasis") %>% pull(A5_ID)
sample_sets <- list(all_samples=all_samples, no_e167_1=no_e167_1,treated=treated, primary=primary, met=met)

estimate_list <- list()
min_rank <- 2
max_rank <- 10
for (current_sample_set in names(sample_sets))
{
  message("Processing estimate for ", current_sample_set, "...")
  mut_mat_current = mut_mat[,sample_sets[[current_sample_set]]]
  mut_mat_md5 <- digest::digest(mut_mat_current,algo="md5")
  checkpoint_file <- paste0("./a5/wgs/quickload_checkpoints/de_novo_nmf_estimate_",mut_mat_md5,".rds")
  if (file.exists(checkpoint_file)) {
    message("Quickloading estimate for ", current_sample_set, ".")
    estimate_list[[current_sample_set]] <- readRDS(checkpoint_file)
  } else {
    estimate_list[[current_sample_set]] <- nmf(mut_mat, rank = min_rank:max_rank, method = "brunet", 
                    nrun = 50, seed = 123456, .opt = "v-p")
    saveRDS(estimate_list[[current_sample_set]], checkpoint_file)
  }
}
#plot(estimate)

min_rank <- 2
max_rank <- 10
sig_ex <- list()
for (current_sample_set in names(sample_sets))
{
  mut_mat_current = mut_mat[,sample_sets[[current_sample_set]]]
  mut_mat_md5 <- digest::digest(mut_mat_current,algo="md5")
  checkpoint_file <- paste0("./a5/wgs/quickload_checkpoints/de_novo_nmf_extract_",mut_mat_md5,".rds")
  if (file.exists(checkpoint_file)) {
    message("Quickloading extraction for ", current_sample_set)
    sig_ex[[current_sample_set]] <- readRDS(checkpoint_file)
  } else {
    sig_ex[[current_sample_set]] <- list()
    
    for (rank in min_rank:max_rank)
    {
      message("Processing extraction for ", current_sample_set, ": Rank", rank, "...")
      sig_ex[[current_sample_set]][[rank]] <- extract_signatures(mut_mat_current, rank = rank, nrun = 200, single_core = FALSE, seed = 123456)
    }
    saveRDS(sig_ex[[current_sample_set]], checkpoint_file)
  }
}



# nmf_res <- rename_nmf_signatures(nmf_res, signatures, cutoff = 0.85)
# colnames(nmf_res$signatures)
# 
# plot_data <- nmf_res_rank6$signatures %>% 
#   as_tibble(rownames="context") %>% 
#   pivot_longer(cols = -context, names_to = "signature", values_to = "contribution") %>% 
#   mutate(Alteration = stringr::str_extract(context,".>."),
#          Alteration = factor(Alteration, levels=c("C>A","C>G","C>T","T>A","T>C","T>G")),
#          context=gsub("\\[.>.\\]",".",context))
# ggplot(plot_data, aes(x=context, y=contribution)) + geom_col() + facet_grid(signature~Alteration)
# 
# current_nmf <- nmf_res_rank6$signatures
# cosine_result <- matrix(nrow=ncol(current_nmf), ncol=ncol(signatures),dimnames = list(colnames(current_nmf), colnames(signatures)) )
# for (i in 1:ncol(current_nmf))
# {
#  for (j in 1:ncol(signatures))
#  {
#    cosine_result[i,j] <- lsa::cosine(current_nmf[,i],  signatures[,j])
#  }
# }
# plot_data <- cosine_result %>% as_tibble(rownames="NMF_sig") %>%  
#   pivot_longer(cols = -NMF_sig, names_to = "SBS_sig", values_to = "cosine_sim") %>% 
#   mutate(SBS_sig = factor(SBS_sig, levels=rev(colnames(signatures))))
# ggplot(plot_data, aes(x=NMF_sig, y=SBS_sig, fill=cosine_sim)) + geom_tile() + scale_fill_gradientn(colours = c("white","white", "orange","red"), values = c(0,0.7,0.8,1))
#               
              
