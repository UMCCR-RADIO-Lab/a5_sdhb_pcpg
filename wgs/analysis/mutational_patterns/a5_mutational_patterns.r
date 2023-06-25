##################################################
# Run Mutational Patterns on A5 somatic variants #
# applied across the cohort                      #
# Author: Aidan Flynn                            #
# Date: 15/10/2022                               #
# Languages: R                                   #
##################################################


library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(purrr)
library(ggplot2)

setwd("/g/data/pq08/projects/ppgl/")

#Load black listed variants
blacklist="/g/data/pq08/projects/ppgl/a5/wgs/analysis/mpileup/blacklists/blacklist_readsupport_gteq3_samplesupport_gteq3.tsv"
blacklisted_variants <- read.delim(blacklist)
blacklisted_variants_gr <- makeGRangesFromDataFrame(blacklisted_variants)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"

vcf_files <- list.files("./a5/wgs/symlinks/somatic_vcf",pattern = "-somatic-PASS.vcf.gz$", full.names = TRUE,recursive = T)

sample_names <- gsub("__.*","", basename(vcf_files))

all_grl <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "all")


#Remove blacklisted SNVs
all_grl_noblacklist <- map(.x = all_grl, .f = function(vcf, blacklist_gr=blacklisted_variants_gr) {
  hits <- findOverlaps(subject=vcf, query=blacklist_gr)
  mcols(vcf)[["FILTER"]][subjectHits(hits)] <- "blacklist"
  return(vcf[mcols(vcf)[["FILTER"]]=="PASS",])
})

##SNVs
grl <- get_mut_type(all_grl_noblacklist, type = "snv")

type_occurrences <- mut_type_occurrences(grl, ref_genome)

ggsave(file="./a5/wgs/analysis/mutational_patterns/snv_type_occurrences.pdf", 
       plot=plot_spectrum(type_occurrences), 
       device=pdf(),
       width=6, height=8, limitsize=F)

signatures = get_known_signatures(genome = "GRCh38", muttype="snv")
mut_mat <- mut_matrix(vcf_list = grl, ref_genome = ref_genome)
fit_res <- fit_to_signatures(mut_mat, signatures)
write.table(fit_res$contribution, "./a5/wgs/analysis/mutational_patterns/A5_blacklistfiltered_sbs_cosmicv3_fitcontribution.tsv", sep="\t", quote=F) 

#DBS
dbs_grl <- get_mut_type(all_grl_noblacklist, type = "dbs")
dbs_grl <- get_dbs_context(dbs_grl)
dbs_counts <- count_dbs_contexts(dbs_grl)
signatures_dbs = get_known_signatures(genome = "GRCh38", muttype = "dbs")
fit_res_dbs <- fit_to_signatures(dbs_counts, signatures_dbs)
write.table(fit_res_dbs$contribution, "./a5/wgs/analysis/mutational_patterns/A5_blacklistfiltered_dbs_cosmicv3_fitcontribution.tsv", sep="\t", quote=F) 

ggplot2::ggsave(file="./a5/wgs/analysis/mutational_patterns/dbs_contexts.pdf", 
                plot=plot_dbs_contexts(dbs_counts, condensed = TRUE), 
                device=pdf(),
                width=30, height=100, limitsize=F)


#INDEL 
indel_grl <- get_mut_type(all_grl_noblacklist, type = "indel")
indel_grl <- get_indel_context(indel_grl, ref_genome)
indel_counts <- count_indel_contexts(indel_grl)

signatures_indel = get_known_signatures(muttype = "indel")
fit_res_indel <- fit_to_signatures(indel_counts, signatures_indel)
write.table(fit_res_indel$contribution, "./a5/wgs/analysis/mutational_patterns/A5_blacklistfiltered_indel_cosmicv3_fitcontribution.tsv", sep="\t", quote=F) 

ggplot2::ggsave(file="./a5/wgs/analysis/mutational_patterns/indel_contexts.pdf", 
                plot=plot_indel_contexts(indel_counts, condensed = TRUE), 
                device=pdf(),
                width=30, height=100, limitsize=F)

save.image("./a5/wgs/analysis/mutational_patterns/A5_blacklistfiltered_mutational_patterns.rworkspace")





# p <- . %>% as_tibble(rownames = "signature") %>% pivot_longer(cols = -signature, names_to="A5_ID", values_to = "Exposure")
# 
# contribution_noblacklist <- read.delim("./a5/wgs/analysis/mutational_patterns/A5_sbs_cosmicv3_fitcontribution.tsv", 
#                                        check.names = F)
# noblacklist_keep <- apply(X = contribution_noblacklist, MARGIN = 1,FUN =  max) > 100
# 
# contribution_blacklist <- read.delim("./a5/wgs/analysis/mutational_patterns/A5_blacklistfiltered_sbs_cosmicv3_fitcontribution.tsv", 
#                                        check.names = F)
# blacklist_keep <- apply(X = contribution_blacklist, MARGIN = 1,FUN =  max) > 100
# 
# all(names(noblacklist_keep) == names(blacklist_keep) )
# 
# contribution_noblacklist <- contribution_noblacklist[noblacklist_keep | blacklist_keep,]
# 
# contribution_blacklist <- contribution_blacklist[noblacklist_keep | blacklist_keep,]
# 
# 
# plot_data <- contribution_noblacklist %>% p %>% 
#   left_join(contribution_blacklist %>% p, 
#             by=c("signature","A5_ID"), 
#             suffix=c("_no_blacklist","_blacklist")) %>% 
#   mutate(across(.cols = c(Exposure_no_blacklist, Exposure_blacklist), 
#                 .fns = \(x) replace_na(x,0))) %>% 
#   mutate(signature = factor(signature, levels=rownames(contribution_blacklist)))
# 
# ggsave(plot = ggplot(plot_data, aes(x=Exposure_no_blacklist, 
#                                     y=Exposure_blacklist, color=signature)) + 
#          geom_point() + 
#          facet_wrap("A5_ID", scales="free"), 
#        filename = "a5/wgs/results/signatures/sigs_post_blacklist.pdf", 
#        height = 40, 
#        width = 40, 
#        units = "in")
# 
# ggsave(plot = ggplot(plot_data, aes(x=Exposure_no_blacklist, 
#                                     y=Exposure_blacklist, color=signature)) + 
#          geom_point() + 
#          ggrepel::geom_text_repel(mapping = aes(label=A5_ID),size=3) + 
#          facet_wrap("signature", scales="free"), 
#        filename = "a5/wgs/results/signatures/sigs_post_blacklist2.pdf", 
#        height = 40, 
#        width = 40, 
#        units = "in")
