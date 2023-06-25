#############################################
# This script creates a figure illustrating #
# the TERT SV from E123 using data from     # 
# PURPLE/AMBER/COBALT                       #
#                                           #
# Author: Aidan Flynn                       #
# Date: 06/06/2023                          #
#############################################

#############
# Load Data #
#############

E123_segs <- read.delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E123/umccrised/E123-1__E123-T01/purple/E123-1__E123-T01.purple.segment.tsv")
E123_baf <- read.delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E123/umccrised/work/E123-1__E123-T01/purple/amber/E123-1__E123-T01.amber.baf.tsv", comment.char = "#")

E123_ratio <- read.delim("/g/data/pq08/projects/ppgl/a5/wgs/analysis/bcbio/E123/umccrised/work/E123-1__E123-T01/purple/cobalt/E123-1__E123-T01.cobalt.ratio.tsv")
E123_ratio <- E123_ratio %>%  filter(referenceReadCount > 0, tumorReadCount >0)
E123_ratio <- 
  E123_ratio %>% 
  mutate(referenceReadCount_norm = (referenceReadCount/sum(referenceReadCount))*10^6, 
         tumorReadCount_norm = (tumorReadCount /sum(tumorReadCount))*10^6) %>% 
  mutate(log_ratio=log2(tumorReadCount_norm/referenceReadCount_norm))
E123_ratio$chromosome <- factor(E123_ratio$chromosome, levels = paste0("chr", c(1:22,"X","Y")))
E123_ratio$position <- as.numeric(E123_ratio$position)

#########################
# Fetch Coordinate data #
#########################

mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "asia")

exon_positions <- biomaRt::getBM(filters= "ensembl_gene_id",
                                 attributes= c("ensembl_gene_id",
                                               "external_gene_name",
                                               "exon_chrom_start",
                                               "exon_chrom_end"),
                                 values=c("ENSG00000164362","ENSG00000072364"), mart = mart)

##############
# SV details #
##############

#From GRIDSS
TERT_AFF4_SV <- c("TERT"=1295378, "AFF4"=132957282)

regions <- 
  list(AFF4 = list(chr="chr5", start=132875395, end=132963634),
       TERT = list(chr="chr5", start=1253147, end=1295068))

##############
# SV details #
##############

probes <- list()
probes[["TERT"]] <- E123_ratio  %>% 
  filter(chromosome == regions$TERT$chr,  position > regions$TERT$start - (0.5*10^5), position < regions$TERT$end + (1*10^5)) 

probes[["AFF4"]] <- E123_ratio  %>% 
  filter(chromosome == regions$AFF4$chr,  position > regions$AFF4$start - (1*10^5), position < regions$AFF4$end + (1*10^5)) 


E123_segs_gr <- GenomicRanges::makeGRangesFromDataFrame(E123_segs, keep.extra.columns = T)
for (g in names(probes)) {

temp <- GenomicRanges::makeGRangesFromDataFrame(df = probes[[g]], 
                                                              start.field = "position",
                                                              end.field = "position",
                                                              keep.extra.columns = T)

hits <- findOverlaps(subject = temp, query = E123_segs_gr)
 
mcols(temp)["purple_cn"] <- NA
mcols(temp)[["purple_cn"]][subjectHits(hits)] <- mcols(E123_segs_gr)[["tumorCopyNumber"]][queryHits(hits)]

probes[[g]] <- bind_cols(probes[[g]], purple_cn=factor(round(mcols(temp)[["purple_cn"]],0), levels=1:5))
probes[[g]] <- probes[[g]] %>% filter(!is.na(purple_cn))

}

probes[["AFF4"]] <- probes[["AFF4"]] %>% mutate(purple_cn=forcats::fct_recode(.f = purple_cn, "5"="4"))

####################################
# Manually assess CN/logR values #
####################################

CN <- c("CN2"=-0.45,"CN3"=-0.1, "CN4"=0.165, "CN5"=0.371)

### Visualise CN
# E123_ratio_subset <- E123_ratio %>% 
#   group_by(chromosome) %>% 
#   slice_sample(prop = 0.01)
# 
# E123_ratio_subset <- GenomicRanges::makeGRangesFromDataFrame(df = E123_ratio_subset, 
#                                                              start.field = "position",
#                                                              end.field = "position",
#                                                              keep.extra.columns = T)
# 
# E123_segs_gr <- GenomicRanges::makeGRangesFromDataFrame(E123_segs, keep.extra.columns = T)
# 
# hits <- findOverlaps(subject = E123_ratio_subset, query = E123_segs_gr)
# 
# mcols(E123_ratio_subset)["purple_cn"] <- NA
# mcols(E123_ratio_subset)[["purple_cn"]][subjectHits(hits)] <- mcols(E123_segs_gr)[["tumorCopyNumber"]][queryHits(hits)]
# 
# E123_ratio_subset <- GenomicRanges::as.data.frame(E123_ratio_subset) 
# 
# ggplot(E123_ratio_subset, 
#        aes(x=start, y=log_ratio, color=factor(round(purple_cn,0), levels=1:5))) + geom_point(size=0.1) + 
#   facet_wrap("seqnames", scales = "free_x", nrow = 1) + 
#   coord_cartesian(ylim=c(-1,1)) +
#   geom_hline(yintercept = CN)

#########
# Plots #
#########

segplots <- list()
for (GOI in c("AFF4","TERT"))
{
  segplots[[GOI]] <- ggplot(probes[[GOI]], 
                            aes(x=position, 
                                y=log_ratio))  +
    geom_point(size=1, mapping = aes(color=purple_cn)) +
    geom_segment(x=regions[[GOI]]$start-(((regions[[GOI]]$end - regions[[GOI]]$start))*0.1), 
                 xend=regions[[GOI]]$end, 
                 y=-0.5,
                 yend=-0.5,
                 arrow = arrow(angle = 30, length = unit(0.2, "cm"),
                               ends = "first", type = "closed")
    ) +
    geom_segment(data=exon_positions %>% filter(external_gene_name==GOI), 
                 aes(x=exon_chrom_start,
                     xend=exon_chrom_end), 
                 y=-0.5,
                 yend=-0.5, 
                 linewidth=2) +
    annotate(geom = "text", 
             x = regions[[GOI]]$start+((regions[[GOI]]$end-regions[[GOI]]$start)/2), 
             y=-0.55, 
             label=GOI) +
    geom_segment(x=TERT_AFF4_SV[GOI], 
                 xend=TERT_AFF4_SV[GOI], 
                 y=CN[["CN2"]], 
                 yend=CN[["CN5"]],
                 linetype=2) +
    geom_curve(x=TERT_AFF4_SV[GOI], 
               xend=ifelse(GOI=="TERT", max(probes[[GOI]]$position), min(probes[[GOI]]$position)), 
               y=CN[["CN5"]], 
               yend=CN[["CN5"]]+0.3, 
               curvature = ifelse(GOI=="TERT", -0.3, 0.25),
               angle = 95,
               linetype=2) +
    theme_bw() +
    #theme(line = element_blank()) +
    #theme(panel.grid.major.y = element_line()) +
    #scale_y_continuous(breaks = CN, minor_breaks = FALSE) +  
    scale_x_continuous(labels = function(x){x/10^6}) +  
    coord_cartesian(ylim=c(-0.55,0.7)) + 
    xlab("Chromosome Position (MB)")
}


segplots[["TERT"]] + segplots[["AFF4"]] + plot_layout(guides = "collect")


