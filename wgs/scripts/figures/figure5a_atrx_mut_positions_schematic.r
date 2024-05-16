##########################################################
##########################################################
##                                                      ##
##       ATRX MUTATION POSITION SCHEMATIC               ##
##                                                      ##
##       Author: Aidan Flynn                            ##
##      Date: 21-01-2023                                ##
##                                                      ##
##      This script generates a visual schematic of     ##
##      the germline SDHB mutations in the A5 cohort    ##
##                                                      ##
##########################################################
##########################################################

setwd("/g/data/pq08/projects/ppgl/a5")

library(googlesheets4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

#################
# Lookup tables #
#################

aa_abrev <- setNames(c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","O","S","U","T","W","Y","V"),
                     c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Pyl","Ser","Sec","Thr","Trp","Tyr","Val"))


#############
# Functions #
#############

load_domains <- function(gs_account, 
                         use_cache=T, 
                         cache_file="/g/data/pq08/projects/ppgl/a5/offline_cache/atrx_domains.tsv") {
  if(!use_cache) {
    gs4_auth(gs_account)
    ATRX_Domains <- read_sheet("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY", 
                               sheet ="ATRX_Domains", 
                               col_types = 'cccdd')
    write.table(x = ATRX_Domains,file = cache_file,quote = F, sep="\t",row.names = F)
  } else {
    ATRX_Domains <- read.delim(cache_file, header=T)
  }
  
  return(ATRX_Domains)
}

load_mutation_events <- function(gs_account, 
                                 use_cache=T, 
                                 cache_file="/g/data/pq08/projects/ppgl/a5/offline_cache/atrx_changes.tsv") {
  if(!use_cache) {
    gs4_auth(gs_account)
    ATRX_Changes <- read_sheet("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY", 
                               sheet ="ATRX_Events", 
                               col_types = 'cddccddcclc') %>%  arrange(Genomic_coord_start)
    
    write.table(x = ATRX_Changes,file = cache_file,quote = F, sep="\t",row.names = F)
  } else {
    ATRX_Changes <- read.delim(cache_file, header=T)
  }
  
  return(ATRX_Changes)
}  

################
# Data Loading #
################

#colors
source("./sample_annotation/scripts/data_loaders/a5_color_scheme.r")


#clinical annotation
source("./sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#SDHB Domain annotation
ATRX_Domains <- load_domains(gs_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

#A5 germline SDHB alterations
ATRX_Changes <- load_mutation_events(gs_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

ATRX_Changes <- ATRX_Changes %>% inner_join(
  a5_anno %>% dplyr::rename(PatientID=`Patient ID`) %>% 
    dplyr::select(PatientID, tumour_metastasised) %>% distinct()
) %>% 
  mutate(tumour_metastasised=factor(ifelse(grepl("Short|Data",tumour_metastasised),
                                           "Short follow-up/No Data" , 
                                           tumour_metastasised), levels=c("No", "Short follow-up/No Data", "Yes"))) %>% 
  filter(PCGR_CONS)
ATRX_Changes <- ATRX_Changes %>% filter(Include)

#########################
# Fetch Coordinate data #
#########################

EnsIds <- c("ATRX"="ENSG00000085224","MAGT1"="ENSG00000102158")


mart <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", mirror = "asia")

exon_positions <- biomaRt::getBM(filters= "ensembl_gene_id",
                                 attributes= c("ensembl_gene_id",
                                               "external_gene_name",
                                               "exon_chrom_start",
                                               "exon_chrom_end"),
                                 values=EnsIds, mart = mart)



############
#          #
# Plotting #
#          #
############

genome_plot_x_min <- min(ATRX_Changes$Genomic_coord_start, na.rm=T) 
genome_plot_x_max <- max(ATRX_Changes$Genomic_coord_end, na.rm=T)
genome_plot_x_min <-  genome_plot_x_min - ((genome_plot_x_max - genome_plot_x_min) * 0.1)
genome_plot_x_max <-  genome_plot_x_max + ((genome_plot_x_max - genome_plot_x_min) * 0.1)

###
# Exons
###

ggExons <- ggplot()  +
  geom_segment(data=exon_positions %>%
                 group_by(external_gene_name) %>%
                 summarise(gene_start=min(exon_chrom_start),
                           gene_end=max(exon_chrom_end)),
               aes(x=gene_start, 
                   xend=gene_end, 
                   y=0,
                   yend=0
               )) +
  geom_segment(data=exon_positions, 
               aes(x=exon_chrom_start,
                   xend=exon_chrom_end), 
               y=0,
               yend=0, 
               linewidth=2) + 
  geom_text(data=exon_positions %>%
              group_by(external_gene_name) %>%
              summarise(gene_start=min(exon_chrom_start),
                        gene_end=max(exon_chrom_end)),
            aes(x=gene_start+((gene_end-gene_start)/2), label=external_gene_name,
                y=-0.1
            )) +
  theme_void() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5), 
        axis.ticks.x = element_line(), 
        axis.ticks.length = unit(x = 2, units = "mm"),
        axis.line.x.bottom = element_line(),
        plot.margin = margin(0,0,10,0)) + 
  scale_x_continuous(expand = c(0,0), trans = "reverse", labels =  function(x){ paste(x/10^6, "MB")}) +
  coord_cartesian(xlim = c(genome_plot_x_max, genome_plot_x_min)) #Min/Max reversed due to trans=reverse in scale

####
# Mutations 
####

#Genome
plot_data_mutation_genome <-
  ATRX_Changes %>% filter(Type != "Deletion") %>% 
  dplyr::select(Genomic_coord_start, cds_annotation) %>%
  distinct() %>% 
  ungroup() %>% 
  mutate(delta = Genomic_coord_start-lag(Genomic_coord_start, default=.$Genomic_coord_start[[1]])) %>% 
  mutate(offset = ifelse((is.na(delta) | delta < 2000), 3, 0)) %>% 
  mutate(offset = ifelse((delta < 2000 & offset==lag(offset, default=.$offset[[1]])), 0, offset)) %>% 
  mutate(offset = ifelse((delta < 2000 & offset==lag(offset, default=.$offset[[1]])), 1.5, offset)) %>% 
  dplyr::select(-delta)

ggATRX_genomic_label <- ggplot(plot_data_mutation_genome, 
                               aes(x=Genomic_coord_start, 
                                   y=offset, 
                                   label=cds_annotation)) + 
  geom_text(angle=90, size=3) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text.y = element_blank(), 
        #axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(),
        plot.margin=ggplot2::margin(1,2,0,2)) + 
  ylab("CDS Change") + 
  coord_cartesian(xlim = c(genome_plot_x_max, genome_plot_x_min), #Min/Max reversed due to trans=reverse in scale
                  ylim=c(-4,4)) +
  scale_x_continuous(expand = c(0,0), trans = "reverse") 

#Protein
plot_data_mutation_protein <-
  ATRX_Changes %>% filter(!(Type %in% c("Splice Region", "Splice Acceptor", "Deletion"))) %>% 
  dplyr::select(AAold, AAnew, AApos_Start) %>%
  distinct() %>%
  mutate(
    label = paste0(AAold, AApos_Start, AAnew),
    label = purrr::modify(.x = label,
                          .f = stringr::str_replace_all, aa_abrev)
  ) %>%
  group_by(AApos_Start) %>%
  summarise(label = paste(label, collapse = "/")) %>% 
  ungroup() %>% 
  mutate(delta = AApos_Start-lag(AApos_Start, default=.$AApos_Start[[1]])) %>% 
  mutate(offset = ifelse((is.na(delta) | delta < 8), 3, 0)) %>% 
  mutate(offset = ifelse((delta < 5 & offset==lag(offset, default=.$offset[[1]])), 0, offset)) %>% 
  dplyr::select(-delta)

ggATRX_AA_label <- ggplot(plot_data_mutation_protein, 
                          aes(x=AApos_Start, 
                              y=offset, 
                              label=label)) + 
  geom_text(angle=90, size=3) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text.y = element_blank(), 
        #axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(),
        plot.margin=ggplot2::margin(1,2,0,2)) + 
  ylab("AA Change") +
  coord_cartesian(xlim = c(0,2492), ylim=c(-1,1)) +
  scale_x_continuous(expand = c(0,0)) 


#####
# Deletions 
#####

ATRX_Changes.deletion <- ATRX_Changes %>% filter(Type=="Deletion")
ATRX_Changes.deletion <- ATRX_Changes.deletion  %>%  
  group_by(Genomic_coord_start, Genomic_coord_end, AApos_Start, AApos_End) %>% summarise(nPatients=n()) %>% 
  mutate(ylevel=1) %>% 
  arrange(AApos_Start, AApos_End) %>% ungroup()


ggPositionDel_genomic <- ggplot(ATRX_Changes.deletion) + 
  geom_segment(aes(x=Genomic_coord_start,y=ylevel, xend=Genomic_coord_end, yend=ylevel), linewidth=6) +
  geom_text(aes(x=(Genomic_coord_start+((Genomic_coord_end-Genomic_coord_start)/2)), y=ylevel, label=paste0("n=",nPatients)), color="white") +
  coord_cartesian(xlim = c(genome_plot_x_max, genome_plot_x_min)) + #Min/Max reversed due to trans=reverse in scale
  scale_x_continuous(expand = c(0,0), trans = "reverse") +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text.y = element_blank(), 
        #axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(),
        plot.margin=ggplot2::margin(1,2,0,2)) + 
  guides(color=F) +
  ylab("Deletions")


ggPositionDel_protein <- ggplot(ATRX_Changes.deletion) + 
  geom_segment(aes(x=AApos_Start,y=ylevel, xend=AApos_End, yend=ylevel), linewidth=6) +
  coord_cartesian(xlim = c(0,2492)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.text.y = element_blank(), 
        #axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(),
        plot.margin=ggplot2::margin(1,2,0,2)) + 
  guides(color=F) +
  ylab("Deletions")

#####
# Domains 
#####

ggDomains <- ggplot(ATRX_Domains) + 
  geom_segment(data=(. %>% filter(Source=="fill")), aes(x=Start,y=0, xend=End, yend=0), size=2) +
  geom_segment(data=(. %>% filter(!is.na(Domain))), aes(x=Start,y=0, xend=End, yend=0), size=6) +
  geom_text(data=(. %>% filter(!is.na(Domain))), aes(x=(Start+((End-Start)/2)), y=0, label=Domain), color="white") +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.x = element_line(),
        plot.margin=ggplot2::margin(0,2,0,2)) + guides(color=F) +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(0,2492)) +
  xlab("Amino Acid Position")

####
# Combine
####

ggATRX_genomic_label + 
  ggPositionDel_genomic + 
  ggExons +
  ggATRX_AA_label + 
  ggPositionDel_protein +
  ggDomains +
  plot_layout(nrow=6)






