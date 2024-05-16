##########################################################
##########################################################
##                                                      ##
##       SDHB GERMLINE MUTATION POSITION SCHEMATIC      ##
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

########
# Lookup tables
########

aa_abrev <- setNames(c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","O","S","U","T","W","Y","V"),
                     c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His","Ile","Leu","Lys","Met","Phe","Pro","Pyl","Ser","Sec","Thr","Trp","Tyr","Val"))


########
# Functions
########

load_domains <- function(gs_account, 
                         use_cache=T, 
                         cache_file="/g/data/pq08/projects/ppgl/a5/offline_cache/sdhb_domains.tsv") {
  if(!use_cache) {
    gs4_auth(gs_account)
    SDHB_Domains <- read_sheet("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY", 
                               sheet ="SDHB_Domains", 
                               col_types = 'cccdd')
    write.table(x = SDHB_Domains,file = cache_file,quote = F, sep="\t",row.names = F)
  } else {
    SDHB_Domains <- read.delim(cache_file, header=T)
  }
  
  return(SDHB_Domains)
}

load_germline_events <- function(gs_account, 
                         use_cache=T, 
                         cache_file="/g/data/pq08/projects/ppgl/a5/offline_cache/sdhb_changes.tsv") {
  if(!use_cache) {
    gs4_auth(gs_account)
    SDHB_Changes <- read_sheet("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY", 
                               sheet ="SDHB_Events", 
                               col_types = 'ccddccdcl')
    write.table(x = SDHB_Changes,file = cache_file,quote = F, sep="\t",row.names = F)
  } else {
    SDHB_Changes <- read.delim(cache_file, header=T)
  }
  
  return(SDHB_Changes)
}  

########
# Data Loading
########

#colors
source("./sample_annotation/scripts/data_loaders/a5_color_scheme.r")

#clinical annotation
source("./sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#SDHB Domain annotation
SDHB_Domains <- load_domains(gs_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

#A5 germline SDHB alterations
SDHB_Changes <- load_germline_events(gs_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

SDHB_Changes <- SDHB_Changes %>%  inner_join(
  a5_anno %>% dplyr::rename(PatientID=`Patient ID`) %>% 
    dplyr::select(PatientID, tumour_metastasised) %>% distinct()
) %>% 
  mutate(tumour_metastasised=factor(ifelse(grepl("Short|Data",tumour_metastasised),
                                             "Short follow-up/No Data" , 
                                             tumour_metastasised), levels=c("No", "Short follow-up/No Data", "Yes")))
SDHB_Changes <- SDHB_Changes %>% filter(Include)

####
# Mutations 
####

clinvar_levels <- c("Pathogenic",
                    "Pathogenic/Likely pathogenic",
                    "Likely pathogenic",  
                    "Benign/Likely benign",
                    "Conflicting",
                    "Uncertain",
                    "No record")

clinvar_colors <- setNames(c(ColorPalette[["DarkRed2"]],
                             ColorPalette[["LightOrange1"]],
                             ColorPalette[["LightBrown1"]],
                             "#ffcc33",
                             "black",
                             ColorPalette[["Purple1"]],
                             ColorPalette[["DarkGrey1"]]), clinvar_levels)

SDHB_Changes.small_variants.clincourse <- 
  SDHB_Changes %>%  
  filter(Type %in% c("Missense","Stop", "Frameshift")) %>%  
  group_by(AAold, AAnew, AApos_Start, AApos_End, Clinvar, tumour_metastasised) %>% summarise(nPatients=n()) %>% 
  mutate(Clinvar = case_when(
    grepl("Uncertain", Clinvar) ~ "Uncertain",
    grepl("Conflicting", Clinvar) ~ "Conflicting",
    TRUE ~ Clinvar
  ))

plot.data <-
  SDHB_Changes.small_variants.clincourse %>%
  dplyr::select(AAold, AAnew, AApos_Start, Clinvar) %>%
  distinct() %>%
  mutate(
    label = paste0(AAold, AApos_Start, AAnew),
    label = purrr::modify(.x = label,
                          .f = stringr::str_replace_all, aa_abrev)
  ) %>%
  group_by(AApos_Start, Clinvar) %>%
  summarise(label = paste(label, collapse = "/")) %>% 
  ungroup() %>% 
    mutate(delta = AApos_Start-lag(AApos_Start, default=.$AApos_Start[[1]])) %>% 
  mutate(offset = ifelse((is.na(delta) | delta < 8), 3, 0)) %>% 
  mutate(offset = ifelse((delta < 5 & offset==lag(offset, default=.$offset[[1]])), 0, offset)) %>% 
  dplyr::select(-delta) %>% 
  mutate(Clinvar=factor(Clinvar, levels=clinvar_levels))

ggSDHB_AA_label <- ggplot(plot.data, 
                          aes(x=AApos_Start, 
                              y=offset, 
                              label=label, color=Clinvar)) + 
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
  coord_cartesian(xlim = c(0,280), ylim=c(-2,5)) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_color_manual(values = clinvar_colors)



ggPositionFreq.clincourse <- ggplot(SDHB_Changes.small_variants.clincourse %>% mutate(AApos_Start=factor(AApos_Start, levels=c(1:280)))) + 
  geom_col(aes(x=AApos_Start, y=nPatients, fill=tumour_metastasised)) +
  #geom_text(data=SDHB_Changes, aes(x=AApos_Start, y=1.6, label=paste(AAold,AAnew, sep="/"))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        #axis.text.y = element_blank(), 
        #axis.title.y = element_blank(), 
        #axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(linetype = 2, colour = ColorPalette[["LightGrey1"]]),
        plot.margin=ggplot2::margin(2,2,1,2)) + guides(color="none") +
  scale_fill_manual(values=c(ColorPalette[["LightGreen1"]],ColorPalette[["DarkBlue1"]],ColorPalette[["LightRed1"]])) +
  scale_x_discrete(drop=F) +
  scale_y_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(0,280))

SDHB_Changes.small_variants.type <- SDHB_Changes %>%  filter(Type %in% c("Missense","Stop", "Frameshift")) %>%  group_by(AAold, AAnew, AApos_Start, AApos_End, Type) %>% summarise(nPatients=n())

ggPositionFreq.type <- ggplot(SDHB_Changes.small_variants.type %>% mutate(AApos_Start=factor(AApos_Start, levels=c(1:280)))) + 
  geom_col(aes(x=AApos_Start, y=nPatients, fill=Type)) +
  #geom_text(data=SDHB_Changes, aes(x=AApos_Start, y=1.6, label=paste(AAold,AAnew, sep="/"))) +
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        #axis.text.y = element_blank(), 
        #axis.title.y = element_blank(), 
        #axis.ticks.y = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        axis.line.y = element_line(),
        axis.line.x = element_line(linetype = 2, colour = ColorPalette[["LightGrey1"]]),
        plot.margin=ggplot2::margin(2,2,1,2)) + guides(color=F) +
  scale_fill_manual(values=c(ColorPalette[["Yellow1"]],ColorPalette[["LightOrange1"]],ColorPalette[["Purple1"]])) +
  scale_x_discrete(drop=F, position = "top") +
  scale_y_reverse(expand = c(0,0)) +
  coord_cartesian(xlim = c(0,280))


####
# Splice 
####

SDHB_Changes.splice.clincourse <- 
  SDHB_Changes %>%  
  filter(Type %in% c("Splice donor" )) %>%  #, "Splice site"
  group_by(AAold, AAnew, AApos_Start, AApos_End, Clinvar, tumour_metastasised) %>% summarise(nPatients=n()) %>% 
  mutate(Clinvar = case_when(
    grepl("Uncertain", Clinvar) ~ "Uncertain",
    grepl("Conflicting", Clinvar) ~ "Conflicting",
    TRUE ~ Clinvar
  ))

plot.data <-
  SDHB_Changes.splice.clincourse %>%
  dplyr::select(AAold, AAnew, AApos_Start, AApos_End, Clinvar) %>%
  distinct() %>%
  mutate(
    label = AAnew #paste0(AAold, AApos_Start, AAnew),
    #label = purrr::modify(.x = label,
    #                      .f = stringr::str_replace_all, aa_abrev)
  ) %>%
  group_by(AApos_Start, AApos_End, Clinvar) %>%
  summarise(label = paste(label, collapse = "/")) %>% 
  ungroup() %>% 
  mutate(Clinvar=factor(Clinvar, levels=clinvar_levels),
         label=gsub("__","\n",label)) %>% 
  pivot_longer(cols = c(AApos_Start,AApos_End),names_to = "Start_End", values_to = "Pos")

ggSDHB_splice <- ggplot(plot.data, 
                          aes(x=Pos, 
                              y=1, 
                              label=label, color=Clinvar,
                              group=label)) + 
  geom_text(size=3, nudge_y = -1) + 
  geom_point() +
  geom_line() +
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
  ylab("Splice Donor Mutations") +
  coord_cartesian(xlim = c(0,280), ylim=c(-1,1)) +
  scale_x_continuous(expand = c(0,0)) + 
  scale_color_manual(values = clinvar_colors)



#####
# Deletions 
#####

SDHB_Changes.deletion <- SDHB_Changes %>% filter(Type=="Deletion")
SDHB_Changes.deletion <- SDHB_Changes.deletion  %>%  
  group_by(AApos_Start, AApos_End, tumour_metastasised) %>% summarise(nPatients=n()) %>% 
  mutate(ylevel=1) %>% 
  arrange(AApos_Start, AApos_End) %>% ungroup()

does_overlap <- function(start_query, end_query, start_reference, end_reference)
{
  if(length(start_query) > 1)
  {
    return_vector <- vector(mode="logical", length = length(start_query))
    for (i in 1:length(start_query))
    {
      return_vector[[i]] <- does_overlap(start_query[[i]], end_query[[i]], start_reference, end_reference)
    }
    return(return_vector)
  }
  else
  {
    s=start_query; S=start_reference; e=end_query; E=end_reference;
    
    if (
      (s<=S & s<=E & e>=S & e<=E) | #Starts Before. Ends During.
      (s>=S & s<=E & e>=S & e>=E) | #Starts During. Ends After
      (s>=S & s<=E & e>=S & e<=E) | #Starts During. Ends During
      (s<=S & s<=E & e>=S & e>=E))  #Starts Before. Ends After
    { return(TRUE)} else 
    {return(FALSE)}
  }
}

for (i in 2:nrow(SDHB_Changes.deletion))
{
  overlaps <- SDHB_Changes.deletion %>% 
    filter(row_number() != i) %>% 
    mutate(does_overlap=does_overlap(
      start_query = AApos_Start, 
      end_query=AApos_End, 
      start_reference=SDHB_Changes.deletion$AApos_Start[[i]], 
      end_reference=SDHB_Changes.deletion$AApos_Start[[i]]
    )) %>% filter(does_overlap)
  
  if(nrow(overlaps)>0)
  {
    new_y_level <- max(overlaps$ylevel) + 1
    SDHB_Changes.deletion$ylevel[[i]] <- new_y_level
  } 
  
}

ggPositionDel <- ggplot(SDHB_Changes.deletion) + 
  geom_segment(aes(x=AApos_Start,y=ylevel, xend=AApos_End, yend=ylevel, color=tumour_metastasised), linewidth=6, linejoin = "round",  lineend = "round") +
  geom_text(aes(x=(AApos_Start+((AApos_End-AApos_Start)/2)), y=ylevel, label=paste0("n=",nPatients))) +
  scale_color_manual(values=c(ColorPalette[["LightGreen1"]],ColorPalette[["DarkBlue1"]],ColorPalette[["LightRed1"]])) +
  coord_cartesian(xlim = c(0,280)) +
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

ggDomains <- ggplot(SDHB_Domains) + 
  geom_segment(data=(. %>% filter(Source=="fill")), aes(x=Start,y=0, xend=End, yend=0), size=2) +
  geom_segment(data=(. %>% filter(!is.na(Domain))), aes(x=Start,y=0, xend=End, yend=0), size=6) +
  geom_text(data=(. %>% filter(!is.na(Domain))), aes(x=(Start+((End-Start)/2)), y=0, label=Domain), color="white") +
  #geom_text(data=SDHB_Changes, aes(x=AApos_Start, y=1.6, label=paste(AAold,AAnew, sep="/"))) +
  theme_bw() + 
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.line.x = element_line(),
        plot.margin=ggplot2::margin(0,2,0,2)) + guides(color=F) +
  scale_x_continuous(expand = c(0,0)) +
  coord_cartesian(xlim = c(0,280)) +
  xlab("Amino Acid Position")

####
# Combine
####
pdf(file = "/g/data/pq08/projects/ppgl/a5/sample_annotation/figures/sdhb_germline_schematic.pdf", width = 12,height = 7)
ggSDHB_AA_label + 
  ggPositionFreq.clincourse + 
  ggSDHB_splice  + 
  ggPositionDel + 
  ggDomains + 
  plot_layout(nrow=5, heights = c(5,8,3,3,2), guides ="collect")
dev.off()
#ggSDHB_AA_label/ggPositionFreq.clincourse/ggPositionDel/ggDomains/ggPositionFreq.type + plot_layout(heights = c(5,8,3,2,8), guides ="collect")
