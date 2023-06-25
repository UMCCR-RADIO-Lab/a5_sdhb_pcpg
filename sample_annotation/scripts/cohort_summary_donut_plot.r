
setwd("/g/data/pq08/projects/ppgl")
renv::activate("./a5")

library(googlesheets4)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(patchwork)


ColorPalette <- c(LightBlue1="#cadae8ff",LightBlue2="#abc4e2ff",LightBlue3="#8aa9d6ff",
                  DarkBlue1="#7884baff",DarkBlue2="#606d9fff",DarkBlue3="#4c5a88ff",
                  Purple1="#765b92ff",Purple2="#9370a5ff",Purple3="#b280a9ff",Purple4="#cc85b1ff",
                  LightRed1="#ee7474ff",LightRed2="#e2696aff",LightRed3="#de5353ff",
                  DarkRed1="#c25858ff",DarkRed2="#c04241ff",DarkOrange1="#c65a44ff",
                  DarkOrange2="#eb7b54ff",LightOrange1="#f18d73ff",LightOrange2="#f5a697ff",
                  LightBrown1="#f5b9b0ff",LightBrown2="#d6a097ff",DarkBrown1="#c17963ff",DarkBrown2="#96665aff",
                  DarkGreen1="#637b62ff",DarkGreen2="#6ea668ff",LightGreen1="#98ca8cff",LightGreen2="#e2eab5ff",
                  Yellow1="#fbf2adff",Yellow2="#fef8c6ff",Yellow3="#f4e764ff",Salmon="#f2c5a7ff",LightSalmon="#f2dec4ff",
                  LemonGrey1="#f7f0e2ff",LemonWhite1="#fefbe5ff",LemonWhite2="#f7f4dcff",LemonGrey2="#efecdcff",
                  LightGrey1="#e2e0d9ff",LightGrey2="#d6d6d4ff",DarkGrey1="#b4b4b4ff",DarkGrey2="#9a9999ff")

ggplot(data.frame(Color=factor(names(ColorPalette), levels=names(ColorPalette)),Hex=ColorPalette), 
       mapping = aes(x=Color,y=1,fill=Hex)) + 
  geom_tile() + 
  scale_fill_identity() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) )


ExcludeSamples <- c()#"E167-1"

##########
# Themes #
##########

nox <-  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
nogrid <- theme(panel.grid = element_blank()) 
noxgrid <- theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
noygrid <- theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

pie_theme <- function () { theme_bw() + theme(axis.text = element_blank(), panel.grid = element_blank(), axis.title.y = element_blank(), panel.border = element_blank(), axis.ticks = element_blank())}


######################
# Load a5_annotation #
######################

#clinical annotation
source("./a5/scripts/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

################
# Primary Site #
################

plot.data.location.patient <- a5_anno %>% dplyr::select(`Patient ID`, Primary_Location_Simplified) %>% distinct() %>% 
  mutate(Primary_Location_Simplified=gsub("Head_neck","Head and neck",Primary_Location_Simplified)) %>% 
  mutate(Primary_Location_Simplified=gsub("_.+","",Primary_Location_Simplified)) %>% 
  group_by(Primary_Location_Simplified) %>%  
  dplyr::count() %>% 
  mutate(Primary_Location_Simplified=factor(as.character(Primary_Location_Simplified), levels=c("Extraadrenal","Adrenal","Head and neck","Unspecified")))

plot.data.location.sample <- a5_anno %>% dplyr::select(`A5_ID`, Primary_Location_Simplified) %>% distinct() %>% 
  mutate(Primary_Location_Simplified=gsub("Head_neck","Head and neck",Primary_Location_Simplified)) %>% 
  mutate(Primary_Location_Simplified=gsub("_.+","",Primary_Location_Simplified)) %>% 
  group_by(Primary_Location_Simplified) %>%  
  dplyr::count() %>% 
  mutate(Primary_Location_Simplified=factor(as.character(Primary_Location_Simplified), levels=c("Extraadrenal","Adrenal","Head and neck","Unspecified")))



###################
# Clinical Course #
###################

plot.data.clinicalcourse.patient <- a5_anno %>% dplyr::select(`Patient ID`, Primary_Location_Simplified, tumour_metastasised) %>% distinct() %>% 
  mutate(Primary_Location_Simplified=gsub("Head_neck","Head and neck",Primary_Location_Simplified)) %>% 
  mutate(Primary_Location_Simplified=gsub("_.+","",Primary_Location_Simplified)) %>% 
  group_by(Primary_Location_Simplified, tumour_metastasised) %>%  
  dplyr::count() %>% 
  mutate(Primary_Location_Simplified=factor(as.character(Primary_Location_Simplified), levels=c("Extraadrenal","Adrenal","Head and neck","Unspecified"))) %>% 
  mutate(tumour_metastasised=factor(as.character(tumour_metastasised), levels=c("Yes","No","Short follow up"))) %>% 
  mutate(tumour_metastasised=forcats::fct_recode(tumour_metastasised, Metastasic="Yes", `No-Metastasis`="No"))

plot.data.clinicalcourse.sample <- a5_anno %>% dplyr::select(`A5_ID`, Primary_Location_Simplified, tumour_metastasised) %>% distinct() %>% 
  mutate(Primary_Location_Simplified=gsub("Head_neck","Head and neck",Primary_Location_Simplified)) %>% 
  mutate(Primary_Location_Simplified=gsub("_.+","",Primary_Location_Simplified)) %>% 
  group_by(Primary_Location_Simplified, tumour_metastasised) %>%  
  dplyr::count() %>% 
  mutate(Primary_Location_Simplified=factor(as.character(Primary_Location_Simplified), levels=c("Extraadrenal","Adrenal","Head and neck","Unspecified"))) %>% 
  mutate(tumour_metastasised=factor(as.character(tumour_metastasised), levels=c("Yes","No","Short follow up"))) %>% 
  mutate(tumour_metastasised=forcats::fct_recode(tumour_metastasised, Metastasic="Yes", `No-Metastasis`="No"))



####################
# Metastatic state #
####################

plot.data.metstate.patient <- a5_anno %>% dplyr::select(`Patient ID`, Primary_Location_Simplified, is_primary_or_met) %>% distinct() %>% 
  mutate(Primary_Location_Simplified=gsub("Head_neck","Head and neck",Primary_Location_Simplified)) %>% 
  mutate(Primary_Location_Simplified=gsub("_.+","",Primary_Location_Simplified)) %>% 
  group_by(Primary_Location_Simplified, is_primary_or_met) %>%  
  dplyr::count() %>% 
  mutate(Primary_Location_Simplified=factor(as.character(Primary_Location_Simplified), levels=c("Extraadrenal","Adrenal","Head and neck","Unspecified"))) %>% 
  mutate(is_primary_or_met=factor(as.character(is_primary_or_met), levels=c("Primary","Recurrent","Metastasis"))) %>% 
  mutate(is_primary_or_met=forcats::fct_recode(is_primary_or_met, Recurrence="Recurrent"))

plot.data.metstate.sample <- a5_anno %>% dplyr::select(`A5_ID`, Primary_Location_Simplified, is_primary_or_met) %>% distinct() %>% 
  mutate(Primary_Location_Simplified=gsub("Head_neck","Head and neck",Primary_Location_Simplified)) %>% 
  mutate(Primary_Location_Simplified=gsub("_.+","",Primary_Location_Simplified)) %>% 
  group_by(Primary_Location_Simplified, is_primary_or_met) %>%  
  dplyr::count() %>% 
  mutate(Primary_Location_Simplified=factor(as.character(Primary_Location_Simplified), levels=c("Extraadrenal","Adrenal","Head and neck","Unspecified"))) %>% 
  mutate(is_primary_or_met=factor(as.character(is_primary_or_met), levels=c("Primary","Recurrent","Metastasis"))) %>% 
  mutate(is_primary_or_met=forcats::fct_recode(is_primary_or_met, Recurrence="Recurrent", Metastasis="Metastasis"))


##############################
# Concentric ring style plot #
##############################

plot.data <- 
  bind_rows(plot.data.location.sample %>% mutate(ring="inner"), 
            plot.data.metstate.sample %>% mutate(ring="middle"), 
            plot.data.clinicalcourse.sample %>% mutate(ring="outer")) %>% ungroup() %>% 
  mutate(category=paste(Primary_Location_Simplified, is_primary_or_met, tumour_metastasised, sep="-"), category=gsub("-NA","",category)) %>% 
  mutate(category=factor(category, levels=c("Adrenal","Extraadrenal","Head and neck","Unspecified",
                                            "Adrenal-Metastasis","Adrenal-Primary","Adrenal-Recurrence",
                                            "Extraadrenal-Metastasis","Extraadrenal-Primary","Extraadrenal-Recurrence",
                                            "Head and neck-Primary",
                                            "Unspecified-Metastasis",
                                            "Adrenal-No-Metastasis","Adrenal-Short follow up","Adrenal-Metastasic",
                                            "Extraadrenal-No-Metastasis","Extraadrenal-Short follow up","Extraadrenal-Metastasic",
                                            "Head and neck-No-Metastasis","Head and neck-Short follow up","Head and neck-Metastasic",
                                            "Unspecified-Metastasic"
  ))) %>% 
  dplyr::select(ring, category, n) %>% 
  mutate(color=case_when(
    category == "Adrenal" ~ ColorPalette[["Purple3"]],
    category == "Extraadrenal" ~ ColorPalette[["DarkBlue3"]],
    category == "Head and neck" ~ ColorPalette[["DarkGreen2"]],
    category == "Unspecified" ~ ColorPalette[["LightGrey2"]],
    grepl("Metastasic", category) ~ ColorPalette[["LightRed1"]],
    grepl("No-Metastasis", category) ~ ColorPalette[["LightBlue2"]],
    grepl("Short follow up", category) ~ ColorPalette[["LightBrown2"]],
    grepl("Metastasis", category) ~ ColorPalette[["DarkOrange2"]],
    grepl("Recurrence", category) ~ ColorPalette[["Yellow1"]],
    grepl("Primary", category) ~ ColorPalette[["LightGreen2"]]
  ))

colors.temp <- plot.data %>% dplyr::select(category, color) %>% distinct()
colors.use <- colors.temp$color
names(colors.use) <- colors.temp$category

ggplot(plot.data, aes(x=ring, y=n, fill=category)) +
  geom_bar(stat = "identity") + #, width = c(rep(1.6,4),rep(0.2,8),rep(0.2,10))) + 
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  coord_polar("y", start=0) + 
  scale_fill_manual(values = colors.use) +
  pie_theme() + 
  ylab("")

########################
# Multi-pie style plot #
########################

gg_AnatomicalLocation <- ggplot(plot.data.location.patient, 
                                aes(x="", y=n, fill=Primary_Location_Simplified)) +
  geom_bar(width = 1, stat = "identity") + 
  geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
  coord_polar("y", start=5.5) + 
  pie_theme() + 
  ylab("Anatomical Site") + 
  #  guides(fill=F) + 
  scale_fill_manual(values=c(Extraadrenal=ColorPalette[["DarkBlue3"]], 
                             Adrenal=ColorPalette[["Purple3"]], 
                             `Head and neck`=ColorPalette[["DarkGreen2"]], 
                             Unspecified=ColorPalette[["LightGrey2"]]))

gg_clinicalcourse <- list()
for (site in plot.data.clinicalcourse.patient$Primary_Location_Simplified)
{
  
  temp <- plot.data.clinicalcourse.patient %>% filter(Primary_Location_Simplified==site)
  
  gg_clinicalcourse[[site]] <-
    ggplot(temp, 
           aes(x="", y=n, fill=tumour_metastasised)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
    coord_polar("y", start=0) + 
    pie_theme() + 
    ylab("") + 
    #guides(fill=F) + 
    scale_fill_manual(values=c(Metastasic=ColorPalette[["LightRed1"]], 
                               `No-Metastasis`=ColorPalette[["LightBlue2"]], 
                               `Short follow up`=ColorPalette[["LightBrown2"]]))
}

gg_metstate <- list()
for (site in plot.data.clinicalcourse.patient$Primary_Location_Simplified)
{
  
  temp <- plot.data.metstate.patient %>% filter(Primary_Location_Simplified==site)
  
  gg_metstate[[site]] <-
    ggplot(temp, 
           aes(x="", y=n, fill=is_primary_or_met)) +
    geom_bar(width = 1, stat = "identity") + 
    geom_text(aes(label = n), position = position_stack(vjust = 0.5)) +
    coord_polar("y", start=0) + 
    pie_theme() + 
    ylab("") + 
    #guides(fill=F) + 
    scale_fill_manual(values=c(Metastasis=ColorPalette[["DarkOrange2"]], 
                               Recurrence=ColorPalette[["Yellow1"]], 
                               Primary=ColorPalette[["LightGreen2"]]))
}



layout <- "
BCDE
FGHI
AAAA
AAAA
"

gg_AnatomicalLocation + theme(plot.margin = margin(0,0,0,0), legend.position = "bottom") +
  gg_clinicalcourse[["Extraadrenal"]] + theme(plot.margin = margin(0,0,0,0), legend.position = "bottom")  + 
  gg_clinicalcourse[["Unspecified"]] + theme(plot.margin = margin(0,0,0,0), legend.position = "bottom")  + 
  gg_clinicalcourse[["Head and neck"]] + theme(plot.margin = margin(0,0,0,0), legend.position = "bottom")  + 
  gg_clinicalcourse[["Adrenal"]] + theme(plot.margin = margin(0,0,0,0), legend.position = "bottom")  + 
  gg_metstate[["Extraadrenal"]] + theme(plot.margin = margin(0,0,0,0), legend.position = "bottom")  + 
  gg_metstate[["Unspecified"]] + theme(plot.margin = margin(0,0,0,0), legend.position = "bottom")  + 
  gg_metstate[["Head and neck"]] + theme(plot.margin = margin(0,0,0,0), legend.position = "bottom")  + 
  gg_metstate[["Adrenal"]] + theme(plot.margin = margin(0,0,0,0), legend.position = "bottom")  + 
  plot_layout(design=layout, guides = "collect")#, widths = c(1,1,10,2), heights = c(5,5,1))




, 
Metastasis=ColorPalette[["DarkOrange2"]], 
Recurrence=ColorPalette[["Yellow1"]], 
Primary=ColorPalette[["LightGreen2"]])



"Adrenal","Extraadrenal","Head and neck","Unspecified",
"Adrenal-Metastasis","Adrenal-Primary","Adrenal-Recurrence",
"Extraadrenal-Metastasis","Extraadrenal-Primary","Extraadrenal-Recurrence",
"Head and neck-Primary",
"Unspecified-Metastasis",
"Adrenal-No-Metastasis","Adrenal-Short follow up","Adrenal-Metastasic",
"Extraadrenal-No-Metastasis","Extraadrenal-Short follow up","Extraadrenal-Metastasic",
"Head and neck-No-Metastasis","Head and neck-Short follow up","Head and neck-Metastasic",
"Unspecified-Metastasic"

"Adrenal","Adrenal-Metastasis","Adrenal-Primary","Adrenal-Recurrence",
"Extraadrenal","Extraadrenal-Primary","Extraadrenal-Recurrence","Extraadrenal-Metastasis",
"Head and neck","Head and neck-Primary",
"Unspecified","Unspecified-Metastasis"