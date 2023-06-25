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


###################
# Load a5_annotation #
###################

gs4_auth("aidan.flynn@umccr-radio-lab.page")
a5_anno <- read_sheet("1hnXdXI29KvvuLxsaTBID1-EbE7mTSk6bG05HcSFfhgo", col_types = "c")
a5_anno <- a5_anno %>% dplyr::rename(`A5_ID`=`A5 ID`)
a5_anno <- a5_anno %>% filter(!is.na(`A5_ID`))

################
# Primary Site #
################

plot.data.location <- a5_anno %>% dplyr::select(`Patient ID`, Primary_Location_Simplified) %>% distinct() %>% 
  mutate(Primary_Location_Simplified=gsub("Head_neck","Head and neck",Primary_Location_Simplified)) %>% 
  mutate(Primary_Location_Simplified=gsub("_.+","",Primary_Location_Simplified)) %>% 
  group_by(Primary_Location_Simplified) %>%  
  dplyr::count() %>% 
  mutate(Primary_Location_Simplified=factor(as.character(Primary_Location_Simplified), levels=c("Extraadrenal","Adrenal","Head and neck","Unspecified")))

gg_AnatomicalLocation <- ggplot(plot.data.location, 
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

###################
# Clinical Course #
###################

plot.data.clinicalcourse <- a5_anno %>% dplyr::select(`Patient ID`, Primary_Location_Simplified, tumour_metastasised) %>% distinct() %>% 
  mutate(Primary_Location_Simplified=gsub("Head_neck","Head and neck",Primary_Location_Simplified)) %>% 
  mutate(Primary_Location_Simplified=gsub("_.+","",Primary_Location_Simplified)) %>% 
  group_by(Primary_Location_Simplified, tumour_metastasised) %>%  
  dplyr::count() %>% 
  mutate(Primary_Location_Simplified=factor(as.character(Primary_Location_Simplified), levels=c("Extraadrenal","Adrenal","Head and neck","Unspecified"))) %>% 
  mutate(tumour_metastasised=factor(as.character(tumour_metastasised), levels=c("Yes","No","Short follow up"))) %>% 
  mutate(tumour_metastasised=forcats::fct_recode(tumour_metastasised, Metastasic="Yes", `No-Metastasis`="No"))

gg_clinicalcourse <- list()
for (site in plot.data.clinicalcourse$Primary_Location_Simplified)
{
  
temp <- plot.data.clinicalcourse %>% filter(Primary_Location_Simplified==site)
    
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

####################
# Metastatic state #
####################

plot.data.metstate <- a5_anno %>% dplyr::select(`Patient ID`, Primary_Location_Simplified, is_primary_or_met) %>% distinct() %>% 
  mutate(Primary_Location_Simplified=gsub("Head_neck","Head and neck",Primary_Location_Simplified)) %>% 
  mutate(Primary_Location_Simplified=gsub("_.+","",Primary_Location_Simplified)) %>% 
  group_by(Primary_Location_Simplified, is_primary_or_met) %>%  
  dplyr::count() %>% 
  mutate(Primary_Location_Simplified=factor(as.character(Primary_Location_Simplified), levels=c("Extraadrenal","Adrenal","Head and neck","Unspecified"))) %>% 
  mutate(is_primary_or_met=factor(as.character(is_primary_or_met), levels=c("Primary","Recurrent","Metastatic"))) %>% 
  mutate(is_primary_or_met=forcats::fct_recode(is_primary_or_met, Recurrence="Recurrent", Metastasis="Metastatic"))

gg_metstate <- list()
for (site in plot.data.clinicalcourse$Primary_Location_Simplified)
{
  
  temp <- plot.data.metstate %>% filter(Primary_Location_Simplified==site)
  
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
