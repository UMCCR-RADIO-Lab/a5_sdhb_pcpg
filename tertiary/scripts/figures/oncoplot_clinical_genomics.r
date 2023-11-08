library(googlesheets4)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(patchwork)

################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_somatic_variants(quickload = T)
# data_loader_germline_variants()
data_loader_gene_cn()
data_loader_cna_segs()
data_loader_sv_gridsslinx()


#Annotate Gene CN as loss/gain
A5_gene_cn_lossgain <- A5_gene_cn %>% mutate(A5_ID=gsub("T0", "", A5_ID)) %>%  
  mutate(loss_gain = case_when(
  !(chromosome %in% c("chrX","chrY")) & minCopyNumber < 1.2 ~ "Loss", 
  !(chromosome %in% c("chrX","chrY")) & maxCopyNumber_ploidyAdj > 5.2 ~ "Gain",
  !(chromosome %in% c("chrX","chrY")) & minCopyNumber > 1.8 & minMinorAlleleCopyNumber < 0.2 ~ "CNLOH",
  chromosome =="chrX" & Gender == "Female" & minCopyNumber < 1.2 ~ "Loss", 
  chromosome =="chrX" & Gender == "Female" & maxCopyNumber_ploidyAdj > 2.2 ~ "Gain",
  chromosome =="chrX" & Gender == "Male" & minCopyNumber < 0.2 ~ "Loss", 
  chromosome =="chrX" & Gender == "Male" & maxCopyNumber_ploidyAdj > 1.2 ~ "Gain",
  TRUE ~ NA_character_)) %>% 
  filter(!is.na(loss_gain), chromosome != "chrY") %>%  
  dplyr::select(A5_ID, gene, loss_gain)

source("/g/data/pq08/projects/ppgl/a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
data_loader_a5_wts_counts(count_file_dir = "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/truseq/gene", 
                          count_file_pattern = ".gene.counts")
# data_loader_arriba(arriba_out_dir = "/g/data/pq08/projects/ppgl/a5/wts/analysis/arriba/truseq/")

wts_log_cpm <- edgeR::cpm(a5_wts_dge_list$qc_ok, log=T) %>% 
  data.frame(check.names = F) %>% 
  tibble::rownames_to_column("ensgid_symbol") %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), sep="_", extra = "merge") %>% 
  pivot_longer(cols = c(-ensgid,-symbol), names_to="A5_ID", values_to="log2cpm") %>% 
  group_by(ensgid) %>% 
  mutate(log2cpm_z=(log2cpm-mean(log2cpm, na.rm=T))/sd(log2cpm, na.rm=T)) 


#source("/g/data/pq08/projects/ppgl/a5/methylation/scripts/data_loaders/a5_methylation_dataloader.r")


##########
# Themes #
##########

nox <-  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
nogrid <- theme(panel.grid = element_blank()) 
noxgrid <- theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
noygrid <- theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())

no_margin <- theme(plot.margin = margin(0,6,0,6))
shrink_legend <- theme(legend.key.size = unit(0.5,"cm"), legend.position = "top")

##########
# Colors #
##########

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")
location_cols[["Extraadrenal (abdominal/thoracic)"]] <-  location_cols[["Extraadrenal"]]
location_cols[["Extraadrenal (aortic)"]] <-  location_cols[["Extraadrenal_aortic"]]
location_cols[["Extraadrenal (bladder)"]] <-  location_cols[["Extraadrenal_bladder"]]

###################
# Sample Ordering #
###################

a5_anno.use <- a5_anno %>% filter(Exclude != "Y") %>% 
  filter(!(A5_ID %in% c("E124-1", "E145-1"))) 

a5_anno.use$is_primary_or_met <- factor(as.character(a5_anno.use$is_primary_or_met), levels=c("Primary","Recurrent","Metastasis"))
a5_anno.use$tumour_metastasised <- factor(as.character(a5_anno.use$tumour_metastasised), levels=c("Yes","No", "Short follow up", "No Data"))
a5_anno.use$`Patient ID` <- factor(as.character(a5_anno.use$`Patient ID`), 
                                   levels=unique(as.character(a5_anno.use$`Patient ID`))[order(as.numeric(gsub("E","",unique(as.character(a5_anno.use$`Patient ID`)))))])
a5_anno.use$TERT_ATRX_Mutation <- factor(as.character(a5_anno.use$TERT_ATRX_Mutation), levels=c("TERT","ATRX", "Unknown"))
# a5_anno.use$A5_ID <- factor(as.character(a5_anno.use$A5_ID), 
#                      levels=unique(as.character(a5_anno.use$`A5_ID`))[order(as.numeric(gsub("E([0-9]{3})-([0-4])","\\1",unique(as.character(a5_anno.use$`A5_ID`)))),
#                                                                      as.numeric(gsub("E([0-9]{3})-([0-4])","\\2",unique(as.character(a5_anno.use$`A5_ID`)))))])

a5_anno.use <- a5_anno.use %>% 
  mutate(Primary_Location_Simplified=dplyr::recode(Primary_Location_Simplified,
                                                   "Extraadrenal_abdominal"="Extraadrenal (abdominal/thoracic)",
                                                   "Extraadrenal_thoracic"="Extraadrenal (abdominal/thoracic)",
                                                   "Extraadrenal_thoracic_cardiac"="Extraadrenal (abdominal/thoracic)",
                                                   "Extraadrenal_bladder"="Extraadrenal (bladder)",
                                                   "Extraadrenal_thoracic_aortic"="Extraadrenal (aortic)",
                                                   "Adrenal_left"="Adrenal",
                                                   "Adrenal_right"="Adrenal",
                                                   "Head_neck"="Head and neck"),
         Primary_Location_Simplified=factor(Primary_Location_Simplified,
                                            levels=c("Extraadrenal (abdominal/thoracic)",
                                                     "Extraadrenal (bladder)",
                                                     "Adrenal",
                                                     "Unspecified",
                                                     "Extraadrenal (aortic)",
                                                     "Head and neck")))

a5_anno.use <- a5_anno.use %>% group_by(`Patient ID`) %>% 
  mutate(max_clinical=case_when(any(differential_group_sampletype_strict=="Metastasis") ~ "Metastasis",
                                any(differential_group_sampletype_strict=="Metastatic primary") ~ "Metastatic primary",
                                any(differential_group_sampletype_strict=="Metastatic local recurrence") ~ "Metastatic local recurrence",
                                any(differential_group_sampletype_strict=="Primary (metastasis present)") ~ "Primary (metastasis present)",
                                any(differential_group_sampletype_strict=="Local recurrence (metastasis present)") ~ "Local recurrence (metastasis present)",
                                any(differential_group_sampletype_strict=="Non-metastatic local recurrence") ~ "Non-metastatic local recurrence",
                                any(differential_group_sampletype_strict=="Primary (short follow up)") ~ "Primary (short follow up)",
                                any(differential_group_sampletype_strict=="Non-metastatic primary") ~ "Non-metastatic primary"),
         max_clinical=factor(max_clinical, levels=c("Metastasis",
                                                    "Metastatic primary",
                                                    "Metastatic local recurrence",
                                                    "Local recurrence (metastasis present)",
                                                    "Primary (metastasis present)",
                                                    "Non-metastatic local recurrence",
                                                    "Non-metastatic primary"))) %>% 
  arrange(desc(is_head_and_neck), max_clinical, `Patient ID`, differential_group_sampletype_strict)
SampleOrder.A5_ID <- as.character(a5_anno.use$A5_ID)
a5_anno.use$A5_ID <- factor(as.character(a5_anno.use$A5_ID), levels = SampleOrder.A5_ID)


##############
##############
## CLINICAL ##
##############
##############


##################
# Disease Course #
##################

# plot.data.metastasis <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, tumour_metastasised)
# 
# 
# ggDiseaseCourse <- ggplot() + 
#   geom_tile(data=plot.data.metastasis, mapping=aes(x=A5_ID, y="Disease Course", fill=tumour_metastasised), width=0.8,height=0.8) +
#   scale_fill_manual(values=c(ColorPalette[["DarkGreen1"]],ColorPalette[["LightGreen1"]],ColorPalette[["LemonGrey1"]],ColorPalette[["LightGrey1"]])) +
#   scale_x_discrete(drop=F) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank())

###############
# Sample Type #
###############

plot.data.sampletype <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, differential_group_sampletype_strict)
#plot.data.linker <- plot.data.sampletype %>% group_by(`Patient ID`) %>% mutate(membercount=n(), idx=as.numeric(A5_ID)) %>%  filter(membercount>1) %>%  group_by(`Patient ID`) %>%  dplyr::slice(which(row_number()==1 | row_number()==n()))

plot.data.sampletype <- 
  plot.data.sampletype %>% 
  left_join(plot.data.sampletype %>% 
              group_by(`Patient ID`) %>% 
              filter(n() >1) %>% 
              group_by(`Patient ID`) %>%  
              mutate(group=cur_group_id()))

ggSampleType <- ggplot() + 
  geom_tile(data=plot.data.sampletype, mapping=aes(x=A5_ID, y="Specimen Type", fill=differential_group_sampletype_strict), width=0.8,height=0.8) +
  geom_line(data=plot.data.sampletype %>% filter(!is.na(group)), mapping = aes(x=A5_ID, group=`Patient ID`, y="Specimen Type"), linewidth=0.8) +
  geom_point(data=plot.data.sampletype %>% filter(!is.na(group)), mapping = aes(x=A5_ID, y="Specimen Type"), size=1.2) +
  #geom_point(data=plot.data.linker, mapping = aes(x=A5_ID, y="Specimen Type"), size=1.2) +
  #geom_point(data=plot.data.sampletype, mapping = aes(x=A5_ID, y="Specimen Type", shape=group), size=1.2) +
  #scale_fill_manual(values=c(ColorPalette[["LightBlue2"]],ColorPalette[["LightOrange2"]],ColorPalette[["LightRed1"]])) +
  scale_shape_identity() + 
  scale_x_discrete(drop=F) +
  scale_fill_manual(name = "Clinical Course", values =sampletype_strict_cols) +
  # scale_fill_manual(
  #   name = "Clinical Course",
  #   values = c("#6ea668", "#6ea6cc", "#d6d6d4", "#f18d73", "#ee7474ff", "#969696", "#c04241ff"),
  #   labels = c("Non-metastatic primary", "Primary (short follow up)", "Non-metastatic local recurrence",
  #              "Primary (metastasis present)", "Metastatic primary","Metastatic local recurrence","Metastasis")
  # ) + 
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank())

# ggSampleType.monotone <- ggplot() + 
#   geom_tile(data=plot.data.sampletype, mapping=aes(x=A5_ID, y="Specimen Type", fill=differential_group_sampletype_strict), width=0.8,height=0.8) +
#   #geom_line(data=plot.data.linker, mapping = aes(x=A5_ID, group=`Patient ID`, y="Specimen Type"), size=0.8) +
#   #geom_point(data=plot.data.linker, mapping = aes(x=A5_ID, y="Specimen Type"), size=1.2) +
#   geom_point(data=plot.data.sampletype, mapping = aes(x=A5_ID, y="Specimen Type", shape=group), size=1.2) +
#   scale_fill_manual(values=c(ColorPalette[["LightGrey1"]],ColorPalette[["DarkGrey2"]])) +
#   scale_shape_identity() + 
#   scale_x_discrete(drop=F) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank())


##################
# Post Treatment #
##################

plot.data.treatment <- a5_anno.use %>% 
  dplyr::select(A5_ID, `Patient ID`, Resection_post_dna_damaging_treatment) %>% 
  mutate(Resection_post_dna_damaging_treatment=recode(Resection_post_dna_damaging_treatment,
                                                      "Unknown" = "Uncertain",
                                                      "Possibly (CVD)" = "Uncertain",
                                                      "Yes (CVD - low)" = "Yes",
                                                      "Yes (CVD)" = "Yes",
                                                      "Yes (CVD, MIBG)" = "Yes",
                                                      "Yes (MIBG)" = "Yes",
                                                      "Yes (MIBG,Carb,Cis, historic)" = "Yes"))

ggTreatment.monotone <- ggplot() + 
  geom_tile(data=plot.data.treatment, 
            mapping=aes(x=A5_ID, y="Resct. post treat.", 
                        fill=Resection_post_dna_damaging_treatment), 
            width=0.8,
            height=0.8) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, 
                                   vjust = 0.5,
                                   hjust=1), 
        panel.grid = element_blank()) + 
  labs(fill="Resct. post treat.") +
  scale_x_discrete(drop=F) +
  scale_fill_manual(values=c(
    "No"=ColorPalette[["DarkGrey2"]],
    "Yes"="black", 
     "Uncertain"= ColorPalette[["LightGrey1"]]))




###############
# Tumour Size #
###############

plot.data.size <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, Largest_primary_category)
plot.data.size$Largest_primary_category <- tidyr::replace_na(plot.data.size$Largest_primary_category, "Not Available")


# ggSize <- ggplot() + 
#   geom_tile(data=plot.data.size, mapping=aes(x=A5_ID, y="Largest Primary", fill=Largest_primary_category), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Largest Primary") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values=as.character(c(ColorPalette["DarkBlue3"],ColorPalette["Purple3"],ColorPalette["LightGrey2"])))

ggSize.monotone <- ggplot() + 
  geom_tile(data=plot.data.size, mapping=aes(x=A5_ID, y="Largest Primary", fill=Largest_primary_category), width=0.8,height=0.8) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Largest Primary") +
  scale_x_discrete(drop=F) +
  scale_fill_manual(values=as.character(c("<5cm"=ColorPalette[["DarkGrey2"]],">5cm"="black", "No data"=ColorPalette[["LightGrey1"]])))


############
# Location #
############

plot.data.location <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, Primary_Location_Simplified) 


ggLocation <- ggplot() + 
  geom_tile(data=plot.data.location, mapping=aes(x=A5_ID, y="Tumour Location", fill=Primary_Location_Simplified), width=0.8,height=0.8) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Tumour Location") +
  scale_x_discrete(drop=F) +
  scale_fill_manual(
    name = "Anatomical Location (Primary)",
    values = location_cols)  

##################
# Head And Neck  #
##################
# 
# plot.data.handn <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, Primary_Location_Simplified) %>%  
#   mutate(TumourLocation=case_when(
#     Primary_Location_Simplified=="Head_neck" ~ "Head and Neck PGL",  
#     Primary_Location_Simplified=="Extraadrenal_thoracic_aortic" ~ "Aortic PGL",  
#     Primary_Location_Simplified=="Unspecified" ~ "No Data",  
#     TRUE ~ "PCC/PGL"))
# plot.data.handn$TumourLocation <- factor(plot.data.handn$TumourLocation, 
#                                          levels = c("PCC/PGL", "Head and Neck PGL", 
#                                                     "Aortic PGL", "No Data"))
# 
# 
# ggHandN <- ggplot() + 
#   geom_tile(data=plot.data.handn, mapping=aes(x=A5_ID, y="Tumour Location", fill=TumourLocation), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Tumour Location") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values=as.character(c("PCC/PGL"=ColorPalette[["LightBrown2"]], 
#                                           "Head and Neck PGL"=ColorPalette[["DarkBrown2"]],
#                                           "Aortic PGL"=ColorPalette[["DarkBrown1"]],
#                                           "No Data"=ColorPalette[["DarkGrey2"]])))
# 
# ggHandN.monotone <- ggplot() + 
#   geom_tile(data=plot.data.handn, mapping=aes(x=A5_ID, y="Tumour Location", fill=TumourLocation), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Tumour Location") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values=as.character(c("PCC/PGL"=ColorPalette[["DarkGrey2"]], 
#                                           "Head and Neck PGL"=ColorPalette[["LightGrey1"]],
#                                           "Aortic PGL"="black",
#                                           "No Data"=ColorPalette[["LemonWhite2"]])))

#########################
# Catecholamine profile #
#########################

plot.data.catecholamine <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, `Catecholamine_profile`) %>% 
  mutate(`Catecholamine_profile`=case_when(
    grepl("Multiple.+| and ", `Catecholamine_profile`) ~"Mixed",
    `Catecholamine_profile`=="Non secreting" ~ "Non-secreting",
    `Catecholamine_profile`=="Norepinephrine*" ~ "Norepinephrine",
    TRUE ~ `Catecholamine_profile`)) %>% 
  mutate(`Catecholamine_profile`=factor(as.character(`Catecholamine_profile`), 
                                        levels=c("Norepinephrine", "Epinephrine", "Dopamine", "Mixed", "Non-secreting", "No data"))) %>% 
  mutate(code=forcats::fct_recode(`Catecholamine_profile`, N="Norepinephrine", E="Epinephrine", D="Dopamine", M="Mixed", "*"="Non-secreting"))

##Colored Tiles
# ggCatecholamine.color <- ggplot() + 
#   geom_tile(data=plot.data.catecholamine, mapping=aes(x=A5_ID, y="Catecholamines", fill=`Catecholamine_profile`), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Catecholamine Profile") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values=as.character(c(ColorPalette["DarkBlue3"],ColorPalette["Purple3"],
#                                           ColorPalette["LightBlue2"],ColorPalette["LightRed1"],
#                                           ColorPalette["LightGreen1"],ColorPalette["LightGrey2"])))


##Text
ggCatecholamine.text <- ggplot() + 
  geom_tile(data=plot.data.catecholamine, mapping=aes(x=A5_ID, y="Catecholamines"), fill=ColorPalette["LightGrey2"], width=0.8,height=0.8) +
  geom_text(data=plot.data.catecholamine %>% filter(`Catecholamine_profile` != "No data"), mapping=aes(x=A5_ID, y="Catecholamines", label=code), size=3)  +
  scale_x_discrete(drop=F) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Catecholamines") 

#############
#############
## GENOMIC ##
#############
#############

#######
# TMB #
#######

plot.data.tmb <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, wgs_tmb)
plot.data.tmb$wgs_tmb <- as.numeric(plot.data.tmb$wgs_tmb)

##Free Axis
# ggTMB.free <- ggplot() + 
#   geom_col(data=plot.data.tmb, mapping=aes(x=A5_ID, y=wgs_tmb)) +
#   scale_x_discrete(drop=F) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + ylab("TMB")

##Log 10
# ggTMB.log10 <- ggplot() + 
#   geom_col(data=plot.data.tmb, mapping=aes(x=A5_ID, y=wgs_tmb)) + scale_y_log10() +
#   scale_x_discrete(drop=F) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + ylab("TMB") #, panel.grid = element_blank()

##Cut Axis
ggTMB.cutoff <- ggplot() +
  geom_col(data=plot.data.tmb, mapping=aes(x=A5_ID, y=wgs_tmb)) +
  coord_cartesian(ylim=c(0,4.5)) +
  geom_text_repel(data=plot.data.tmb %>% filter(wgs_tmb>5),
                  aes(x=A5_ID, y=wgs_tmb,label=wgs_tmb), color=ColorPalette["LightOrange1"]) +
  scale_x_discrete(drop=F) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + ylab("TMB")

##Heatmap
# ggTMB.heatmap <- ggplot() + 
#   geom_tile(data=plot.data.tmb, 
#             mapping=aes(x=A5_ID, 
#                         y="TMB", 
#                         fill=wgs_tmb),
#             width=0.8,
#             height=0.8) +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradient(low="white", high="blue", limits=c(0,5), na.value="red") +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1),
#         axis.title.y = element_blank(),
#         panel.grid = element_blank())

############
# SV count #
############

plot.data.sv <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, structural_variant_count)
plot.data.sv$structural_variant_count <- as.numeric(plot.data.sv$`structural_variant_count`)
plot.data.sv <- plot.data.sv %>% ungroup() %>% 
  mutate(structural_variant_count_z=
           ((structural_variant_count-mean(structural_variant_count))/
              sd(structural_variant_count)))

##Free Axis
# ggSV <- ggplot() + 
#   geom_col(data=plot.data.sv, mapping=aes(x=A5_ID, y=structural_variant_count)) +
#   scale_x_discrete(drop=F) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + ylab("SV count")
# 
# ggSV.cutoff <- ggplot() + 
#   geom_col(data=plot.data.sv, mapping=aes(x=A5_ID, y=structural_variant_count)) +
#   coord_cartesian(ylim=c(0,500)) +
#   geom_text_repel(data=plot.data.sv %>% filter(`structural_variant_count`>500), 
#                   aes(x=A5_ID, y=`structural_variant_count`,label=`structural_variant_count`), color=ColorPalette["LightOrange1"]) +
#   scale_x_discrete(drop=F) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + ylab("nSV")

##Heatmap
ggSV.heatmap <- ggplot() + 
  geom_tile(data=plot.data.sv, 
            mapping=aes(x=A5_ID, 
                        y="nSV", 
                        fill=structural_variant_count),
            width=0.8,
            height=0.8) +
  scale_x_discrete(drop=F) +
  scale_fill_gradientn(colours = c("white", "blue", "green", "orange", "red"), values = c(0,0.05,0.1,0.5,1)) +
  #scale_fill_gradient(low="white", high="blue", limits=c(0,100), na.value="red") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1),
        axis.title.y = element_blank(),
        panel.grid = element_blank())

# ggSV.heatmap_z <- ggplot() + 
#   geom_tile(data=plot.data.sv, mapping=aes(x=A5_ID, y="nSV", fill=structural_variant_count_z)) +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradient(low="white", high="blue") +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1))

###################
# Telomere Length #
###################


# plot.data.telomere <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, telhunter_log2_telcontentratio) %>% 
#   mutate(telhunter_log2_telcontentratio= as.numeric(telhunter_log2_telcontentratio),
#          Elongated_Telomeres=ifelse(telhunter_log2_telcontentratio>0.5,"log2(T/C) > 0.5", "log2(T/C) < 0.5"),
#          Elongated_Telomeres=factor(Elongated_Telomeres, levels=c("log2(T/C) < 0.5", "log2(T/C) > 0.5")))
# 
# #Heat Map
# ggTelomere.heat <- ggplot() + 
#   geom_tile(data=plot.data.telomere, mapping=aes(x=A5_ID, y="Telomere Ratio", fill=telhunter_log2_telcontentratio), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="log2 Telomere Ratio") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradient2(low = ColorPalette["DarkGreen1"], mid = "white",high = ColorPalette["Purple1"], midpoint = 0) 
# 
# #col graph
# ggTelomere.col <- ggplot() + 
#   geom_col(data=plot.data.telomere, mapping=aes(x=A5_ID, y=telhunter_log2_telcontentratio)) +
#   scale_x_discrete(drop=F) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + ylab("log2 Telomere Ratio")
# 
# #Binary
# ggTelomere.monotone <- ggplot() + 
#   geom_tile(data=plot.data.telomere %>% filter(Elongated_Telomeres=="log2(T/C) > 0.5"), mapping=aes(x=A5_ID, y="Telomere Ratio", fill=Elongated_Telomeres), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="log2 Telomere Ratio") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values= c(as.character(ColorPalette["DarkGrey1"]))) 


#######
# WGD #
#######


plot.data.wgd <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, `sample_ploidy`) %>% mutate(WGD=ifelse(as.numeric(`sample_ploidy`) > 3.2,"Yes","No"))
plot.data.wgd$WGD <- factor(as.character(plot.data.wgd$WGD), levels=c("Yes","No"))

#Heat Map
# ggWGD <- ggplot() + 
#   geom_tile(data=plot.data.wgd, mapping=aes(x=A5_ID, y="WGD", fill=WGD), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Whole Genome Doubling") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values = c("Black", ColorPalette[["LightGrey1"]]))

ggWGD.monotone <- ggplot() + 
  geom_tile(data=plot.data.wgd, mapping=aes(x=A5_ID, y="WGD", fill=WGD), width=0.8,height=0.8) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Whole Genome Doubling") +
  scale_x_discrete(drop=F) +
  scale_fill_manual(values = c("Yes"="black", "No"=ColorPalette[["LightGrey1"]]))


########################
# CHROMOTHRIPSIS EVENT #
########################
plot.data.chromothripsis <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, `chromothriptic_event`) %>% 
  mutate(`chromothriptic_event`=ifelse(`chromothriptic_event`=="NA", NA_character_,`chromothriptic_event`),
         `chromothriptic_event`=gsub("chr|p|q","",`chromothriptic_event`),
         `chromothriptic_event`=gsub("[(}]Exchange[)]","",`chromothriptic_event`),
         `chromothriptic_event`=gsub("[&]","\n",`chromothriptic_event`),
         dummy=T,
         chromothriptic_event_binary=ifelse(is.na(chromothriptic_event), "No","Yes"))  


ggChromothripsis <- ggplot(data=plot.data.chromothripsis, mapping=aes(x=A5_ID, y="Chromothripsis", fill=dummy, label=`chromothriptic_event`)) + 
  geom_tile(width=0.8, height=0.8, color=as.character(ColorPalette[["LightGrey1"]])) +
  geom_text( size=1.5) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + 
  #labs(fill="Chromothripsis") +
  scale_x_discrete(drop=F) +  
  scale_fill_manual(values = c(ColorPalette[["LightGrey1"]])) + guides(fill="none")

ggChromothripsis.binary <- ggplot(data=plot.data.chromothripsis, mapping=aes(x=A5_ID, y="Chromothripsis", fill=chromothriptic_event_binary)) + 
  geom_tile(width=0.8, height=0.8, color=as.character(ColorPalette[["LightGrey1"]])) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + 
  #labs(fill="Chromothripsis") +
  scale_x_discrete(drop=F) +  
  scale_fill_manual(values = c("Yes"="black","No"=ColorPalette[["LightGrey1"]])) + guides(fill="none")


############
# FEATURES #
############
# offline = F 
# if(offline) { 
#   features <- read.delim("/g/data/pq08/projects/ppgl/a5/offline_cache/sample_features.tsv")
# } else {
#   feature_sheet <- as_sheets_id("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY")
#   gs4_auth("aidan.flynn@umccr-radio-lab.page")
#   features <- read_sheet(ss = feature_sheet,sheet = "SampleFeatures", col_types = "c")
# }

keep_genes_curated <- a5_somatic_variants_keep_curated %>% 
  ungroup() %>%  
  filter(A5_ID != "E167-1") %>% 
  left_join(a5_anno %>% dplyr::select(A5_ID,`Patient ID`)) %>% 
  dplyr::select(`Patient ID`, PCGR_SYMBOL) %>% 
  distinct() %>% 
  group_by(PCGR_SYMBOL) %>% 
  dplyr::count() %>% 
  filter(n>1) %>% 
  arrange(desc(n)) %>%
  pull(PCGR_SYMBOL)

keep_genes_pcawg <- a5_somatic_variants_keep %>% 
  ungroup() %>%
  filter(A5_ID != "E167-1") %>% 
  left_join(a5_anno %>% dplyr::select(A5_ID,`Patient ID`)) %>% 
  dplyr::select(`Patient ID`, PCGR_SYMBOL) %>% 
  distinct() %>% 
  group_by(PCGR_SYMBOL) %>% 
  dplyr::count() %>% 
  filter(n>1) %>% 
  arrange(desc(n)) %>%  
  filter(PCGR_SYMBOL %in% pcawg_drivers$Gene) %>% 
  pull(PCGR_SYMBOL) 


tcga_seen_gteq_3 <- c("HRAS", "NF1", "TTN", "EPAS1", "ATRX", "MUC16", "RET", "SCRIB", "ABCA13", "MUC17", 
                    "FLG", "SPHKAP", "ADGRG7", "SETD2", "HSPG2", "USP9X", "CSDE1", "PARG", "AHNAK", 
                    "HYDIN", "RYR1", "KCNH5", "KMT2C", "MLLT6", "UNC79", "TGDS", "RP1L1", "MUC5B", 
                    "PCDHB4", "LYST", "PRUNE2", "WDR11", "FKBP9", "COG7", "TENM1", "TRHDE", "BPTF", 
                    "BIRC6", "ANK3", "BRCA2", "HUWE1", "LAMA3", "VHL", "VIT", "TRRAP", "COL12A1", 
                    "ABCA12", "ALMS1", "ZFC3H1", "PDE4DIP")

keep_genes <- unique(c(keep_genes_curated,keep_genes_pcawg,tcga_seen_gt_3))

features_expr <- wts_log_cpm %>%  
  ungroup() %>% 
  filter(symbol %in% keep_genes, abs(log2cpm_z) > 2) %>% 
  mutate(Event=ifelse(log2cpm_z > 2, "High expression", "Low expression")) %>% 
  dplyr::select(Gene=symbol, A5_ID, Event)

features_sv <- A5_gridss_keep %>% mutate(A5_ID=gsub("T0", "", A5_ID)) %>% 
  filter(GeneStartDisrupted | GeneEndDisrupted, 
         (GeneStartName %in% keep_genes | GeneEndName %in% keep_genes))

features_sv <- bind_rows(
  features_sv %>% dplyr::select(A5_ID, Gene=GeneStartName, disrupted=GeneStartDisrupted),
  features_sv %>% dplyr::select(A5_ID, Gene=GeneEndName, disrupted=GeneEndDisrupted)) %>% 
  mutate(Event="Structural variant") %>% 
  left_join(features_expr %>%  dplyr::rename("Expr_Outlier"=Event)) %>% 
  filter(as.logical(disrupted) | !is.na(Expr_Outlier), Gene %in% keep_genes) %>%  distinct() %>% 
  dplyr::select(-Expr_Outlier, -disrupted)


features_mutation <- a5_somatic_variants_keep %>%  ungroup() %>%  filter(PCGR_SYMBOL %in% keep_genes) %>% dplyr::select(A5_ID, PCGR_SYMBOL, PCGR_CONSEQUENCE) %>% 
  mutate(PCGR_CONSEQUENCE=ifelse(PCGR_SYMBOL == "TERT" & PCGR_CONSEQUENCE== "upstream_gene_variant", "Promoter mutation", PCGR_CONSEQUENCE )) %>% 
  mutate(PCGR_CONSEQUENCE=dplyr::recode(PCGR_CONSEQUENCE, 
                                        "missense_variant"="Missense", 
                                        "frameshift_variant"="Frameshift", 
                                        "stop_gained"="Stop gained",
                                        "splice_acceptor_variant"="Splice acceptor",
                                        "splice_donor_variant"="Splice donor",
                                        "splice_region_variant|_intron_variant"="Splice region",
                                        "missense_variant|_splice_region_variant"= "Missense",
                                        "inframe_deletion"="In-frame deletion"
                                        )) %>% 
  dplyr::rename("Event"="PCGR_CONSEQUENCE", "Gene"="PCGR_SYMBOL") 

features_cn <- A5_gene_cn_lossgain %>%
  dplyr::select(A5_ID, gene, loss_gain) %>% 
  dplyr::rename("Gene"="gene", Event=loss_gain) %>% 
  inner_join(features_mutation %>%  dplyr::select(A5_ID,Gene)) 

features <- bind_rows(
  features_mutation %>%  mutate(Source="Mutation"),
  features_expr %>% mutate(Source="Expression"),
  features_sv %>% mutate(Source="SV"),
  features_cn %>% mutate(Source="CN")
  ) %>% 
  distinct()

features <- features %>%  filter(A5_ID %in% SampleOrder.A5_ID, Gene %in% c(features_mutation$Gene[features_mutation$A5_ID != "E167-1"], features_sv$Gene))
#features$Gene <- factor(as.character(features$Gene), rev(keep_genes))
features$Gene <- factor(as.character(features$Gene), features %>% filter(Source %in% c("Mutation", "SV")) %>% group_by(Gene) %>% dplyr::count() %>% arrange(n) %>% pull(Gene))
features$Event <- factor(as.character(features$Event), levels = c(names(genomic_alteration_cols), "High expression", "Low expression", "Gain", "CNLOH", "Loss"))
features$Source <- factor(as.character(features$Source), levels = c("Mutation","SV","Expression", "CN"))
features$A5_ID <- factor(as.character(features$A5_ID), levels = SampleOrder.A5_ID)
features <- features %>% filter(!is.na(Gene)) # Remove genes without mutation


ggFeatures <- ggplot() + 
  geom_tile(data=features %>% filter(Source %in% c("Mutation","SV","Expression")) %>% group_by(A5_ID,Gene) %>% slice_min(n= 1, order_by = Event), mapping=aes(x=A5_ID, y=Gene, fill=Event), width=0.8,height=0.8) + #, alpha=Clonal
  geom_point(data=features %>% filter(Source %in% c("CN")), 
             mapping=aes(x=A5_ID, y=Gene, shape=Event), color="black") + #, alpha=Clonal
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_line(colour = "#f7f7f7ff")) + 
  scale_fill_manual(values = c(genomic_alteration_cols, "Low expression"="#0000ff", "High expression"="#ff0000")) +
  #scale_alpha_manual(values = c(1,0.4)) +
  scale_shape_manual(values=c("Gain"=24,"Loss"=25,"CNLOH"=23)) +
  scale_x_discrete(drop=F) +
  labs(fill="Mutation Type") 

# ggFeatures.TERTATRX <- ggplot() + 
#   geom_tile(data=features %>% filter(Gene %in% c("TERT","ATRX")), mapping=aes(x=A5_ID, y=Gene, fill=Event), width=0.8,height=0.8) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_line(colour = "#f7f7f7ff")) + 
#   scale_fill_manual(values = genomic_alteration_cols) +
#   scale_x_discrete(drop=F) +
#   labs(fill="Mutation Type") 

################
################
## EXPRESSION ##
################
################

#######################
# Catecholamine Genes #
#######################

# MarkerGenes <- c("SLC18A1","ARG2","CHGB","CHGA","TRIB3","DBH","DDC","LEF1","CARTPT","SLC18A2","TH","PENK","ALK","ASCL1","ATF3","CDH9","ERBB4","FOS","GATA3","HAND2","ISL1","JUN","JUNB","NPY","NTRK3","PHOX2A","PHOX2B","PNMT","RET","SLC6A2","SOX11","TFAP2B")
# TERT_OverExpr_Genes <- c("FAM83D","ESPL1","UBE2C","KIF18B","DLGAP5","TOP2A","CDCA2","HJURP","AURKB","BIRC5","MKI67","CENPF","NEK2")
# 
# plot.data.markerexp <- wts_log_cpm %>% filter(symbol %in% MarkerGenes)
# plot.data.markerexp <- crossing(A5_ID=SampleOrder.A5_ID, symbol=MarkerGenes) %>%  left_join(plot.data.markerexp)  
# plot.data.markerexp <- plot.data.markerexp %>%  mutate(zExpClamped=case_when(
#   log2cpm_z > 2.5 ~ 2.5, 
#   log2cpm_z < -2.5 ~ -2.5, 
#   TRUE ~ log2cpm_z
# ),
# zExpOutOfRange=case_when(
#   log2cpm_z > 2.5 ~ TRUE, 
#   log2cpm_z < -2.5 ~ TRUE, 
#   TRUE ~ NA
# ))
# plot.data.markerexp$A5_ID <- factor(as.character(plot.data.markerexp$A5_ID), levels = SampleOrder.A5_ID)
# 
# g <- ggplot() + 
#   geom_tile(data=plot.data.markerexp, mapping=aes(x=A5_ID, y=symbol, fill=log2cpm_z), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="log2-CPM (z-scaled)") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradient2(low = ColorPalette["DarkBlue1"], mid = "white",high = ColorPalette["DarkRed1"], midpoint = 0, na.value = "grey") 
# 
# ggMarker.exp.clamped <- ggplot() + 
#   geom_tile(data=plot.data.markerexp, mapping=aes(x=A5_ID, y=symbol, fill=zExpClamped), width=0.8,height=0.8) +
#   geom_point(data=plot.data.markerexp %>%  filter(zExpOutOfRange), mapping=aes(x=A5_ID, y=symbol), size=0.4) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="log2-CPM (z-scaled)") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradient2(low = ColorPalette["DarkBlue1"], mid = "white",high = ColorPalette["DarkRed1"], midpoint = 0, na.value = "grey") 


#################
# TERT/ATRX Exp #
#################

#ADD  

# plot.data.taexp <- wts_log_cpm %>% filter(symbol %in% c("TERT", "ATRX"))
# plot.data.taexp$A5_ID <- factor(as.character(plot.data.taexp$A5_ID), levels = SampleOrder.A5_ID)
# 
# 
# ggTERTATRXExp <- ggplot() + 
#   geom_tile(data=plot.data.taexp, mapping=aes(x=A5_ID, y=symbol, fill=log2cpm_z), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="log2-CPM (z-scaled)") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradient2(low = ColorPalette["DarkBlue1"], mid = "white",high = ColorPalette["DarkRed1"], midpoint = 0) 
# 
# 
# ## TERT Copy Number
# plot.data.tertcn <- A5_gene_cn %>% mutate(A5_ID=gsub("-T0","-",A5_ID)) %>%  filter(gene=="TERT") %>%  
#   dplyr::select(A5_ID, maxCopyNumber) %>% 
#   dplyr::rename(TERT_max_CN=maxCopyNumber) %>% mutate(TERT_max_CN=factor(as.character(round(TERT_max_CN,0)), levels=as.character(0:10))) %>% 
#   filter(A5_ID %in% SampleOrder.A5_ID)
# 
# plot.data.tertcn$A5_ID <- factor(as.character(plot.data.tertcn$A5_ID), levels = SampleOrder.A5_ID)
# 
# ggTERTCN <- ggplot() + 
#   geom_tile(data=plot.data.tertcn, mapping=aes(x=A5_ID, y="TERT CN", fill=TERT_max_CN), width=0.8,height=0.8) +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values=as.character(c(ColorPalette["LightBlue1"], ColorPalette["DarkBlue1"], 
#                                           ColorPalette["Purple1"],ColorPalette["Purple3"],
#                                           ColorPalette["DarkRed1"])),
#                     guide = guide_legend(direction = "horizontal", title.position = "right",
#                                          label.position="bottom", label.hjust = 0.5, label.vjust = 0.5)) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + 
#   labs(fill="TERT Max. CN") 



########
# KI67 #
########

####
# Staining
####

plot.data.ki67stain <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, Ki67_Staining) %>% 
  mutate(Ki67_Staining_numeric = as.numeric(gsub("<","",Ki67_Staining)), 
         Ki67_Staining_categorical = case_when(
           Ki67_Staining == "<0.5" ~ "<0.5",
           Ki67_Staining %in% c("0.5","<1","1.5", as.character(1:2)) ~ "0.5-2",
           Ki67_Staining %in% as.character(3:5) ~ "3-5",
           Ki67_Staining %in% as.character(6:8) ~ "6-8",
           Ki67_Staining %in% as.character(8:10) ~ "8-10",
           Ki67_Staining %in% as.character(10:100) ~ "10+",
           TRUE ~ "No data"),
         Ki67_Staining_categorical = factor(Ki67_Staining_categorical, levels=c("<0.5", "0.5-2",  "3-5", "6-8", "8-10", "10+",  "No Data")),
         Ki67_Staining_binary=case_when(
           Ki67_Staining_numeric>=3 ~ ">3%",
           Ki67_Staining_numeric<3 ~ "<3%",
           TRUE ~ "No data"))


#numeric
# ggKi67stain.numeric <- ggplot() + 
#   geom_tile(data=plot.data.ki67stain, mapping=aes(x=A5_ID, y="Ki67 % (IHC)", fill=Ki67_Staining_numeric), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Ki67 IHC(%)") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradientn(colours = c("LightBlue1", ColorPalette["LightBlue2"], ColorPalette["LightOrange1"], ColorPalette["DarkRed1"]), values=c(0,0.1,0.3,0.6,1)) 
# 
# #categorical
# ggKi67stain.category <- ggplot() + 
#   geom_tile(data=plot.data.ki67stain, mapping=aes(x=A5_ID, y="Ki67 % (IHC)", fill=Ki67_Staining_categorical), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Ki67 IHC(%)") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values = c(ColorPalette[["LightBlue1"]], ColorPalette[["DarkBlue1"]], ColorPalette[["Yellow1"]], 
#                                ColorPalette[["LightOrange1"]], ColorPalette[["DarkRed1"]], ColorPalette[["LightGrey1"]])) 

#binary
ggKi67stain.binary <- ggplot() + 
  geom_tile(data=plot.data.ki67stain, mapping=aes(x=A5_ID, y="Ki67 (IHC)", fill=Ki67_Staining_binary), width=0.8,height=0.8) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Ki67 IHC(%)") +
  scale_x_discrete(drop=F) +
  scale_fill_manual(values = c("<3%"=ColorPalette[["DarkGrey2"]], ">3%"="Black", "No data"=ColorPalette[["LightGrey1"]])) 


##Ki67 Expression
plot.data.ki67expr <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, MKI67_log2_cpm) %>%  mutate(MKI67_log2_cpm=as.numeric(MKI67_log2_cpm)) 

ggKi67.expr <- ggplot() + 
  geom_tile(data=plot.data.ki67expr, mapping=aes(x=A5_ID, y="MKi67 (log2-CPM)", fill=MKI67_log2_cpm), width=0.8,height=0.8) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="MKi67 Expr (log2-CPM)") +
  scale_x_discrete(drop=F) +
  scale_fill_gradientn(colours = c("white", ColorPalette["LightGreen1"], ColorPalette["LightOrange1"], ColorPalette["DarkRed1"]), values=c(0,0.4,0.6,0.9,1)) 


############
# SDHB IHC #
############

# plot.data.sdhbihc <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, SDHB_Staining) %>% 
#   dplyr::mutate(SDHB_Staining=case_when(
#     SDHB_Staining=="neg" ~ "Negative",
#     SDHB_Staining=="pos" ~ "Positive",
#     TRUE ~ "No Data"))
# 
# plot.data.sdhbihc$SDHB_Staining <- factor(as.character(plot.data.sdhbihc$SDHB_Staining), levels=c("Negative","Positive", "No Data"))
# 
# ggsdhbihc <- ggplot() + 
#   geom_tile(data=plot.data.sdhbihc, mapping=aes(x=A5_ID, y="SDHB IHC", fill=SDHB_Staining), width=0.8,height=0.8) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + 
#   labs(fill="SDHB IHC") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values = c(Negative=ColorPalette[["LightBrown1"]],  Positive=ColorPalette[["LightGreen1"]], "No Data"=ColorPalette[["LightGrey1"]]))
# 
# ggsdhbihc.monotone <- ggplot() + 
#   geom_tile(data=plot.data.sdhbihc, mapping=aes(x=A5_ID, y="SDHB IHC", fill=SDHB_Staining), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="SDHB IHC") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values = c(Negative=ColorPalette[["DarkGrey1"]],  Positive="black", "No Data"=ColorPalette[["LightGrey1"]]))
# 

###################
# Chr1p36.13_Loss #
###################
 

# plot.data.chr1ploss <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, Chr1p36.13_Loss) 
# 
# plot.data.chr1ploss$Chr1p36.13_Loss <- factor(as.character(plot.data.chr1ploss$Chr1p36.13_Loss), levels=c("Loss","CNLOH","Diploid"))
# 
# ggchr1ploss <- ggplot() + 
#   geom_tile(data=plot.data.chr1ploss, mapping=aes(x=A5_ID, y="Chr1p36.13 Loss", fill=Chr1p36.13_Loss), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Chr1p36.13 Loss") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values = c("Loss"=ColorPalette[["LightRed2"]],  
#                                "CNLOH"=ColorPalette[["Purple3"]], 
#                                "Diploid"=ColorPalette[["LightGreen1"]]))
# 
# ggchr1ploss.monotone <- ggplot() + 
#   geom_tile(data=plot.data.chr1ploss, mapping=aes(x=A5_ID, y="Chr1p36.13 Loss", fill=Chr1p36.13_Loss), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Chr1p36.13 Loss") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values = c("Loss"="black",  
#                                "CNLOH"=ColorPalette[["DarkGrey1"]], 
#                                "Diploid"=ColorPalette[["LightGreen1"]]))

#################
# Tumour Purity #
#################

# plot.data.purity <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, `sample_purity`) %>% 
#   dplyr::mutate(`sample_purity`=as.numeric(`sample_purity`)) %>% dplyr::rename(Tumour_Purity=`sample_purity`)
# 
# ggPurity.heat <-  ggplot() + 
#   geom_tile(data=plot.data.purity, mapping=aes(x=A5_ID, y="Tumour Purity", fill=Tumour_Purity), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Tumour_Purity") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradient2(low = ColorPalette[["DarkRed1"]], mid = ColorPalette[["LightBlue1"]], high = ColorPalette[["DarkBlue1"]], midpoint = 0.4)
# 
# ggPurity.col <-  ggplot() + 
#   geom_col(data=plot.data.purity, mapping=aes(x=A5_ID, y=Tumour_Purity), width=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Tumour_Purity") +
#   scale_x_discrete(drop=F) + coord_cartesian(ylim=c(0,1)) + scale_y_continuous(labels = c("0","","","","1"))

###################
# STAR Expression #
###################

# plot.data.star <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, STAR_log2cpm) %>% 
#   dplyr::mutate(STAR_log2cpm=as.numeric(STAR_log2cpm)) 
# 
# ggStar <-  ggplot() + 
#   geom_tile(data=plot.data.star, mapping=aes(x=A5_ID, y="STAR log2-CPM", fill=STAR_log2cpm), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="STAR (log2-CPM)") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradient2(low = "white", mid = "white", high = ColorPalette[["DarkRed1"]], midpoint = 1)

##################
# SDHB Genotypes #
##################


# plot.data.genotype <- a5_anno.use %>% dplyr::select(A5_ID, `Patient ID`, `germline_SDHB_mutation`, `germline_SDHB_protein_change`) %>% 
#   dplyr::mutate(`Germline_SDHB_mutation`=case_when(
#     `germline_SDHB_mutation` %in% c("c.418G>T","c.136C>T","c.268C>T","c.380T>G","c.137G>A") ~`germline_SDHB_mutation`,
#     grepl("c.[0-9]+[ATCG]>[ATCG]", `germline_SDHB_mutation`) ~ "Other Misense",
#     grepl("c.[0-9]+[+][0-9]",`germline_SDHB_mutation`) ~ "Splice",
#     grepl("[Dd]el",`germline_SDHB_mutation`)  ~ "Deletion"
#   )) %>% 
# dplyr::mutate(`protein change`=case_when(
#   `germline_SDHB_protein_change` %in% c("p.Val140Phe","p.Arg46*","p.Arg90*","p.Arg46Gln","p.Ile127Ser") ~`germline_SDHB_protein_change`,
#   grepl("p.[A-Za-z]{3}[0-9]+[A-Za-z]{3}", `germline_SDHB_protein_change`) ~ "Other misense",
#   grepl("p.[A-Za-z]{3}[0-9]+(fs)?[*]", `germline_SDHB_protein_change`) ~ "Other frameshift/stop",
#   grepl("Splice",`germline_SDHB_protein_change`) ~ "Splice",
#   grepl("[Dd]el",`germline_SDHB_protein_change`)  ~ "Deletion"
# ))
# 
# plot.data.genotype$`protein change` <- factor(as.character(plot.data.genotype$`protein change`), levels= c("p.Val140Phe","p.Arg46Gln","p.Arg46*","p.Arg90*","p.Ile127Ser","Other misense","Deletion","Splice","Other frameshift/stop"))
# 
# ggGenotype <- ggplot() + 
#   geom_tile(data=plot.data.genotype, mapping=aes(x=A5_ID, y="Germline SDHB", fill=`protein change`), width=0.8,height=0.8) +
#   theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) + labs(fill="Germline SDHB") +
#   scale_x_discrete(drop=F) +
#   scale_fill_manual(values = c(ColorPalette[["LightBlue1"]], ColorPalette[["DarkBlue1"]], ColorPalette[["DarkGreen1"]], 
#                                ColorPalette[["LightOrange1"]], ColorPalette[["LightRed1"]], ColorPalette[["LightGreen2"]],
#                                ColorPalette[["DarkBrown2"]], ColorPalette[["Purple1"]], ColorPalette[["Yellow1"]])) 

##############
# Signatures #
##############

#The analysis relies on data previously generated by:
# "/g/data/pq08/projects/ppgl/a5/wgs/analysis/mutational_patterns/A5_mutational_patterns.r"
load("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mutational_patterns/A5_mutational_patterns.rworkspace")


sig_plot_prep = . %>% tibble::rownames_to_column("Signature") %>%
  tidyr::pivot_longer(cols = -Signature,
                      names_to = "A5_ID",
                      values_to = "Contribution") %>%
  mutate(`Patient ID`=gsub("-.","",A5_ID)) %>%
  mutate(Signature=factor(Signature,
                          levels = colnames(current_signature_set)))  %>%
  group_by(A5_ID) %>%
  mutate(Total=sum(Contribution),
         Prop_Contribution=Contribution/Total) %>%
  dplyr::select(-Total) %>%
  mutate(A5_ID=gsub("[.]","-", A5_ID)) %>%
  group_by(Signature) %>%
  mutate(Contribution_z = (Contribution - mean(Contribution))/sd(Contribution),
         Prop_Contribution_z = (Prop_Contribution - mean(Prop_Contribution))/sd(Prop_Contribution)) %>%
  filter(A5_ID %in% SampleOrder.A5_ID)

# find_scaled_zero <- function(plot.data.sigs, feature_column) {
#   scaled_zero_value = tibble(raw=plot.data.sigs[[feature_column]],
#                              scaled=scales::rescale(x = plot.data.sigs[[feature_column]], to = c(0, 1))) %>%
#     arrange(abs(raw)) %>%
#     slice_head(n=1) %>%
#     pull(scaled)
#   return(scaled_zero_value)
# }

#######
# SBS #
#######

current_signature_set <- signatures

a5_sbs <- 
  data.frame(fit_res$contribution, 
             check.names = F) %>% 
  sig_plot_prep

topsigs_sbs <- a5_sbs %>% group_by(Signature) %>% 
  summarise(max_contrib=max(Contribution), max_prop=max(Prop_Contribution)) %>% 
  filter(max_prop > 0.15, max_contrib > 500) %>% 
  pull(Signature)
sig_order_sbs <- a5_sbs %>% group_by(Signature) %>% summarise(mean_prop=mean(Prop_Contribution)) %>% arrange(mean_prop) %>% pull(Signature) 



plot.data.sigs_sbs <- a5_sbs %>% filter(Signature %in% topsigs_sbs)
plot.data.sigs_sbs$Signature <- factor(as.character(plot.data.sigs_sbs$Signature), levels=sig_order_sbs)

plot.data.sigs_sbs$A5_ID <- factor(as.character(plot.data.sigs_sbs$A5_ID), levels = SampleOrder.A5_ID)

# ggSignatures_sbs <- ggplot() + 
#   geom_tile(data=plot.data.sigs_sbs, mapping=aes(x=A5_ID, y=Signature, fill=Prop_Contribution), width=0.8,height=0.8) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
#         panel.grid = element_blank()) + 
#   labs(fill="Signature Proportion") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradient(low = "white", high = ColorPalette[["Purple2"]])
# 
# 
# scaled_zero_value = tibble(raw=plot.data.sigs_sbs$Contribution_z,
#                            scaled=scales::rescale(x = plot.data.sigs_sbs$Contribution_z, to = c(0, 1))) %>% 
#   arrange(abs(raw)) %>% 
#   slice_head(n=1) %>% 
#   pull(scaled)


ggSignatures_sbs_contrib <- ggplot() + 
  geom_tile(data=plot.data.sigs_sbs, mapping=aes(x=A5_ID, y=Signature, fill=Contribution), width=0.8,height=0.8) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
        panel.grid = element_blank()) + 
  labs(fill="Signature Proportion") +
  scale_x_discrete(drop=F) +
  scale_fill_gradientn(values= c(0, 0.05, 0.1, 0.2, 0.5, 1), colors = c("white",ColorPalette[["DarkBlue1"]],ColorPalette[["LightGreen1"]],ColorPalette[["LightOrange1"]]), na.value = "red", limits = c(0,2500)) 

# ggSignatures_sbs_contrib_z <- ggplot() + 
#   geom_tile(data=plot.data.sigs_sbs, mapping=aes(x=A5_ID, y=Signature, fill=Contribution_z), width=0.8,height=0.8) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
#         panel.grid = element_blank()) + 
#   labs(fill="Signature Proportion") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradientn(values = c(0, find_scaled_zero(plot.data.sigs_sbs, "Contribution_z"), 0.5,1), 
#                        colors = c(ColorPalette[["Purple2"]],
#                                   "white",
#                                   ColorPalette[["LightRed2"]],
#                                   "Red"))
# 
# ggSignatures_sbs_contrib_prop_z <- ggplot() + 
#   geom_tile(data=plot.data.sigs_sbs, mapping=aes(x=A5_ID, y=Signature, fill=Prop_Contribution_z), width=0.8,height=0.8) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
#         panel.grid = element_blank()) + 
#   labs(fill="Signature Proportion") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradientn(values = c(0, find_scaled_zero(plot.data.sigs_sbs, "Prop_Contribution_z"), 0.5,1), 
#                        colors = c(ColorPalette[["Purple2"]],
#                                   "white",
#                                   ColorPalette[["LightRed2"]],
#                                   "Red"))


#######
# DBS #
#######

current_signature_set <- signatures_dbs

a5_dbs <- 
  data.frame(fit_res_dbs$contribution, 
             check.names = F) %>% 
  sig_plot_prep

topsigs_dbs <- a5_dbs %>% group_by(Signature) %>% 
  summarise(max_contrib=max(Contribution), max_prop=max(Prop_Contribution)) %>% 
  filter(max_prop > 0.15, max_contrib > 10) %>% 
  pull(Signature)
sig_order_dbs <- a5_dbs %>% group_by(Signature) %>% summarise(mean_prop=mean(Prop_Contribution)) %>% arrange(mean_prop) %>% pull(Signature) 

plot.data.sigs_dbs <- a5_dbs %>% filter(Signature %in% topsigs_dbs)
plot.data.sigs_dbs$Signature <- factor(as.character(plot.data.sigs_dbs$Signature), levels=sig_order_dbs)

plot.data.sigs_dbs$A5_ID <- factor(as.character(plot.data.sigs_dbs$A5_ID), levels = SampleOrder.A5_ID)


ggSignatures_dbs_contrib <- ggplot() +
  geom_tile(data=plot.data.sigs_dbs, mapping=aes(x=A5_ID, y=Signature, fill=Contribution), width=0.8,height=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1),
        panel.grid = element_blank()) +
  labs(fill="Signature Proportion") +
  scale_x_discrete(drop=F) +
  scale_fill_gradientn(values= c(0, 0.05, 0.1, 0.2, 0.5, 1), colors = c("white",ColorPalette[["DarkBlue1"]],ColorPalette[["LightGreen1"]],ColorPalette[["LightOrange1"]])) 

# ggSignatures_dbs_contrib_z <- ggplot() + 
#   geom_tile(data=plot.data.sigs_dbs, mapping=aes(x=A5_ID, y=Signature, fill=Contribution_z), width=0.8,height=0.8) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
#         panel.grid = element_blank()) + 
#   labs(fill="Signature Proportion") +
#   scale_x_discrete(drop=F) + 
#   scale_fill_gradientn(values = c(0, find_scaled_zero(plot.data.sigs_dbs, "Contribution_z"), 0.5,1), 
#                        colors = c(ColorPalette[["Purple2"]],
#                                   "white",
#                                   ColorPalette[["LightRed2"]],
#                                   "Red"))
# 
# ggSignatures_dbs_contrib_prop_z <- ggplot() + 
#   geom_tile(data=plot.data.sigs_dbs, mapping=aes(x=A5_ID, y=Signature, fill=Prop_Contribution_z), width=0.8,height=0.8) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
#         panel.grid = element_blank()) + 
#   labs(fill="Signature Proportion") +
#   scale_x_discrete(drop=F) + 
#   scale_fill_gradientn(values = c(0, find_scaled_zero(plot.data.sigs_dbs, "Prop_Contribution_z"), 0.5,1), 
#                        colors = c(ColorPalette[["Purple2"]],
#                                   "white",
#                                   ColorPalette[["LightRed2"]],
#                                   "Red"))

#########
# INDEL #
#########

current_signature_set <- signatures_indel

a5_indel <- 
  data.frame(fit_res_indel$contribution, 
             check.names = F) %>% 
  sig_plot_prep

topsigs_indel <- a5_indel %>% group_by(Signature) %>% 
  summarise(max_contrib=max(Contribution), max_prop=max(Prop_Contribution)) %>% 
  filter(max_prop > 0.15, max_contrib > 50) %>% 
  pull(Signature)

sig_order_indel <- a5_indel %>% 
  group_by(Signature) %>% 
  summarise(mean_prop=mean(Prop_Contribution)) %>% 
  arrange(mean_prop) %>% 
  pull(Signature) 

plot.data.sigs_indel <- a5_indel %>% filter(Signature %in% topsigs_indel)
plot.data.sigs_indel$Signature <- factor(as.character(plot.data.sigs_indel$Signature), levels=sig_order_indel)

plot.data.sigs_indel$A5_ID <- factor(as.character(plot.data.sigs_indel$A5_ID), levels = SampleOrder.A5_ID)

ggSignatures_indel_contrib <- ggplot() +
  geom_tile(data=plot.data.sigs_indel, mapping=aes(x=A5_ID, y=Signature, fill=Contribution), width=0.8,height=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1),
        panel.grid = element_blank()) +
  labs(fill="Signature Proportion") +
  scale_x_discrete(drop=F) +
  scale_fill_gradientn(values= c(0, 0.05, 0.1, 0.2, 0.5, 1), 
                       colors = c("white",ColorPalette[["DarkBlue1"]],ColorPalette[["LightGreen1"]],
                                  ColorPalette[["LightOrange1"]])) 

# ggSignatures_indel_contrib_z <- ggplot() + 
#   geom_tile(data=plot.data.sigs_indel, mapping=aes(x=A5_ID, y=Signature, fill=Contribution_z), width=0.8,height=0.8) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
#         panel.grid = element_blank()) + 
#   labs(fill="Signature Proportion") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradientn(values = c(0, find_scaled_zero(plot.data.sigs_indel, "Contribution_z"), 0.5,1), 
#                        colors = c(ColorPalette[["Purple2"]],
#                                   "white",
#                                   ColorPalette[["LightRed2"]],
#                                   "Red"))
# 
# ggSignatures_indel_contrib_prop_z <- ggplot() + 
#   geom_tile(data=plot.data.sigs_indel, mapping=aes(x=A5_ID, y=Signature, fill=Prop_Contribution_z), width=0.8,height=0.8) +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), 
#         panel.grid = element_blank()) + 
#   labs(fill="Signature Proportion") +
#   scale_x_discrete(drop=F) +
#   scale_fill_gradientn(values = c(0, find_scaled_zero(plot.data.sigs_indel, "Prop_Contribution_z"), 0.5,1), 
#                        colors = c(ColorPalette[["Purple2"]],
#                                   "white",
#                                   ColorPalette[["LightRed2"]],
#                                   "Red"))



####
# CNA
####

ggCNVSeg <- ggplot(A5_seg_keep %>% filter(A5_ID %in% SampleOrder.A5_ID) %>% 
                            mutate(Class=forcats:::fct_recode(.f=Class,  
                                                              Loss="Loss + Subclonal CNLOH",
                                                              Gain="WGD+Gain",
                                                              Gain="Gain+LOH",
                                                              None="Minor Subclonal Loss", 
                                                              None="Diploid/Haploid-X"),
                                   A5_ID = factor(A5_ID, levels=SampleOrder.A5_ID)), 
                          aes(x=A5_ID, xend=A5_ID, y=start_offset, yend=end_offset, color=Class)) + 
  geom_segment(linewidth=2.5) +
  geom_hline(data = chr_offsets, mapping=aes(yintercept=offset), linewidth=0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1), panel.grid = element_blank()) +
  scale_y_reverse(breaks = chr_offsets$offset, labels = chr_offsets$chromosome, expand=c(0,0)) +
  scale_color_manual(values=cna_palette) + 
  ylab("Copy Number")


###########
# COMBINE #
###########

gg_oncoplot <- ggLocation + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSampleType + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggTreatment.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSize.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggCatecholamine.text + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggKi67stain.binary + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  #ggKi67stain.category + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggKi67.expr + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  # ggsdhbihc.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  # ggchr1ploss.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  #ggGenotype + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  #ggPurity.heat + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  # ggPurity.col + nox + no_margin + shrink_legend  + xlab("") + theme(axis.title.y = element_text(angle = 0)) +
  #ggStar + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggTMB.cutoff + nox + no_margin + shrink_legend + xlab("") +
  ggSV.heatmap + nox + no_margin + shrink_legend + xlab("") +
  ggChromothripsis.binary + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggWGD.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSignatures_sbs_contrib  + nox  + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSignatures_dbs_contrib  + nox  + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSignatures_indel_contrib  + nox  + no_margin + shrink_legend + ylab("") + xlab("") +
  # ggTelomere.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  # ggTERTATRXExp + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggFeatures + nox  + no_margin + shrink_legend + ylab("") + xlab("") +
  ggCNVSeg + no_margin + shrink_legend + ylab("") + 
  plot_layout(ncol = 1, heights = 
                c(0.1, #ggLocation
                  0.1, #ggDiseaseCourse
                  0.1, #ggDiseaseCourse
                  # 0.1, #ggSampleType
                  0.1, #ggSize
                  0.1, #ggCatecholamine
                  0.1, #ggKi67stain.numeric
                  #0.1, #ggKi67stain.category
                  0.1, #ggKi67.expr
                  # 0.1, #ggsdhbihc
                  # 0.1, #ggchr1ploss
                  #0.1, #ggGenotype
                  #0.1, #ggPurity.heat
                  # 0.1, #ggPurity.col
                  #0.1, #ggStar
                  0.4, #ggTMB
                  0.1, #ggSV
                  0.1, #ggChromothripsis
                  0.1, #ggWGD
                  0.75, #gg_SBS_Signatures
                  0.45, #gg_DBS_Signatures
                  0.45, #gg_ID_Signatures
                  # 0.1, #ggTelomere
                  # 0.2, #ggTERTATRX.exp
                  0.9, #ggFeatures
                  2 #ggCNA
                ), guides = "collect") & theme(legend.position = "bottom")

ggsave(filename = "/g/data/pq08/projects/ppgl/a5/tertiary/results/figures/oncoplot_clinical_genomics.pdf",
       plot = gg_oncoplot, 
       device = pdf(),
       width = 7,
       height= 6.5,
       units = "in",
       scale = 2)
dev.off()

ggMarker.exp +
  ggMarker.exp.clamped +
  plot_layout(ncol = 1, guides = "collect")




###############
# Clinical Features
###############

ggHandN + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggDiseaseCourse + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSampleType.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSize.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggCatecholamine.text + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggKi67stain.binary + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  #ggKi67stain.category + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggKi67.expr + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggsdhbihc.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggchr1ploss.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggGenotype + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  #ggPurity.heat + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggPurity.col + nox + no_margin + shrink_legend  + xlab("") + theme(axis.title.y = element_text(angle = 0)) +
  ggStar + no_margin + shrink_legend + ylab("") +
  
  plot_layout(ncol = 1, heights = 
                c(0.1, #ggHandN
                  0.1, #ggDiseaseCourse
                  0.1, #ggSampleType
                  0.1, #ggSize
                  0.1, #ggCatecholamine
                  0.1, #ggKi67stain.numeric
                  #0.1, #ggKi67stain.category
                  0.1, #ggKi67.expr
                  0.1, #ggsdhbihc
                  0.1, #ggchr1ploss
                  0.1, #ggGenotype
                  #0.1, #ggPurity.heat
                  0.1, #ggPurity.col
                  0.1 #ggStar
                  
                ), guides = "collect")


###############
# Basic
###############

ggLocation + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSampleType + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSize.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggCatecholamine.text + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggKi67stain.numeric + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  #ggKi67stain.category + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggsdhbihc.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggTMB.cutoff + nox + no_margin + shrink_legend + xlab("") +
  ggFeatures +
  plot_layout(ncol = 1, heights = 
                c(0.1, #ggLocation
                  0.1, #ggSampleType
                  0.1, #ggSize
                  0.1, #ggCatecholamine
                  0.1, #ggKi67stain.numeric
                  0.1, #ggsdhbihc
                  0.5, #ggFeatures
                  0.3 #ggFeatures
                ), guides = "collect")

ggLocation + shrink_legend + ylab("") + xlab("") +
  ggSampleType + shrink_legend + ylab("") + xlab("") +
  ggSize.monotone + shrink_legend + ylab("") + xlab("") +
  ggCatecholamine.text + shrink_legend + ylab("") + xlab("") +
  ggKi67stain.numeric + shrink_legend + ylab("") + xlab("") +
  #ggKi67stain.category + shrink_legend + ylab("") + xlab("") +
  ggsdhbihc.monotone + shrink_legend + ylab("") + xlab("") +
  ggFeatures +
  plot_layout(ncol = 1, heights = 
                c(0.1, #ggLocation
                  0.1, #ggSampleType
                  0.1, #ggSize
                  0.1, #ggCatecholamine
                  0.1, #ggKi67stain.numeric
                  0.1, #ggsdhbihc
                  0.5 #ggFeatures
                ), guides = "collect")

## Genomic Features
ggHandN + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggDiseaseCourse + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSampleType.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggTMB.cutoff + nox + no_margin + shrink_legend + xlab("") +
  ggSV.cutoff + nox + no_margin + shrink_legend + xlab("") +
  ggWGD.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggChromothripsis + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSignatures  + nox  + no_margin + shrink_legend + ylab("") + xlab("") +
  ggCNVSeg + nox + no_margin + shrink_legend  + xlab("") +
  ggFeatures + no_margin + shrink_legend + ylab("") + 
  plot_layout(ncol = 1, heights = 
                c(0.1, #ggHandN
                  0.1, #ggDiseaseCourse
                  0.1, #ggSampleType
                  0.4, #ggTMB
                  0.4, #ggSV
                  0.1, #ggWGD
                  0.15, #ggChromothripsis
                  1, #ggSignatures
                  1.3, #ggCNVSeg
                  2 #ggFeatures
                ), guides = "collect")

ggMarker.exp +
  ggMarker.exp.clamped +
  plot_layout(ncol = 1, guides = "collect")

###########################
# SDHB MUTATION POSITIONS #
###########################

#googlesheets4::gs4_auth()

SDHB_Domains <- read_sheet("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY", sheet ="SDHB_Domains", range = "A1:E6", col_types = 'cccdd')
SDHB_Changes <- read_sheet("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY", sheet ="SDHB_Events", range = "A1:H86", col_types = 'ccddccdl')

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

SDHB_Changes.small_variants.clincourse <- SDHB_Changes %>%  filter(Type %in% c("Missense","Stop", "Frameshift")) %>%  group_by(AAold, AAnew, AApos_Start, AApos_End, tumour_metastasised) %>% summarise(nPatients=n())

ggSDHB_AA_label <- ggplot(SDHB_Changes.small_variants.clincourse %>% select(AAold, AAnew,AApos_Start) %>% distinct(), 
                          aes(x=AApos_Start, y=1, label=paste0(AAold, AApos_Start, AAnew))) + 
  geom_text(angle=90) + 
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
  coord_cartesian(xlim = c(0,280)) +
  scale_x_continuous(expand = c(0,0))

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
        plot.margin=ggplot2::margin(2,2,1,2)) + guides(color=F) +
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
  geom_segment(aes(x=AApos_Start,y=ylevel, xend=AApos_End, yend=ylevel, color=tumour_metastasised), size=6, linejoin = "round",  lineend = "round") +
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

ggSDHB_AA_label/ggPositionFreq.clincourse/ggPositionDel/ggDomains/ggPositionFreq.type + plot_layout(heights = c(5,8,3,2,8), guides ="collect")

#################
# Age@Dx vs Pos #
#################

ggSDHB_Onset <- ggplot(SDHB_Changes %>%  filter(Type %in% c("Missense","Stop", "Frameshift")) %>% select(Age_dx,AApos_Start) %>% 
                         mutate(AApos_Start_raw=factor(AApos_Start, levels=1:280),
                                AApos_Start_50aa=cut(AApos_Start, breaks = c(0,50,100,150,200,250,300)),
                                AApos_Start_100aa=cut(AApos_Start, breaks = c(0,100,200,300)),
                                AApos_Start_150aa=cut(AApos_Start, breaks = c(0,150,300))) %>% 
                         pivot_longer(cols = c(-Age_dx, -AApos_Start), names_to="break_class", values_to="AA_pos") %>% 
                         mutate(break_class=factor(as.character(break_class), levels=c("AApos_Start_raw","AApos_Start_50aa","AApos_Start_100aa","AApos_Start_150aa"))), 
                       aes(x=AA_pos, y=Age_dx)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_point() +
  theme_bw() +
  facet_wrap("break_class", ncol=1, scale="free_x") +
  ylab("Age at Dx") 

#############
# TERT/ATRX #
#############


plot.data.tertatrxexp <- a5_anno.use %>%  select(A5_ID, TERT_ATRX_Mutation, telhunter_log2_telcontentratio, telhunter_RNA_telcontent)  %>% inner_join(
  A5_exp %>% ungroup() %>%  filter(Symbol %in% c("TERT","ATRX")) %>% select(A5_ID, Symbol,log2cpm) %>% 
    pivot_wider(id_cols=A5_ID, names_from = "Symbol", values_from="log2cpm") ) %>% 
  mutate(Mutation_Status=ifelse(TERT_ATRX_Mutation=="Unknown", 
                                "WT", 
                                paste(TERT_ATRX_Mutation, "Mutant")),
         telhunter_log2_telcontentratio=as.numeric(telhunter_log2_telcontentratio),
         telhunter_RNA_telcontent=as.numeric(telhunter_RNA_telcontent))

ggATRXExpr.boxplot <- ggplot(plot.data.tertatrxexp, aes(x=Mutation_Status, y=ATRX, fill=Mutation_Status)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.05) + 
  scale_x_discrete(labels=function(x){ stringr::str_wrap(x,width = 10)}) +
  scale_fill_manual(values=c(`ATRX Mutant`=ColorPalette[["LightBlue3"]],
                             `TERT Mutant`=ColorPalette[["LightOrange1"]], 
                             WT=ColorPalette[["DarkGrey1"]])) + 
  theme_bw() + guides(fill="none") + ylab("ATRX Expression\n(log2-CPM)") + xlab("")

ggTERTExpr.boxplot <- ggplot(plot.data.tertatrxexp, aes(x=Mutation_Status, y=TERT, fill=Mutation_Status)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.05) + 
  scale_x_discrete(labels=function(x){ stringr::str_wrap(x,width = 10)}) +
  scale_fill_manual(values=c(`ATRX Mutant`=ColorPalette[["LightBlue3"]],
                             `TERT Mutant`=ColorPalette[["LightOrange1"]], 
                             WT=ColorPalette[["DarkGrey1"]])) + 
  theme_bw() + guides(fill="none") + ylab("TERT Expression\n(log2-CPM)") + xlab("")

ggTelLength.boxplot <- ggplot(plot.data.tertatrxexp, aes(x=Mutation_Status, y=telhunter_log2_telcontentratio, fill=Mutation_Status)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.05) + 
  scale_x_discrete(labels=function(x){ stringr::str_wrap(x,width = 10)}) +
  scale_fill_manual(values=c(`ATRX Mutant`=ColorPalette[["LightBlue3"]],
                             `TERT Mutant`=ColorPalette[["LightOrange1"]], 
                             WT=ColorPalette[["DarkGrey1"]])) + 
  theme_bw() + guides(fill="none") + ylab("Telomere Content\n(log2-Tumour/Normal)") + xlab("")

ggTelLength.boxplot <- ggplot(plot.data.tertatrxexp, aes(x=Mutation_Status, y=telhunter_log2_telcontentratio, fill=Mutation_Status)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.05) + 
  scale_x_discrete(labels=function(x){ stringr::str_wrap(x,width = 10)}) +
  scale_fill_manual(values=c(`ATRX Mutant`=ColorPalette[["LightBlue3"]],
                             `TERT Mutant`=ColorPalette[["LightOrange1"]], 
                             WT=ColorPalette[["DarkGrey1"]])) + 
  theme_bw() + guides(fill="none") + ylab("Telomere Content\n(log2-Tumour/Normal)") + xlab("")

ggRNATelContent.boxplot <- ggplot(plot.data.tertatrxexp, aes(x=Mutation_Status, y=telhunter_RNA_telcontent, fill=Mutation_Status)) + 
  geom_boxplot(outlier.alpha = 0) +
  geom_jitter(width=0.05) + 
  scale_x_discrete(labels=function(x){ stringr::str_wrap(x,width = 10)}) +
  scale_fill_manual(values=c(`ATRX Mutant`=ColorPalette[["LightBlue3"]],
                             `TERT Mutant`=ColorPalette[["LightOrange1"]], 
                             WT=ColorPalette[["DarkGrey1"]])) + 
  theme_bw() + guides(fill="none") + ylab("RNA Telomere\nRead Content") + xlab("")
# ggRNATelContent.boxplot + 
#   inset_element(left=0.1,top=0.9, bottom=0.3, right=0.6, 
#                 ggRNATelContent.boxplot + 
#                   coord_cartesian(ylim=c(0,250)) + 
#                   theme(axis.title.y = element_blank(),
#                         axis.text.x = element_blank(),
#                         axis.ticks.x = element_blank()))

plot.data.tert_atrx_types <- features %>% mutate(Patient=gsub("-.", "" , A5_ID)) %>% 
  dplyr::select(Patient, Gene, Event) %>% distinct() %>% 
  filter(Gene %in% c("ATRX","TERT")) %>% 
  filter(!(Gene=="TERT" & Event=="Missense")) %>% 
  group_by(Gene, Event) %>%  dplyr::count(name = "nPatients")

ggTERTATRXEvents <- ggplot(plot.data.tert_atrx_types, aes(x=Gene, y=nPatients, fill=Event)) + geom_col(position = position_stack()) +
  scale_fill_manual(values = c(Missense=ColorPalette[["DarkBlue1"]], `Stop Gained`=ColorPalette[["DarkRed1"]],
                               `Stop Lost`=ColorPalette[["DarkOrange2"]], Frameshift=ColorPalette[["LightBrown2"]], 
                               `Promoter Mutation`=ColorPalette[["Purple3"]], `Structural Variant`= ColorPalette[["LightGreen1"]], 
                               `Splice Acceptor`=ColorPalette[["DarkGrey2"]], `Splice Donor`= ColorPalette[["DarkGrey1"]],
                               `Splice Region`=ColorPalette[["LightGrey1"]], `Homozyg. Del.`= ColorPalette[["DarkGreen1"]])) +
  theme_bw() + theme(axis.title.x = element_blank())



ggHandN + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggDiseaseCourse + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggSampleType.monotone + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggFeatures.TERTATRX + nox + no_margin + shrink_legend + ylab("") + 
  ggTERTCN + nox + no_margin + shrink_legend + ylab("") + 
  ggTERTATRXExp + nox + no_margin + shrink_legend + ylab("") + xlab("") +
  ggTelomere.heat + no_margin + shrink_legend + ylab("") + xlab("") +
  ggATRXExpr.boxplot +
  ggTERTExpr.boxplot +
  ggTelLength.boxplot +
  ggRNATelContent.boxplot +
  ggTERTATRXEvents +
  plot_layout(ncol = 1, heights = 
                c(1, #ggHandN
                  1, #ggDiseaseCourse
                  1, #ggSampleType
                  2, #ggFeatures
                  1, #ggTERTCN
                  2, #ggTERTATRX.exp
                  1, #ggTelomere
                  6
                  
                  
                ), guides = "collect", design = "AAAAA\nBBBBB\nCCCCC\nDDDDD\nEEEEE\nFFFFF\nGGGGG\nHIJKL")


#############
# c-circles #
#############

plot.data <- a5_anno.use %>% 
  mutate(c_circle_result=factor(c_circle_result, levels=c("No data","Negative","Positive - Low","Positive")), 
         clinical_course=case_when(is_primary_or_met=="Primary" & tumour_metastasised == "Yes" ~ "Metastatic Primary", 
                                   is_primary_or_met=="Primary" & tumour_metastasised == "No" ~ "Non-Metastatic Primary", 
                                   is_primary_or_met=="Primary" ~ "Primary - Short Follow-up", 
                                   TRUE ~ is_primary_or_met))
ggplot(plot.data, aes(x=c_circle_result, y=as.numeric(telhunter_log2_telcontentratio), color=TERT_ATRX_Mutation, shape=clinical_course)) + 
  geom_jitter(width = 0.15) + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle=45, hjust=1)) + 
  ylab("log2 T/N telomere ratio")
