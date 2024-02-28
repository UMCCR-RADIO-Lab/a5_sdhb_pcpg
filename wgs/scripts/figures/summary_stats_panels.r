library(ggplot2)
library(patchwork)
################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

# source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
# data_loader_cna_segs()
# 
# #source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/pcnt_genome_altered/percent_genome_altered.r")

a5_anno.use <- a5_anno %>% filter(Exclude=="N") %>% 
  mutate(differential_group_anatomy=case_when(
    A5_ID == "E185-1" ~ "Head_neck",
    A5_ID == "E129-1" ~ "Abdominal_Thoracic",
    A5_ID == "E128-1" ~ "Abdominal_Thoracic",
    A5_ID %in% c("E148-1","E155-1") ~ "Head_neck",
    TRUE ~ differential_group_anatomy)) %>% 
  mutate(differential_group_anatomy=
           dplyr::case_match(differential_group_anatomy, 
                         "Abdominal_Thoracic" ~ "Abdominal/Thoracic", 
                         "Head_neck" ~ "Head and Neck/Aortic")) %>% 
  mutate(MetastaticLineage=case_when(
    differential_group_sampletype_strict %in% c("Metastatic primary", "Metastasis") ~ "Confirmed",
    differential_group_sampletype_strict %in% c("Primary (metastasis reported)") ~ "Unconfirmed",
    TRUE ~ "Not applicable")) %>% 
  mutate(TERT_ATRX_Mutation = dplyr::case_when(A5_ID == "E171-1" ~ "TERT",
                                               TRUE ~ TERT_ATRX_Mutation)) %>% 
  mutate(wgs_tmb = as.numeric(wgs_tmb),
         structural_variant_count = as.numeric(structural_variant_count),
         Largest_primary_dimensions_cm=as.numeric(Largest_primary_dimensions_cm))
  

#######
# Percent Genome Altered
#######



plot_data <- a5_anno.use %>% 
               dplyr::select(A5_ID, differential_group_sampletype,
                             differential_group_sampletype_strict, 
                             differential_group_anatomy, MetastaticLineage, 
                             proportion_genome_diploid, TERT_ATRX_Mutation) 

#Clinical

p1 <- ggplot(plot_data %>%  filter(differential_group_sampletype_strict != "Primary (short follow up)"), 
       aes(x=differential_group_sampletype, 
           y=proportion_genome_diploid)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(mapping = aes(shape=MetastaticLineage), width=0.2) + 
  scale_shape_manual(values=c(Confirmed=8, Unconfirmed=1, `Not applicable`=19)) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + 
  #geom_text(mapping=aes(label=A5_ID)) + 
  facet_grid(.~differential_group_anatomy, scales="free_x", space="free") +
  coord_cartesian(y=c(0.5,1))

p2 <- ggplot(plot_data %>%  filter(differential_group_sampletype_strict != "Primary (short follow up)"), 
       aes(x=TERT_ATRX_Mutation, 
           y=proportion_genome_diploid)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_jitter(mapping = aes(shape=MetastaticLineage, color=differential_group_sampletype), width=0.2) + 
  scale_shape_manual(values=c(Confirmed=8, Unconfirmed=1, `Not applicable`=19)) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + 
  #geom_text(mapping=aes(label=A5_ID)) + 
  facet_grid(.~differential_group_anatomy, scales="free_x", space="free") +
  scale_color_manual(values=specimen_type_cols) +
  coord_cartesian(y=c(0.5,1))
 

#######
# TERT/ATRX/WT - TMB
#######



jp=position_jitter(width = 0.3, seed=99)
p3 <- ggplot(a5_anno.use %>%  filter(differential_group_sampletype_strict != "Primary (short follow up)"), 
       aes(x=TERT_ATRX_Mutation, 
           y=wgs_tmb)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(mapping = aes(shape=MetastaticLineage, color=differential_group_sampletype), position = jp) + 
  geom_path(mapping = aes(group=`Patient ID`), color="grey", position = jp) +
  scale_shape_manual(values=c(Confirmed=8, Unconfirmed=1, `Not applicable`=19)) +
  scale_color_manual(values=specimen_type_cols) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + 
  facet_grid(.~differential_group_anatomy, scales="free_x", space="free") + coord_cartesian(ylim=c(0,5))

jp=position_jitter(width = 0.3, seed=99)
p4 <- ggplot(a5_anno.use %>%  filter(differential_group_sampletype_strict != "Primary (short follow up)") %>% 
         mutate(treatment_exposed = ifelse(Resection_post_dna_damaging_treatment %in% c("No","Uncertain"),"No", "Yes")), 
       aes(x=TERT_ATRX_Mutation, 
           y=wgs_tmb)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(mapping = aes(shape=treatment_exposed, color=differential_group_sampletype), position = jp) + 
  geom_path(mapping = aes(group=`Patient ID`), color="grey", position = jp) +
  scale_shape_manual(values=c(Yes=8, No=19)) +
  scale_color_manual(values=specimen_type_cols) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + 
  facet_grid(.~differential_group_anatomy, scales="free_x", space="free") + coord_cartesian(ylim=c(0,5))


#######
# Largest Primary - TMB
#######

plot_data <- a5_anno.use %>% 
  filter(differential_group_sampletype_strict != "Primary (short follow up)") %>% 
  filter(!is.na(Largest_primary_dimensions_cm)) %>% 
  mutate(Largest_primary_dimensions_cm_category = cut(Largest_primary_dimensions_cm, breaks=c(0,5,50), labels=c("0-5cm", "5+cm"))) 

jp=position_jitter(width = 0.3, seed=99)
p5 <-ggplot(plot_data, 
       aes(
         x=Largest_primary_dimensions_cm_category, 
         #x=paste(Largest_primary_dimensions_cm_category, TERT_ATRX_Mutation), 
         y=wgs_tmb)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(mapping = aes(shape=MetastaticLineage, color=differential_group_sampletype), position = jp) + 
  geom_path(mapping = aes(group=`Patient ID`), color="grey", position = jp) +
  scale_shape_manual(values=c(Confirmed=8, Unconfirmed=1, `Not applicable`=19)) +
  scale_color_manual(values=specimen_type_cols) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + 
  facet_grid(.~differential_group_anatomy, scales="free_x", space="free") + coord_cartesian(ylim=c(0,3.5))


#######
# Largest Primary - TERT/ATRX
#######

plot_data <- a5_anno.use %>% 
  filter(differential_group_sampletype_strict != "Primary (short follow up)") %>% 
  filter(!is.na(Largest_primary_dimensions_cm)) %>% 
  dplyr::select(`Patient ID`, Largest_primary_dimensions_cm, TERT_ATRX_Mutation, differential_group_anatomy) %>% 
  distinct()

jp=position_jitter(width = 0.2, seed=99)
p6 <- ggplot(plot_data, 
       aes(
         y=Largest_primary_dimensions_cm, 
         #x=paste(Largest_primary_dimensions_cm_category, TERT_ATRX_Mutation), 
         x=TERT_ATRX_Mutation)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = jp) + 
  geom_path(mapping = aes(group=`Patient ID`), color="grey", position = jp) +
  scale_shape_manual(values=c(Confirmed=8, Unconfirmed=1, `Not applicable`=19)) +
  scale_color_manual(values=specimen_type_cols) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + 
  facet_grid(.~differential_group_anatomy, scales="free_x", space="free") #+ coord_cartesian(ylim=c(0,3.5))

#######
# TERT/ATRX/WT - SV count
#######

plot_data <- a5_anno.use %>% 
  filter(differential_group_sampletype_strict != "Primary (short follow up)") %>% 
  mutate(structural_variant_count= as.numeric(structural_variant_count))

jp=position_jitter(width = 0.3, seed=99)
p7 <- ggplot(plot_data, 
       aes(x=TERT_ATRX_Mutation, 
           y=structural_variant_count)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(mapping = aes(shape=MetastaticLineage, color=differential_group_sampletype), position = jp) + 
  geom_path(mapping = aes(group=`Patient ID`), color="grey", position = jp) +
  scale_shape_manual(values=c(Confirmed=8, Unconfirmed=1, `Not applicable`=19)) +
  scale_color_manual(values=specimen_type_cols) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + 
  facet_grid(.~differential_group_anatomy, scales="free_x", space="free") + 
  scale_y_log10() #+
  #coord_cartesian(ylim=c(0,110))

jp=position_jitter(width = 0.3, seed=99)
p8 <- ggplot(plot_data, 
       aes(x=differential_group_sampletype, 
           y=structural_variant_count)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(mapping = aes(shape=MetastaticLineage, color=TERT_ATRX_Mutation), position = jp) + 
  geom_path(mapping = aes(group=`Patient ID`), color="grey", position = jp) +
  scale_shape_manual(values=c(Confirmed=8, Unconfirmed=1, `Not applicable`=19)) +
  scale_color_manual(values=driver_cols) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + 
  facet_grid(.~differential_group_anatomy, scales="free_x", space="free") + 
  scale_y_log10() #+
  #coord_cartesian(ylim=c(0,110))


jp=position_jitter(width = 0.3, seed=99)
p9 <- ggplot(plot_data %>% filter(is.na(chromothriptic_event)), 
       aes(x=differential_group_sampletype, 
           y=structural_variant_count)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(mapping = aes(shape=MetastaticLineage, color=TERT_ATRX_Mutation), position = jp) + 
  geom_path(mapping = aes(group=`Patient ID`), color="grey", position = jp) +
  scale_shape_manual(values=c(Confirmed=8, Unconfirmed=1, `Not applicable`=19)) +
  scale_color_manual(values=driver_cols) +
  theme_bw() + theme(axis.text.x = element_text(angle=90, vjust = 0.5,hjust=1)) + 
  facet_grid(.~differential_group_anatomy, scales="free_x", space="free") #+ 
  #scale_y_log10() #+
#coord_cartesian(ylim=c(0,110))

pdf(file = "/g/data/pq08/projects/ppgl/a5/wgs/results/figures/genome_summary_plots.pdf", width = 14, height = 5, onefile = T)
p1 + p2  + plot_annotation(title ="% genome diploid")
  p3 + p4   + plot_annotation(title ="TMB")
  p5 + p6   + plot_annotation(title ="largest primary") 
  p7 + p8 + ylab("SV count (log 10)") + p9 + ylab("SV count (chrtps ex.)")  + plot_annotation(title ="Strutural variant count") 
dev.off()

