#################################################
# Generate an upset plot of technologies        #
# applied across the cohort                     #
# Author: Aidan Flynn                           #
# Date: 30/05/2023                              #
# Languages: R                                  #
#################################################

library(dplyr)
library(tidyr)
library(ggupset)
library(grid)

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")

if(!exists("a5_anno")) {
  data_loader_a5_clinical_anno(use_cache = T) }

gs4_auth("aidan.flynn@umccr-radio-lab.page")
a5_assays <- read_sheet("1hnXdXI29KvvuLxsaTBID1-EbE7mTSk6bG05HcSFfhgo", 
                      sheet="AssaysPerformed", 
                      col_types = "c")


upset_data <- a5_assays %>%
  filter(Excluded=="N") %>% 
  dplyr::select(A5_ID, WGS, `RNA-Seq`, Methylation, `small-RNA`, snATAC, `snRNA-Seq`, `C-circle`) %>%  
  pivot_longer(cols = -A5_ID,names_to = "Assay", values_to = "Completed") %>%
  filter(Completed=="1") %>% 
  group_by(A5_ID) %>% 
  summarise(Assay=list(Assay))
  
gg_upset <- ggplot(upset_data, aes(x=Assay)) +
  geom_bar() +
  #geom_text(mapping=aes(y=1, label=after_stat(count))) +
  scale_x_upset() + ylab("Samples (n)")
  
assay_groups <- upset_data %>% 
  rowwise() %>% 
  mutate(Assay=toString(Assay)) %>% 
  group_by(Assay) %>% 
  mutate(n_assay_combo=n())

assay_combo_by_anatomy <- a5_anno %>% filter(Exclude=="N") %>% 
  dplyr::select(A5_ID, Primary_Location_Simplified) %>% 
  mutate(Primary_Location_Simplified=
           dplyr::recode(Primary_Location_Simplified, "Adrenal_left"="Adrenal",
                         "Adrenal_right"="Adrenal",
                         "Head_neck"="Head and neck",
                         "Extraadrenal_abdominal"="Extraadrenal",
                         "Extraadrenal_thoracic"="Extraadrenal",
                         "Extraadrenal_abdominal"="Extraadrenal",
                         "Extraadrenal_bladder"="Extraadrenal",
                         "Extraadrenal_thoracic_cardiac"="Extraadrenal")) %>% 
  inner_join(assay_groups) %>% 
  group_by(Assay, Primary_Location_Simplified, n_assay_combo) %>% dplyr::count(name="n_samples")

assay_combo_by_anatomy <- assay_combo_by_anatomy %>% arrange(-n_assay_combo, Assay, Primary_Location_Simplified) %>% 
  mutate(Assay=factor(Assay, levels=unique(.$Assay)))

gg_assay_by_anatomy <- ggplot(assay_combo_by_anatomy, aes(x=Assay, y=n_samples, fill=Primary_Location_Simplified)) + 
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = location_cols, name="Anatomical Location\n(Primary)") + 
  theme_bw() +
  scale_y_continuous(breaks = seq(0,60,5))

gg_upset_gtable <- ggplotGrob(gg_upset)
gg_assay_by_anatomy_gtable <- ggplotGrob(gg_assay_by_anatomy)

gg_assay_by_anatomy_gtable_panel <-  gtable::gtable_filter(gg_assay_by_anatomy_gtable, "panel")
gg_assay_by_anatomy_gtable_axis_l <-  gtable::gtable_filter(gg_assay_by_anatomy_gtable, "axis-l")

gg_upset_gtable$grobs[[6]] <- gg_assay_by_anatomy_gtable_panel$grobs[[1]]
gg_upset_gtable$grobs[[3]] <- gg_assay_by_anatomy_gtable_axis_l$grobs[[1]]

gg_assay_by_anatomy_gtable_guide_box <-  gtable::gtable_filter(gg_assay_by_anatomy_gtable, "guide-box")
gg_upset_gtable <- gtable::gtable_add_cols(gg_upset_gtable, gg_assay_by_anatomy_gtable$widths[[9]], 9)
gg_upset_gtable <- gtable::gtable_add_rows(gg_upset_gtable, heights = gg_assay_by_anatomy_gtable$heights[[9]], 15)

gg_upset_gtable <- gtable::gtable_add_grob(x = gg_upset_gtable,
                                           grobs = gg_assay_by_anatomy_gtable_guide_box,
                                           t = 7,
                                           l = 9,
                                           b = 9,
                                           r = 10,
                                           z = 14,
                                           name="guide-box")

grid.newpage()
grid.draw(gg_upset_gtable)
