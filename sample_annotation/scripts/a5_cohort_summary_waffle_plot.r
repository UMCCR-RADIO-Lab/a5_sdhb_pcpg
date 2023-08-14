#################################################
# Generate a waffle plot of cohort sample types #
# Author: Aidan Flynn                           #
# Date: 03/04/2023                              #
# Languages: R                                  #
#################################################

library(dplyr)
library(tidyr)

swap_position <- function(Sample_ID, position_data, new_x, new_y)
{
  swap_reciever_differential_group_sampletype_strict <- position_data %>% 
    filter(A5_ID==Sample_ID) %>% 
    pull(differential_group_sampletype_strict) %>% 
    as.character()
  
  swap_reciever_primarylocation <- position_data %>% 
    filter(A5_ID==Sample_ID) %>% 
    pull(Primary_Location_Simplified) %>% 
    as.character()
  
  swap_donor_sample_type <- position_data %>% 
    filter(x_pos == new_x, y_pos == new_y) %>% 
    pull(differential_group_sampletype_strict) %>% 
    as.character()
  
  swap_donor_sampleid <- position_data %>% 
    filter(x_pos == new_x, y_pos == new_y) %>% 
    pull(A5_ID) %>% 
    as.character()
  
  if(swap_donor_sample_type != swap_reciever_sample_type) {
    stop("Donor/Reciever sample type do not match")
  }
  
  swap_receiver_x <- position_data$x_pos[position_data$A5_ID==Sample_ID]
  swap_receiver_y <- position_data$y_pos[position_data$A5_ID==Sample_ID]
  
  position_data$x_pos[position_data$A5_ID==Sample_ID] <- new_x
  position_data$y_pos[position_data$A5_ID==Sample_ID] <- new_y
  
  position_data$x_pos[position_data$A5_ID==swap_donor_sampleid] <- swap_receiver_x
  position_data$y_pos[position_data$A5_ID==swap_donor_sampleid] <- swap_receiver_y
  
  return(position_data)
  
}


if(!require(waffle)){
  devtools::install_github("hrbrmstr/waffle")
  library(waffle)
}

if(!exists("a5_anno"))
{
  source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
  data_loader_a5_clinical_anno("aidan.flynn@umccr-radio-lab.page", use_cache = T)
}


##########
# Colors #
##########

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")
location_cols[["Extraadrenal (abdominal/thoracic)"]] <-  location_cols[["Extraadrenal"]]
location_cols[["Extraadrenal (aortic)"]] <-  location_cols[["Extraadrenal_aortic"]]
location_cols[["Extraadrenal (bladder)"]] <-  location_cols[["Extraadrenal_bladder"]]

################################
# Prepare data for waffle plot #
################################

waffle_data <- a5_anno %>% filter(Exclude != "Y") %>%  
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
                                                     "Extraadrenal (aortic)",
                                                     "Adrenal",
                                                     "Head and neck",
                                                     "Unspecified"))) %>% 
  dplyr::select(A5_ID, `Patient ID`, differential_group_sampletype_strict, Primary_Location_Simplified) %>% 
  group_by(`Patient ID`) %>% 
  mutate(group=ifelse(n()>1, `Patient ID`, NA))

###########
# Prepare squares for clinical outcome
###########

waffle_squares <- waffle_data %>% group_by(differential_group_sampletype_strict) %>% 
  dplyr::count()

###########
# Prepare dots for primary location
###########

n_samples=(a5_anno %>% filter(Exclude=="N") %>% nrow())
waffle_dots <- waffle_data %>% 
  ungroup() %>% 
  arrange(differential_group_sampletype_strict, Primary_Location_Simplified) %>% 
  mutate(y_pos=rep(1:10,10)[1:n_samples], 
         x_pos=rep(1:10,each=10)[1:n_samples]) 


#Assign linked/paired samples to set positions for joining by lines
#E128
waffle_dots <- swap_position("E128-1",waffle_dots, new_x = 1, new_y = 10)
waffle_dots <- swap_position("E128-2",waffle_dots, new_x = 2, new_y = 10)
#E229
waffle_dots <- swap_position("E229-1",waffle_dots, new_x = 1, new_y = 9)
waffle_dots <- swap_position("E229-2",waffle_dots, new_x = 2, new_y = 9)
#E122
waffle_dots <- swap_position("E122-1",waffle_dots, new_x = 8, new_y = 7)
waffle_dots <- swap_position("E122-2",waffle_dots, new_x = 8, new_y = 6)
#E132
waffle_dots <- swap_position("E132-2",waffle_dots, new_x = 8, new_y = 3)
waffle_dots <- swap_position("E132-1",waffle_dots, new_x = 9, new_y = 3)
#E143
waffle_dots <- swap_position("E143-3",waffle_dots, new_x = 7, new_y = 10)
waffle_dots <- swap_position("E143-1",waffle_dots, new_x = 8, new_y = 10)
waffle_dots <- swap_position("E143-2",waffle_dots, new_x = 9, new_y = 10)
#E146
waffle_dots <- swap_position("E146-1",waffle_dots, new_x = 8, new_y = 1)
waffle_dots <- swap_position("E146-2",waffle_dots, new_x = 9, new_y = 1)
#E158
waffle_dots <- swap_position("E158-1",waffle_dots, new_x = 8, new_y = 2)
waffle_dots <- swap_position("E158-2",waffle_dots, new_x = 9, new_y = 2)
#E159
waffle_dots <- swap_position("E159-2",waffle_dots, new_x = 4, new_y = 9)
waffle_dots <- swap_position("E159-3",waffle_dots, new_x = 7, new_y = 9)
waffle_dots <- swap_position("E159-1",waffle_dots, new_x = 8, new_y = 9)
waffle_dots <- swap_position("E159-4",waffle_dots, new_x = 9, new_y = 9)
#E166
waffle_dots <- swap_position("E166-1",waffle_dots, new_x = 9, new_y = 6)
waffle_dots <- swap_position("E166-2",waffle_dots, new_x = 10, new_y = 6)
#E167
waffle_dots <- swap_position("E167-1",waffle_dots, new_x = 9, new_y = 5)
waffle_dots <- swap_position("E167-2",waffle_dots, new_x = 10, new_y = 5)
#E169
waffle_dots <- swap_position("E169-1",waffle_dots, new_x = 9, new_y = 4)
waffle_dots <- swap_position("E169-2",waffle_dots, new_x = 10, new_y = 4)
#E225
waffle_dots <- swap_position("E225-1",waffle_dots, new_x = 7, new_y = 8)
waffle_dots <- swap_position("E225-2",waffle_dots, new_x = 8, new_y = 8)


#Re-sort to bring similar combinations back into blocks after required position swaps  
#Filter rows without a "required" location, isolate values and sort 
waffle_dots_unpinned_values <- 
  waffle_dots %>% 
    filter(is.na(group)) %>% 
    dplyr::select(-x_pos,-y_pos) %>% 
    arrange(differential_group_sampletype_strict, Primary_Location_Simplified)

#Filter rows without a "required" location, isolate locations and sort  
waffle_dots_unpinned_locations <- 
  waffle_dots %>% 
    filter(is.na(group)) %>% 
    dplyr::select(x_pos,y_pos) %>% 
    arrange(x_pos,y_pos)

#Re-attach locations to sorted values
waffle_dots_unpinned <- 
  bind_cols(waffle_dots_unpinned_values, 
            waffle_dots_unpinned_locations)

#Filter rows with a "required" location
waffle_dots_pinned <- 
  waffle_dots %>% 
  filter(!is.na(group))

#Recombine resorted free-moving samples with pinned samples
waffle_dots <- 
  bind_rows(waffle_dots_unpinned, 
            waffle_dots_pinned) %>% 
  arrange(x_pos,y_pos)

#################
# Generate plot #
#################

ggplot( ) +
  geom_waffle(data = waffle_squares, mapping = aes(fill = differential_group_sampletype_strict, values = n), n_rows = 10, size = 0.33, colour = "white", flip = FALSE, radius = unit(4, "pt")) +
  geom_point(data=waffle_dots, mapping = aes(x=x_pos,y=y_pos), color="black", size=4) +
  geom_point(data=waffle_dots, mapping = aes(x=x_pos,y=y_pos, color=Primary_Location_Simplified), size=3) +
  geom_line(data=waffle_dots %>% filter(A5_ID != "E159-2"), mapping = aes(x=x_pos,y=y_pos, group=`Patient ID`)) +
  geom_line(data=waffle_dots %>% filter(A5_ID %in% c("E159-2","E159-3")), 
            mapping = aes(x=x_pos,y=y_pos, group=`Patient ID`),
            linetype=2) +
  scale_fill_manual(
    name = "Clinical Course",
    values = c("#6ea668", "#6ea6cc", "#d6d6d4", "#f18d73", "#ee7474ff", "#969696", "#c04241ff"),
    labels = c("Non-metastatic primary", "Primary (short follow up)", "Non-metastatic local recurrence",
               "Primary (metastasis present)", "Metastatic primary","Metastatic local recurrence","Metastasis")
  ) +
  scale_color_manual(
    name = "Anatomical Location (Primary)",
    values = location_cols)  +
  coord_equal() +
  theme_enhance_waffle() +
  theme(axis.ticks = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())

