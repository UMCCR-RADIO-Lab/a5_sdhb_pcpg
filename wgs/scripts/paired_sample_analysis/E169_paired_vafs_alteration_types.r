#####################################
# Analysis of VAF densities of      #
# Transition/Transversion mutations #
# in paired samples from E169       #
#                                   #
# Author: Aidan Flynn               #
# Date: 19/06/2023                  #
#####################################

################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/wgs_dataloaders.r")
data_loader_somatic_variants(quickload = T)

load("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mutational_patterns/A5_mutational_patterns.rworkspace")

################
# Isolate E169 data and rev/comp REF/ALT to match SNV96
################

comp <- c("A","T","C","G")
names(comp) <- c("T","A","G","C")

E169_VAF <- a5_somatic_variants %>% 
  ungroup() %>% 
  rowwise() %>% 
  filter(A5_ID %in% c("E169-1","E169-2"),
         REF %in% c("A","T","C","G"),
         ALT %in% c("A","T","C","G")) %>% 
  mutate(ALT2=case_when(
    REF %in% c("A","G") ~ comp[[ALT]],
    TRUE ~ ALT),
    REF2=case_when(
      REF == "A" ~ "T",
      REF == "G" ~ "C",
      TRUE ~ REF
    )
  ) %>% dplyr::select(A5_ID, REF2,ALT2, Tumour_AF) %>% 
  mutate(Alteration= paste0(REF2,">",ALT2)) %>% 
  ungroup()
E169_VAF$Tumour_AF <- as.numeric(E169_VAF$Tumour_AF)

####################
# Plot VAF density #
####################

ggplot(E169_VAF, aes(x=Tumour_AF, after_stat(count),color=Alteration)) + 
  geom_density() + 
  scale_color_manual(values=c("#03bdee","#000000","#e52a25","#cdc9ca","#a3ce62","#edc5c5")) +
  facet_wrap("A5_ID") + 
  theme_bw()




plot_data <- fit_res$contribution %>% 
  as_tibble(rownames = "Signature") %>% 
  pivot_longer(cols = -Signature, names_to = "A5_ID", values_to = "Contribution") %>% 
  mutate(Signature=factor(Signature, levels=rownames(fit_res$contribution))) %>% 
  filter(A5_ID %in% c("E169-1","E169-2"), Contribution > 100) 


ggplot(plot_data, aes(x=Signature,y=Contribution)) + geom_col() + facet_wrap("A5_ID") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) 

plot_data <- fit_res_indel$contribution %>% 
  as_tibble(rownames = "Signature") %>% 
  pivot_longer(cols = -Signature, names_to = "A5_ID", values_to = "Contribution") %>% 
  mutate(Signature=factor(Signature, levels=rownames(fit_res_indel$contribution))) %>% 
  filter(A5_ID %in% c("E169-1","E169-2")) 


ggplot(plot_data , aes(x=Signature,y=Contribution)) + geom_col() + facet_wrap("A5_ID") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) 

ggplot(a5_anno %>%  
         filter(Exclude != "Y", wgs_tmb <500) %>% 
         mutate(color=ifelse(A5_ID %in% c("E167-1","E167-2"), "red","black")), 
       aes(x=reorder(PublicationID, wgs_tmb), y=wgs_tmb, fill=color)) +
  geom_col() +
  scale_fill_identity() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) + 
  ylab("TMB") + 
  xlab("")
