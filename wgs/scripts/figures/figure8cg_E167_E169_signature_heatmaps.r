#######################################
# Heatmaps of signature contributions #
# contributions mutations in paired   #
# samples from E167 and E169          #
#                                     #
# Author: Aidan Flynn                 #
# Date: 19/06/2023                    #
#######################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(MutationalPatterns)
library(patchwork)

################
# Data Loaders #
################

load("/g/data/pq08/projects/ppgl/a5/wgs/analysis/mutational_patterns/A5_mutational_patterns.rworkspace")

###########################
# Plot Signature Heatmaps #
###########################

gg_sbs <- list()
gg_indel <- list()

for (patient in c("E167", "E169"))
{
                        
plot_data_sbs <- fit_res$contribution %>% 
  as_tibble(rownames = "Signature") %>% 
  pivot_longer(cols = -Signature, names_to = "A5_ID", values_to = "Contribution") %>% 
  mutate(Signature=factor(Signature, levels=rev(rownames(fit_res$contribution)))) %>% 
  filter(grepl(patient, A5_ID)) %>% 
  group_by(Signature) %>% 
  filter(any(Contribution > 100)) %>% 
  mutate(A5_ID= dplyr::recode(A5_ID, !!!setNames(a5_anno$PublicationID, a5_anno$A5_ID)))

gg_sbs[[patient]] <- ggplot(plot_data_sbs, aes(y=Signature,fill=Contribution, x=A5_ID)) + geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5))

scale_colors <- c("white","#7884baff","#98ca8cff","#eb7b54ff","#c04241ff")  
if(patient=="E167")
{
  gg_sbs[[patient]] <- 
    gg_sbs[[patient]] + 
    scale_fill_gradientn(values = c(0,0.01,0.05,0.1,0.15,0.3,1),
                         colors = scale_colors) 
} else {
  gg_sbs[[patient]] <- gg_sbs[[patient]] + 
    scale_fill_gradientn(colors = scale_colors) 
  
}

gg_sbs[[patient]] <- gg_sbs[[patient]] + theme(panel.grid = element_blank(), axis.title.y = element_blank())

plot_data_indel <- fit_res_indel$contribution %>% 
  as_tibble(rownames = "Signature") %>% 
  pivot_longer(cols = -Signature, 
               names_to = "A5_ID", 
               values_to = "Contribution") %>% 
  mutate(Signature=factor(Signature, 
                          levels=rev(rownames(fit_res_indel$contribution)))) %>% 
  filter(grepl(patient, A5_ID)) %>% 
  group_by(Signature) %>% 
  filter(any(Contribution > 10)) %>% 
  mutate(A5_ID= dplyr::recode(A5_ID, !!!setNames(a5_anno$PublicationID, a5_anno$A5_ID)))

gg_indel[[patient]] <- 
  ggplot(plot_data_indel, 
         aes(y=Signature,fill=Contribution, x=A5_ID)) + 
  geom_tile() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5)) +
  scale_fill_gradientn(colors = scale_colors) +
  theme(panel.grid = element_blank(), 
        axis.title.y = element_blank())

}


gg_sbs[["E167"]] / gg_indel[["E167"]] + plot_layout(heights=c(3,1))

gg_sbs[["E169"]] / gg_indel[["E169"]] + plot_layout(ncol=2)

gg_sbs[["E169"]] / gg_indel[["E169"]] + gg_mgmt_expr$E169 + gg_mmr_expr + plot_layout(design = "ABC\nDDD", widths=c(1,1,2))
