###########################################
###########################################
##                                       ##
##   ___________A5 UMAPS______________   ##
##                                       ##
##   Script for producing UMAPs from     ##
##   combined published datasets and     ##
##   A5 WTS/Methylation and small-RNA    ##
##   datasets                            ##
##                                       ##
###########################################
###########################################


setwd("/g/data/pq08/projects/ppgl")
#renv::activate("./a5")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(googlesheets4)
library(patchwork)
library(umap)

##################
# Color Definitions
##################

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

#Color pallets and themes
source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")


############################################
# Import data loaders and run data mergers #
############################################

#clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#Methylation
source("./a5/tertiary/scripts/data_mergers/combine_tcga_comete_a5_methylation_data.r")

#WTS
source("./a5/tertiary/scripts/data_mergers/combine_tcga_flynn_a5_wts_data.r")

#small-rna
source("./a5/tertiary/scripts/data_mergers/combine_tcga_comete_a5_smallrna_data.r")

#############
# Functions #
#############

#Perform UMAP and annotate plotting coordinates with sample annotation
make_umap_plot_data <- function (current_data_set_name, current_annotation_set_name, n_neighbors=10, umap_seed=10) { 
  message("UMAP processing ", current_data_set_name)
  
  #pull inputs from global environment
  data_matrix <- get(current_data_set_name, envir = globalenv())
  anno <- get(current_annotation_set_name, envir = globalenv())
  
  
  if("A5_ID" %in% colnames(anno)){ 
    anno <- anno %>% dplyr::rename(Sample=A5_ID)  
  }
  
  #UMAP config
  set.seed(umap_seed)
  umap_config <- umap.defaults
  umap_config$n_neighbors=n_neighbors
  
  #run UMAP
  umap_result <- umap(t(data_matrix), config = umap_config)
  
  umap_coord <- data.frame(umap_result$layout)
  
  umap_layout_annotated <-  
    tibble(UMAP1 = umap_coord$X1,
           UMAP2 = umap_coord$X2,
           Sample = rownames(umap_coord)) %>% 
    left_join(anno, by="Sample") 
  
  if("differential_group" %in% colnames(umap_layout_annotated)){   
    umap_layout_annotated <-  
      umap_layout_annotated %>% 
      tidyr::separate(differential_group, 
                      into=c("SampleType"), 
                      extra="drop", 
                      sep="_", 
                      remove=F)
  } 
  
  if("Primary_Location_Simplified" %in% colnames(umap_layout_annotated)){   
    umap_layout_annotated <- 
      umap_layout_annotated %>% 
      mutate(Primary_Location_Simplified = 
             gsub("_abdominal|_thoracic|_bladder|_[Ll]eft|_[Rr]ight",
                  "",
                  Primary_Location_Simplified)) 
  }
  
  return(umap_layout_annotated)
}

####
# Prepare datasets
####

a5_wts_logcpm <- cpm(a5_wts_dge_list$qc_ok,log = T)

a5_smallrna_logcpm <- cpm(a5_smallrna_dge_list$qc_ok, log = T)

rownames(a5_methylation_mval) <- a5_methylation_mval$Probe
a5_methylation_mval <- a5_methylation_mval %>%  dplyr::select(-Probe)

data_sets <- list(wts_all="wts_tcga_flynn_a5_lcpm.batch_removed", wts_a5="a5_wts_logcpm",
                  meth_all="meth_tcga_comete_a5_27k_mval.batch_removed", meth_a5="a5_methylation_mval",
                  smallrna_all="smallrna_tcga_comete_a5_lcpm.batch_removed", smallrna_a5="a5_smallrna_logcpm") 

annotation_sets <- list(wts_all="wts_tcga_flynn_a5_anno", wts_a5="a5_anno",
                  meth_all="meth_tcga_comete_a5_27k_anno", meth_a5="a5_anno",
                  smallrna_all="smallrna_tcga_comete_a5_anno", smallrna_a5="a5_anno")



#######
# Generate and annotate UMAPs
#######

umap_coord_annotated <- purrr::map2(.x = data_sets, .y = annotation_sets, .f = make_umap_plot_data)

########
# Plots
########

plot_data__umap__tcga_flynn_a5_comete <-  
  bind_rows(
    umap_coord_annotated[["wts_all"]] %>% mutate(Platform = "WTS"),
    umap_coord_annotated[["meth_all"]] %>% mutate(Platform = "Methylation"),
    umap_coord_annotated[["smallrna_all"]] %>% mutate(Platform = "Small-RNA")
  ) %>%  mutate(Platform=factor(Platform, levels = c("WTS", "Methylation", "Small-RNA")))

plot_data__umap__tcga_flynn_a5_comete <- plot_data__umap__tcga_flynn_a5_comete %>% mutate(Dataset = case_when(
  Dataset == "Flynn" ~ "Flynn et al.",
  Dataset == "E-MTAB-733" ~ "Loriot et al.",
  TRUE ~ Dataset
))
  
plot_data__umap__a5 <-  
  bind_rows(
    umap_coord_annotated[["wts_a5"]] %>% mutate(Platform = "WTS"),
    umap_coord_annotated[["meth_a5"]] %>% mutate(Platform = "Methylation"),
    umap_coord_annotated[["smallrna_a5"]] %>% mutate(Platform = "Small-RNA")
  ) %>%  mutate(Platform=factor(Platform, levels = c("WTS", "Methylation", "Small-RNA")))



####
# Colour = Cluster
####

#TCGA + COMETE + A5 Merged
gg_tcga_a5_umap <- ggplot(mapping= aes(x=UMAP1, 
                                       y=UMAP2, 
                                       color=new_naming, 
                                       shape=Dataset)) + 
  geom_point(data = plot_data__umap__tcga_flynn_a5_comete %>% filter(Dataset == "A5"), alpha = 0.8) +
  geom_point(data = plot_data__umap__tcga_flynn_a5_comete %>% filter(Dataset != "A5")) + 
  #Annotate sample count
  geom_text(data=plot_data__umap__tcga_flynn_a5_comete %>% 
              group_by(Platform) %>% 
              summarise(n=n(), x_max=max(UMAP1), y_max=max(UMAP2), 
                        x_min=min(UMAP1), y_min=min(UMAP2),
                        x_range=x_max-x_min, y_range = y_max-y_min) %>% 
              mutate(Label=paste0("n=",n), x_pos=x_max-(x_range*0.08), y_pos=y_max-(y_range*0.01)),
            aes(label=Label, x=x_pos, y=y_pos, shape=NULL), color="black") + 
  #label E124/E145
  geom_text_repel(data=plot_data__umap__tcga_flynn_a5_comete %>% 
                    filter(Sample %in% c("E124-1", "E145-1"), 
                           Platform=="WTS"),
            aes(label=Sample, shape=NULL), color="black", nudge_y =3) + 
  xlab("UMAP1") + ylab("UMAP2") + 
  labs(color="Subtype") + 
  scale_color_manual(values =subtype_cols) +
  scale_shape_manual(values = c(A5=19, TCGA=0, "Loriot et al."=2, "Flynn et al."=5)) +
  theme_bw() + 
  theme(panel.grid = element_line(colour = "#f6f6f6ff"), aspect.ratio = 1) + 
  facet_wrap("Platform", scale="free")

# A5 Only
gg_a5_umap <- ggplot(plot_data__umap__a5, aes(x=UMAP1, y=UMAP2, color=Primary_Location_Simplified)) + 
  geom_point() + 
  geom_text(data=plot_data__umap__a5 %>% 
              group_by(Platform) %>% 
              summarise(n=n(), x_max=max(UMAP1), y_max=max(UMAP2), 
                        x_min=min(UMAP1), y_min=min(UMAP2),
                        x_range=x_max-x_min, y_range = y_max-y_min) %>% 
              mutate(Label=paste0("n=",n), x_pos=x_max-(x_range*0.08), y_pos=y_max-(y_range*0.01)),
            aes(label=Label, x=x_pos, y=y_pos, shape=NULL), color="black") + 
  geom_text_repel(data=plot_data__umap__a5 %>% 
                    filter(Sample %in% c("E185-1")),
                  aes(label=Sample, shape=NULL), color="black") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  labs(color="Primary site") + 
  scale_color_manual(values = location_cols) +
  #scale_shape_manual() +
  theme_bw() + 
  theme(panel.grid = element_line(colour = "#f6f6f6ff"), aspect.ratio = 1) + 
  facet_wrap("Platform", scale="free", nrow=1)

gg_subtype_umap <- gg_tcga_a5_umap + gg_a5_umap + plot_layout(nrow = 2)

ggsave(filename = "a5/tertiary/results/plots/umap/meth_wts_smallrna__tcga_comete_flynn_a5_umap__seed10_nn10.pdf",
       plot = gg_subtype_umap,
       width = 14, height = 9, 
       units = "in")

####
# Colour = Clinical Outcome
####

#TCGA + COMETE + A5 Merged
gg_tcga_a5_umap <- ggplot(mapping= aes(x=UMAP1, 
                                       y=UMAP2, 
                                       color=Malignancy, 
                                       shape=Dataset)) + 
  geom_point(data = plot_data__umap__tcga_flynn_a5_comete %>% filter(Dataset == "A5"), alpha = 0.8) +
  geom_point(data = plot_data__umap__tcga_flynn_a5_comete %>% filter(Dataset != "A5")) + 
  #Annotate sample count
  geom_text(data=plot_data__umap__tcga_flynn_a5_comete %>% 
              group_by(Platform) %>% 
              summarise(n=n(), x_max=max(UMAP1), y_max=max(UMAP2), 
                        x_min=min(UMAP1), y_min=min(UMAP2),
                        x_range=x_max-x_min, y_range = y_max-y_min) %>% 
              mutate(Label=paste0("n=",n), x_pos=x_max-(x_range*0.08), y_pos=y_max-(y_range*0.01)),
            aes(label=Label, x=x_pos, y=y_pos, shape=NULL), color="black") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  labs(color="Clinical Course") + 
  scale_color_manual(values =c("Metastatic"=specimen_type_cols[["Metastasis"]],
                               "Non-Metastatic"=specimen_type_cols[["NonMetastaticPrimary"]])) +
  scale_shape_manual(values = c(A5=19, TCGA=0, "Loriot et al."=2, "Flynn et al."=5)) +
  theme_bw() + 
  theme(panel.grid = element_line(colour = "#f6f6f6ff"), aspect.ratio = 1) + 
  facet_wrap("Platform", scale="free")

# A5 Only
gg_a5_umap <- ggplot(plot_data__umap__a5, aes(x=UMAP1, y=UMAP2, color=SampleType)) + 
  geom_point() + 
  geom_text(data=plot_data__umap__a5 %>% 
              group_by(Platform) %>% 
              summarise(n=n(), x_max=max(UMAP1), y_max=max(UMAP2), 
                        x_min=min(UMAP1), y_min=min(UMAP2),
                        x_range=x_max-x_min, y_range = y_max-y_min) %>% 
              mutate(Label=paste0("n=",n), x_pos=x_max-(x_range*0.08), y_pos=y_max-(y_range*0.01)),
            aes(label=Label, x=x_pos, y=y_pos, shape=NULL), color="black") + 
  # geom_text_repel(data=plot_data__umap__a5 %>% 
  #                   filter(Sample %in% c("E185-1")),
  #                 aes(label=Sample, shape=NULL), color="black") + 
  xlab("UMAP1") + ylab("UMAP2") + 
  labs(color="Primary site") + 
  scale_color_manual(values = specimen_type_cols) +
  #scale_shape_manual() +
  theme_bw() + 
  theme(panel.grid = element_line(colour = "#f6f6f6ff"), aspect.ratio = 1) + 
  facet_wrap("Platform", scale="free", nrow=1)

gg_clinicalcourse_umap <- gg_tcga_a5_umap + gg_a5_umap + plot_layout(nrow = 2)

ggsave(filename = "a5/tertiary/results/plots/umap/meth_wts_smallrna__tcga_comete_flynn_a5_umap_clinicalcolours__seed10_nn10.pdf",
       plot = gg_clinicalcourse_umap,
       width = 14, height = 9, 
       units = "in")
