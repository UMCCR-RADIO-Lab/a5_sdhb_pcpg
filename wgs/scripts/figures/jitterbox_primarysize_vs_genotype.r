###########################################
# This script generates a jitter-box plot #
# of Tumour Mutation Burden versus        #
# clinical outcome/category               #
#                                         #
# Author: Aidan Flynn                     #
# Date: 17/08/2023                        #
###########################################

library(ggplot2)

################
# Data Loaders #
################

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account = "aidan.flynn@umccr-radio-lab.page", use_cache = T)

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

############
# Plotting #
############

plot_data <- a5_anno %>%
  filter(!is.na(Largest_primary_dimensions_cm), differential_group_sampletype_strict != "Primary (short follow up)") %>% 
  filter(differential_group_anatomy != "Ambiguous", A5_ID != "E146-1") %>% 
  dplyr::select(`Patient ID`, Largest_primary_dimensions_cm, 
                TERT_ATRX_Mutation, differential_group_anatomy, 
                cell_of_origin, tumour_metastasised) %>% 
  filter() %>% 
  mutate(
    class = paste(TERT_ATRX_Mutation, ifelse(tumour_metastasised == "Yes", "(metastasis reported)", "(metastasis not reported)")),
    group = paste(cell_of_origin, class,sep="#")) %>% 
  distinct() %>% 
  mutate(class = factor(class, levels = c("TERT (metastasis reported)", "ATRX (metastasis reported)", 
                                          "WT (metastasis reported)", "WT (metastasis not reported)")),
         cell_of_origin=factor(cell_of_origin, levels=c( "Chromaffin", "Non_chromaffin")))

##t-test

comp_list <- combn(x = unique(plot_data$group), m = 2) %>% 
  t() %>% 
  data.frame() %>% 
  setNames(c("g1", "g2")) %>% 
  filter(grepl("^Non_chromaffin",g1) & grepl("^Non_chromaffin",g2) |
           grepl("^Chromaffin",g2) & grepl("^Chromaffin",g1))

t_results <- purrr::pmap(comp_list, .f = \(g1, g2) {
    t_test_data <- plot_data %>% filter(group %in% c(g1, g2))

    t_results[[paste(g1,g2,sep="_vs_")]] <- t.test(Largest_primary_dimensions_cm~group,  
                                     t_test_data, 
                                     alternative = c("two.sided"))
    
  }) %>% 
  setNames(paste(comp_list$g1, comp_list$g2, sep="_vs_"))


sig_tests <- purrr::map2(.x = t_results, 
            .y = names(t_results), 
            .f = \(t_result, comparison) {
              comp_members <- stringr::str_split_1(string = comparison, pattern = "_vs_")
              df_result <- data.frame(group1=comp_members[[1]], group2 = comp_members[[2]], pval = t_result$p.value)
              return(df_result)}) %>% bind_rows() %>% 
  filter(pval < 0.05)

sig_tests <- sig_tests %>% 
  separate_wider_delim(cols = c(group1, group2), names = c("cell_of_origin", "class"), delim = "#", names_sep = "_") %>% 
  dplyr::select(-group2_cell_of_origin) %>% 
  dplyr::rename(cell_of_origin=group1_cell_of_origin)

sig_tests$group1_class <- factor(sig_tests$group1_class, levels=levels(plot_data$class))
sig_tests$group2_class <- factor(sig_tests$group2_class, levels=levels(plot_data$class))
sig_tests$cell_of_origin <- factor(sig_tests$cell_of_origin,  levels =  levels(plot_data$cell_of_origin))

pvalue_lines_df <- tibble(y=seq(17, (17+(nrow(sig_tests)-1)), 1),
                          yend=seq(17, (17+(nrow(sig_tests)-1)), 1),
                          x=as.numeric(sig_tests$group1_class), 
                          xend=as.numeric(sig_tests$group2_class),
                          cell_of_origin = sig_tests$cell_of_origin)

pvalue_text_df <- tibble(y=seq(17.5, (17.5+(nrow(sig_tests)-1)), 1),
       x=as.numeric(sig_tests$group1_class)+((as.numeric(sig_tests$group2_class)-as.numeric(sig_tests$group1_class))/2), 
       label=paste0("p=",(round(sig_tests$pval,4))),
       cell_of_origin =  factor(sig_tests$cell_of_origin, levels=levels(plot_data$cell_of_origin)))

jp = position_jitter(width = 0.1, seed=25)
gg_primarysize <- ggplot(plot_data, 
                    aes(x=class, y=Largest_primary_dimensions_cm)) + 
  geom_boxplot(outlier.alpha = 0) + 
  geom_point(position = jp) + 
  geom_segment(data=pvalue_lines_df,
               mapping=aes(x=x,y=y,xend=xend,yend=yend, color=NULL,fill=NULL,group=NULL)) +
  geom_text(data=pvalue_text_df,
            mapping=aes(x=x,y=y,label=label, color=NULL,fill=NULL,group=NULL)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, vjust = 1,hjust=1)) +
  ylab("Largest primary dimensions (cm)") +
  xlab("") + 
  facet_grid(~cell_of_origin, space="free_x", scales = "free_x")

