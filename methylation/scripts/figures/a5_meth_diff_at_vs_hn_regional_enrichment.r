##########################################
# This script generates plots            #
# focused on co-located clusters of      #
# genes found to be DM between           #
# abdo-thoracic and head and neck PPGL   #
#                                        #
# Author: Aidan Flynn                    #
# Date: 29/08/2023                       #
#                                        #
##########################################

library(ggplot2)
library(ggrepel)
library(patchwork)

setwd("/g/data/pq08/projects/ppgl")


#################
# Run DM script #
#################

#Performs DM
# creates globals:
# - epic_array_annotation_hg38
# - diff_meth_result
source("./a5/methylation/scripts/a5_methylation_analysis_v2.r")

########################################
# Annotate DM TopTable with probe data #
########################################

top_table <- diff_meth_result[["hn"]][["Chromaffin_vs_Non_chromaffin"]][["top_control2_adj"]]

top_table <- top_table %>% as_tibble(rownames = "Name") %>% 
  left_join(epic_array_annotation_hg38 %>% 
              as_tibble() %>% 
              mutate(POS_hg38=as.numeric(Start_hg38)+1) %>% 
              dplyr::select(Name, CHR_hg38, POS_hg38, UCSC_RefGene_Name) )
  

##################################
# Extract significant gene peaks #
##################################

region_probes <- top_table %>% 
  filter(CHR_hg38 %in% c("chr7", "chr12", "chr17", "chr17", "chr20"))  %>% 
  dplyr::select(Name, F.p.BH, CHR_hg38, POS_hg38, UCSC_RefGene_Name)
 
peak_probes <- top_table %>% 
  filter(
    (CHR_hg38=="chr7"  & POS_hg38 > 2.65*10^7  & POS_hg38 < 2.72*10^7) |
      (CHR_hg38=="chr12"  & POS_hg38 > 5.38*10^7 & POS_hg38 < 5.41*10^7) | 
      (CHR_hg38=="chr17" & POS_hg38 > 3.355*10^7 & POS_hg38 < 3.5*10^7) | 
      (CHR_hg38=="chr17" & POS_hg38 > 4.12*10^7 & POS_hg38 < 4.18*10^7) | 
      (CHR_hg38=="chr17" & POS_hg38 > 4.82*10^7 & POS_hg38 < 4.9*10^7) |
      (CHR_hg38=="chr20" & POS_hg38 > 6.1*10^7 & POS_hg38 < 6.2*10^7))  %>% 
  mutate(cytoband=case_when(
    CHR_hg38=="chr7" ~ "chr7p15.2", 
    CHR_hg38=="chr12" ~ "chr12q13.13", 
    CHR_hg38=="chr17" & POS_hg38 > 4.68 *10^7 ~ "chr17q21.32", 
    CHR_hg38=="chr17" & POS_hg38 > 4.02 *10^7 ~ "chr17q21.2", 
    CHR_hg38=="chr17" & POS_hg38 > 3.35*10^7 ~ "chr17q12", 
    CHR_hg38=="chr20" ~ "chr20q13.33")) %>% 
  dplyr::select(Name, F.p.BH, CHR_hg38, POS_hg38, UCSC_RefGene_Name, cytoband)



#####################
# Prepare plot data #
#####################

########
# Join per sample expression 
########

plot_data <- region_probes %>% inner_join(b_vals[region_probes$Name,] %>% 
                                         as_tibble(rownames = "Name") %>% 
                                         pivot_longer(cols=-Name, 
                                                      names_to = "A5_ID", 
                                                      values_to = "beta"), 
                                       by=c("Name"))

########
# Join clinical annotation
########

plot_data <- plot_data %>% inner_join(a5_anno %>% dplyr::select(A5_ID, differential_group_anatomy))

########
# Compute delta beta
########

beta_delta <- plot_data %>% 
  filter(differential_group_anatomy %in% c("Abdominal_Thoracic", "Head_neck")) %>% 
  group_by(Name, differential_group_anatomy) %>% 
  summarise(mean_beta=mean(beta)) %>% 
  pivot_wider(id_cols = "Name", names_from = differential_group_anatomy, values_from = mean_beta) %>% 
  mutate(beta_delta = Abdominal_Thoracic - Head_neck) %>% 
  dplyr::select(Name, beta_delta)


########
# Make Plot friendly names for clinical values
########

plot_data <- plot_data %>% 
  mutate(`Tumour Location`=dplyr::recode(differential_group_anatomy, 
                                         "Abdominal_Thoracic"="Abdominal/Thoracic",
                                         "Head_neck"="Head and neck"))

########
# Finalise plot data
########

plot_data_manhatten <- plot_data %>% 
  dplyr::select(Name, F.p.BH, CHR_hg38, POS_hg38, UCSC_RefGene_Name) %>%  
  distinct() %>% 
  inner_join(beta_delta) %>% mutate(CHR_hg38 = factor(CHR_hg38, levels=c("chr7", "chr12", "chr17", "chr20")))

plot_data_regions <- plot_data %>% 
  inner_join(peak_probes %>% dplyr::select(Name, cytoband))

########
# Factorise probes in chromosome positional order
########

plot_data_regions <- plot_data_regions %>% 
  arrange(CHR_hg38, POS_hg38) %>% 
  mutate(Name=factor(Name, levels=unique(.$Name)))

 
##################
# Generate Plots #
##################

########
# Gene expression by anatomical type for each region
########

gg_meth_region <- list()
for (current_cytoband in c("chr7p15.2", "chr12q13.13", "chr17q12", "chr17q21.2", "chr17q21.32", "chr20q13.33"))
{
  current_plot_data <- plot_data_regions  %>% 
    filter(cytoband == current_cytoband, 
           `Tumour Location` %in% c("Abdominal/Thoracic","Head and neck")) %>% 
    mutate(Significant = ifelse(F.p.BH < 0.01, "p < 0.01", "Not significant")) %>% 
    arrange(POS_hg38) %>% 
    group_by(Name, CHR_hg38, POS_hg38, UCSC_RefGene_Name, `Tumour Location`, cytoband, Significant) %>% 
    summarise(median=median(beta), q1=quantile(beta,0.25), q3=quantile(beta,0.75))
  
  probe_region <- paste0(unique(current_plot_data$CHR_hg38),":", min(current_plot_data$POS_hg38),"-",max(current_plot_data$POS_hg38))
  
  polygon_data <- bind_rows(current_plot_data %>% arrange(Name, POS_hg38,`Tumour Location`) %>% dplyr::select(-q3) %>% dplyr::rename(beta=q1),
                            current_plot_data %>% arrange(desc(Name),`Tumour Location`) %>% dplyr::select(-q1) %>% dplyr::rename(beta=q3))
  
  gene_labels <- current_plot_data %>% filter(UCSC_RefGene_Name != "") %>% 
    separate_rows(UCSC_RefGene_Name, sep=";") %>% distinct() %>% 
    group_by(UCSC_RefGene_Name) %>% 
    mutate(region_type=case_when(row_number() == 1 ~ "start",
                                 row_number() == floor(n()/2) ~ "middle",
                                 row_number() == n() ~ "end",
    TRUE ~ NA_character_)) %>% 
    filter(!is.na(region_type)) %>% 
    dplyr::select(Name, POS_hg38, UCSC_RefGene_Name, cytoband, region_type) %>% 
    pivot_wider(id_cols = c(UCSC_RefGene_Name, cytoband),
                names_from = region_type, 
                values_from = c(Name,POS_hg38)) %>% 
    mutate(Name_middle=ifelse(is.na(Name_middle), as.character(Name_start), as.character(Name_middle)),
           POS_hg38_middle=ifelse(is.na(POS_hg38_middle), POS_hg38_start, POS_hg38_middle),
           y=sample(seq(-0.5,0,0.005), 1))
                            
  
  gg_meth_region[[current_cytoband]] <- 
    ggplot(data=current_plot_data) + 
    geom_point(mapping=aes(x=Name, 
                           y=median, 
                           # fill=`Tumour Location`,
                           color=`Tumour Location`,
                           shape= Significant),
               size=0.8) +
    # geom_line(mapping=aes(x=Name,
    #                      y=median,
    #                      # fill=`Tumour Location`,
    #                      color=`Tumour Location`,
    #                      group=`Tumour Location`)) +
    geom_segment(data = gene_labels,
      mapping=aes(x=Name_start, 
                  y=y,
                  yend=y,
                  xend=Name_end),
              ) +
    geom_text(data = gene_labels,
                 mapping=aes(x=Name_middle,
                             y=y,
                             label=UCSC_RefGene_Name)) +
    geom_polygon(data=polygon_data,
      mapping=aes(x=Name, 
                          y=beta, 
                          fill=`Tumour Location`,
                          #color=`Tumour Location`,
                          group=`Tumour Location`,
                  ), alpha=0.3) + 
    scale_x_discrete(labels=\(labs) { paste0(labs, rep(c("-------------------","--------------------------------------",""), times=ceiling(length(labs)/3))[1:length(labs)]) }) +
    # scale_fill_manual(values=c(`Abdominal/Thoracic`=location_cols[["Head_neck"]],
    #                            `Head and neck`=location_cols[["Extraadrenal"]])) +
    scale_color_manual(values=c(`Abdominal/Thoracic`=location_cols[["Head_neck"]],
                                `Head and neck`=location_cols[["Extraadrenal"]])) +
    scale_shape_manual(values=c(`p < 0.01`=8, "Not significant" = 20)) +
    theme_bw() +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          #axis.text.x = element_text(angle = 90, vjust=0.5,hjust=1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) + 
    ylab(bquote('beta')) +
    facet_wrap("cytoband", scales="free_x", ncol=1) + 
    xlab(paste0("Probes (",probe_region,")")) 
  
}


########
# Manhattan style plot (chr. position by -log10(p)) from limma top-table
########

find_scaled_zero <- function(pdata, feature_column) {
  scaled_zero_value = tibble(raw=pdata[[feature_column]],
                             scaled=scales::rescale(x = pdata[[feature_column]], to = c(0, 1))) %>%
    arrange(abs(raw)) %>%
    slice_head(n=1) %>%
    pull(scaled)
  return(scaled_zero_value)
}

gg_manhattan <- ggplot(plot_data_manhatten %>% filter(CHR_hg38 %in% c("chr7","chr12","chr17")), aes(x=POS_hg38, y=-log10(F.p.BH), color=beta_delta)) + 
  geom_point(size=1) + 
  geom_hline(yintercept = -log10(0.01), linetype=2, color="red") +
  facet_wrap("CHR_hg38", scales="free_x") +
  scale_x_continuous(labels = \(x) paste(x/10^6, "Mb")) +
  scale_color_gradientn(colours = c("blue", "grey", "orange", "red"), 
                        values = c(0, find_scaled_zero(plot_data_manhatten, "beta_delta"),0.6, 1),
                        name= "delta mean-beta") +
  ylab(bquote('-'~log[10]~'(adjusted p-value)')) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust=1,
                                   hjust=1),
        axis.title.x = element_blank()) 

################
# Format plots #
################

#layout_design <- "AABBB\nCCCCC\nDDDEE"

gg_manhattan +
gg_meth_region[["chr7p15.2"]] +
  gg_meth_region[["chr12q13.13"]] +
  gg_meth_region[["chr17q12"]] + 
  gg_meth_region[["chr17q21.2"]] + 
  gg_meth_region[["chr17q21.32"]] +
  plot_layout(guides="collect", 
              #design = layout_design, 
              nrow=6,
              #heights = c(1,1.5,1)
              ) 


