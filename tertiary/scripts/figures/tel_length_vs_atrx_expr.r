################################
# Plot of telomere length vs   #
# atrx expression highlighting #
# c-circle and TERT/ATRX       #
# mutation status              #
#                              #  
# Author: Aidan Flynn          #
# Date: 31/07/2023             #
################################

#################
# Clinical Data #
#################

if(!exists("a5_anno"))
{
  source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
  data_loader_a5_clinical_anno("aidan.flynn@umccr-radio-lab.page", use_cache = T)
}

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")


#######
# WTS #
#######

source("/g/data/pq08/projects/ppgl/a5/wts/scripts/data_loaders/a5_wts_dataloader.r")

data_loader_a5_wts_counts(count_file_dir = "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/truseq/gene/")

#####################
# Prepare Expr Data #
#####################

log2cpm <- edgeR::cpm(a5_wts_dge_list[["SDHB"]], log=T) %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), sep="_", extra = "merge") %>% 
  pivot_longer(cols = c(-ensgid,-symbol), names_to = "A5_ID", values_to = "log2_cpm") %>% 
  group_by(symbol) %>% mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm))/sd(log2_cpm))

############
# Plotting #
############

plot_data <- log2cpm  %>% 
  filter(symbol=="ATRX") %>%  
  left_join(a5_anno %>%  
              dplyr::select(A5_ID, `Patient ID`, TERT_ATRX_Mutation, 
                            telhunter_log2_telcontentratio, Primary_Location_Simplified , 
                            differential_group_sampletype_strict, c_circle_result) %>% 
              mutate(c_circle_result = recode(c_circle_result, "Positive - Low"="Positive"))) %>% 
  arrange(c_circle_result, differential_group_sampletype_strict, TERT_ATRX_Mutation)
  
  
ggplot(plot_data ,aes(x=as.numeric(telhunter_log2_telcontentratio), 
             y=log2_cpm_z, 
             color=c_circle_result,
             fill=differential_group_sampletype_strict, 
             shape=TERT_ATRX_Mutation, 
             group=`Patient ID`)) + 
  geom_point(size=3, stroke=1) + 
  geom_line( alpha=0.3, size=0.5, color="black") +
  #geom_text(aes(label=A5_ID)) + 
  scale_fill_manual(values = sampletype_strict_cols) +
  scale_color_manual(values = c("Positive"="Red", 
                                "Negative"="Black", 
                                "No data"="Grey")) +
  scale_shape_manual(values = c("WT"=21,
                                "TERT"=24,
                                "ATRX"=25)) +
  #scale_shape_manual(values = c(8,18,15,17,19,4)) +
  # scale_size_manual(values=c(1.3,2.2,2.2)) +
  theme_bw() +
  theme(legend.position = "bottom", aspect.ratio = 1) +
  guides(fill = guide_legend(override.aes = list(shape=21)),
         color = guide_legend(override.aes = list(shape=1, linetype=c(0,0,0)))) +
  ylab("ATRX Expression (z-scaled)") +
  xlab("Tumour/Normal Telomere Content Ratio (log2)")

