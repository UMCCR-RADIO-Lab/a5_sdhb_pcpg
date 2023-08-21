################################
# Plot of telomere length vs   #
# atrx expression highlighting #
# c-circle and TERT/ATRX       #
# mutation status              #
#                              #  
# Author: Aidan Flynn          #
# Date: 31/07/2023             #
################################

library(ggplot2)
library(dplyr)
library(tidyr)

#################
# Clinical Data #
#################

if(!exists("a5_anno"))
{
  source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
  data_loader_a5_clinical_anno("aidan.flynn@umccr-radio-lab.page", use_cache = T)
}

source("/g/data/pq08/projects/ppgl/a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")


# #######
# # WTS #
# #######
# 
# source("/g/data/pq08/projects/ppgl/a5/wts/scripts/data_loaders/a5_wts_dataloader.r")
# 
# data_loader_a5_wts_counts(count_file_dir = "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/truseq/gene/")

#################################
# ATRX mutation type annotation #
#################################

offline = F 
if(offline) { 
  features <- read.delim("/g/data/pq08/projects/ppgl/a5/offline_cache/sample_features.tsv")
} else {
  feature_sheet <- as_sheets_id("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY")
  gs4_auth("aidan.flynn@umccr-radio-lab.page")
  features <- read_sheet(ss = feature_sheet,sheet = "SampleFeatures", col_types = "c")
}

features <- features %>%  filter(as.logical(Include))
features$Gene <- factor(as.character(features$Gene), features %>% group_by(Gene) %>% dplyr::count() %>% arrange(n) %>% pull(Gene))
features$Event <- factor(as.character(features$Event), levels = c("Missense", "Stop Gained", "Stop Lost", "Frameshift", "Promotor Mutation", "Structural Variant",
                                                                  "Splice Acceptor", "Splice Region", "Homozyg. Del."))
if (!("TERT_ATRX_Mutation_Event" %in% colnames(a5_anno)))
{
  a5_anno <- a5_anno %>%  left_join(features %>% 
                                      filter(Gene %in% c("TERT","ATRX")) %>%  
                                      #mutate(Event=paste(Gene, Event, sep =" - ")) %>% 
                                      dplyr::select(A5_ID,Gene, Event) %>% 
                                      dplyr::rename(TERT_ATRX_Mutation_Event=Event, TERT_ATRX_Mutation=Gene)) %>% 
    mutate(TERT_ATRX_Mutation_Event=ifelse(TERT_ATRX_Mutation=="WT","None",as.character(TERT_ATRX_Mutation_Event)))
  
  a5_anno$TERT_ATRX_Mutation_Event <- factor(as.character(a5_anno$TERT_ATRX_Mutation_Event), levels=c("Missense", "Splice Region", "Splice Acceptor", "Stop Gained", "Frameshift", "Structural Variant", "Promotor Mutation", "None"))
  levels(a5_anno$TERT_ATRX_Mutation_Event) <- c("Missense", "Splice", "Splice", "Stop/Frameshift", "Stop/Frameshift", "Structural Variant", "Promotor Mutation", "None")
  
  a5_anno$TERT_ATRX_Mutation <- factor(as.character(a5_anno$TERT_ATRX_Mutation), levels=c("ATRX", "TERT", "WT"))
}

#####################
# Prepare Expr Data #
#####################
# 
# log2cpm <- edgeR::cpm(a5_wts_dge_list[["SDHB"]], log=T) %>% 
#   as_tibble(rownames = "ensgid_symbol") %>% 
#   separate(ensgid_symbol, into=c("ensgid","symbol"), sep="_", extra = "merge") %>% 
#   pivot_longer(cols = c(-ensgid,-symbol), names_to = "A5_ID", values_to = "log2_cpm") %>% 
#   group_by(symbol) %>% mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm))/sd(log2_cpm))

############
# Plotting #
############

plot_data <- a5_anno %>%  
              dplyr::select(A5_ID, `Patient ID`, TERT_ATRX_Mutation, TERT_ATRX_Mutation_Event,
                            telhunter_log2_telcontentratio, Primary_Location_Simplified , 
                            differential_group_sampletype_strict, c_circle_result) %>% 
              mutate(c_circle_result = recode(c_circle_result, "Positive - Low"="Weakly Positive"),
                     c_circle_result=factor(c_circle_result, levels=c("No data","Negative", "Weakly Positive", "Positive"))) %>% 
  arrange(TERT_ATRX_Mutation, c_circle_result, differential_group_sampletype_strict)
  
pj=position_jitter(width=0.2, seed=10)  
ggplot(plot_data ,aes(y=as.numeric(telhunter_log2_telcontentratio), 
                      x=TERT_ATRX_Mutation, 
                      color=c_circle_result,
                      fill=differential_group_sampletype_strict, 
                      group=`Patient ID`)) + 
  geom_point(position=pj, size=3, stroke=1, shape=21) + 
  geom_path(position=pj, alpha=0.3, size=0.5, color="black", linetype=2) +
  #geom_text(aes(label=A5_ID)) + 
  scale_fill_manual(values = sampletype_strict_cols) +
  scale_color_manual(values = c("Positive"="red", 
                                "Weakly Positive"="orange",
                                "Negative"="black", 
                                "No data"="Grey")) +
  # scale_shape_manual(values = c("WT"=21,
  #                               "TERT"=24,
  #                               "ATRX"=25)) +
  #scale_shape_manual(values = c(8,18,15,17,19,4)) +
  # scale_size_manual(values=c(1.3,2.2,2.2)) +
  theme_bw() +
  theme(legend.position = "bottom", aspect.ratio = 1) +
  guides(fill = guide_legend(override.aes = list(shape=21)),
         color = guide_legend(override.aes = list(shape=1, linetype=c(0,0,0,0)))) +
  ylab("Tumour/Normal Telomere Content Ratio (log2)") +
  xlab("")

