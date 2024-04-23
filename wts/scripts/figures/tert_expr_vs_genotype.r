################################
# Plot of telomere length vs   #
# tert expression highlighting #
# tert mutation status         #
#                              #  
# Author: Aidan Flynn          #
# Date: 08/08/2023             #
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

data_loader_a5_wts_counts(count_file_dir = "/g/data/pq08/projects/ppgl/a5/wts/analysis/htseq/neb/gene/")

#################################
# TERT mutation type annotation #
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

log2cpm <- edgeR::cpm(a5_wts_dge_list[["SDHB"]], log=T) %>% 
  as_tibble(rownames = "ensgid_symbol") %>% 
  separate(ensgid_symbol, into=c("ensgid","symbol"), sep="_", extra = "merge") %>% 
  pivot_longer(cols = c(-ensgid,-symbol), names_to = "A5_ID", values_to = "log2_cpm") %>% 
  group_by(symbol) %>% mutate(log2_cpm_z=(log2_cpm-mean(log2_cpm))/sd(log2_cpm))

############
# Plotting #
############

plot_data <- log2cpm  %>% 
  filter(symbol=="TERT") %>%  
  left_join(a5_anno %>%  
              dplyr::select(A5_ID, `Patient ID`, TERT_ATRX_Mutation, TERT_ATRX_Mutation_Event,
                            telhunter_log2_telcontentratio, Primary_Location_Simplified , 
                            differential_group_sampletype_strict, c_circle_result) %>% 
              mutate(c_circle_result = recode(c_circle_result, "Positive - Low"="Positive"))) %>% 
  mutate(TERT_Mutation=recode(TERT_ATRX_Mutation, "ATRX"="WT")) %>% 
  mutate(TERT_ATRX_Mutation_Event=ifelse(TERT_ATRX_Mutation == "TERT", as.character(TERT_ATRX_Mutation_Event), "None")) %>% 
  arrange(TERT_ATRX_Mutation, log2_cpm_z, `Patient ID`, differential_group_sampletype_strict) 

jp=position_jitter(width = 0.2, seed=80)
ggplot(plot_data ,aes(x=TERT_ATRX_Mutation, 
                      y=log2_cpm, 
                      color=differential_group_sampletype_strict, 
                      shape=TERT_ATRX_Mutation_Event, 
                      group=`Patient ID`)) + 
  geom_point(position = jp, size=3) + 
  geom_path(position = jp, alpha=0.3, linewidth=0.5, color="black") +
  #geom_text(aes(label=A5_ID)) + 
  scale_color_manual(values = sampletype_strict_cols) +
  scale_shape_manual(values = c(`Structural Variant`=8,
                                `Promotor Mutation`=4,
                                `Stop/Frameshift`=15,
                                Missense=17,
                                None=19,
                                Splice=18),
                     name="TERT Mutation Type") +
  
  theme_bw() +
  ylab("TERT Expression (log2 CPM)") +
  xlab("")

