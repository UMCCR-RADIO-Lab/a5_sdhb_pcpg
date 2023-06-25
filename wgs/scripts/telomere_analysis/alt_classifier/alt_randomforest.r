

#FEATURES
#telomere content tumor/control log2 ratio, 
#number of telomere insertions, 
#number of breakpoints, 
#distance of TGAGGG, TCAGGG, TTGGGG, TTCGGG, and TTTGGG singletons 

gs4_auth("aidan.flynn@umccr-radio-lab.page")
A5_anno <- read_sheet("1hnXdXI29KvvuLxsaTBID1-EbE7mTSk6bG05HcSFfhgo", col_types = "c")
A5_anno <- A5_anno %>% dplyr::rename(`A5_ID`=`A5 ID`)
A5_anno <- A5_anno %>% filter(!is.na(`A5_ID`))

A5_anno <- A5_anno %>% dplyr::select(-TERT_ATRX, -TERT_ATRX_Event)



feature_sheet <- as_sheets_id("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY")
gs4_auth("aidan.flynn@umccr-radio-lab.page")
features <- read_sheet(ss = feature_sheet,sheet = "SampleFeatures", col_types = "c")
features <- features %>%  filter(as.logical(Include))
features$Gene <- factor(as.character(features$Gene), features %>% group_by(Gene) %>% dplyr::count() %>% arrange(n) %>% pull(Gene))
features$Event <- factor(as.character(features$Event), levels = c("Missense", "Stop Gained", "Stop Lost", "Frameshift", "Promotor Mutation", "Structural Variant",
                                                                  "Splice Acceptor", "Splice Region", "Homozyg. Del."))
A5_anno <- A5_anno %>%  left_join(features %>% 
                              filter(Gene %in% c("TERT","ATRX")) %>%  
                              #mutate(Event=paste(Gene, Event, sep =" - ")) %>% 
                              dplyr::select(A5_ID,Gene, Event) %>% 
                              dplyr::rename(TERT_ATRX_Event=Event, TERT_ATRX=Gene)) %>% 
  mutate(TERT_ATRX=ifelse(is.na(TERT_ATRX),"WT", as.character(TERT_ATRX)), TERT_ATRX_Event=ifelse(is.na(TERT_ATRX_Event),"None",as.character(TERT_ATRX_Event)))

A5_anno$TERT_ATRX_Event <- factor(as.character(A5_anno$TERT_ATRX_Event), levels=c("Missense", "Splice Region", "Splice Acceptor", "Stop Gained", "Frameshift", "Structural Variant", "Promotor Mutation", "None"))
levels(A5_anno$TERT_ATRX_Event) <- c("Missense", "Splice", "Splice", "Stop/Frameshift", "Stop/Frameshift", "Structural Variant", "Promotor Mutation", "None")

A5_anno$TERT_ATRX <- factor(as.character(A5_anno$TERT_ATRX), levels=c("ATRX", "TERT", "WT"))


set.seed(123)

# Read Data
setwd("C:/Users/AFFLY/OneDrive - The University of Melbourne/ResearchData/RADIO/A5/Scripts/ALT_Classifier")
sieverling <- read.delim("Sieverling_altrf_input.tsv")
sieverling.input <- sieverling %>%  dplyr::select(icgc_specimen_id, TMM_associated_mut, TMM_associated_mut_summary, tel_content_log2, telomere_insertions, 
                                             breakpoints, TGAGGG_singleton_dist, TCAGGG_singleton_dist, 
                                             TTGGGG_singleton_dist, TTCGGG_singleton_dist, TTTGGG_singleton_dist)

#Remove cases with missing input values
sieverling.input <- sieverling.input %>% filter(!apply(sieverling.input,1,anyNA))

ExcludeFeatures <- c("TMM_associated_mut")

#Exlude non-training features
sieverling.input <- sieverling.input[,-which(colnames(sieverling.input) %in% ExcludeFeatures)] 

#Subsample data to max 100 cases
sieverling.input <- sieverling.input %>% group_by(TMM_associated_mut_summary) %>% slice_sample(n = 100,replace = F) %>%  ungroup()

#Factorise outcome variable
sieverling.input <- sieverling.input %>% filter(TMM_associated_mut_summary %in% c("ATRX_DAXX_trunc", "TERT_mod", "Other")) 
sieverling.input$TMM_associated_mut_summary <- factor(as.character(sieverling.input$TMM_associated_mut_summary), levels = c("ATRX_DAXX_trunc", "TERT_mod", "Other"))
sieverling.input <- sieverling.input %>%  
  dplyr::rename(PID=icgc_specimen_id) %>% 
  mutate(ALTExpected=factor(ifelse(TMM_associated_mut_summary=="ATRX_DAXX_trunc", "ATRX_DAXX_trunc", "Other"), levels=c("ATRX_DAXX_trunc", "Other"))) %>% 
  dplyr::select(-TMM_associated_mut_summary)

# Data Partition
ind <- sample(2, nrow(sieverling.input), replace = TRUE, prob = c(0.7, 0.3))
train <- sieverling.input[ind==1,-1]

test <- sieverling.input[ind==2,-1]

# Random Forest
library(randomForest)

rf <- randomForest(ALTExpected~., data=train,
                   ntree = 100,
                   mtry = 3,
                   importance = TRUE,
                   proximity = TRUE)
print(rf)

# Prediction & Confusion Matrix - train data
library(caret)
p1 <- predict(rf, train)
confusionMatrix(p1, train$ALTExpected)

# # Prediction & Confusion Matrix - test data
p2 <- predict(rf, test)
confusionMatrix(p2, test$ALTExpected)

# Error rate of Random Forest
plot(rf)

# Tune mtry
t <- tuneRF(train[,-which(colnames(train)=="ALTExpected")], train$ALTExpected,
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 100,
            trace = TRUE,
            improve = 0.05)

# No. of nodes for the trees
hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")

# Variable Importance
varImpPlot(rf,
           sort = T,
           n.var = 8,
           main = "Top - Variable Importance")
importance(rf)
varUsed(rf)

# Partial Dependence Plot
partialPlot(rf, train, telomere_insertions, "ATRX_DAXX_trunc")

# Extract Single Tree
getTree(rf, 1, labelVar = TRUE)

# Multi-dimensional Scaling Plot of Proximity Matrix
MDSplot(rf, factor(train$ALTExpected),k = 3, cex=1.4, palette = c("red","blue","grey"))

a5th.summary <- read.delim("c:/ResearchData/RADIO/A5/Data/WGS/hg38/Telomerehunter/telomerehunter_a5_summary.tsv")
a5th.summary.tumour <- a5th.summary %>% filter(sample=="tumor")
a5th.summary.control <- a5th.summary %>% filter(sample=="control")
a5th.summary.ratio <- a5th.summary %>% filter(sample=="log2(tumor/control)")

a5th.summary.ratio %>% 
  inner_join(A5_anno %>% 
               dplyr::select(A5_ID, `Assumed driver of metastasis`) %>% 
               dplyr::rename(PID=A5_ID, ALT_Driver=`Assumed driver of metastasis`) %>% 
               mutate(PID=gsub("-(.)","-T0\\1",PID))) %>% 
  ggplot(aes(x=ALT_Driver,y=tel_content)) + geom_point() + geom_text(aes(label=PID))
  

a5th.summary.ratio %>% select(-total_reads, -read_length, -repeat_threshold_set, 
                              -repeat_threshold_used, -tel_reads, -intratel_reads, 
                              -gc_bins_for_correction, -total_reads_with_tel_gc, -sample) %>% 
  pivot_longer(cols = -PID, names_to="Feature", values_to="Value") %>% filter(!is.na(Value)) %>% 
  inner_join(A5_anno %>% 
               select(A5_ID, `Assumed driver of metastasis`) %>% 
               rename(PID=A5_ID, ALT_Driver=`Assumed driver of metastasis`) %>% 
               mutate(PID=gsub("-(.)","-T0\\1",PID))) %>%  
  ggplot(aes(x=ALT_Driver, y=Value, color=ALT_Driver)) + geom_boxplot() + geom_jitter(width=0.1, height=0) + facet_wrap("Feature", scales = "free_y") + guides(color=F)


a5th.singletons <- read.delim("c:/ResearchData/RADIO/A5/Data/WGS/hg38/Telomerehunter/telomerehunter_a5_singletons.tsv")

a5th.singletons.dist <- a5th.singletons %>%  
  dplyr::select(PID, pattern, distance_to_expected_singleton_log2_ratio) %>% 
  pivot_wider(id_cols = PID, names_from=pattern, values_from=distance_to_expected_singleton_log2_ratio) %>% 
  rename_with(~ paste0(.x, "_singleton_dist"), .cols=-PID) 

a5th.norm_tvr <- read.delim("c:/ResearchData/RADIO/A5/Data/WGS/hg38/Telomerehunter/telomerehunter_a5_normalized_TVR_counts.tsv")

a5th.top_tvr <- read.delim("c:/ResearchData/RADIO/A5/Data/WGS/hg38/Telomerehunter/telomerehunter_a5_TVR_top_contexts.tsv")

a5th.telosv <- read.delim("c:/ResearchData/RADIO/A5/Data/WGS/hg38/Telomerehunter/telomere_rep_inserts_count.tsv")

a5.test <- a5th.singletons.dist %>%  dplyr::select(PID, TGAGGG_singleton_dist, TCAGGG_singleton_dist, 
           TTGGGG_singleton_dist, TTCGGG_singleton_dist, TTTGGG_singleton_dist) %>% 
  full_join(a5th.singletons %>% dplyr::select(PID, tel_content_log2_ratio) %>% dplyr::rename(tel_content_log2=tel_content_log2_ratio) %>% distinct()) %>% 
  full_join(a5th.telosv %>% dplyr::rename(PID=A5_ID,telomere_insertions=n_telomere_insert_events, breakpoints=total_svs))


p3 <- predict(rf, a5.test)
a5.test$ALTExpected <- predict(rf, a5.test)
p3 <- predict(rf, a5.test, type = "prob")
a5.test <- cbind(a5.test, data.frame(p3)) %>% 
  inner_join(A5_anno %>%  dplyr::select(A5_ID, `Assumed driver of metastasis`) %>%  dplyr::rename(PID=A5_ID, ALT_Driver=`Assumed driver of metastasis`) %>% mutate(PID=gsub("-(.)","-T0\\1",PID)) )

ggplot(a5.test, aes(y=ATRX_DAXX_trunc, x=ALT_Driver, color=ALT_Driver)) + geom_point()

confusionMatrix(p3, a5.test$ALTExpected)

a5.test %>%  dplyr::select(-Other,-ATRX_DAXX_trunc, -ALTExpected) %>% pivot_longer(cols = c(-PID,-ALT_Driver), names_to="Feature", values_to="Value") %>% 
  ggplot(aes(x=ALT_Driver, y=Value, color=ALT_Driver)) + geom_boxplot() + geom_jitter(width=0.1, height=0) + facet_wrap("Feature", scales = "free_y") + guides(color=F)


a5th.summary.rna <- read.delim("c:/ResearchData/RADIO/A5/Data/WGS/hg38/Telomerehunter/telomerehunter_a5_rna_summary.tsv")
a5th.summary.rna %>% 
  inner_join(A5_anno %>% dplyr::select(A5_ID, TERT_ATRX, TERT_ATRX_Event) %>% dplyr::rename(PID=A5_ID)) %>% 
  ggplot(aes(x=TERT_ATRX, y=tel_content)) + geom_jitter(width = 0.2, height = 0) 

p1 <- a5th.summary.rna %>% 
  inner_join(A5_anno %>% dplyr::select(A5_ID, `Patient ID`, TERT_ATRX, TERT_ATRX_Event, is_primary_or_met, tumour_metastasised) %>% 
               mutate(isMet_didMet=paste(is_primary_or_met, tumour_metastasised, sep="-"), 
                      isMet_didMet=gsub("Recurrent-.+","Recurrent",isMet_didMet)) %>% 
               dplyr::rename(PID=A5_ID)) %>% 
  inner_join(a5th.singletons %>% dplyr::select(PID, tel_content_log2_ratio) %>% mutate(PID=gsub("T0","",PID)) %>% distinct()) %>% 
  ggplot(aes(x=tel_content_log2_ratio, y=tel_content, color=TERT_ATRX, shape=is_primary_or_met, group=`Patient ID`)) + 
  geom_point() + 
  geom_line(alpha=0.3) +
  #geom_text_repel(aes(label=PID)) + 
  ylab("RNA Telomere Content") + 
  xlab("log2 Tumour/Normal Telomere Ratio") +
  theme_bw()



plot.data <- a5th.top_tvr   %>% pivot_wider(id_cols = c(PID, pattern), 
                               names_from = Sample, 
                               values_from=c(Bases, Percent)) %>% 
  mutate(PID=gsub("T0","",PID)) %>% 
  inner_join(A5_anno %>% 
               dplyr::select(A5_ID, `Patient ID`, TERT_ATRX, TERT_ATRX_Event, `Assumed driver of metastasis`, `Tumour location simplified` ,is_primary_or_met, tumour_metastasised) %>% 
               dplyr::rename(PID=A5_ID, 
                             location_simple=`Tumour location simplified`,
                             Primary_Met=is_primary_or_met, 
                             Met_Case=tumour_metastasised)) %>% 
               
  inner_join(a5th.summary.rna %>% dplyr::select(PID, tel_content) %>% dplyr::rename(TERRA_content=tel_content))

  ggplot(plot.data, aes(x=Percent_tumor, y=Percent_control, label=PID, color=TERT_ATRX, shape=Primary_Met, group=`Patient ID`)) + 
  #geom_text(size=3) + 
  geom_line(alpha=0.3) +   
  geom_point() +
  facet_wrap("pattern", scales="free")
  
  ggplot(plot.data, aes(y=Percent_tumor/Percent_control, x=pattern, color=TERT_ATRX)) + 
    #geom_text(size=3) + 
    #geom_line(alpha=0.3) +   
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position=position_jitterdodge(jitter.width = 0.1)) 
    #facet_wrap("pattern", scales="free")
  

  
############
#  Expression vs Telomere content
###########  
  
  
  
plot.list <- list()
for (g in c("ATRX", "TERT"))
{  
  plot.list[[g]] <-  A5_zexp %>% 
    filter(SYMBOL==g) %>%  
    left_join(A5_anno %>%  
                dplyr::select(A5_ID, `Patient ID`, TERT_ATRX, TERT_ATRX_Event, 
                              telhunter_log2_telcontentratio, `Assumed driver of metastasis`, 
                              `Tumour location simplified` ,is_primary_or_met, 
                              tumour_metastasised) %>% 
                dplyr::rename(location_simple=`Tumour location simplified`,
                              Primary_Met=is_primary_or_met)) %>% 
    mutate(Primary_Met=factor(as.character(Primary_Met), levels=c("Primary","Recurrent", "Metastatic")))  %>% 
    
    ggplot(aes(x=as.numeric(telhunter_log2_telcontentratio), y= zExp, color=TERT_ATRX, shape=tumour_metastasised, group=`Patient ID`, size=Primary_Met)) + 
    geom_point() + 
    geom_line( alpha=0.3, size=0.5) +
    #geom_text(aes(label=A5_ID)) + 
    scale_color_manual(values = c(ColorPalette[["DarkBlue1"]], 
                                  ColorPalette[["LightOrange1"]], 
                                  ColorPalette[["LightGreen1"]])) +
    scale_shape_manual(values = c(16,8,17)) +
    #scale_shape_manual(values = c(8,18,15,17,19,4)) +
    scale_size_manual(values=c(1.3,2.2,2.2)) +
    theme_bw() +
    ylab(paste(g," Expression (z-scaled)")) +
    xlab("Tumour/Normal Telomere Content Ratio (log2)")
  
}
  plot.list[[1]]   +  plot.list[[2]] + geom_text(aes(label=A5_ID), size=5) + plot_layout(guides = "collect")
    
  A5_zexp %>% 
    filter(SYMBOL=="ATRX") %>%  
    left_join(A5_anno %>%  
                dplyr::select(A5_ID, `Patient ID`, TERT_ATRX, TERT_ATRX_Event, telhunter_log2_telcontentratio, `Assumed driver of metastasis`, `Tumour location simplified` ,is_primary_or_met, tumour_metastasised) %>% 
                dplyr::rename(location_simple=`Tumour location simplified`,
                              Primary_Met=is_primary_or_met)) %>% mutate(Primary_Met=factor(as.character(Primary_Met), levels=c("Primary","Recurrent", "Metastatic")))  %>% 
    ggplot(aes(x=as.numeric(telhunter_log2_telcontentratio), y= zExp, color=TERT_ATRX, shape=TERT_ATRX_Event, size=Primary_Met)) + 
    geom_point() + 
    #geom_text(aes(label=A5_ID)) + 
    scale_color_manual(values = c(ColorPalette[["DarkBlue1"]], 
                                  ColorPalette[["LightOrange1"]], 
                                  ColorPalette[["LightGreen1"]])) +
    #scale_shape_manual(values = c(16,8,17)) +
    scale_shape_manual(values = c(8,18,15,17,19,4)) +
    scale_size_manual(values=c(1.3,2.2,2.2)) +
    theme_bw() +
    ylab("ATRX Expression (z-scaled)") +
    xlab("Tumour/Normal Telomere Content Ratio (log2)")
 
  
########
# TVRs normal by all reads
########
  
  anova.data <- a5th.norm_tvr %>%  mutate(PID=gsub("T0","",PID)) %>%  
    mutate(log2_ratio_count_norm_by_all_reads=log2(Count_norm_by_all_reads_T/Count_norm_by_all_reads_C)) %>% 
    dplyr::select(PID, Pattern, log2_ratio_count_norm_by_all_reads) %>% 
    inner_join(A5_anno %>% 
                 dplyr::select(A5_ID, TERT_ATRX) %>% 
                 dplyr::rename(PID=A5_ID)) %>% dplyr::select(-PID) 
  anova.data$Pattern <- as.factor(anova.data$Pattern)
  anova.data$TERT_ATRX <- as.factor(anova.data$TERT_ATRX)
  summary(anova.data)
  
  sig_patterns <- c()
  for (p in levels(anova.data$Pattern))
  {
    
    ooa <- aov(formula = log2_ratio_count_norm_by_all_reads ~ TERT_ATRX, data = anova.data %>% filter(Pattern==p, !is.infinite(log2_ratio_count_norm_by_all_reads)))
    if(summary(ooa)[[1]][['Pr(>F)']][[1]] < 0.05)
    {
      sig_patterns <- c(sig_patterns, p)  
    }
  }

  plot.data <-  a5th.norm_tvr %>%  filter(Pattern %in% sig_patterns) %>%  mutate(PID=gsub("T0","",PID)) %>%  
    mutate(log2_ratio_count_norm_by_all_reads=log2(Count_norm_by_all_reads_T/Count_norm_by_all_reads_C)) %>%
    inner_join(A5_anno %>% 
                 dplyr::select(A5_ID, `Patient ID`, TERT_ATRX, TERT_ATRX_Event, `Assumed driver of metastasis`, `Tumour location simplified` ,is_primary_or_met, tumour_metastasised) %>% 
                 dplyr::rename(PID=A5_ID, 
                               location_simple=`Tumour location simplified`,
                               Primary_Met=is_primary_or_met, 
                               Met_Case=tumour_metastasised)) # %>% filter(Pattern %in% sig_patterns)
    ggplot(plot.data, aes(y=log2_ratio_count_norm_by_all_reads, x=Pattern, color=TERT_ATRX)) + 
    #geom_text(size=3) + 
    #geom_line(alpha=0.3) +   
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position=position_jitterdodge(jitter.width = 0.1), alpha=0.4) + 
      theme_bw()  + 
      theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5)) +
      coord_cartesian(ylim=c(-3,3)) + 
    scale_color_manual(values = c(ColorPalette[["DarkBlue1"]], 
                                  ColorPalette[["LightOrange1"]], 
                                  ColorPalette[["DarkGrey1"]])) + 
      ylab("Tumour/Normal Ratio (log2)") + xlab("")
      
  #facet_wrap("pattern", scales="free")
  
    


    a5th.singletons %>%  filter(!is.na(distance_to_expected_singleton_log2_ratio )) %>% 
       dplyr::select(PID, pattern, distance_to_expected_singleton_log2_ratio) %>%  mutate(PID=gsub("T0","",PID)) %>% inner_join(A5_anno %>% 
                                      dplyr::select(A5_ID, `Patient ID`, TERT_ATRX, TERT_ATRX_Event, `Assumed driver of metastasis`, `Tumour location simplified` ,is_primary_or_met, tumour_metastasised) %>% 
                                      dplyr::rename(PID=A5_ID, 
                                                    location_simple=`Tumour location simplified`,
                                                    Primary_Met=is_primary_or_met, 
                                                    Met_Case=tumour_metastasised)) %>% 
    
    
    ggplot(aes(y=distance_to_expected_singleton_log2_ratio, x=pattern, color=TERT_ATRX)) + 
      #geom_text(size=3) + 
      #geom_line(alpha=0.3) +   
      geom_boxplot(outlier.alpha = 0) +
      geom_point(position=position_jitterdodge(jitter.width = 0.1), alpha=0.4) + 
      theme_bw()  + 
      theme(axis.text.x = element_text(angle=90, hjust=1,vjust=0.5)) +
      coord_cartesian(ylim=c(-3,3)) + 
      scale_color_manual(values = c(ColorPalette[["DarkBlue1"]], 
                                    ColorPalette[["LightOrange1"]], 
                                    ColorPalette[["LightGreen1"]])) 
    
    summary_sheet.tvr <- a5th.norm_tvr %>%  filter(Pattern %in% sig_patterns) %>%  
      mutate(log2_ratio_count_norm_by_all_reads=log2(Count_norm_by_all_reads_T/Count_norm_by_all_reads_C)) %>% 
      select(PID, Pattern, log2_ratio_count_norm_by_intratel_reads, log2_ratio_count_norm_by_all_reads) %>%  
      pivot_longer(cols = c(-PID, -Pattern)) %>%  
      pivot_wider(id_cols = PID, names_from=c(Pattern, name), values_from=value, names_sep="_")

    summary_sheet <- A5_anno %>% mutate(PID=gsub("-","-T0",A5_ID)) %>% 
      select(PID, is_primary_or_met, 
             tumour_metastasised, `Assumed driver of metastasis`,
             `is_head_and_neck`) %>% 
      full_join(a5th.summary.ratio %>% 
      select(-sample, -total_reads, -read_length, -repeat_threshold_set, 
             -repeat_threshold_used, -tel_reads, -intratel_reads, 
             -gc_bins_for_correction, -total_reads_with_tel_gc) %>% 
      rename(dna_tel_content_tumour_control_log2ratio=tel_content)) %>% 
      full_join(a5th.summary.rna %>% 
                  select(PID, tel_content) %>% 
                  rename(terra_rna_tel_content=tel_content) %>% 
                  mutate(PID=gsub("-","-T0",PID))) %>% 
      full_join(a5th.singletons.dist) %>% full_join(summary_sheet.tvr) %>% 
      relocate(terra_rna_tel_content, .after=dna_tel_content_tumour_control_log2ratio) %>% 
      full_join(a5th.telosv %>% rename(PID=A5_ID)) 
    summary_sheet
write_sheet(summary_sheet %>% mutate(across(.cols = everything(),.fns=function(x){ ifelse(is.infinite(x), NA, x) })), ss = "15aaq-0BD3ZVcZDqK2SwLmb5YntWO42xWIRmIoayl8HI", sheet = "Combined_summary")     
