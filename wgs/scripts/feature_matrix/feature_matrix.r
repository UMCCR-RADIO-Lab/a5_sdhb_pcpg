library(googlesheets4)
source("C:/Users/AFFLY/OneDrive - The University of Melbourne/ResearchData/RADIO/A5/Scripts/dataloaders/load_All_Data.R")

feature_sheet <- as_sheets_id("1IOczHu91crprHvjxTKDRpqLfpeuYLaGwNXZVsWmXIKY")

Xgenes <- pcawg_drivers_hgnc_syn  %>% filter(chromosome_name=="chrX") %>% pull(Symbol) %>% sort()

##########################
# Feature Matrix Overall #
##########################

BiAllelicInactivation.nonX <- bind_rows(
  A5_gene_cn.keep.pcawg.brief %>% ungroup() %>% filter(!(Gene %in% Xgenes),grepl("CNDel\\(|CNGeneBreak\\([01]", Annotation)) %>% inner_join(A5_vcfs.keep.pcawg.brief.anno %>% filter(A5_ID != "E167-1"), by=c("A5_ID", "Gene")),
  A5_gene_cn.keep.pcawg.brief %>% ungroup() %>% filter(!(Gene %in% Xgenes),grepl("CNDel\\(|CNGeneBreak\\([01]", Annotation)) %>% inner_join(A5_germline_pcgr.keep.pcawg.brief %>% filter(grepl("Pathogenic",Annotation)), by=c("A5_ID", "Gene")),
  A5_gene_cn.keep.pcawg.brief %>% ungroup() %>% filter(!(Gene %in% Xgenes), grepl("HomoDel\\(", Annotation)),
  A5_gridss.keep.pcawg.brief %>% ungroup() %>% filter(!(Gene %in% Xgenes), grepl("SVInterrupt", Annotation))  %>% inner_join(A5_vcfs.keep.pcawg.brief.anno %>% filter(A5_ID != "E167-1"), by=c("A5_ID", "Gene"))
) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID) %>%  distinct() %>% group_by(Gene) %>% dplyr::count(name="BiAllelicInactivation") %>% 
  mutate(BiAllelicInactivation=as.character(BiAllelicInactivation))

BiAllelicInactivation.X <- bind_rows(
  A5_vcfs.keep.pcawg.brief %>% filter(Gene %in% Xgenes, A5_ID != "E167-1"),
  A5_germline_pcgr.keep.pcawg.brief %>% filter(Gene %in% Xgenes, grepl("Pathogenic",Annotation)),
  A5_gene_cn.keep.pcawg.brief %>% ungroup() %>% filter(Gene %in% Xgenes, grepl("HomoDel\\(", Annotation)),
  A5_gridss.keep.pcawg.brief %>% ungroup() %>% filter(Gene %in% Xgenes, grepl("SVInterrupt", Annotation))
) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% left_join(Anno %>% dplyr::select(A5_ID, Gender)) %>% 
  ungroup() %>% dplyr::select(Gene, Gender, A5_PatientID) %>%  distinct() %>% 
  group_by(Gene,Gender) %>% 
  dplyr::count(name="BiAllelicInactivation") %>% pivot_wider(names_from="Gender", values_from='BiAllelicInactivation') %>% 
  mutate(across(.cols = c(male,female), replace_na, replace = 0)) %>% 
  mutate(BiAllelicInactivation=ifelse(female>0,paste0(as.character(male),"(Female:+",female,")"), as.character(male))) %>% 
  dplyr::select(-male,-female)

BiAllelicInactivation <- bind_rows(BiAllelicInactivation.nonX, BiAllelicInactivation.X)

Mutated <- A5_vcfs.keep.pcawg.brief %>% ungroup() %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="Mutated")

GermlineMutated <- A5_germline_pcgr.keep.pcawg.brief %>% ungroup() %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="GermlineMutated")

GermlineMutated.PathOnly <- A5_germline_pcgr.keep.pcawg.brief %>% filter(grepl("Pathogenic",Annotation)) %>% ungroup() %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="GermlinePathogenic")

OverExpressed <- A5_zexp.keep.pcawg.brief %>% ungroup() %>% filter(grepl("High-Expression", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="ExprOutlier_Up")
UnderExpressed <- A5_zexp.keep.pcawg.brief %>% ungroup() %>% filter(grepl("Low-Expression", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="ExprOutlier_Down")

OverMeth <- A5_zMeth.keep.pcawg.brief %>% ungroup() %>% filter(grepl("High-Methylation", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="MethOutlier_Up")
UnderMeth <- A5_zMeth.keep.pcawg.brief %>% ungroup() %>% filter(grepl("Low-Methylation", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="MethOutlier_Down")


HomoDel <- A5_gene_cn.keep.pcawg.brief %>% ungroup() %>% filter(grepl("HomoDel\\(", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="HomoDel")
CNGain <- A5_gene_cn.keep.pcawg.brief %>% ungroup() %>% filter(grepl("CNGain\\(", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="CNGain")
CNDelOneCopyLeft <- A5_gene_cn.keep.pcawg.brief %>% filter(grepl("CNDel\\([01]", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="CNDel_OneCopyLeft")
CNDelOneCopyLeft_Focal <- A5_gene_cn.keep.pcawg.brief %>% ungroup() %>% filter(grepl("CNDelFocal\\([01]", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="CNDelFocal_OneCopyLeft")

SVInterrupt <- A5_gridss.keep.pcawg.brief %>% ungroup()  %>% filter(grepl("SVInterrupt", Annotation))  %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>% group_by(Gene) %>% dplyr::count(name="SVInterrupt")
NearSV <- A5_gridss.keep.pcawg.brief %>% ungroup() %>% filter(grepl("NearSV", Annotation))  %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>% group_by(Gene) %>% dplyr::count(name="NearSV")

RNAFusion <- A5_arriba.keep.pcawg.brief %>% ungroup() %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID) %>% distinct() %>% group_by(Gene) %>% dplyr::count(name="RNA_Fusion")


AllFeatures <- bind_rows(
  list(A5_vcfs.keep.pcawg.brief,
       A5_germline_pcgr.keep.pcawg.brief %>% filter(grepl("Pathogenic",Annotation)),
       A5_zexp.keep.pcawg.brief,
       A5_zMeth.keep.pcawg.brief,
       A5_gridss.keep.pcawg.brief, 
       A5_gene_cn.keep.pcawg.brief %>% filter(!grepl("X/Y-Male", Annotation)), 
       A5_arriba.keep.pcawg.brief
  )
) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  group_by(A5_ID, Gene) %>% distinct() %>% 
  summarise(Annotation=paste(Annotation, sep = ";", collapse = ";")) %>% filter(!grepl("^CNDel\\(.+\\)$",Annotation))

AffectedAny <- AllFeatures %>% mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(A5_ID, Gene) %>% distinct()  %>% group_by(Gene) %>% dplyr::count(name="AnyAlteration")

SummaryMatrix <- purrr::reduce(.x = list(Mutated, OverExpressed, UnderExpressed, UnderMeth, OverMeth, HomoDel, CNDelOneCopyLeft, CNDelOneCopyLeft_Focal, CNGain, SVInterrupt, NearSV, RNAFusion, BiAllelicInactivation, GermlineMutated, GermlineMutated.PathOnly, AffectedAny), .f = full_join)

SummaryMatrix <- SummaryMatrix %>% mutate(across(2:ncol(SummaryMatrix)-1, replace_na, replace = 0))

# Annotation_Header <- Anno %>% dplyr::select(`A5_ID`, `Assumed driver of metastais`, `is_head_and_neck`, 
#                                             tumour_metastasised, is_primary_or_met) %>% 
#   pivot_longer(cols=-A5_ID, names_to="Gene", values_to="ClinicalFeature") %>% pivot_wider(id_cols=Gene, names_from=A5_ID, values_from=ClinicalFeature)

FeatureMatrix <- SummaryMatrix %>% full_join(
  AllFeatures  %>% filter(A5_ID %in% Anno$A5_ID) %>% 
    pivot_wider(id_cols = Gene,names_from=A5_ID, values_from=Annotation) %>% mutate(across(starts_with("E"),replace_na, replace="")), by="Gene") %>% 
  filter(!is.na(Gene), Gene %in% pcawg_drivers_hgnc_syn$Symbol) %>% 
  left_join(pcawg_drivers_hgnc_syn %>% dplyr::select(Symbol, chromosome_name) %>% dplyr::rename(chromosome=chromosome_name, Gene=Symbol)) %>% relocate(chromosome, .after=Gene) 

# SampleOrder <- Anno %>% arrange(`Assumed driver of metastais`, desc(tumour_metastasised),is_primary_or_met) %>% pull(`A5_ID`)

# FeatureMatrix <- bind_rows(Annotation_Header, FeatureMatrix) %>% relocate(!matches("E[0-9]{3}-[0-9]"))
# FeatureMatrix <- FeatureMatrix[,c(colnames(FeatureMatrix)[1:13],SampleOrder)]

FeatureMatrix_t <- AllFeatures  %>% ungroup() %>% filter(!is.na(Gene), Gene %in% pcawg_drivers_hgnc_syn$Symbol) %>% 
  pivot_wider(id_cols = A5_ID,names_from=Gene, values_from=Annotation) %>% mutate(across(2:ncol(.)-1,replace_na, replace="")) %>% 
  left_join(Anno %>% dplyr::select(`A5_ID`, `Assumed driver of metastais`, `is_head_and_neck`, 
                                   tumour_metastasised, is_primary_or_met, Gender)) %>% 
  relocate(c(`Assumed driver of metastais`, `is_head_and_neck`, 
             tumour_metastasised, is_primary_or_met, Gender),.after=A5_ID)
FeatureMatrix_t <- FeatureMatrix_t[,c(colnames(FeatureMatrix_t)[1:6],sort(colnames(FeatureMatrix_t)[7:ncol(FeatureMatrix_t)]))]  
FeatureMatrix_t <- FeatureMatrix_t %>% filter(!is.na(`Assumed driver of metastais`))

# write.table(FeatureMatrix, "C:/Users/AFFLY/OneDrive - The University of Melbourne/ResearchData/RADIO/A5/Tables/FeatureSummary_GeneCentric.tsv", sep="\t", row.names=F, quote=F)
# write.table(FeatureMatrix_t, "C:/Users/AFFLY/OneDrive - The University of Melbourne/ResearchData/RADIO/A5/Tables/FeatureSummary_SampleCentric.tsv", sep="\t", row.names=F, quote=F)

write_sheet(FeatureMatrix, ss=feature_sheet, sheet = "GeneCentric")
write_sheet(FeatureMatrix_t, ss=feature_sheet, sheet = "SampleCentric")

################################
# Feature Matrix By Annotation #
################################

Anno.subset <- Anno %>% dplyr::select(`A5_ID`, `Assumed driver of metastais`, `is_head_and_neck`, 
                                      tumour_metastasised, is_primary_or_met)
Anno.subset <- Anno.subset %>% mutate(`Assumed driver of metastais`=case_when(
  tumour_metastasised == "Yes" & `Assumed driver of metastais` == "Unknown" ~ "Met - No TERT/ATRX",
  tumour_metastasised == "No" & `Assumed driver of metastais` == "Unknown" ~ "Benign",
  tumour_metastasised == "Short follow up" & `Assumed driver of metastais` == "Unknown" ~ "SFU",
  tumour_metastasised == "Short follow up" & `Assumed driver of metastais` != "Unknown" ~ paste("SFU", `Assumed driver of metastais`),
  tumour_metastasised == "No" & `Assumed driver of metastais` %in% c("TERT", "ATRX") ~ paste("Benign", `Assumed driver of metastais`),
  TRUE ~ `Assumed driver of metastais`
))

sheetfriendlynames <-  setNames(object = c("SummaryCounts_MetastaticDriver", "SummaryCounts_PrimaryMet", "SummaryCounts_MetastaticYesNo"),nm = 
                                  c("Assumed driver of metastais", "Tumour metastasised", "is_primary_or_met"))

A5_vcfs.keep.pcawg.brief.anno <- A5_vcfs.keep.pcawg.brief %>% left_join(Anno.subset) 
A5_germline_pcgr.keep.pcawg.brief.anno <- A5_germline_pcgr.keep.pcawg.brief %>% left_join(Anno.subset) 
A5_zexp.keep.pcawg.brief.anno <- A5_zexp.keep.pcawg.brief %>% left_join(Anno.subset)
A5_zMeth.keep.pcawg.brief.anno <- A5_zMeth.keep.pcawg.brief %>% left_join(Anno.subset)
A5_gridss.keep.pcawg.brief.anno <- A5_gridss.keep.pcawg.brief %>% left_join(Anno.subset) 
A5_gene_cn.keep.pcawg.brief.anno <- A5_gene_cn.keep.pcawg.brief %>% left_join(Anno.subset) 
A5_arriba.keep.pcawg.brief.anno <- A5_arriba.keep.pcawg.brief %>% left_join(Anno.subset)

BiAllelicInactivation <- list()
Mutated <- list()
OverExpressed <- list()
UnderExpressed <- list()
OverMeth <- list()
UnderMeth <- list()
HomoDel <- list()
CNGain <- list()
CNDelOneCopyLeft <- list()
CNDelOneCopyLeft_Focal <- list()
SVInterrupt <- list()
NearSV <- list()
RNAFusion <- list()
AffectedAny <- list()
SummaryMatrix <- list()
AllFeatures <- list()
for (AnnoCol in c("Assumed driver of metastais","Tumour metastasised", "is_primary_or_met"))
{

BiAllelicInactivation.nonX <- bind_rows(
  A5_gene_cn.keep.pcawg.brief.anno %>% ungroup() %>% filter(!(Gene %in% Xgenes),grepl("CNDel\\(|CNGeneBreak\\([01]", Annotation)) %>% inner_join(A5_vcfs.keep.pcawg.brief.anno %>% filter(A5_ID != "E167-1"), by=c("A5_ID", "Gene", AnnoCol)),
  A5_gene_cn.keep.pcawg.brief.anno %>% ungroup() %>% filter(!(Gene %in% Xgenes),grepl("CNDel\\(|CNGeneBreak\\([01]", Annotation)) %>% inner_join(A5_germline_pcgr.keep.pcawg.brief.anno %>% filter(grepl("Pathogenic",Annotation)), by=c("A5_ID", "Gene")),
  A5_gene_cn.keep.pcawg.brief.anno %>% ungroup() %>% filter(!(Gene %in% Xgenes), grepl("HomoDel\\(", Annotation)),
  A5_gridss.keep.pcawg.brief.anno %>% ungroup() %>% filter(!(Gene %in% Xgenes), grepl("SVInterrupt", Annotation))  %>% inner_join(A5_vcfs.keep.pcawg.brief.anno %>% filter(A5_ID != "E167-1"), by=c("A5_ID", "Gene", AnnoCol))
) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>% group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="BiAllelicInactivation") %>% 
  mutate(BiAllelicInactivation=as.character(BiAllelicInactivation))

BiAllelicInactivation.X <- bind_rows(
  A5_vcfs.keep.pcawg.brief.anno %>% filter(Gene %in% Xgenes, A5_ID != "E167-1"),
  A5_germline_pcgr.keep.pcawg.brief.anno %>% filter(Gene %in% Xgenes, grepl("Pathogenic",Annotation)),
  A5_gene_cn.keep.pcawg.brief.anno %>% ungroup() %>% filter(Gene %in% Xgenes, grepl("HomoDel\\(", Annotation)),
  A5_gridss.keep.pcawg.brief.anno %>% ungroup() %>% filter(Gene %in% Xgenes, grepl("SVInterrupt", Annotation))
) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% left_join(Anno %>% dplyr::select(A5_ID, Gender)) %>% 
  ungroup() %>% dplyr::select(Gene, Gender, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>% 
  group_by(Gene,Gender, !!rlang::sym(AnnoCol)) %>% 
  dplyr::count(name="BiAllelicInactivation") %>% pivot_wider(names_from="Gender", values_from='BiAllelicInactivation') %>% 
  mutate(across(.cols = c(male,female), replace_na, replace = 0)) %>% 
  mutate(BiAllelicInactivation=ifelse(female>0,paste0(as.character(male),"(Female:+",female,")"), as.character(male))) %>% 
  dplyr::select(-male,-female)


BiAllelicInactivation[[AnnoCol]] <- bind_rows(BiAllelicInactivation.nonX, BiAllelicInactivation.X)

Mutated[[AnnoCol]] <- A5_vcfs.keep.pcawg.brief.anno %>% ungroup() %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>%  group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="Mutated")

OverExpressed[[AnnoCol]] <- A5_zexp.keep.pcawg.brief.anno %>% ungroup() %>% filter(grepl("High-Expression", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct()%>%  group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="ExprOutlier_Up")
UnderExpressed[[AnnoCol]] <- A5_zexp.keep.pcawg.brief.anno %>% ungroup() %>% filter(grepl("Low-Expression", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>%  group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="ExprOutlier_Down")

OverMeth[[AnnoCol]] <- A5_zMeth.keep.pcawg.brief.anno %>% ungroup() %>% filter(grepl("High-Methylation", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>%  group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="MethOutlier_Up")
UnderMeth[[AnnoCol]] <- A5_zMeth.keep.pcawg.brief.anno %>% ungroup() %>% filter(grepl("Low-Methylation", Annotation))%>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>% distinct() %>%  group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="MethOutlier_Down")


HomoDel[[AnnoCol]] <- A5_gene_cn.keep.pcawg.brief.anno %>% ungroup() %>% filter(grepl("HomoDel\\(", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>%  group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="HomoDel")
CNGain[[AnnoCol]] <- A5_gene_cn.keep.pcawg.brief.anno %>% ungroup() %>% filter(grepl("CNGain\\(", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>%  group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="CNGain")
CNDelOneCopyLeft[[AnnoCol]] <- A5_gene_cn.keep.pcawg.brief.anno %>% ungroup() %>% filter(grepl("CNDel\\([01]", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>%  group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="CNDel_OneCopyLeft")
CNDelOneCopyLeft_Focal[[AnnoCol]] <- A5_gene_cn.keep.pcawg.brief %>% ungroup() %>% filter(grepl("CNDelFocal\\([01]", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID)  %>%  distinct() %>%  group_by(Gene) %>% dplyr::count(name="CNDelFocal_OneCopyLeft")

SVInterrupt[[AnnoCol]] <- A5_gridss.keep.pcawg.brief.anno %>% ungroup() %>% filter(grepl("SVInterrupt", Annotation)) %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>% group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="SVInterrupt")
NearSV[[AnnoCol]] <- A5_gridss.keep.pcawg.brief.anno %>% ungroup() %>% filter(grepl("NearSV", Annotation))  %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>% group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="NearSV")

RNAFusion[[AnnoCol]] <- A5_arriba.keep.pcawg.brief.anno %>% ungroup() %>% 
  mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% 
  dplyr::select(Gene, A5_PatientID, all_of(AnnoCol)) %>%  distinct() %>% group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="RNA_Fusion")



AllFeatures[[AnnoCol]] <- bind_rows(
  list(A5_vcfs.keep.pcawg.brief.anno, 
       A5_zexp.keep.pcawg.brief.anno,
       A5_zMeth.keep.pcawg.brief.anno,
       A5_gridss.keep.pcawg.brief.anno, 
       A5_gene_cn.keep.pcawg.brief.anno %>% filter(!grepl("X/Y-Male", Annotation)), 
       A5_arriba.keep.pcawg.brief.anno
  )
) %>% 
  group_by(A5_ID, Gene, !!rlang::sym(AnnoCol)) %>% 
  summarise(Annotation=paste(Annotation, sep = ";", collapse = ";")) %>% filter(!grepl("^CNDel\\(.+\\)$",Annotation))

AffectedAny[[AnnoCol]] <- AllFeatures[[AnnoCol]] %>% ungroup() %>% mutate(A5_PatientID=gsub("-.$","",A5_ID)) %>% dplyr::select(A5_PatientID, Gene, !!rlang::sym(AnnoCol)) %>% distinct()  %>% group_by(Gene, !!rlang::sym(AnnoCol)) %>% dplyr::count(name="AnyAlteration")

SummaryMatrix[[AnnoCol]] <- purrr::reduce(.x = list(Mutated[[AnnoCol]], OverExpressed[[AnnoCol]], UnderExpressed[[AnnoCol]], 
                                         UnderMeth[[AnnoCol]], OverMeth[[AnnoCol]], HomoDel[[AnnoCol]], 
                                         CNDelOneCopyLeft[[AnnoCol]], CNGain[[AnnoCol]], SVInterrupt[[AnnoCol]], 
                                         NearSV[[AnnoCol]], RNAFusion[[AnnoCol]], BiAllelicInactivation[[AnnoCol]], 
                                         AffectedAny[[AnnoCol]]), .f = full_join)

extract_XBIA <- function (xinfo, nSamples) {
  returnvector <- vector(length = length(xinfo))
  for (i in 1:length(xinfo))
  {  
    if(!grepl("Female", xinfo[[i]]))
    {
      returnvector[[i]] <- paste0(round((as.numeric(xinfo[[i]])/nSamples[[i]])*100,1),"%(", xinfo[[i]], ")")
    } else {
      MaleCount <- as.numeric(gsub("([0-9]+)\\(Female:\\+([0-9]+)\\)","\\1",xinfo[[i]]))
      FemaleCount <- as.numeric(gsub("([0-9]+)\\(Female:\\+([0-9]+)\\)","\\2",xinfo[[i]]))
      returnvector[[i]] <- paste0(round((MaleCount/nSamples[[i]])*100,1),"%(",MaleCount,")(Female:+", round((FemaleCount/nSamples[[i]])*100,1),"%(",FemaleCount,"))")
    }
  }
  return(returnvector)
}

SummaryMatrix[[AnnoCol]] <- SummaryMatrix[[AnnoCol]] %>% mutate(across(3:ncol(SummaryMatrix[[AnnoCol]])-2, replace_na, replace = 0)) %>% left_join(Anno.subset %>% group_by(!!rlang::sym(AnnoCol)) %>% dplyr::count(name="nSamples")) %>% 
  ungroup() %>% 
  mutate(across(.cols = c(-Gene,-BiAllelicInactivation,-!!rlang::sym(AnnoCol), -nSamples), 
                ~ paste0(round((.x/nSamples)*100,1),"%(", .x, ")"))) %>% 
  mutate(BiAllelicInactivation=extract_XBIA(BiAllelicInactivation,nSamples)) %>% dplyr::select(-nSamples) 


# write.table(SummaryMatrix[[AnnoCol]], paste("C:/Users/AFFLY/OneDrive - The University of Melbourne/ResearchData/RADIO/A5/Tables/SummaryMatrix_",AnnoCol,".tsv"), sep="\t", row.names=F, quote=F)

write_sheet(SummaryMatrix[[AnnoCol]], ss=feature_sheet, sheet = sheetfriendlynames[AnnoCol])

}

for (AnnoCol in c("Assumed driver of metastais", "is_head_and_neck","Tumour metastasised", "is_primary_or_met"))
{
  
write.table(SummaryMatrix[[AnnoCol]], paste("C:/Users/AFFLY/OneDrive - The University of Melbourne/ResearchData/RADIO/A5/Tables/SummaryMatrix_",AnnoCol,".tsv"), sep="\t", row.names=F, quote=F)
}

