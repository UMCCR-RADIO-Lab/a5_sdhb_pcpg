library(StructuralVariantAnnotation)
library(VariantAnnotation)
library(dplyr)
library(tidyr)


merge_and_fill <- function(compA, compB) {
  compA %>% full_join(compB) %>% ungroup() %>% 
    mutate(across(.cols=ends_with("_Called"), .fns=~replace_na(.x,"N"))) %>% 
    mutate(across(.cols=ends_with("VAF"), .fns = ~replace_na(.x,0)))
}

make_pairs <- function(input)
{
  return_frame=data.frame(Sample1=vector(mode='character'),Sample2=vector('character'))
  for (i in 1:length(input))
  {
    for (j in i+1:length(input))
    {
      if (j>length(input)) {next}
      
      a = input[i]
      b = input[j]
      return_frame <- rbind(return_frame, data.frame(Sample1=a,Sample2=b))
    }
  }
  return(return_frame)
}

compare_small_variant_vcfs <- function (vcf1, vcf2, normal_name_vcf1, tumour_name_vcf1, normal_name_vcf2, tumour_name_vcf2, genome, collapse_annotation=F, Keep_Only_Highest_Consequence=F, Keep_Only_One_Annotation=F, drop_info=F) {
  
  if(!((class(vcf1)=="character" | class(vcf1)=="CollapsedVCF") & 
       (class(vcf2)=="character" | class(vcf2)=="CollapsedVCF")))
  {
    stop("VCF1 and VCF2 must be supplied as a path to a VCF file or a Collapsed VCF loaded by the Variant Annotation package")  
  }
  
  if(class(vcf1)=="character")
  {
    message("Reading VCF1...",vcf1)
    vcf1 <- readVcf(vcf1, genome = genome)
  } 
  
  if(class(vcf2)=="character")
  {
    message("Reading VCF2:",vcf2)
    vcf2 <- readVcf(vcf2, genome = genome)
  } 
  
  rownames(colData(vcf1))[which(rownames(colData(vcf1))==normal_name_vcf1)] <- "Normal"
  rownames(colData(vcf1))[which(rownames(colData(vcf1))==tumour_name_vcf1)] <- "Tumour"
  rownames(colData(vcf2))[which(rownames(colData(vcf2))==normal_name_vcf2)] <- "Normal"
  rownames(colData(vcf2))[which(rownames(colData(vcf2))==tumour_name_vcf2)] <- "Tumour"
  
  message("Finding shared/unique variants ...")
  common <- GenomicRanges::findOverlaps(granges(vcf1), granges(vcf2))
  unique_to_vcf1 <- setdiff(1:length(vcf1), queryHits(common))
  unique_to_vcf2 <- setdiff(1:length(vcf2), subjectHits(common))
  
  vcf1_to_append <- vcf1[unique_to_vcf1,]
  
  all <- VariantAnnotation::rbind(vcf2, vcf1_to_append)
  
  in_vcf1 <- GenomicRanges::findOverlaps(granges(vcf1), granges(all))
  in_vcf2 <- GenomicRanges::findOverlaps(granges(vcf2), granges(all))
  
  mcols(all)[[paste0(tumour_name_vcf1,"_Called")]] <- "N"
  mcols(all)[[paste0(tumour_name_vcf2,"_Called")]] <- "N"
  
  mcols(all)[[paste0(tumour_name_vcf1,"_Called")]][subjectHits(in_vcf1)] <- "Y"
  mcols(all)[[paste0(tumour_name_vcf2,"_Called")]][subjectHits(in_vcf2)] <- "Y"
  
  message("Annotating VAFs ...")
  mcols(all)[[paste0(normal_name_vcf1,"_VAF")]] <- 0
  vcf1_nvafs <- geno(vcf1)[["AF"]][queryHits(in_vcf1),"Normal"]
  vcf1_nvafs[lapply(vcf1_nvafs, length)==0] <- 0 #fill missing VAF values
  mcols(all)[[paste0(normal_name_vcf1,"_VAF")]][subjectHits(in_vcf1)] <- unlist(lapply(vcf1_nvafs,"[[",1)) #fetch first VAF in rare case of multi-allele
  
  if(normal_name_vcf1 != normal_name_vcf2)
  {
    mcols(all)[[paste0(normal_name_vcf2,"_VAF")]] <- 0
  }  
  
  vcf2_nvafs <- geno(vcf2)[["AF"]][queryHits(in_vcf2),"Normal"]
  vcf2_nvafs[lapply(vcf2_nvafs, length)==0] <- 0 #fill missing VAF values
  mcols(all)[[paste0(normal_name_vcf2,"_VAF")]][subjectHits(in_vcf2)] <- unlist(lapply(vcf2_nvafs,"[[",1)) #fetch first VAF in rare case of multi-allele
  
  mcols(all)[[paste0(tumour_name_vcf1,"_VAF")]] <- 0
  vcf1_vafs <- geno(vcf1)[["AF"]][queryHits(in_vcf1),"Tumour"]
  vcf1_vafs[lapply(vcf1_vafs, length)==0] <- 0 #fill missing VAF values
  mcols(all)[[paste0(tumour_name_vcf1,"_VAF")]][subjectHits(in_vcf1)] <- unlist(lapply(vcf1_vafs,"[[",1)) #fetch first VAF in rare case of multi-allele
  
  mcols(all)[[paste0(tumour_name_vcf2,"_VAF")]] <- 0
  vcf2_vafs <- geno(vcf2)[["AF"]][queryHits(in_vcf2),"Tumour"]
  vcf2_vafs[lapply(vcf2_vafs, length)==0] <- 0 
  mcols(all)[[paste0(tumour_name_vcf2,"_VAF")]][subjectHits(in_vcf2)] <- unlist(lapply(vcf2_vafs,"[[",1)) 
  
  if(drop_info) {
    message("Dropping VCF INFO and returning annotated ranges ... ")
    info(all)[names(info(all))] <- NULL
    all.df <- as.data.frame(granges(all),row.names = NULL)
    all.df$ALT <- unstrsplit(CharacterList(all.df$ALT), sep = ",")
    
    return(all.df)
  }
  
  message("Processing annotation ... ")
  max_entries <- 30000
  if (length(all) > max_entries){
    message("Combined VCF has ", length(all), " rows which greater than the cutoff of ", max_entries,". Annotation handling will be skipped.")
    all.df <- data.frame(seqnames= seqnames(all), 
                         start = start(all), 
                         end = end(all), 
                         REF = mcols(all)[["REF"]], 
                         ALT = unstrsplit(CharacterList(mcols(all)[["ALT"]]), sep = ","))
    all.df[[paste0(tumour_name_vcf1,"_Called")]] <- mcols(all)[[paste0(tumour_name_vcf1,"_Called")]]
    all.df[[paste0(tumour_name_vcf2,"_Called")]] <- mcols(all)[[paste0(tumour_name_vcf2,"_Called")]]
    all.df[[paste0(tumour_name_vcf1,"_VAF")]] <- mcols(all)[[paste0(tumour_name_vcf1,"_VAF")]]
    all.df[[paste0(tumour_name_vcf2,"_VAF")]] <- mcols(all)[[paste0(tumour_name_vcf2,"_VAF")]]
    all.df[[paste0(normal_name_vcf1,"_VAF")]] <- mcols(all)[[paste0(normal_name_vcf1,"_VAF")]]
    
    if("ANN" %in% names(info(all)))
    {
      all.df[["Annotation"]] <- unstrsplit(info(all)[["ANN"]], sep = ",")
    }
    
    if("CSQ" %in% names(info(all)))
    {
      all.df[["Annotation"]] <- unstrsplit(info(all)[["CSQ"]], sep = ",")
    }
    
  } else {
    #HMF-5 SAGE/PURPLE VCF
    if("ANN" %in% names(info(all)))
    {
      mcols(all)[["Annotation"]] <- unlist(lapply(info(all)[["ANN"]],paste, collapse=";"))
      AnnoColNames <- c("Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID","Feature_Type",
                        "Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos_length",
                        "CDS.pos_length","AA.pos_length","Distance","ERRORS_WARNINGS_INFO")
      impact_col="Annotation_Impact"
    }
    
    #UMCCRISED VCF
    if("CSQ" %in% names(info(all)))
    {
      mcols(all)[["Annotation"]] <- unlist(lapply(info(all)[["CSQ"]],paste, collapse=";"))
      AnnoColNames <- c("Allele","Consequence","IMPACT","SYMBOL","Gene","Feature_type","Feature",
                        "BIOTYPE","EXON","INTRON","HGVSc","HGVSp","cDNA_position","CDS_position",
                        "Protein_position","Amino_acids","Codons","Existing_variation", "ALLELE_NUM",
                        "DISTANCE","STRAND","FLAGS","PICK","VARIANT_CLASS","SYMBOL_SOURCE","HGNC_ID",
                        "CANONICAL","MANE","TSL","APPRIS","CCDS","ENSP","SWISSPROT","TREMBL","UNIPARC","RefSeq",
                        "DOMAINS","HGVS_OFFSET","AF","AFR_AF","AMR_AF","EAS_AF","EUR_AF","SAS_AF",
                        "gnomAD_AF","gnomAD_AFR_AF","gnomAD_AMR_AF","gnomAD_ASJ_AF","gnomAD_EAS_AF",
                        "gnomAD_FIN_AF","gnomAD_NFE_AF","gnomAD_OTH_AF","gnomAD_SAS_AF",
                        "CLIN_SIG","SOMATIC","PHENO","CHECK_REF","NearestExonJB")
      PCGR_col_idx <- grep("PCGR",names(info(all)))
      PCGR_col_names <- names(info(all))[PCGR_col_idx]
      
      mcols(all)[PCGR_col_names] <- info(all)[PCGR_col_idx]
      impact_col="IMPACT"
    }
    
    
    
    all.df <- as.data.frame(granges(all),row.names = NULL, check.names=F)
    #all.df$ALT <- unlist(lapply(lapply(all.df$ALT,as.character),"paste",sep=",", collapse=",")) #coverts alt allele from array of DNAstring to character (slow)
    all.df$ALT <- unstrsplit(CharacterList(all.df$ALT), sep = ",")
    all.df <- all.df %>% tidyr::separate_rows(Annotation, sep=";")
    all.df <- all.df %>% tidyr::separate(col = Annotation, sep="\\|", into=AnnoColNames)
    
    if ("QUAL" %in% colnames(all.df)) 
    {
      all.df <- all.df %>% dplyr::select(-QUAL)
    }
    
    if ("FILTER" %in% colnames(all.df)) 
    {
      all.df <- all.df %>% dplyr::select(-FILTER)
    }
    
    all.df <- all.df %>% dplyr::distinct() #Handle duplicate variant entries
    
    all.df <- all.df %>% mutate("{impact_col}":=factor(!!sym(impact_col), levels=c("MODIFIER", "LOW", "MODERATE", "HIGH" ))) 
    
    
    if(Keep_Only_Highest_Consequence | Keep_Only_One_Annotation)
    {
      all.df <- all.df %>% group_by(seqnames, start, end, Allele) %>% 
        arrange(!!sym(impact_col)) %>% 
        {
          if(Keep_Only_One_Annotation)
          {
            mutate(.data=., keep=ifelse(((as.numeric(!!sym(impact_col)) == max(as.numeric(!!sym(impact_col)))) & (row_number()==1)), T, F))  #Keep if MODIFIER and first entry
          } else if(Keep_Only_Highest_Consequence){
            mutate(.data=., keep=(as.numeric(!!sym(impact_col)) == max(as.numeric(!!sym(impact_col))))) #Keep only annotations with the highest IMPACT for variant
          }  
        } %>%  
        relocate(keep, .after = !!sym(impact_col)) %>% 
        filter(keep) %>% dplyr::select(-keep)
    }
    
    if(!Keep_Only_One_Annotation){
      options(dplyr.summarise.inform = FALSE)
      all.df <- all.df %>% group_by(across(setdiff(colnames(all.df),AnnoColNames))) %>% 
        { 
          if(collapse_annotation)
          {
            summarise(.data = .,across(all_of(AnnoColNames),~paste(unique(.x),sep=";", collapse=";")))
          } else {
            summarise(.data = .,across(all_of(AnnoColNames),~paste(.x,sep=";", collapse=";")))
          }
        }
      options(dplyr.summarise.inform = TRUE)
    }
  }
  
  return(all.df)
}


multi_compare_small_variant_vcf <- function(vcfs, normal_names, tumour_names, genome,  collapse_annotation=F, Keep_Only_Highest_Consequence=F, Keep_Only_One_Annotation=F, drop_info=F)
{
  pairs <- make_pairs(tumour_names)
  
  pair_comps <- list()
  for (i in 1:nrow(pairs))
  {
    t1_idx = which(tumour_names == pairs$Sample1[i])
    t2_idx = which(tumour_names == pairs$Sample2[i])
    
    normal_name_vcf1 = normal_names[[t1_idx]]
    normal_name_vcf2 = normal_names[[t2_idx]]
    
    tumour_name_vcf1 = tumour_names[[t1_idx]]
    tumour_name_vcf2 = tumour_names[[t2_idx]]
    
    vcf1=vcfs[[t1_idx]]
    vcf2=vcfs[[t2_idx]]
    
    pair_comps[[paste(tumour_name_vcf1,tumour_name_vcf2,sep="_")]] <- compare_small_variant_vcfs(vcf1, 
                                                                                   vcf2, 
                                                                                   normal_name_vcf1 = normal_name_vcf1, 
                                                                                   tumour_name_vcf1 = tumour_name_vcf1, 
                                                                                   normal_name_vcf2 = normal_name_vcf2, 
                                                                                   tumour_name_vcf2 = tumour_name_vcf2, 
                                                                                   genome= genome, 
                                                                                   collapse_annotation = collapse_annotation,
                                                                                   Keep_Only_Highest_Consequence=Keep_Only_Highest_Consequence, 
                                                                                   Keep_Only_One_Annotation=Keep_Only_One_Annotation,
                                                                                   drop_info = drop_info)
  }
  
  if (nrow(pairs) == 1)
  {
    return(pair_comps[[1]])
  } else {
    return(purrr:::reduce(pair_comps, merge_and_fill) %>% 
      arrange(seqnames, start,end) %>%  
      relocate(ends_with("_Called"), ends_with("VAF"), .after = ALT)) 
  }
}
  
compare_sv_vcfs <- function (vcf1, vcf2, 
                             normal_name_vcf1=NULL, tumour_name_vcf1, 
                             tumour_name_vcf1_override=NULL, 
                             normal_name_vcf2=NULL, tumour_name_vcf2, 
                             tumour_name_vcf2_override=NULL, 
                             genome, maxgap=0, sizemargin=0,
                             info_filters) {
  
  if(!((class(vcf1)=="character" | class(vcf1)=="CollapsedVCF") & 
       (class(vcf2)=="character" | class(vcf2)=="CollapsedVCF")))
  {
    stop("VCF1 and VCF2 must be supplied as a path to a VCF file or a Collapsed VCF loaded by the Variant Annotation package")  
  }
  
  if(class(vcf1)=="character")
  {
    vcf1 <- readVcf(vcf1, genome = genome)
  } 
  
  if(class(vcf2)=="character")
  {
    vcf2 <- readVcf(vcf2, genome = genome)
  } 
  
  if (!is.null(normal_name_vcf1) && !(normal_name_vcf1 %in% rownames(colData(vcf1))))
  {
    stop(paste("Sample name", normal_name_vcf1, "not found in VCF supplied as vcf1"))
  }
  
  if (!(tumour_name_vcf1 %in% rownames(colData(vcf1))))
  {
    stop(paste("Sample name", tumour_name_vcf1, "not found in VCF supplied as vcf1"))
  }
  
  if (!is.null(normal_name_vcf2) && !(normal_name_vcf2 %in% rownames(colData(vcf2))))
  {
    stop(paste("Sample name", normal_name_vcf2, "not found in VCF supplied as vcf2"))
  }
  
  if (!(tumour_name_vcf2 %in% rownames(colData(vcf2))))
  {
    stop(paste("Sample name", tumour_name_vcf2, "not found in VCF supplied as vcf2"))
  }
  
  if (!is.null(normal_name_vcf1))
  {
    rownames(colData(vcf1))[which(rownames(colData(vcf1))==normal_name_vcf1)] <- "Normal"
  }
  rownames(colData(vcf1))[which(rownames(colData(vcf1))==tumour_name_vcf1)] <- "Tumour"
  
  if (!is.null(normal_name_vcf2))
  {
    rownames(colData(vcf2))[which(rownames(colData(vcf2))==normal_name_vcf2)] <- "Normal"
  }
  rownames(colData(vcf2))[which(rownames(colData(vcf2))==tumour_name_vcf2)] <- "Tumour"
  
  if(!is.null(tumour_name_vcf1_override)) { tumour_name_vcf1 = tumour_name_vcf1_override }
  if(!is.null(tumour_name_vcf2_override)) { tumour_name_vcf2 = tumour_name_vcf2_override }
  
  purple_info_fields= c("REFPAIR","RP","SR", "PURPLE_AF","PURPLE_CN",  "PURPLE_CN_CHANGE" )
  
 

 breakformat = list()
 breakformat[[tumour_name_vcf1]] <- list()
 breakformat[[tumour_name_vcf1]][["breakpoint"]] <- breakpointRanges(vcf1,info_columns=purple_info_fields)
 breakformat[[tumour_name_vcf1]][["breakend"]] <- breakendRanges(vcf1)
 
 breakformat[[tumour_name_vcf2]] <- list()
 breakformat[[tumour_name_vcf2]][["breakpoint"]] <- breakpointRanges(vcf2,info_columns=purple_info_fields)
 breakformat[[tumour_name_vcf2]][["breakend"]] <- breakendRanges(vcf2)

 all.breakpoint <- sv_combine(set_a = breakformat[[tumour_name_vcf1]][["breakpoint"]], 
            set_a_name = tumour_name_vcf1, 
            set_b = breakformat[[tumour_name_vcf2]][["breakpoint"]], 
            set_b_name = tumour_name_vcf2, 
            mode = "breakpoint", 
            maxgap= maxgap, 
            sizemargin = sizemargin, markcommon = T, remove_dup = T)
 
 all.breakend <- sv_combine(set_a = breakformat[[tumour_name_vcf1]][["breakend"]], 
                              set_a_name = tumour_name_vcf1, 
                              set_b = breakformat[[tumour_name_vcf2]][["breakend"]], 
                              set_b_name = tumour_name_vcf2, 
                              mode = "breakend", 
                              maxgap= maxgap, 
                              sizemargin = sizemargin, markcommon = T, remove_dup = T)
 
 all <- GenomicRanges::bindROWS(x=all.breakpoint, objects = list(all.breakend)) 
 
  return(GenomicRanges::as.data.frame(all, row.names = NULL))
}


mark_common_sv <- function(mode, query, subject, queryname, subjectname, maxgap, type="equal", sizemargin)
{
  mcols(query)[[paste0(queryname,"_Called")]] <- vector(length = length(query))
  mcols(subject)[[paste0(subjectname,"_Called")]] <- vector(length = length(subject))
  
  mcols(query)[[paste0(subjectname,"_Called")]] <- vector(length = length(query))
  mcols(subject)[[paste0(queryname,"_Called")]] <- vector(length = length(subject))
  
 
  
  if(length(query) > 0 & 
     length(subject) > 0)
  {
    mcols(subject)[[paste0(subjectname,"_Called")]] <- TRUE
    mcols(query)[[paste0(queryname,"_Called")]] <- TRUE
    
    if(maxgap > 0 & sizemargin > 0)
    {
      if (mode=="breakend"){
        hits <- findOverlaps(query = query, subject = subject,maxgap = maxgap, type=type)
      } else if (mode=="breakpoint"){
        hits <- findBreakpointOverlaps(query = query,
                                       subject = subject,
                                       maxgap = maxgap, 
                                       sizemargin = sizemargin)
      }
      print(paste0(queryname," vs ",subjectname," - ",mode,": ",length(hits)," of ", "(", length(query),"/", length(subject),")"))
      
      mcols(query)[[paste0(subjectname,"_Called")]] <- FALSE
      mcols(subject)[[paste0(queryname,"_Called")]] <- FALSE
      
      if(length(hits)>0)
      {
        #Mark observed hits
        mcols(query)[[paste0(subjectname,"_Called")]][queryHits(hits)] <- TRUE
        mcols(subject)[[paste0(queryname,"_Called")]][subjectHits(hits)] <- TRUE
        
        
        #Cross annotate hit names from vcf row names
        query_hit_to_subject_name_map <-  data.frame(queryHits=queryHits(hits), 
                                                     subjectNames=names(subject)[subjectHits(hits)]) %>% 
          group_by(queryHits) %>% 
          summarise(subjectNames=paste(subjectNames, collapse = ",", sep=","))
          
        subject_hit_to_query_name_map <-  data.frame(subjectHits=subjectHits(hits), 
                                                     queryNames=names(query)[queryHits(hits)]) %>% 
          group_by(subjectHits) %>% 
          summarise(queryNames=paste(queryNames, collapse = ",", sep=","))
        
        mcols(query)[[paste0(queryname,"_svID")]] <- names(query)
        mcols(subject)[[paste0(subjectname,"_svID")]] <- names(subject)
        
        mcols(query)[[paste0(subjectname,"_svID")]] <- NA_character_
        mcols(subject)[[paste0(queryname,"_svID")]] <- NA_character_
        
        mcols(query)[[paste0(subjectname,"_svID")]][query_hit_to_subject_name_map$queryHits] <- query_hit_to_subject_name_map$subjectNames
        mcols(subject)[[paste0(queryname,"_svID")]][subject_hit_to_query_name_map$subjectHits] <- subject_hit_to_query_name_map$queryNames
        
          
      } 
      
    } else {
      queryKey <- paste(seqnames(query),start(query),end(query), sep = "-")
      subjectKey <- paste(seqnames(subject),start(subject),end(subject), sep = "-")

      mcols(query)[[paste0(subjectname,"_Called")]] <- queryKey %in% subjectKey
      mcols(subject)[[paste0(queryname,"_Called")]] <-subjectKey %in% queryKey
      
      mcols(query)[[paste0(queryname,"_svID")]] <- names(query)
      mcols(subject)[[paste0(subjectname,"_svID")]] <- names(subject)
      
      mcols(query)[[paste0(subjectname,"_svID")]] <- NA_character_
      mcols(subject)[[paste0(queryname,"_svID")]] <- NA_character_
      
      query_hit_to_subject_name_map <- (data.frame(queryHits=1:length(queryKey), subjectNames=names(subject)[match(queryKey, subjectKey)]))
      query_hit_to_subject_name_map <- query_hit_to_subject_name_map[!is.na(query_hit_to_subject_name_map$subjectNames),]
      
      subject_hit_to_query_name_map <- (data.frame(subjectHits=1:length(subjectKey), queryNames=names(query)[match(subjectKey, queryKey)]))
      subject_hit_to_query_name_map <- subject_hit_to_query_name_map[!is.na(subject_hit_to_query_name_map$queryNames),]
      
      mcols(query)[[paste0(subjectname,"_svID")]][query_hit_to_subject_name_map$queryHits] <- query_hit_to_subject_name_map$subjectNames
      mcols(subject)[[paste0(queryname,"_svID")]][subject_hit_to_query_name_map$subjectHits] <- subject_hit_to_query_name_map$queryNames
      
      print(paste0(queryname," vs ",subjectname," - ",mode,": ",sum(queryKey %in% subjectKey)," of ", "(", length(query),"/", length(subject),")"))
    }
  }
  else {
    if(length(query)>0)
    {
      mcols(query)[[paste0(queryname,"_Called")]] <- TRUE
      mcols(query)[[paste0(subjectname,"_Called")]] <- FALSE
    }
    
    if(length(subject)>0)
    {
      mcols(subject)[[paste0(queryname,"_Called")]] <- FALSE
      mcols(subject)[[paste0(subjectname,"_Called")]] <- TRUE
    }
    print(paste0(queryname," vs ",subjectname," - ",mode,": ",0," of ", "(", length(query),"/", length(subject),")"))
  }
  returnval = list()
  returnval[[queryname]] <- query
  returnval[[subjectname]] <- subject
  return(returnval)
  
}

sv_combine <- function(set_a, set_b, set_a_name=deparse(substitute(set_a)), set_b_name=deparse(substitute(set_b)), 
                       mode, markcommon=TRUE, maxgap=NULL, type = "equal", sizemargin = NULL, remove_dup=FALSE)
{
  if(markcommon)
  {
    if (!is.null(maxgap) & !is.null(type) & !is.null(sizemargin))
    {
      marked <- mark_common_sv(mode = mode, query = set_a, subject = set_b, 
                               queryname = set_a_name, subjectname = set_b_name, 
                               maxgap = maxgap, type = type, sizemargin = sizemargin)
      mcols(marked[[set_a_name]])[["Present"]] <- ifelse(mcols(marked[[set_a_name]])[[paste0(set_b_name,"_Called")]], "Common", paste(set_a_name,"only",sep="_"))
      mcols(marked[[set_b_name]])[["Present"]] <- ifelse(mcols(marked[[set_b_name]])[[paste0(set_a_name,"_Called")]], "Common", paste(set_b_name,"only",sep="_"))
      # mcols(marked[[set_a_name]])[[paste0(set_b_name,"_Called")]] <- NULL
      # mcols(marked[[set_b_name]])[[paste0(set_a_name,"_Called")]] <- NULL
      set_a <- marked[[set_a_name]]
      set_b <- marked[[set_b_name]]
      
      if (remove_dup){
        set_b <- set_b[mcols(set_b)[["Present"]] != "Common",]
      }
      else
      {
        mcols(set_a)[mcols(set_a)[["Present"]]=="Common","Present"] <- paste("Common", set_a_name, sep="_")
        mcols(set_b)[mcols(set_b)[["Present"]]=="Common","Present"] <- paste("Common", set_b_name, sep="_")
      }
    }
    else {
      stop("Maxgap, type, and sizemargin are required by sv_combine() when markcommon=TRUE")
    }
  }
  
  return(GenomicRanges::bindROWS(x=set_a, objects = list(set_b)))
}
