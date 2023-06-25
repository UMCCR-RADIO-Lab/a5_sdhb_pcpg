#basedir <- "C:/ResearchData/RADIO/A5/Data"
library(dplyr)
library(tidyr)
library(VariantAnnotation)

# load_mafs <- function(maf_dir, combine=T) {
#   
#   A5_mafs_bcbio.files <- list.files(maf_dir,full.names = T, pattern = "ensemble-annotated.maf")
#   names(A5_mafs_bcbio.files) <- gsub("^.+/(E.+)-ensemble-annotated.maf","\\1",A5_mafs_bcbio.files)
#   
#   A5_mafs_umccr.files <- list.files(maf_dir,full.names = T, pattern="somatic.grch37.maf")
#   names(A5_mafs_umccr.files) <- gsub("^.+/((FL|E).+)__.+-somatic.grch37.maf","\\1",A5_mafs_umccr.files)
#   names(A5_mafs_umccr.files) <- gsub("FL-2","E233",names(A5_mafs_umccr.files))
#   
#   A5_mafs_bcbio <- lapply(A5_mafs_bcbio.files, read.delim, comment.char = '#')
#   
#   
#   A5_mafs_umccr <- lapply(A5_mafs_umccr.files,read.delim, comment.char = '#')
#   
#   annotate_UMMCR <- function(SampleName, umccr_mafs, bcbio_mafs)
#   {
#     umccr_maf <- umccr_mafs[[SampleName]] 
#     bcbio_maf <- bcbio_mafs[[SampleName]]
#     message(SampleName)
#     
#     
#     bcbio_maf.passfilter <- bcbio_maf %>% mutate(TUMOR_AF = as.numeric(TUMOR_AF)) %>%
#       filter(TUMOR_AF >= 0.1) %>%
#       # Remove anything labelled as a comon variant by vcf2maf
#       filter(FILTER != "common_variant")%>%
#       # Remove stuff with a higher gnomad freq
#       mutate(gnomAD_AF = as.numeric(gnomAD_AF)) %>%
#       mutate(gnomAD_AF = replace(gnomAD_AF, is.na(gnomAD_AF),0)) %>%
#       filter(gnomAD_AF == 0) %>%
#       filter(ExAC_FILTER == "PASS" | is.na(ExAC_FILTER) | ExAC_FILTER=="") %>%
#       mutate(TUMOR_MQ = as.numeric(TUMOR_MQ)) %>%
#       mutate(NORMAL_MQ = as.numeric(NORMAL_MQ)) %>%
#       mutate(TUMOR_VD = as.numeric(TUMOR_VD)) %>%
#       # Only keep low AF things if they have good mapping quality or a lot of tumour reads
#       filter(TUMOR_AF >= 0.2 | TUMOR_MQ >= 59 & NORMAL_MQ >= 59) %>% select(Chromosome, Start_Position, End_Position, Tumor_Seq_Allele2)
#     
#     return_df <- bcbio_maf %>% mutate(Tumor_Sample_Barcode=SampleName) %>%
#              left_join(umccr_maf %>% 
#                          select(Chromosome, Start_Position, End_Position, Tumor_Seq_Allele2, variant_id, variant_qual) %>% 
#                          mutate(umccrise_present=TRUE), by=c("Chromosome", "Start_Position", "End_Position", "Tumor_Seq_Allele2")) %>% 
#              left_join(bcbio_maf.passfilter %>% mutate(Andrew_PassFilter=TRUE)) %>% 
#              mutate(umccrise_present=replace_na(umccrise_present,FALSE),Andrew_PassFilter=replace_na(Andrew_PassFilter,FALSE)) %>% 
#              mutate(Present=ifelse(umccrise_present & Andrew_PassFilter,"Both", ifelse(Andrew_PassFilter, "AndrewOnly",ifelse(umccrise_present,"umccronly","Neither"))))
#     
#     return_df$PUBMED <- as.character(return_df$PUBMED)
#     
#     
#     return(return_df)
#   }
#   
#   A5_mafs_bcbio_ummcr_matched <- lapply(names(A5_mafs_bcbio), annotate_UMMCR,  umccr_mafs=A5_mafs_umccr, bcbio_mafs=A5_mafs_bcbio)
#   if (combine)
#   {
#     return(do.call("bind_rows", A5_mafs_bcbio_ummcr_matched))
#   } else {
#     return(A5_mafs_bcbio_ummcr_matched)
#   }
# }

restructure_partial_vector <- function(x)
{
  returnvector <- unlist(x)
  if(length(returnvector) != length(x)) # & any(isEmpty(x)))
  {
    vectorclass <- unique(unlist(lapply(x,class)))
    returnvector.placeholder <- vector(mode=vectorclass, length = length(x))
    returnvector.placeholder[!vapply(x, isEmpty, FUN.VALUE = F)] <- returnvector
    returnvector.placeholder[vapply(x, isEmpty, FUN.VALUE = F)] <- NA
    return(returnvector.placeholder)
  } else {
    return(unlist(x))
  }
  
}

collapse_nested_lists <- function(df) {
  data.frame(lapply(
    df,
    function(x) { 
      if(any(unlist(lapply(x, length) > 1))) 
      {
        return(as.character(unlist(lapply(x, paste ,collapse = ",", sep="," ))))
      } else
      {
          if("AsIs" %in% class(x) )
          {
            class(x) <- "list"
            x <- restructure_partial_vector(x)
            
          } else if (any(vapply(x, isEmpty, FUN.VALUE = F)))
          {
            x <- restructure_partial_vector(x)
          }
          else {
            return(unlist(x))}
      }
    }))
}

genofield_to_df <- function(genofield, vcf)
{
  if(length(dim(geno(vcf)[[genofield]]))==2)
  {
    
    gfield_df <- data.frame(geno(vcf)[[genofield]]) %>% tibble::rownames_to_column("VariantID")
    gfield_df <- collapse_nested_lists(gfield_df)
    
    Nidx <- which(grepl("[.]B0[0-9]",colnames(gfield_df)))
    colnames(gfield_df)[Nidx] <- paste("Normal", genofield, sep="_")
    
    if (ncol(gfield_df)==3) #Somatic VCF
    {
      Tidx <- which(grepl("[.]T0[0-9]",colnames(gfield_df)))
      colnames(gfield_df)[Tidx] <- paste("Tumour", genofield, sep="_")
    }
    
    return(gfield_df)
    
  } else if (length(dim(geno(vcf)[[genofield]]))==3){
    
    dftemp <- list()
    for (i in 1:dim(geno(vcf)[[genofield]])[3])
    {
      gfield_df <- data.frame(geno(vcf)[[genofield]][,,i,drop=F]) %>% tibble::rownames_to_column("VariantID")
      gfield_df <- collapse_nested_lists(gfield_df)
      
      Nidx <- which(grepl("[.]B0[0-9]",colnames(gfield_df)))
      colnames(gfield_df)[Nidx] <- paste("Normal", genofield, i, sep="_")
      
      if (ncol(gfield_df)==3) #Somatic VCF
      {
        Tidx <- which(grepl("[.]T0[0-9]",colnames(gfield_df)))
        colnames(gfield_df)[Tidx] <- paste("Tumour", genofield, sep="_")
      }
      
      dftemp[[i]] <- gfield_df
      
    }
    return(purrr::reduce(.x=dftemp,.f=full_join, by="VariantID"))
  } else {
    stop("Unhandled number of dimensions in genofield_to_df")
  }
  
}

extract_geno_fields <- function(vcf)
{
  geno_data_list <- list()
  for (gfield in names(geno(vcf)))
  {
    gfield_df <- genofield_to_df(gfield, vcf)
    geno_data_list[[gfield]] <- gfield_df
  }
  purrr::reduce(.x = geno_data_list, .f = full_join, by="VariantID")
}

vcf_to_dataframe <- function(vcf)
{
  message("Processing ", paste(samples(header(vcf)), sep="/", collapse = "/"), "...")
  fixed_df <- data.frame(rowRanges(vcf))
  fixed_df$ALT <- unstrsplit(CharacterList(fixed_df$ALT), sep = ",") #convert DNA string to character and collapse multiple alleles
  fixed_df$VariantID <- names(rowRanges(vcf))
  
  info_df <- collapse_nested_lists(data.frame(info(vcf)))
  
  geno_df <- extract_geno_fields(vcf)
  
  return_frame <- cbind(fixed_df, info_df)
  colnames(return_frame) <- make.unique(colnames(return_frame))
  return_frame <- return_frame %>% left_join(geno_df, by="VariantID")
  
  return(return_frame)
  
}

load_somatic_vcfs <- function(somatic_vcf_dir, combine=T, threads=1) {
  
  
  A5_vcfs_umccr.files <- list.files(somatic_vcf_dir, 
                                    pattern="-somatic-PASS.vcf.gz$", 
                                    full.names = T, 
                                    recursive = T)
  
  names(A5_vcfs_umccr.files) <- gsub("E.+__(.+)-somatic-PASS.vcf.gz",
                                     "\\1",
                                     basename(A5_vcfs_umccr.files))
  
  if(threads==1){
    message("Using single thread processing ...")
    
    A5_vcfs_umccr <- lapply(A5_vcfs_umccr.files,readVcf)
    A5_vcfs_umccr <- lapply(A5_vcfs_umccr, vcf_to_dataframe)
    
  } else {
    
    library(parallel)
    system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))
    cl <- makeCluster(threads, outfile="")
    clusterEvalQ(cl, library("dplyr"))
    clusterEvalQ(cl, library("purrr"))
    clusterEvalQ(cl, library("VariantAnnotation"))
    clusterExport(cl, c("collapse_nested_lists", "genofield_to_df", "extract_geno_fields", "vcf_to_dataframe", "restructure_partial_vector"))
    A5_vcfs_umccr <- parLapply(cl, A5_vcfs_umccr.files,readVcf)
    A5_vcfs_umccr <- parLapply(cl, A5_vcfs_umccr,vcf_to_dataframe)
    stopCluster(cl)
    rm(cl)
  }
  
  if (combine)
  {
    for (i in 1:length(A5_vcfs_umccr))
    {
      for (col in colnames(A5_vcfs_umccr[[i]]))
      {
        
        if(class(A5_vcfs_umccr[[i]][,col])=="AsIs")
        {
          A5_vcfs_umccr[[i]][,col]  <- as.list(A5_vcfs_umccr[[i]][,col])
        }
        A5_vcfs_umccr[[i]][,col] <- as.character(A5_vcfs_umccr[[i]][,col])
      }
    }
    
    A5_vcfs_umccr <- bind_rows(A5_vcfs_umccr,.id = "A5_ID")
    A5_vcfs_umccr <- A5_vcfs_umccr[,!(unlist(lapply(A5_vcfs_umccr, function (x) {all(is.na(x))} )))]
  }
  
  return(A5_vcfs_umccr)
}


load_germline_vcfs <- function(germline_cpsr_vcf_dir, combine=T, threads=1) {
  
  
  A5_vcfs_umccr.files <- list.files(germline_cpsr_vcf_dir, pattern="germline.predispose_genes.vcf.gz$", full.names = T, recursive = T)
  names(A5_vcfs_umccr.files) <- gsub("^.+/E.+__(.+)-germline.predispose_genes.vcf.gz","\\1",A5_vcfs_umccr.files)
  
  if(threads==1){
    A5_vcfs_umccr <- lapply(A5_vcfs_umccr.files,readVcf)
    A5_vcfs_umccr <- lapply(A5_vcfs_umccr, vcf_to_dataframe)
  } else {
    library(parallel)
    system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))
    cl <- makeCluster(threads, outfile="")
    clusterEvalQ(cl, library("dplyr"))
    clusterEvalQ(cl, library("purrr"))
    clusterEvalQ(cl, library("VariantAnnotation"))
    clusterExport(cl, c("collapse_nested_lists", "genofield_to_df", "extract_geno_fields", "vcf_to_dataframe", "restructure_partial_vector"))
    A5_vcfs_umccr <- parLapply(cl, A5_vcfs_umccr.files,readVcf)
    A5_vcfs_umccr <- parLapply(cl, A5_vcfs_umccr,vcf_to_dataframe)
    stopCluster(cl)
    rm(cl)
  }
  
  if (combine)
  {
    A5_vcfs_umccr <- lapply(A5_vcfs_umccr, function (vcf_df) { data.frame(lapply(vcf_df, as.character)) })
    
    # classbase <- unlist(lapply(A5_vcfs_umccr[[1]],class))
    
    #Edge case fix for one VCF that ends up with a character vector instead of a list
    # for (i in 1:length(A5_vcfs_umccr))
    # {
    #   classcurrent <- unlist(lapply(A5_vcfs_umccr[[i]],class))
    #   if (any(classcurrent != classbase))
    #   {
    #     mismatch <- which(classcurrent != classbase)
    #     if(classbase[[mismatch]]=="AsIs")
    #     {
    #       A5_vcfs_umccr[[i]][[mismatch]]  <- as.list(A5_vcfs_umccr[[i]][,mismatch])
    #     }
    #   }
    # }
    #A5_vcfs_umccr
    A5_vcfs_umccr <- bind_rows(A5_vcfs_umccr,.id = "A5_ID")
    
    A5_vcfs_umccr <- A5_vcfs_umccr[,!(unlist(lapply(A5_vcfs_umccr, function (x) {all(is.na(x))} )))]
    
    AnnoColNames <- c("Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID","Feature_Type",
                      "Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos_length",
                      "CDS.pos_length","AA.pos_length","Distance","ERRORS_WARNINGS_INFO")
    
    if ("ANN" %in% colnames(A5_vcfs_umccr))
      {
      A5_vcfs_umccr <- A5_vcfs_umccr %>% tidyr::separate(col = ANN, sep="\\|", into=AnnoColNames)
    }
    return(A5_vcfs_umccr)
  } else {
    return(A5_vcfs_umccr)
  }
}

load_germline_cpsr_tsv <- function(germline_cpsr_tsv_dir, exclude_clinvar_classes=NULL) {
  
  A5_germline_cpsr.files <- list.files(germline_cpsr_tsv_dir, pattern="-normal.cpsr.snvs_indels.tiers.tsv", full.names = T, recursive = T)
  names(A5_germline_cpsr.files) <- gsub("^.+/(E.+)-.__.+-normal.cpsr.snvs_indels.tiers.tsv","\\1",A5_germline_cpsr.files)
  A5_germline_cpsr.files <- A5_germline_cpsr.files[!duplicated(names(A5_germline_cpsr.files))]
  
  A5_germline_cpsr <- lapply(A5_germline_cpsr.files,read.delim)
  
  returnval <- bind_rows(A5_germline_cpsr,.id = "Patient_ID")
  
  returnval <- returnval %>% mutate(A5_ID=paste0(Patient_ID, "-B01")) %>% dplyr::relocate(A5_ID, .before = Patient_ID)
  
  if(!is.null(exclude_clinvar_classes))
  {
    returnval <- returnval %>% filter(!(CPSR_CLASSIFICATION %in% exclude_clinvar_classes))
  }
  
  return(returnval)
}

  