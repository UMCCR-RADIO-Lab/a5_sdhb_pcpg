library(dplyr)
library(purrr)
library(StructuralVariantAnnotation)

#####################
# CLASS DEFINITIONS #
#####################

class_definitions_file <- "/g/data/pq08/projects/ppgl/a5/wgs/scripts/data_loaders/linx_class_definitions.r" 
source(class_definitions_file)

#############
# FUNCTIONS #
#############


load_linx_all_samples <- function(linx_base_dir, purple_base_dir, threads=1)
{
  
  suffixes <- c(breakend=".linx.breakend.tsv$", clusters=".linx.clusters.tsv$", 
                fusions=".linx.fusion.tsv$", svs=".linx.svs.tsv$", 
                purple_sv_vcf=".purple.sv.vcf.gz$")
  
  file_names <- list()
  for (s in names(suffixes)[1:4])
  {
    file_names[[s]] <- list.files(linx_base_dir, pattern = suffixes[[s]], recursive = T, full.names = T)
  }
  
  s=names(suffixes)[5]
  file_names[[s]] <- list.files(purple_base_dir, pattern = suffixes[[s]], recursive = T, full.names = T)
  
  if(!all(sapply(file_names, length)==sum(sapply(file_names, length))/length(file_names)))
  {
    stop("Sanity check failed: Linx/purple file list lengths are not all identical. Files may be missing.")
  }
  
  file_names <- purrr::map2(.x = file_names, .y=names(suffixes), function (filename, suffix_name, suffixes) { 
    return_frame <- tibble(A5_ID=gsub(suffixes[[suffix_name]],"", basename(filename)), "{suffix_name}" := filename) }, 
    suffixes=suffixes) %>% purrr::reduce(.f = left_join, by="A5_ID")
  
  
  
  if(threads == 1){
  return_svs <- list()
  
  for (i in 1:nrow(file_names))
  {
    message(paste("Processing", file_names$Sample[[i]], "..."))
    return_svs[[file_names$Sample[[i]]]] <- load_linxsv_from_output(sv_file = file_names$svs[[i]], 
                                                                    breakend_file = file_names$breakend[[i]], 
                                                                    fusion_file = file_names$fusions[[i]],
                                                                    cluster_file = file_names$clusters[[i]],
                                                                    purple_sv_vcf = file_names$purple_sv_vcf[[i]])
  } 
  } else {
    
    library(parallel)
    system(sprintf("taskset -p 0xffffffff %d", Sys.getpid()))
    cl <- makeCluster(threads, outfile="")
    #clusterEvalQ(cl, library("dplyr"))
    clusterEvalQ(cl, library("StructuralVariantAnnotation"))
    clusterExport(cl, c("class_definitions_file")) 
    clusterEvalQ(cl, source(class_definitions_file))
    clusterExport(cl, c("parse_svs", "parse_breakends","parse_clusters","parse_fusions","get_sv_from_vcf",
                        "load_linxsv_from_output"))
    
    
    return_svs <- parLapply(cl, split(file_names, file_names$A5_ID), function (sample_files) {
      
      func_args <- as.list(sample_files)
      
      if(!all(lapply(func_args, length) == 1)) { stop("Sanity check failed when splitting samples for processing")  }
      
      return(load_linxsv_from_output(sv_file = func_args$svs, 
                              breakend_file = func_args$breakend, 
                              fusion_file = func_args$fusions,
                              cluster_file = func_args$clusters,
                              purple_sv_vcf = func_args$purple_sv_vcf))
    })
    
    stopCluster(cl)
  }
  
  return(return_svs)
}

load_linxsv_from_output <- function(sv_file=NULL, breakend_file=NULL, fusion_file=NULL, cluster_file=NULL, purple_sv_vcf=NULL)
{
  if(!is.null(sv_file) & !is.null(breakend_file) & !is.null(fusion_file) & !is.null(cluster_file) & !is.null(purple_sv_vcf))
  {
    svs_input <- read.delim(sv_file)
    breakend_input <- read.delim(breakend_file)
    fusion_input <- read.delim(fusion_file)
    cluster_input <- read.delim(cluster_file)
      
    purple_vcf <- VariantAnnotation::readVcf(purple_sv_vcf)
    
    svs <- parse_svs(svs_input, breakend_input, fusion_input, cluster_input, purple_vcf)
    
    return(svs)
    
  }
  else(
    stop("All linx sv files must be supplied")
  )
  
}

get_sv_from_vcf <- function (sv_vcfId, purple_vcf, purple_vcf.breakpointRanges, purple_vcf.breakendRanges)
{
  
  if(sv_vcfId %in% names(purple_vcf))
  {
    if(grepl("(o|h)$",sv_vcfId))
    {
      mate_sv_vcfId <- unlist(info(purple_vcf[sv_vcfId])[["MATEID"]])
      sv_rows <- purple_vcf[names(purple_vcf) %in% c(sv_vcfId, mate_sv_vcfId),]
      bedpe <- data.frame(chrom1 = GenomeInfoDb::seqnames(sv_rows[sv_vcfId,]), 
                          start1 = start(sv_rows[sv_vcfId,]), end1 = end(sv_rows[sv_vcfId,]), chrom2 = GenomeInfoDb::seqnames(sv_rows[mate_sv_vcfId,]), 
                          start2 = start(sv_rows[mate_sv_vcfId,]), end2 = end(sv_rows[mate_sv_vcfId,]), 
                          RP = info(sv_rows[sv_vcfId,])[["RP"]], SR = info(sv_rows[sv_vcfId,])[["SR"]], BVF = info(sv_rows[sv_vcfId,])[["BVF"]])
    }
    else {
      sv_rows <- purple_vcf[names(purple_vcf) == sv_vcfId,]
      
      bedpe <- data.frame(chrom1 = GenomeInfoDb::seqnames(sv_rows[sv_vcfId,]), 
                          start1 = start(sv_rows[sv_vcfId,]), end1 = end(sv_rows[sv_vcfId,]), chrom2 = NA, 
                          start2 =  NA, end2 = NA, 
                          RP = info(sv_rows[sv_vcfId,])[["RP"]], SR = info(sv_rows[sv_vcfId,])[["SR"]], BVF = info(sv_rows[sv_vcfId,])[["BVF"]])
      
    }
  } else { stop(paste("sv_vcfId not found in Purple VCF:", sv_vcfId))}
  
  return(bedpe)
 
}

parse_svs <- function(svs_input=NULL, breakend_input=NULL, fusion_input=NULL, cluster_input=NULL, purple_vcf=NULL)
{
  return_sv <- list()
  total_sv <- nrow(svs_input)
  for (i in 1:total_sv)
  {
    
    return_sv[[i]] <- LINX_SV( 
              vcfId=(svs_input$vcfId[i]),
              svId=as.numeric(svs_input$svId[i]),
              clusterId=as.numeric(svs_input$clusterId[i]),
              clusterReason=as.character(svs_input$clusterReason[i]),
              fragileSiteStart=as.logical(svs_input$fragileSiteStart[i]),
              fragileSiteEnd=as.logical(svs_input$fragileSiteEnd[i]),
              isFoldback=as.logical(svs_input$isFoldback[i]),
              lineTypeStart=as.character(svs_input$lineTypeStart[i]),
              lineTypeEnd=as.character(svs_input$lineTypeEnd[i]),
              junctionCopyNumberMin=svs_input$junctionCopyNumberMin[i],
              junctionCopyNumberMax=svs_input$junctionCopyNumberMax[i],
              geneStart=svs_input$geneStart[i],
              geneEnd=as.character(svs_input$geneEnd[i]),
              replicationTimingStart=svs_input$replicationTimingStart[i],
              replicationTimingEnd=svs_input$replicationTimingEnd[i],
              localTopologyIdStart=svs_input$localTopologyIdStart[i],
              localTopologyIdEnd=svs_input$localTopologyIdEnd[i],
              localTopologyStart=svs_input$localTopologyStart[i],
              localTopologyEnd=svs_input$localTopologyEnd[i],
              localTICountStart=svs_input$localTICountStart[i],
              localTICountEnd=svs_input$localTICountEnd[i],
              breakends=parse_breakends(as.numeric(svs_input$svId[i]), breakend_input, fusion_input),
              clusters=parse_clusters(as.numeric(svs_input$clusterId[i]), cluster_input),
              vcf_info=get_sv_from_vcf(svs_input$vcfId[i], purple_vcf))
    message("Processed ", i, " of ", total_sv)
    
  }
  
  return(return_sv)
}

parse_breakends <- function(svId, breakend_input=NULL, fusion_input=NULL)
{
  BOI <- breakend_input[breakend_input$svId==svId,]
  return_breaks <- list()
  
  if(nrow(BOI) > 0 )
  {
    for (i in 1:nrow(BOI))
    {
      return_breaks[[i]] <- LINX_BREAKEND(
                                id=BOI$id[[i]],
                                svId=BOI$svId[[i]],
                                isStart=as.logical(BOI$isStart[[i]]),
                                gene=BOI$gene[[i]],
                                transcriptId=BOI$transcriptId[[i]],
                                canonical=as.logical(BOI$canonical[[i]]),
                                geneOrientation=BOI$geneOrientation[[i]],
                                disruptive=as.logical(BOI$disruptive[[i]]),
                                reportedDisruption=as.logical(BOI$reportedDisruption[[i]]),
                                undisruptedCopyNumber=BOI$undisruptedCopyNumber[[i]],
                                regionType=BOI$regionType[[i]],
                                codingContext=BOI$codingContext[[i]],
                                biotype=BOI$biotype[[i]],
                                exonicBasePhase=BOI$exonicBasePhase[[i]],
                                nextSpliceExonRank=BOI$nextSpliceExonRank[[i]],
                                nextSpliceExonPhase=BOI$nextSpliceExonPhase[[i]],
                                nextSpliceDistance=BOI$nextSpliceDistance[[i]],
                                totalExonCount=BOI$totalExonCount[[i]],
                                type=BOI$type[[i]],
                                chromosome=BOI$chromosome[[i]],
                                orientation=BOI$orientation[[i]],
                                strand=BOI$strand[[i]],
                                chrBand=BOI$chrBand[[i]],
                                exonUp=BOI$exonUp[[i]],
                                exonDown=BOI$exonDown[[i]],
                                junctionCopyNumber=BOI$junctionCopyNumber[[i]],
                                fusions=parse_fusions(BOI$id[[i]], fusion_input))
    }
  }
  
  return(return_breaks)
}

parse_fusions  <- function(breakendId, fusion_input=NULL)
{
  FOI <- fusion_input[fusion_input$fivePrimeBreakendId==breakendId | fusion_input$threePrimeBreakendId==breakendId,]
  return_fusions <- list()
  
  if(nrow(FOI) > 0 )
  {
    for (i in 1:nrow(FOI))
    {
      return_fusions[[i]] <- LINX_FUSION(fivePrimeBreakendId=FOI$fivePrimeBreakendId[[i]],
                                 threePrimeBreakendId=FOI$threePrimeBreakendId[[i]],
                                 name=FOI$name[[i]],
                                 reported=as.logical(FOI$reported[[i]]),
                                 reportedType=FOI$reportedType[[i]],
                                 phased=FOI$phased[[i]],
                                 likelihood=as.character(FOI$likelihood[[i]]),
                                 chainLength=FOI$chainLength[[i]],
                                 chainLinks=FOI$chainLinks[[i]],
                                 chainTerminated=as.logical(FOI$chainTerminated[[i]]),
                                 domainsKept=ifelse(is.na(FOI$domainsKept[[i]]), "", FOI$domainsKept[[i]]),
                                 domainsLost=ifelse(is.na(FOI$domainsLost[[i]]), "", FOI$domainsLost[[i]]),
                                 skippedExonsUp=FOI$skippedExonsUp[[i]],
                                 skippedExonsDown=FOI$skippedExonsDown[[i]],
                                 fusedExonUp=FOI$fusedExonUp[[i]],
                                 fusedExonDown=FOI$fusedExonDown[[i]],
                                 geneStart=FOI$geneStart[[i]],
                                 geneContextStart=FOI$geneContextStart[[i]],
                                 transcriptStart=FOI$transcriptStart[[i]],
                                 geneEnd=FOI$geneEnd[[i]],
                                 geneContextEnd=FOI$geneContextEnd[[i]],
                                 transcriptEnd=FOI$transcriptEnd[[i]],
                                 junctionCopyNumber=FOI$junctionCopyNumber[[i]])
    }
  }
  
  return(return_fusions)
}


parse_clusters  <- function(clusterId, cluster_input=NULL)
{
  
  COI <- cluster_input[cluster_input$clusterId==clusterId,]
  return_clusters <- list()
  
  if(nrow(COI) > 0 )
  {
    for (i in 1:nrow(COI))
    {
      return_clusters[[i]] <- LINX_CLUSTER(
                                  clusterId=as.numeric(COI$clusterId[[i]]),
                                  category=COI$category[[i]],
                                  synthetic=as.logical(COI$synthetic[[i]]),
                                  resolvedType=COI$resolvedType[[i]],
                                  clusterCount=as.numeric(COI$clusterCount[[i]]),
                                  clusterDesc=COI$clusterDesc[[i]])
    }
  }
  
  return(return_clusters)
}

is_fusion <- function(sv) {
  if (class(sv) != "LINX_SV") { stop("Expected object of class LINX_SV") }
  if (length(sv@breakends) > 0)
  {
    if (length(sv@breakends[[1]]@fusions) > 0) {return (TRUE)}
  }
  return (FALSE)
}

summarise_sv_as_row <- function(sv)
{
  GeneStartInfo <- list()
  if(!is.na(sv@geneStart) & sv@geneStart != "")
  {
    if (grepl(";", sv@geneStart)) {
      gene_list <- unlist(strsplit(sv@geneStart, ";"))
      
    } else {
      gene_list <- c(sv@geneStart)
    }
    
    for (gene in gene_list)
    {
      GeneStartInfo[[gene]] <- get_gene_info(gene, sv@breakends, TRUE)
      names(GeneStartInfo[[gene]]) <- gsub("^","GeneStart", names(GeneStartInfo[[gene]]))
    }
  }
  
  
  GeneEndInfo <- list()
  if(!is.na(sv@geneEnd) & sv@geneEnd != "")
  {
    if (grepl(";", sv@geneEnd)) {
      gene_list <- unlist(strsplit(sv@geneEnd, ";"))
      
    } else {
      gene_list <- c(sv@geneEnd)
    }
    
    for (gene in gene_list)
    {
      GeneEndInfo[[gene]] <- get_gene_info(gene, sv@breakends, FALSE)
      names(GeneEndInfo[[gene]]) <- gsub("^","GeneEnd", names(GeneEndInfo[[gene]]))
    }
  }
  
  sv_info <- list()
  sv_info[["isFusion"]] <- is_fusion(sv)
  sv_info[["fragileSite"]] <- sv@fragileSiteEnd |  sv@fragileSiteEnd
  sv_info[["Topology"]] <- ifelse(sv@localTopologyStart == sv@localTopologyEnd, sv@localTopologyStart, paste(sv@localTopologyStart,"/", sv@localTopologyEnd))
  
  vcf_info <- sv@vcf_info
  vcf_info <- data.frame(lapply(vcf_info, as.character))
  
  if (length(sv@clusters) > 0)
  {
    sv_info[["resolvedClusterType"]] <- sv@clusters[[1]]@resolvedType
    sv_info[["clusterId"]] <- sv@clusters[[1]]@clusterId
    sv_info[["clusterCount"]] <- sv@clusters[[1]]@clusterCount
  }
  
  sv_info[["vcfId"]] <- sv@vcfId
  
  if(length(GeneStartInfo) > 0)
  {
    GeneStartInfo <- dplyr::bind_rows(GeneStartInfo)
  }
  else
  {
    GeneStartInfo <- data.frame(GeneStartName=NA, GeneStartDisrupted=NA, GeneStartType=NA, GeneStartReportable=NA, GeneStartBiotype=NA, GeneStartRegion=NA)
  }
  
  if(length(GeneEndInfo) > 0)
  {
    GeneEndInfo <- dplyr::bind_rows(GeneEndInfo)
  }
  else
  {
    GeneEndInfo <- data.frame(GeneEndName=NA, GeneEndDisrupted=NA, GeneEndType=NA, GeneEndReportable=NA, GeneEndBiotype=NA, GeneEndRegion=NA)
  }
  
  return(tidyr::crossing(GeneStartInfo, dplyr::bind_rows(GeneEndInfo), data.frame(sv_info), vcf_info))
  
}


get_gene_info <- function(gene_name, breakend_list, start)
{
  info_list <- list()
  info_list[["Name"]] <- gene_name
  for (be in breakend_list)
  {
    if(be@isStart==start & be@gene==gene_name)
    {
      info_list[["Disrupted"]] <- be@disruptive
      info_list[["Type"]] <- be@type
      info_list[["Reportable"]] <- be@reportedDisruption
      info_list[["Biotype"]] <- be@biotype
      
      if(be@regionType == "EXONIC")
      {
        info_list[["Region"]] <- paste0("Exon-",be@exonUp)
      } else if(be@regionType == "INTRONIC")
      {
        info_list[["Region"]] <- paste0("Intronic-(",be@exonUp,"/",be@exonDown,")")
      } else
      {
        info_list[["Region"]] <- be@regionType
      }
    }
  }
  return(info_list)
}
