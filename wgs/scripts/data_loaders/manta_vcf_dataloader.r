library(StructuralVariantAnnotation)
library(dplyr)
library(tidyr)

fetchPurpleMantaFileNames <- function(basedir)
{
  A5_manta.files <- list.files(paste(basedir,"WGS/hg38/sv/manta/purple", sep="/"), pattern = ".purple.sv.vcf.gz", full.names = T, recursive = T)
  names(A5_manta.files) <- gsub("__.+$","",basename(A5_manta.files))
  return(A5_manta.files)
}

readPurpleMantaVCFs <- function(vcfFileNames, genome="hg38")
{
  return_vcfs <- list()
  
  for (i in 1:length(vcfFileNames))
  {
    if(file.size(vcfFileNames[[i]]) < 10000000)
    {
      print(vcfFileNames[[i]])
      return_vcfs[[names(vcfFileNames)[i]]] <- readVcf(vcfFileNames[[i]], genome = genome)
    }
    else
    { message("SV file ", vcfFileNames[i], " larger than 10MB, not processing") }
  }
  
  return(return_vcfs)

}


PurpleMantaVCFtoDataFrame <- function(PurpleMantaVCF, annotation_mode="full", removeGenelessAnnotation=T, removeHGVSp=T, removeLOF=T)
{

AnnoColNames <- c("Allele","Annotation","Annotation_Impact","Gene_Name","Gene_ID","Feature_Type",
              "Feature_ID","Transcript_BioType","Rank","HGVS.c","HGVS.p","cDNA.pos_length",
              "CDS.pos_length","AA.pos_length","Distance","ERRORS_WARNINGS_INFO")
SimpleAnnoColNames <- c("SVTYPE","EFFECT","GENE(s)","TRANSCRIPT","PRIORITY","PRIORITY_ORDINAL")

fixedcols <- data.frame(rowRanges(PurpleMantaVCF)) %>% dplyr::select(-strand, -paramRangeID, -QUAL)

fixedcols$ALT <- unlist(lapply(X = fixedcols$ALT, FUN = paste, sep=",", collapse=","))


if(annotation_mode=="full")
{
anno_df <- lapply(info(PurpleMantaVCF)[,"ANN"], function(x) {
  x <- unlist(x)
  x <- x[!grepl("interaction",x)]
  returnframe <- tibble(Anno=unlist(x)) %>% tidyr::separate(col = Anno, sep="\\|", into=AnnoColNames)
  
  if (removeHGVSp) { 
    returnframe <- returnframe[,-which(colnames(returnframe)=="HGVS.p")] 
    returnframe <- returnframe[,-which(colnames(returnframe)=="Rank")] 
    returnframe <- unique(returnframe)
  }
  if (removeGenelessAnnotation) { returnframe <- returnframe[returnframe$Gene_Name != "",] }
  if (removeLOF  & ("LOF" %in% colnames(returnframe))) { returnframe <- returnframe[,-which(colnames(returnframe)=="LOF")]  }
  if(nrow(returnframe)>0)
  {return(returnframe)} else
  { return (NULL)}
})

for (i in 1:nrow(fixedcols))
{
  if(!is.null(anno_df[[i]])){
    anno_df[[i]] <- bind_cols(fixedcols[i,c("seqnames","start","end","ALT")], anno_df[i])
  }
}

anno_df <- do.call("rbind",anno_df)
anno_df$HGVS.c <- gsub("%3B",";",anno_df$HGVS.c)  
if(!removeHGVSp) { anno_df$HGVS.p <- gsub("%3B",";",anno_df$HGVS.p) }
anno_df$cDNA.pos_length <- gsub("%3B",";",anno_df$cDNA.pos_length)
anno_df$CDS.pos_length <- gsub("%3B",";",anno_df$CDS.pos_length)
anno_df$AA.pos_length <- gsub("%3B",";",anno_df$AA.pos_length)
anno_df <- anno_df[,-which(colnames(anno_df)=="Allele")] 
#return_df <- fixedcols %>% left_join(anno_df, by=c("seqnames","start","end","ALT")) %>% left_join(inf, by=c("seqnames","start","end","ALT"))

} else
{
  
  anno_df <- lapply(info(PurpleMantaVCF)[,"SIMPLE_ANN"], function(x) {
    data.frame(Anno=unlist(x)) %>% mutate(Anno=gsub("\\|(?=\\||$)", "|-", Anno, perl = T)) %>% tidyr::separate(col = Anno, sep="\\|", into=SimpleAnnoColNames)
  })
  
  for (i in 1:nrow(fixedcols))
  {
    anno_df[[i]] <- bind_cols(fixedcols[i,c("seqnames","start","end","ALT")], anno_df[i])
  }
  
  anno_df <- do.call("rbind",anno_df)
  
}

inf <- cbind(fixedcols[,c("seqnames","start","end","ALT")], info(PurpleMantaVCF)[,-which(colnames(info(PurpleMantaVCF)) %in% c("ANN", colnames(anno_df)))])



return_df <- fixedcols %>% left_join(anno_df, by=c("seqnames","start","end","ALT")) %>% left_join(inf, by=c("seqnames","start","end","ALT"))
#return_df <- return_df %>% mutate(A5_ID=A5_ID) %>% relocate(A5_ID, .after = NULL)
 
return_df[,"BPI_AF"] <- unlist(lapply(return_df[,"BPI_AF"], paste, collapse=","))
return_df[,"CIEND"]  <- unlist(lapply(return_df[,"CIEND"], paste, collapse=","))
return_df[,"CIPOS"]  <- unlist(lapply(return_df[,"CIPOS"], paste, collapse=","))
return_df[,"LOF"]  <- unlist(lapply(return_df[,"LOF"], paste, collapse=","))
return_df[,"PURPLE_AF"]  <- unlist(lapply(return_df[,"PURPLE_AF"], paste, collapse=","))
return_df[,"PURPLE_CN"] <- unlist(lapply(return_df[,"PURPLE_CN"], paste, collapse=","))
return_df[,"PURPLE_CN_CHANGE"] <- unlist(lapply(return_df[,"PURPLE_CN_CHANGE"], paste, collapse=","))

return_df <- return_df[,-which(colnames(return_df)=="SIMPLE_ANN")] 
return_df <- return_df[,-which(colnames(return_df)=="CONTIG")] 

return(return_df)

}
