gometh <- function(sig.cpg, all.cpg=NULL, collection=c("GO","KEGG"), 
                   array.type = c("450K","EPIC"), plot.bias=FALSE, 
                   prior.prob=TRUE, anno=NULL, equiv.cpg = TRUE, 
                   fract.counts = TRUE, 
                   genomic.features = c("ALL", "TSS200","TSS1500","Body",
                                        "1stExon","3'UTR","5'UTR","ExonBnd"),
                   sig.genes = FALSE,
                   offline=F, 
                   offline_cache_dir=NULL)
  # Gene ontology testing or KEGG pathway analysis for Illumina methylation 
  # arrays based on goseq
  # Takes into account probability of differential methylation based on
  # numbers of probes on array per gene
  # Belinda Phipson
  # 28 January 2015. Last updated 1 September 2020.
  # EPIC functionality contributed by Andrew Y.F. Li Yim
{
  array.type <- match.arg(toupper(array.type), c("450K","EPIC"))    
  collection <- match.arg(toupper(collection), c("GO","KEGG"))
  genomic.features <- match.arg(genomic.features, c("ALL", "TSS200","TSS1500",
                                                    "Body", "1stExon","3'UTR",
                                                    "5'UTR","ExonBnd"), 
                                several.ok = TRUE)
  
  if(length(genomic.features) > 1 & any(grepl("ALL", genomic.features))){
    message("All input CpGs are used for testing.") 
    genomic.features <- "ALL"
  } 
  
  if(array.type == "450K" & any(grepl("ExonBnd", genomic.features))){
    stop("'ExonBnd' is not an annotated feature on 450K arrays,\n
           please remove it from your genomic.feature parameter\n
           specification.") 
  }
  
  if(collection == "GO"){
    go <- .getGO()
    result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=go$idList, 
                      array.type=array.type, plot.bias=plot.bias, 
                      prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
                      fract.counts=fract.counts, 
                      genomic.features = genomic.features,
                      sig.genes = sig.genes)
    result <- merge(go$idTable,result,by.x="GOID",by.y="row.names")
    rownames(result) <- result$GOID
    
  } else if(collection == "KEGG"){
    kegg <- .getKEGG(offline = offline, offline_cache_dir = offline_cache_dir)
    result <- gsameth(sig.cpg=sig.cpg, all.cpg=all.cpg, collection=kegg$idList, 
                      array.type=array.type, plot.bias=plot.bias, 
                      prior.prob=prior.prob, anno=anno, equiv.cpg=equiv.cpg,
                      fract.counts=fract.counts, 
                      genomic.features = genomic.features,
                      sig.genes = sig.genes)
    result <- merge(kegg$idTable,result,by.x="PathwayID",by.y="row.names")
    rownames(result) <- result$PathwayID
  }
  
  result[,-1]
}  

.getGO <- function(){
  if(!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("org.Hs.eg.db package required but not installed.")
  egGO2ALLEGS <- utils::getFromNamespace("org.Hs.egGO2ALLEGS", "org.Hs.eg.db")
  GeneID.PathID <- AnnotationDbi::toTable(egGO2ALLEGS)[,c("gene_id", "go_id", "Ontology")]
  d <- !duplicated(GeneID.PathID[, c("gene_id", "go_id")])
  GeneID.PathID <- GeneID.PathID[d, ]
  GOID.TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db, 
                                                      keys=unique(GeneID.PathID$go_id), 
                                                      columns=c("GOID","ONTOLOGY","TERM"), 
                                                      keytype="GOID"))
  go <- tapply(GeneID.PathID$gene_id, GeneID.PathID$go_id, list)
  
  list(idList=go, idTable=GOID.TERM)
}

.getKEGG <- function(offline=F, offline_cache_dir=NULL){
  if(offline & is.null(offline_cache_dir)) { stop("You must provide a cache directory for KEGG offline mode") }
  
  if(offline & !is.null(offline_cache_dir)) { 
    GeneID.PathID <- read.delim(file=paste(offline_cache_dir,"kegg_geneid_pathid.tsv", sep="/"),
                                sep="\t")
  }
  
  if(!offline) { 
    GeneID.PathID <- limma::getGeneKEGGLinks(species.KEGG = "hsa", convert = TRUE)  
    if (!is.null(offline_cache_dir)) {
      write.table(x = GeneID.PathID, 
                  file=paste(offline_cache_dir,"kegg_geneid_pathid.tsv", sep="/"),
                  sep="\t",
                  row.names=F)
    }
  }
  
  isna <- rowSums(is.na(GeneID.PathID[, 1:2])) > 0.5
  GeneID.PathID <- GeneID.PathID[!isna, ]
  ID.ID <- paste(GeneID.PathID[, 1], GeneID.PathID[, 2], sep = ".")
  d <- !duplicated(ID.ID)
  GeneID.PathID <- GeneID.PathID[d, ]
  
  if(offline & !is.null(offline_cache_dir)) { 
    PathID.PathName <- read.delim(file=paste(offline_cache_dir,"kegg_pathid_pathname.tsv", sep="/"),
                                  sep="\t")
  }
  
  if(!offline) { 
    PathID.PathName <- limma::getKEGGPathwayNames(species.KEGG = "hsa", 
                                                  remove.qualifier = TRUE)
    if (!is.null(offline_cache_dir)) {
      write.table(x = PathID.PathName, 
                  file=paste(offline_cache_dir,"kegg_pathid_pathname.tsv", sep="/"),
                  sep="\t",
                  row.names=F)
    }
  }
  
  PathID.PathName$PathwayID <- paste0("path:", PathID.PathName$PathwayID)
  GeneID.PathID <- merge(GeneID.PathID, PathID.PathName, by="PathwayID")
  kegg <- tapply(GeneID.PathID$GeneID, GeneID.PathID$PathwayID, list)
  
  list(idList = kegg, idTable = PathID.PathName)
}


.plotBias <- function(D,bias)
  # Plotting function to show gene level CpG density bias
  # Belinda Phipson
  # 5 March 2015
{
  o <- order(bias)
  splitf <- rep(1:100,each=200)[1:length(bias)]
  avgbias <- tapply(bias[o],factor(splitf),mean)
  sumDM <- tapply(D[o],factor(splitf),sum)
  propDM <- sumDM/table(splitf)
  graphics::par(mar=c(5,5,2,2))
  graphics::plot(avgbias,as.vector(propDM),
                 xlab="Number of CpGs per gene (binned)",
                 ylab="Proportion Differential Methylation",cex.lab=1.5,
                 cex.axis=1.2)
  graphics::lines(stats::lowess(avgbias,propDM),col=4,lwd=2)
}

.estimatePWF <- function(D,bias)
  # An alternative to goseq function nullp, which is transformation invariant
  # Belinda Phipson and Gordon Smyth
  # 6 March 2015
{
  prior.prob <- bias
  o <- order(bias)
  prior.prob[o] <- limma::tricubeMovingAverage(D[o],span=0.5)
  prior.prob
}

.getFlatAnnotation <- function(array.type=c("450K","EPIC"),anno=NULL)
  # flatten 450k or EPIC array annotation
  # Jovana Maksimovic
  # 18 September 2018
  # Updated 18 September 2018
  # Modified version of Belida Phipson's .flattenAnn code
{
  if(is.null(anno)){
    if(array.type=="450K"){
      anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
    } else {
      anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
  }
  
  # get rid of the non-CpG sites
  ann.keep<-anno[grepl("^cg",anno$Name),]
  
  # get rid of CpGs that are not annotated
  missing<-ann.keep$UCSC_RefGene_Name==""
  ann.keep<-ann.keep[!missing,]
  
  # get individual gene names for each CpG
  geneslist<-strsplit(ann.keep$UCSC_RefGene_Name,split=";")
  names(geneslist)<-rownames(ann.keep)
  
  grouplist<-strsplit(ann.keep$UCSC_RefGene_Group,split=";")
  names(grouplist)<-rownames(ann.keep)
  
  flat<-data.frame(symbol=unlist(geneslist),group=unlist(grouplist))
  flat$symbol<-as.character(flat$symbol)
  flat$group <- as.character(flat$group)
  
  flat$cpg<- substr(rownames(flat),1,10)
  
  #flat$cpg <- rownames(flat)
  flat$alias <- suppressWarnings(limma::alias2SymbolTable(flat$symbol))
  
  #eg <- toTable(org.Hs.egSYMBOL2EG)
  eg <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, 
                                               keys=AnnotationDbi::keys(org.Hs.eg.db), 
                                               columns=c("ENTREZID","SYMBOL"), 
                                               keytype="ENTREZID"))
  colnames(eg) <- c("gene_id","symbol")
  
  m <- match(flat$alias,eg$symbol)
  flat$entrezid <- eg$gene_id[m]
  flat <- flat[!is.na(flat$entrezid),]
  
  # keep unique cpg by gene name annotation
  id<-paste(flat$cpg,flat$entrezid,sep=".")
  d <- duplicated(id)
  flat.u <- flat[!d,]
  flat.u
  # This randomly samples only 1 gene ID for multimapping CpGs
  #.reduceMultiMap(flat.u)
}




