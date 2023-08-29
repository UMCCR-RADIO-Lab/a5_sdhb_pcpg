setwd("/g/data/pq08/projects/ppgl/a5/")

library(tidyverse)
library(ChAMP)
library(furrr)

################
# Data Loading #
################

base_dir="/g/data/pq08/projects/ppgl"
data_dir <- paste0(base_dir,"/a5/methylation/raw_data/ILMLEPIC-16614/idat_symlinks")

myLoad <- champ.load(data_dir, arraytype = "EPIC")

control_samples <- c("E158_T01_D_1", "E158_T01_D_2", "E154_T01_D")

data(probe.features.epic)

#################
# CNA Detection #
#################

#The following code is a modified version of the champ.CNA function from the ChAMP package

intensity <- myLoad$intensity
pheno <- ifelse(colnames(intensity) %in% control_samples, "lowcna", "testgroup" )
controlGroup="lowcna"

#Extracts names of samples 
names <- colnames(intensity)
#Quantile normalises intensities	
intsqn <- preprocessCore::normalize.quantiles(as.matrix(intensity))
colnames(intsqn)<-names
#Calculates Log2
intsqnlog<-log2(intsqn)
  
message("<< Calculate mean value difference between each sample to mean control samples >>")
intsqnlogratio <- apply(intsqnlog[,which(!pheno %in% controlGroup)],2,function(x) x - rowMeans(as.data.frame(intsqnlog[,which(pheno %in% controlGroup)])))

ints <- data.frame(intensity,probe.features[rownames(intensity),c("MAPINFO","CHR")])
ints$MAPINFO <- as.numeric(ints$MAPINFO)

#Replaces Chr X and Y with 23 and 24
levels(ints$CHR)[levels(ints$CHR)=='X'] <- '23'
levels(ints$CHR)[levels(ints$CHR)=='Y'] <- '24'

message("<< Generate CHR and MAPINFO information >>")
CHR <- ints$CHR
MAPINFO <- ints$MAPINFO
  
  
#Runs CNA and generates individual DNA Copy Number profiles
threads=6
options(future.debug = TRUE)
options(future.globals.maxSize=10^10) 
plan(strategy="multisession", workers=threads) 
sampleResultFull <- 
  furrr::future_map(.x = 1:ncol(intsqnlogratio), 
                    .f= \(i) {
                      CNA.object <- DNAcopy::CNA(cbind(intsqnlogratio[,i]), CHR, MAPINFO ,data.type = "logratio", sampleid = paste(colnames(intsqnlogratio)[i],"qn"))
                      smoothed.CNA.object <- DNAcopy::smooth.CNA(CNA.object)
                      segment.smoothed.CNA.object <- DNAcopy::segment(smoothed.CNA.object, verbose = 1,alpha=0.001, undo.splits="sdundo", undo.SD=2)
                      return(segment.smoothed.CNA.object)
             
})


############
# Plotting #
############

gg_methcn <- list()
for (i in 1:length(sampleResultFull))
{
  plot_data_probe <- data.frame(sampleResultFull[[i]]$data)
  plot_data_seg <- data.frame(sampleResultFull[[i]]$output)
  offsets <- plot_data_probe %>%
    group_by(chrom) %>% 
    summarise(size=max(maploc)) %>% 
    arrange(as.numeric(chrom)) %>% 
    mutate(offset=cumsum(size),
           offset=offset-size) %>% 
    dplyr::select(-size)
  plot_data_probe <- plot_data_probe %>% inner_join(offsets) %>% mutate(maploc_offset=maploc+offset)
  
  sample=colnames(plot_data_probe)[[3]]
  colnames(plot_data_probe)[[3]] <- "log_ratio"
  plot_data_seg <- plot_data_seg %>% inner_join(offsets) %>% mutate(loc.start_offset=loc.start + offset, loc.end_offset=loc.end + offset)
  
  gg_methcn[[sample]] <- ggplot() + 
    geom_point(data = plot_data_probe %>% group_by(chrom) %>% slice_sample(prop = 0.1), 
               mapping = aes(x=maploc_offset, y=log_ratio, color=chrom), size=0.3) +
    geom_segment(data = plot_data_seg, 
                 mapping = aes(x=loc.start_offset, 
                               xend=loc.end_offset, 
                               y=seg.mean, 
                               yend=seg.mean), 
                 color="red")  + 
    theme_bw() + 
    coord_cartesian(ylim=c(-0.5,0.5)) +
    ggtitle(sample)
}
