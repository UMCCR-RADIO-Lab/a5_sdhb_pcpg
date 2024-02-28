########################
# Consensus Clustering #
########################

setwd("/g/data/pq08/projects/ppgl")
renv::activate("./a5")

library(ConsensusClusterPlus)
library(networkD3)
library(htmlwidgets)

############################################
# Import data loaders and run data mergers #
############################################

#clinical annotation
source("./a5/sample_annotation/scripts/data_loaders/a5_clinical_annotation_dataloader.r")
data_loader_a5_clinical_anno(google_account="aidan.flynn@umccr-radio-lab.page", use_cache=T)

#Methylation
source("./a5/tertiary/scripts/data_mergers/combine_tcga_comete_a5_methylation_data.r")

#WTS
source("./a5/tertiary/scripts/data_mergers/combine_tcga_flynn_a5_wts_data.r")

#small-rna
source("./a5/tertiary/scripts/data_mergers/combine_tcga_comete_a5_smallrna_data.r")

#Color pallets and themes
source("./a5/sample_annotation/scripts/data_loaders/a5_color_scheme.r")

########################
# Function Definitions #
########################

cv_feature_select <- function(data_matrix, cv_cutoff=NULL, keep_proportion=NULL) {
  
  calc_cv <- function(x) { sd(x)/mean(x) }
  
  if(is.null(cv_cutoff) & is.null(keep_proportion)){
    stop("You must specify either coefficient of variation cutoff or a proportion of data to retain")
  }
  
  if(!is.null(cv_cutoff) & !is.null(keep_proportion)){
    message("You have specified both a coefficient of variation cutoff or a proportion of data to retain. CV cutoff will be used.")
  }
  
  data_cv <- apply(data_matrix, 1, calc_cv)
  
  if(!is.null(cv_cutoff)){
     return(data_cv > cv_cutoff)
  }
  
  if(!is.null(keep_proportion)){
     n_total <- length(data_cv)
     n_keep <- floor(n_total*keep_proportion)
     cv_cutoff <- sort(data_cv, decreasing = T)[n_keep]
     return(data_cv > cv_cutoff)
  }
}

triangle = function(m,mode=1){
  #mode=1 for CDF, vector of lower triangle.
  #mode==3 for full matrix.
  #mode==2 for calcICL; nonredundant half matrix coun
  #mode!=1 for summary 
  n=dim(m)[1]
  nm = matrix(0,ncol=n,nrow=n)
  fm = m
  
  
  nm[upper.tri(nm)] = m[upper.tri(m)] #only upper half
  
  fm = t(nm)+nm
  diag(fm) = diag(m)
  
  nm=fm
  nm[upper.tri(nm)] = NA
  diag(nm) = NA
  vm = m[lower.tri(nm)]
  
  if(mode==1){
    return(vm) #vector 		
  }else if(mode==3){
    return(fm) #return full matrix
  }else if(mode == 2){
    return(nm) #returns lower triangle and no diagonal. no double counts.
  }
  
}

compute_delta_cdf=function(cc_out_object,breaks=100){
  #plot CDF distribution
  ml <- lapply(cc_out_object[2:length(cc_out_object)],"[[","ml")
  k=length(ml)
  areaK = c()
  for (i in 1:(length(ml))){
    v=triangle(ml[[i]],mode=1)
    
    #empirical CDF distribution. default number of breaks is 100    
    h = hist(v, plot=FALSE, breaks=seq(0,1,by=1/breaks))
    h$counts = cumsum(h$counts)/sum(h$counts)
    
    #calculate area under CDF curve, by histogram method.
    thisArea=0
    for (bi in 1:(length(h$breaks)-1)){
      thisArea = thisArea + h$counts[bi]*(h$breaks[bi+1]-h$breaks[bi]) #increment by height by width
      bi = bi + 1
    }
    areaK = c(areaK,thisArea)
  }
  
  #plot area under CDF change.
  deltaK=areaK[1] #initial auc at k=2
  for(i in 1:(length(areaK))){
    #proportional increase relative to prior K.
    deltaK = c(deltaK,( areaK[i] - areaK[i-1])/areaK[i-1])
  }
  return(deltaK)
}

run_con_clust <- function(dataset, top_cv_proportion, analysis_name=NULL, cohort_name=NULL, exclude_samples=NULL, scale=F, centre=F)
{
  if("character" %in% class(dataset)){ 
    dataset = get(dataset) 
    
    if(is.null(analysis_name)){ analysis_name <- dataset}
  } else  {
    #use the name of the variable used in function call
    if(is.null(analysis_name)){ analysis_name=deparse(substitute(dataset)) } 
  }
  
  if(!is.null(exclude_samples)){
    dataset <- dataset[,!(colnames(dataset) %in% exclude_samples)]
  }
  
  analysis_suffix <- ifelse(top_cv_proportion==1, "allfeatures", paste0("top",top_cv_proportion*100, "pccv"))
  analysis_prefix <- ifelse(grepl(cohort_name, analysis_name), analysis_name, paste(cohort_name, analysis_name, sep="_"))
  analysis_label=paste(analysis_prefix, analysis_suffix,sep="_")
  
  use_features <- cv_feature_select(dataset, keep_proportion=top_cv_proportion)
  dataset <- as.matrix(dataset[use_features,])
  
  cc_result <- list()
  
  if(scale || centre)
  {
    dataset <- t(scale(x = t(dataset), scale=scale, center=centre))
  }
  
  cc_result[["cc"]] <- ConsensusClusterPlus(dataset,
                                            pItem=0.7, 
                                            clusterAlg='hc', 
                                            distance='pearson', 
                                            innerLinkage='ward.D2', 
                                            finalLinkage='ward.D2', 
                                            reps=1000, 
                                            maxK=10, 
                                            #corUse='pairwise.complete.obs', 
                                            seed=1,           
                                            title=analysis_label,
                                            plot="png")
  
  cc_result[["icl"]]=calcICL(cc_result[["cc"]],title=analysis_label,plot="png")
  
  return(cc_result)
  
}

#Funtion to find the maximum value of k 
# at which all subsequent reductions in 
# the are under the CDF are below a threshold 
find_elbow <- function(delta_cdf, threshold=0.01) {
  names(delta_cdf) <- as.character(2:(length(delta_cdf)+1))
  above_threshold = 1000
  below_threshold = 1
  k=1
  while(max(above_threshold) > min(below_threshold)) {
    k <- k + 1
    delta_delta <- abs(delta_cdf[k:length(delta_cdf)] - delta_cdf[(k-1):(length(delta_cdf)-1)])
    names(delta_delta) <- paste(as.numeric(names(delta_delta))-1,names(delta_delta),sep="-")
    below_threshold <- which(delta_delta < threshold)
    above_threshold <- which(delta_delta > threshold)
    if(length(above_threshold)>0) {
    message("Found deltas above threshold at intervals:", toString(names(delta_delta)[above_threshold]))}
  }
  
  selected_k <- min(below_threshold) + k - 1 
  message ("Selected k=", selected_k)
  return(min(below_threshold) + k - 1 ) #-1 because k starts at 2
}

make_sankey_input <- function(input, flow) {
  if (!inherits(input, "data.frame")){
    stop("Input data must be in the form of a dataframe")
  }
  
  if(class(flow) != "character" & length(flow) < 2)
  {
    stop("Flow must be a character vector of column of length 2 or greater")
  }
  
  if(!all(flow %in% colnames(input)))
  {
    stop("Flow should contain only valid column names from input data")
  }
  
  links <- list()
  for (i in 2:length(flow)){
    f1 <- flow[i-1]
    f2 <- flow[i]
    
    links[[paste(f1,f2,sep="_")]] <-
      input %>% group_by(!!sym(f1),	!!sym(f2))  %>%
      dplyr::count(name="value") %>% 
      dplyr::rename(source=!!sym(f1),	target=!!sym(f2))
  }
  
  # A connection data frame is a list of flows with intensity for each flow
  links <- bind_rows(links)
  
  # From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(name=unique(unlist(links[c(1,2)])))
  
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  return(list(links=links,nodes=nodes))
}

###################
# Prepare A5 data #
###################

a5_wts_lcpm <- cpm(a5_wts_dge_list$NoEx, log = T)
a5_smallrna_lcpm <- cpm(a5_smallrna_dge, log = T)
rownames(a5_methylation_mval) <- a5_methylation_mval$Probe
a5_methylation_mval <- a5_methylation_mval[,-1]

##########################
##########################
# Clustering - All Data ##
##########################
##########################

rerun=T

old_wd <- getwd()
setwd("/g/data/pq08/projects/ppgl/a5/tertiary/results/consensus_clustering")

genotype_subtype_anno <- bind_rows(wts_tcga_flynn_a5_anno, 
                                   meth_tcga_comete_a5_27k_anno,
                                   meth_tcga_comete_a5_450k_anno,
                                   smallrna_tcga_comete_a5_anno) %>% 
  dplyr::select(Sample, Genotype, new_naming) %>% 
  distinct() %>%  dplyr::rename(subtype=new_naming)

#STAR expressing cases to exclude
star <- wts_tcga_flynn_a5_lcpm["ENSG00000147465",]
star.excl <- names(star)[(star - mean(star))/sd(star) >2]

#Genotypes other than SDHB
not_sdhb.excl <- genotype_subtype_anno$Sample[genotype_subtype_anno$Genotype!="SDHB"]


#########################
# Dataset configuration #
#########################

cohorts <- list()

#All cases
cohorts[["all_samples"]] <- list(dataset=c(smrna="smallrna_tcga_comete_a5_lcpm.batch_removed",
                                           wts="wts_tcga_flynn_a5_lcpm.batch_removed",
                                           meth="meth_tcga_comete_a5_27k_mval.batch_removed"),
                                 top_cv_proportion=c(smrna=1,
                                                     wts=0.2,
                                                     meth=1),
                                 exclude_samples=star.excl)

#SDHB Samples Only
cohorts[["all_sdhb"]] <- list(dataset=c(smrna="smallrna_tcga_comete_a5_lcpm.batch_removed",
                                        wts="wts_tcga_flynn_a5_lcpm.batch_removed",
                                        meth="meth_tcga_comete_a5_27k_mval.batch_removed"),
                              top_cv_proportion=c(smrna=1,
                                                  wts=0.2,
                                                  meth=1),
                              exclude_samples=c(star.excl, not_sdhb.excl))
#A5 Only
cohorts[["a5"]]<-  list(dataset=c(smrna="a5_smallrna_lcpm",
                                  wts="a5_wts_lcpm",
                                  meth="a5_methylation_mval"),
                        top_cv_proportion=c(smrna=1,
                                            wts=0.2,
                                            meth=1),
                        excluded_sample=NULL)


####################
# Run Conc. Clust. #
####################

checkpoint_file <- "wts_meth_smallrna_conclustplus_output.rds"

con_clust_output <- list()

if(rerun | !file.exists(checkpoint_file)) {
  for (cohort_name in names(cohorts))
  {
    template_df <- data.frame(cohort_name=cohort_name, cohorts[[cohort_name]][-3])
    template_df$analysis_name <- gsub("_(lcpm|mval)([.]batch_removed)?", "", template_df$dataset)
    rownames(template_df) <- template_df$analysis_name
      
    exclude_samples <- cohorts[[cohort_name]]$exclude_samples
    
    con_clust_output[[cohort_name]] <- purrr::pmap(.l=template_df, .f = run_con_clust, exclude_samples=exclude_samples)
    names(con_clust_output[[cohort_name]]) <- rownames(template_df)
  }
  
  saveRDS(object = con_clust_output, file = checkpoint_file)
  
  } else {
  con_clust_output <- readRDS(file = checkpoint_file)
}



for (cohort_name in names(con_clust_output))
{
  datasets <- con_clust_output[[cohort_name]]
  cohorts[[cohort_name]][["optimal_k"]] <- vector(mode = "numeric", length = length(datasets))
  names(cohorts[[cohort_name]][["optimal_k"]]) <- names(datasets)
  for (dataset in names(datasets))
  {
    delta_cdf <- compute_delta_cdf(con_clust_output[[cohort_name]][[dataset]]$cc, breaks = 100)
    cohorts[[cohort_name]][["optimal_k"]][dataset] <- find_elbow(delta_cdf, threshold = 0.015)
    # par(mfrow=c(1,2))
    # ConsensusClusterPlus:::CDF(c("spacer",lapply(con_clust_output[[dataset]]$cc[-1],"[[","ml")), breaks = 100)
  }
}

#Optimal Ks from CDF
# datatypes["smallrna_tcga_comete_a5_lcpm.batch_removed","optimal_k"] <- 7
# datatypes["wts_tcga_flynn_a5_lcpm.batch_removed","optimal_k"] <- 8
# datatypes["meth_tcga_comete_a5_27k_mval.batch_removed","optimal_k"] <- 6
# datatypes["meth_tcga_comete_a5_450k_mval.batch_removed","optimal_k"] <- 6

################
# Sankey Plots #
################

slice_highest_affinity <- function(input, optimal_k){
  best_cluster <- input$icl$itemConsensus %>% 
    filter(k==optimal_k) %>% 
    group_by(item) %>% 
    slice_max(n=1, order_by=itemConsensus)
  return(best_cluster) 
}

cluster_membership <- list()
sankey_input_raw <- list()
for (cohort_name in names(con_clust_output))
{
  cluster_membership[[cohort_name]] <- list()
  sankey_input_raw[[cohort_name]] <- list()
  for (dataset in names(con_clust_output[[cohort_name]]))
  {
    cluster_membership[[cohort_name]][["optimalK"]] <- map2(.x = con_clust_output[[cohort_name]], 
                                              .y= cohorts[[cohort_name]]$optimal_k, 
                                              .f = slice_highest_affinity) %>% 
      bind_rows(.id="dataset") %>% 
      pivot_wider(id_cols = item, names_from = dataset, values_from = cluster)
    
    cluster_membership[[cohort_name]][["k2"]] <- map2(.x = con_clust_output[[cohort_name]], 
                                                            .y= rep(2,length(con_clust_output[[cohort_name]])), 
                                                            .f = slice_highest_affinity) %>% 
      bind_rows(.id="dataset") %>% 
      pivot_wider(id_cols = item, names_from = dataset, values_from = cluster)
    

    
    for (k in names(cluster_membership[[cohort_name]]))
    {
      
      column_names <- colnames(cluster_membership[[cohort_name]][[k]])
      smallrna_col <- column_names[grepl("smallrna",column_names)]
      wts_col <- column_names[grepl("wts",column_names)]
      meth_col <- column_names[grepl("meth",column_names)]
      
    sankey_input_raw[[cohort_name]][[k]] <- cluster_membership[[cohort_name]][[k]] %>% 
      left_join(genotype_subtype_anno, by=c("item"="Sample")) %>% 
      mutate(across(everything(), as.character)) %>% 
      filter(!if_any(matches("small|wts|meth"), is.na)) %>%  
      mutate("{smallrna_col}":=paste("smRNA", !!sym(smallrna_col), sep="_"),
             "{wts_col}":=paste("WTS", !!sym(wts_col), sep="_"),
             "{meth_col}":=paste("meth27k", !!sym(meth_col), sep="_"))
    }
  }
}






#Define flows for sankey plots
flows <- list()
flows[["all_samples"]] <- list(
  c("Genotype", "subtype", "wts_tcga_flynn_a5", "meth_tcga_comete_a5_27k", "smallrna_tcga_comete_a5"),
  c("Genotype", "wts_tcga_flynn_a5"),
  c("Genotype", "meth_tcga_comete_a5_27k"),
  c("Genotype", "smallrna_tcga_comete_a5"),
  c("Genotype", "wts_tcga_flynn_a5", "meth_tcga_comete_a5_27k", "smallrna_tcga_comete_a5"),
  c("Genotype", "wts_tcga_flynn_a5", "smallrna_tcga_comete_a5"),
  c("subtype", "wts_tcga_flynn_a5")
)
flows[["all_sdhb"]] <- list(
  c("subtype", "wts_tcga_flynn_a5", "meth_tcga_comete_a5_27k", "smallrna_tcga_comete_a5"),
  c("subtype", "wts_tcga_flynn_a5"),
  c("subtype", "meth_tcga_comete_a5_27k"),
  c("subtype", "smallrna_tcga_comete_a5"),
  c("subtype", "wts_tcga_flynn_a5", "meth_tcga_comete_a5_27k", "smallrna_tcga_comete_a5"),
  c("subtype", "wts_tcga_flynn_a5", "smallrna_tcga_comete_a5")
)
flows[["a5"]] <- list(
  c("subtype", "a5_wts", "a5_methylation", "a5_smallrna"),
  c("subtype", "a5_methylation"),
  c("subtype", "a5_smallrna"),
  c("subtype", "a5_wts")
)


#Make sankey networks
sankey_networks <- list()
for (cohort_name in names(sankey_input_raw))
{
  for (k in names(sankey_input_raw[[cohort_name]]))
  {
    for(flow in flows[[cohort_name]]){
      
      sankey_input_ready <- make_sankey_input(input = sankey_input_raw[[cohort_name]][[k]], flow = flow)
      flow_label <- gsub("_tcga|_flynn|_?a5_?|_comete","",paste(flow,collapse="-"))
      
      sankey_label <- paste(cohort_name,k,flow_label,sep="_")
      
      sankey_networks[[sankey_label]] <- 
        sankeyNetwork(Links = sankey_input_ready$links, Nodes = sankey_input_ready$nodes,
                      Source = "IDsource", Target = "IDtarget",
                      Value = "value", NodeID = "name", 
                      sinksRight=FALSE,fontSize = 12,nodePadding = 2, 
                      width = 1400, height = 1000
                      )
    }
  }
}
#Write out hmtl for sankey plots
setwd("/g/data/pq08/projects/ppgl/a5/tertiary/results/consensus_clustering/sankey_widgets")
purrr::walk2(.x=sankey_networks, .y=paste0(names(sankey_networks),".html"), saveWidget)

setwd(old_wd)
