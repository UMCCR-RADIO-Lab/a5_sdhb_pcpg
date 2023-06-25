library(networkD3)
library(htmlwidgets)
library(dplyr)

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

#Reference Set
genotype_subtype_anno <- bind_rows(wts_tcga_flynn_a5_anno %>% dplyr::select(Sample, Genotype, new_naming), 
                             meth_tcga_comete_a5_27k_anno %>% dplyr::select(Sample, Genotype, new_naming),
                             meth_tcga_comete_a5_450k_anno %>% dplyr::select(Sample, Genotype, new_naming),
                             small_rna_tcga_comete_a5_anno %>% dplyr::select(Sample, Genotype, new_naming)) %>% 
  distinct() %>%  dplyr::rename(subtype=new_naming)

sankey_input_raw <- cluster_membership %>% left_join(genotype_subtype_anno, by=c("item"="Sample")) %>% 
  mutate(across(everything(), as.character)) %>% filter(!if_any(matches("small_|wts_|meth_"), is.na)) %>%  
  mutate(small_rna_tcga_comete_a5=paste("smRNA", small_rna_tcga_comete_a5, sep="_"),
         wts_tcga_flynn_a5=paste("WTS", wts_tcga_flynn_a5, sep="_"),
         meth_tcga_comete_a5_27k=paste("meth27k", meth_tcga_comete_a5_27k, sep="_"),
         meth_tcga_comete_a5_450k=paste("meth450k", meth_tcga_comete_a5_450k, sep="_"))

flows <- list(
  c("Genotype", "subtype", "wts_tcga_flynn_a5", "meth_tcga_comete_a5_27k", "small_rna_tcga_comete_a5"),
  c("Genotype", "wts_tcga_flynn_a5"),
  c("Genotype", "meth_tcga_comete_a5_27k"),
  c("Genotype", "small_rna_tcga_comete_a5"),
  c("Genotype", "wts_tcga_flynn_a5", "meth_tcga_comete_a5_27k", "small_rna_tcga_comete_a5"),
  c("Genotype", "wts_tcga_flynn_a5", "small_rna_tcga_comete_a5"),
  c("subtype", "wts_tcga_flynn_a5")
)

# links$group = ifelse(links$target %in% c("Metastatic","Recurrent", "Yes"),"1","2")
# nodes$group <- as.character(c(rep(1,84),2,3,4,5,6,7,8))

sankey_networks <- list()
for(flow in flows){
  sankey_input_ready <- make_sankey_input(input = sankey_input_raw, flow = flow)
  flow_label <- gsub("_tcga|_flynn|_a5|_comete","",paste(flow,collapse="-"))
  sankey_networks[[flow_label]] <- sankeyNetwork(Links = sankey_input_ready$links, Nodes = sankey_input_ready$nodes,
                                                 Source = "IDsource", Target = "IDtarget",
                                                 Value = "value", NodeID = "name", 
                                                 sinksRight=FALSE,fontSize = 12,nodePadding = 2, width = 900, height = 700, 
                                                 #                   NodeGroup = "group", 
                                                 #                  colourScale = 'd3.scaleOrdinal(["#000000","#c25858ff","#7884baff", "#eb7b54ff","#ee7474ff", "#abc4e2ff", "#efecdcff", "#b4b4b4ff"])'
  )
}

setwd("/g/data/pq08/projects/ppgl/a5/tertiary/results/consensus_clustering/sankey_widgets")
purrr::map2(.x=sankey_networks, .y=paste0(names(sankey_networks),".html"), saveWidget)
saveWidget(p, file="C:/ResearchData/RADIO/temp/A5_sankey_label.html")