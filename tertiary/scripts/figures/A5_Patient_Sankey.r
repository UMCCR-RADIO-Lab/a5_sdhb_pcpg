library(networkD3)
library(dplyr)

#Reference Set

SankeyInput <- A5_anno %>% select(`Patient ID`, A5_ID, is_primary_or_met, tumour_metastasised)

# A connection data frame is a list of flows with intensity for each flow
links <-   bind_rows(
  # SankeyInput %>% group_by(`Patient ID`, A5_ID)  %>% 
  #                      dplyr::count(name="value") %>% 
  #                      dplyr::rename(source=`Patient ID`,	target=A5_ID),
                     SankeyInput %>% group_by(`Patient ID`,	is_primary_or_met)  %>%
                       dplyr::count(name="value") %>% 
                       dplyr::rename(source=`Patient ID`,	target=is_primary_or_met),
                     SankeyInput %>% group_by(is_primary_or_met,	tumour_metastasised) %>% 
                       dplyr::count(name="value")  %>% 
                       dplyr::rename(source=is_primary_or_met,	target=tumour_metastasised))

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

links$group = ifelse(links$target %in% c("Metastatic","Recurrent", "Yes"),"1","2")
nodes$group <- as.character(c(rep(1,84),2,3,4,5,6,7,8))

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE,fontSize = 0,nodePadding = 2, width = 400, height = 1000, 
                   NodeGroup = "group", 
                   colourScale = 'd3.scaleOrdinal(["#000000","#c25858ff","#7884baff", "#eb7b54ff","#ee7474ff", "#abc4e2ff", "#efecdcff", "#b4b4b4ff"])')
p
networkD3::col

library(htmlwidgets)
saveWidget(p, file="C:/ResearchData/RADIO/temp/A5_sankey_label.html")
