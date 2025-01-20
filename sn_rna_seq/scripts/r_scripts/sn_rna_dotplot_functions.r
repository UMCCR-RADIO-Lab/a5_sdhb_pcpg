# functions 
library(scales)
library(circlize)
library(ComplexHeatmap)
library(tidyr)
library(dplyr)
#library(ggsci)
#library(pals)

#----
# Dotplot Functions
#----

HeatmapDotPlot = function(colour, size,colour.label="Scaled exp", size.label="Fraction expressing", cell.size = 0.5, scale=TRUE, cluster_columns=FALSE, cluster_rows=FALSE, col=c("blue","red"), show_heatmap_legend=TRUE, row_title_rot = 0, ...){
  size = as.matrix(size)/max(size,na.rm = T)
  size = sqrt(size)
  if(scale){
    colour = colour %>% t() %>% scale() %>% t()
  }
  if(is.character(col)){
    if(scale){
      col_fun = colorRamp2(seq(from=-3, to=3, length.out=length(col)), col)
    } else {
      col_fun = colorRamp2(seq(from=min(colour), to=max(colour), length.out=length(col)), col)
    }
  } else {
    col_fun = col
  }
  hm = Heatmap(matrix = colour,
               row_gap = unit(2, "mm"), 
               column_gap = unit(2, "mm"), 
               cell_fun = function(j, i, x, y, width, height, fill){
                 grid.rect(x = x, y = y, width = width, height = height,
                           gp = gpar(col = "grey", fill = ifelse(is.na(size[i, j]) | is.nan(size[i, j]), "grey", NA)))
                 grid.circle(x = x, y = y, r = 0.5*(size[i,j])*max(unit.c(width, height)), 
                             gp = gpar(fill = col_fun(colour[i, j]), col = NA))
                 
                 # grid.rect(x = x, y = y, width = width, height = height,
                 #           gp = gpar(col = "grey", fill = NA))
                 # grid.circle(x = x, y = y, r = 0.5*(size[i,j])*max(unit.c(width, height)), 
                 #             gp = gpar(fill = col_fun(colour[i, j]), col = NA))
                 
                 # grid.rect(x = x, y = y, width = width, height = height,
                 #           gp = gpar(col = "grey", fill = NA))
                 # grid.circle(x = x, y = y, r = ifelse(is.na(size[i, j]), 0.5*1*max(unit.c(width, height)), 0.5*(size[i,j])*max(unit.c(width, height))), 
                 #             gp = gpar(fill = ifelse(is.na(colour[i, j]), "grey", col_fun(colour[i, j])), col = NA))
               },
               rect_gp = gpar(type="none"), 
               height = unit(nrow(colour)*cell.size, "cm"), 
               width = unit(ncol(colour)*cell.size, "cm"),
               cluster_columns=cluster_columns, 
               cluster_rows=cluster_rows, 
               show_heatmap_legend = TRUE, 
               row_title_rot = row_title_rot, ...)
  hm.legend = list()
  if(!is.null(size.label)){
    hm.legend = c(hm.legend, list(Legend(title = size.label,
                                         labels = c(0.25, 0.50, 0.75, 1.00) %>% as.character,
                                         size=unit.c(unit(sqrt(0.25)*cell.size, "cm"),
                                                     unit(sqrt(0.5)*cell.size, "cm"),
                                                     unit(sqrt(0.75)*cell.size , "cm"),
                                                     unit(sqrt(1.0)*cell.size, "cm")),
                                         type = "points",
                                         grid_height = unit(cell.size,"cm"),
                                         grid_width=unit(cell.size,"cm"),
                                         legend_height=unit(4*cell.size*2, "cm"),
                                         background = NA)))
  }
  if(!is.null(colour.label)){
    hm.legend = c(hm.legend, list(Legend(title=colour.label,
                                         col_fun = col_fun)))
  }
  return(hm) 
}

HeatmapDotPlot.Seurat = function(object, cells = NULL, features = NULL, assay = "RNA", slot="data", aggr.by = "orig.ident", split.by = NULL, gene_grouping = NULL, annot.columns=NULL, annot.colours = NULL, aggr.fun = mean, cell.size=0.5, show_legend=TRUE,show_annotation_name=FALSE, annotation_labels = NULL, legend_title = NULL, add_missing_genes=T, keep_aggr_names = FALSE, ...){
  if(class(aggr.fun) == "character"){aggr.fun = get(aggr.fun)}
  data = GetAssayData(object, slot=slot, assay=assay)
  if(length(features) > 0){
    data = data[intersect(features, rownames(data)), ,drop=F]
    missing.genes = setdiff(features, rownames(data))
    if(length(missing.genes) > 1) {cat(paste("Features not found in specified assay:", paste(missing.genes, collapse =", "), "\n"))}
  }
  aggr = object@meta.data[,aggr.by,drop=F]
  aggr = apply(aggr, 1, paste, collapse="___")
  aggr.levels = levels(factor(object@meta.data[,aggr.by[1]]))
  if(length(aggr.by) > 1){
    for(a in aggr.by[-1]){
      aggr.levels = unlist(lapply(aggr.levels, function(aggr.level){paste(aggr.level, levels(factor(object@meta.data[,a])), sep='___')}))
    }
  }
  aggr = factor(aggr, levels = aggr.levels)
  meta.data = object@meta.data
  
  all_samples <- unique(meta.data$orig.ident)
  all_cell_types <- unique(meta.data$cell_type)
  all_genotypes <- unique(meta.data$TERT_ATRX_Mutation)
  
  if(length(cells) > 0){
    if(is.character(cells)[1]){cells = match(cells, colnames(data))}
    data = data[,cells]
    aggr = aggr[cells]
    meta.data = meta.data[cells,]
  }
  
  aggr = factor(aggr, levels=intersect(levels(aggr), unique(as.character(aggr))))
  data = lapply(levels(aggr), function(x){data[,aggr==x,drop=F]}) %>% setNames(levels(aggr)) 
  colour = do.call(cbind, lapply(data, function(data){apply(data,1,aggr.fun)})) # average expression
  size = do.call(cbind, lapply(data, function(data){apply(as.matrix(data) > min(data),1, mean)})) # % expressed, assuming min(data) corresponds to 0 counts
  meta.data = meta.data[match(names(data), aggr),,drop=F]
  
  if(add_missing_genes)
  {
    missing_data <- matrix(data=NA, ncol=ncol(colour), nrow=length(missing.genes), dimnames=list(missing.genes, colnames(colour)))
    colour <- rbind(colour, missing_data)
    size <- rbind(size, missing_data)
  }
  
  colour <- colour[features,]
  size <- size[features,]
  
  if(!keep_aggr_names) {
  colnames(colour) = gsub("___.*$", "", colnames(colour))
  colnames(size) = gsub("___.*$", "", colnames(size))
  }
  
  grouping = NULL
  if(length(split.by) == 1){
    grouping = factor(meta.data[,split.by], levels = c(all_samples, 
                                                       all_cell_types, 
                                                       all_genotypes))
  }
  if(length(annot.columns) > 0){
    if(length(annot.colours) > 0){
      annot = HeatmapAnnotation(df = meta.data[,annot.columns,drop=F], col = annot.colours, which="column", show_legend = show_legend, show_annotation_name = show_annotation_name)
    } else {
      annot = HeatmapAnnotation(df = meta.data[,annot.columns,drop=F], which="column", show_legend=show_legend, show_annotation_name = show_annotation_name)
    }
    hm = HeatmapDotPlot(colour=colour, size=size, cell.size=cell.size, top_annotation = annot, show_heatmap_legend=show_legend, column_split=grouping, row_split = gene_grouping, ...)
  } else {
    hm = HeatmapDotPlot(colour=colour, size=size, cell.size=cell.size, column_split=grouping, row_split = gene_grouping, ...)
  }
  return(hm)
}


cell.size = 0.5
label="Scaled exp"
size.label="Fraction expressing"
hm.legend = list()
hm.legend = c(hm.legend, list(Legend(title = size.label,
                                     labels = c(0.25, 0.50, 0.75, 1.00) %>% as.character,
                                     size=unit.c(unit(sqrt(0.25)*cell.size, "cm"),
                                                 unit(sqrt(0.5)*cell.size, "cm"),
                                                 unit(sqrt(0.75)*cell.size , "cm"),
                                                 unit(sqrt(1.0)*cell.size, "cm")),
                                     type = "points",
                                     grid_height = unit(cell.size,"cm"),
                                     grid_width=unit(cell.size,"cm"),
                                     legend_height=unit(4*cell.size*2, "cm"),
                                     background = NA)))

# ----
# Heatmap Legend Parameters
# ----

hm_legend_params <- list(title = "Scaled exp",
                         legend_width = unit(3, "cm"),
                         title_position = "topcenter",
                         direction = "horizontal")

