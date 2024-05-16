require(ggplot2)
require(rlang)
library(grid)

geom_split_tile <- function(mapping = NULL, data = NULL,
                            stat = "identity", position = "identity",
                            ...,
                            na.rm = FALSE,
                            show.legend = NA,
                            inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSplitTile,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list2(
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomSplitTile <- ggproto("GeomSplitTile", Geom,
                         draw_panel = function(self, data, panel_params, coord,
                                               lineend = "butt", linejoin = "round", linemitre = 10) {
                           data <- check_linewidth(data, snake_class(self))
                           n <- nrow(data)
                           if (n == 1) return(zeroGrob())
                           
                           if(class(data$side) != "factor") {
                             stop("Side column must be factors")
                           }
                           
                           #Check and force factors here <<---
                           assign(x = "input_data", value = data)
                           data <- make_triangles(data)
                           
                           munched <- coord_munch(coord, data, panel_params)
                           
                           # Sort by group to make sure that colors, fill, etc. come in same order
                           munched <- munched[order(munched$id),]
                           
                           # For gpar(), there is one entry per polygon (not one entry per point).
                           # We'll pull the first value from each group, and assume all these values
                           # are the same within each group.
                           first_idx <- !duplicated(munched$id)
                           first_rows <- munched[first_idx, ]
                           assign(x = "plot_data", value = data, envir = global_env())
                           assign(x = "plot_data_first_rows", value = data, envir = global_env())
                           ggplot2:::ggname(
                             "geom_split_tile",
                             polygonGrob(
                               munched$x, munched$y, default.units = "native",
                               id = munched$id,
                               gp = gpar(
                                 col = first_rows$colour,
                                 fill = alpha(first_rows$fill, first_rows$alpha),
                                 lwd = first_rows$linewidth * .pt,
                                 lty = first_rows$linetype,
                                 lineend = lineend,
                                 linejoin = linejoin,
                                 linemitre = linemitre
                               )
                             )
                           )
                           
                         },
                         
                         default_aes = aes(colour = NA, fill = "grey20", linewidth = 0.5, linetype = 1,
                                           alpha = NA, subgroup = NULL, gap=0.02, size=0.9),
                         
                         handle_na = function(data, params) {
                           data
                         },
                         
                         required_aes = c("x", "y", "side"),
                         
                         draw_key = draw_key_polygon,
                         
                         rename_size = TRUE
)

make_triangles <- function(input_data) {
  if(!inherits(input_data, "data.frame"))
  {
    stop("Input data must be a dataframe")
  }
  
  # pf = data.frame(input_data,
  #                 x_integer=as.numeric(input_data[["x"]]),
  #                 y_integer=as.numeric(input_data[["y"]]), 
  #                 side_integer=as.numeric(input_data$side))
  
  
  
  .functionarray <- function(operators, position, size) 
  {
    return_vec = vector(mode="numeric", length = length(operators))
    for (i in 1:length(operators))
    {
      return_vec[i] <- get(operators[i])(position, size)
      
    }
    return(return_vec)
  }
  
  .make_corners <- function (x, y, side, size, gap, ...) {
    
    
    if (as.numeric(side)==1)
    {
      x_coords = .functionarray(c("-","+","-"), x, size/2)
      y_coords = .functionarray(c("+","+","-"), y, size/2)
      x_coords = x_coords - c(0,gap,0)
      y_coords = y_coords + c(0,0,gap)
    } else {
      x_coords = .functionarray(c("+","-","+"), x, size/2)
      y_coords = .functionarray(c("-","-","+"), y, size/2)
      x_coords = x_coords + c(0,gap,0)
      y_coords = y_coords - c(0,0,gap)
    }
    
    ret_df <- data.frame(x=x_coords,
                         y=y_coords,
                         side=side,
                         id=paste(x,
                                  y,
                                  side, 
                                  sep="_"),
                         row.names = NULL)
    
    ret_df <- data.frame(ret_df, ...)
    
    return(ret_df)
  }
  
  rf <- dplyr::bind_rows(purrr::pmap(input_data, .make_corners))
  rf$id <- factor(rf$id)
  
  # rf[["x_label"]] <- factor(as.character(rf[["x_label"]]), levels= levels(input_data[,"x"]))
  # rf[["y_label"]] <- factor(as.character(rf[["y_label"]]), levels= levels(input_data[,"y"]))
  # rf[["fill"]]=factor(as.character(rf[,"fill"]), levels= levels(input_data[,"fill"]))
  
  
  for (col in intersect(colnames(input_data),colnames(rf)))
  {
    if (inherits(input_data[,col], what = "factor") & !inherits(rf[,col], what = "factor"))
    {
      rf[,col] <- factor(rf[,col], levels= levels(input_data[,col]))
    } else if (inherits(input_data[,col], what = "factor") == "factor" & inherits(rf[,col], what = "factor")){
      levels(rf[,col]) <- levels(input_data[,col])
    }
  }
  
  return(rf)
  
}

check_linewidth <- function(data, name) {
  if (is.null(data$linewidth) && !is.null(data$size)) {
    deprecate_soft0("3.4.0", I(paste0("Using the `size` aesthetic with ", name)), I("the `linewidth` aesthetic"))
    data$linewidth <- data$size
  }
  data
}
