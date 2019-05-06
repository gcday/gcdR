#' Feature expression heatmap
#'
#' Draws a heatmap of single cell feature expression.
#'
#' @param object Seurat object
#' @param features A vector of features to plot, defaults to \code{VariableFeatures(object = object)}
#' @param cells A vector of cells to plot
#' @param disp.min Minimum display value (all values below are clipped)
#' @param disp.max Maximum display value (all values above are clipped); defaults to 2.5
#' if \code{slot} is 'scale.data', 6 otherwise
#' @param group.by A vector of variables to group cells by; pass 'ident' to group by cell identity classes
#' @param group.bar Add a color bar showing group status for cells
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param assay Assay to pull from
# @param check.plot Check that plotting will finish in a reasonable amount of time
#' @param label Label the cell identies above the color bar
#' @param size Size of text above color bar
#' @param cols Colors to use for color bar
#' @param hjust Horizontal justification of text above color bar
#' @param angle Angle of text above color bar
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple dimensions
#'
#' @return A ggplot object
#'
#' @importFrom stats median
#' @importFrom scales hue_pal
#' @importFrom ggplot2 annotation_raster coord_cartesian ggplot_build aes_string
#' @export
#'
#' @examples
#' DoHeatmap(object = pbmc_small)
#'
gcdDoHeatmap <- function(
  object,
  features = NULL,
  cells = NULL,
  group.by = 'ident',
  group.bar = TRUE,
  disp.min = -2.5,
  disp.max = NULL,
  slot = 'scale.data',
  assay = NULL,
  label = TRUE,
  size = 5.5,
  cols = NULL,
  hjust = 0,
  angle = 45,
  combine = TRUE
) {
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(
    test = slot == 'scale.data',
    yes = 2.5,
    no = 6
  )
  # make sure features are present 
  possible.features <- rownames(x = GetAssayData(object = object, slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if(length(x = features) == 0) {
      stop("No requested features found in the ", slot, " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", slot, 
            " slot for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
  }
  data <- as.data.frame(x = t(as.matrix(x = GetAssayData(
    object = object, 
    slot = slot)[features, cells, drop = FALSE])))
  # data <- as.data.frame(x = t(FetchData(object = object,
  #   vars = features, cells = cells, slot = slot)))
  object <- suppressMessages(expr = StashIdent(object = object, save.name = 'ident'))
  group.by <- group.by %||% 'ident'
  groups.use <- object[[group.by]][cells, , drop = FALSE]
  # group.use <- switch(
  #   EXPR = group.by,
  #   'ident' = Idents(object = object),
  #   object[[group.by, drop = TRUE]]
  # )
  # group.use <- factor(x = group.use[cells])
  plots <- vector(mode = 'list', length = ncol(x = groups.use))
  for (i in 1:ncol(x = groups.use)) {
    group.use <- groups.use[, i, drop = TRUE]
    group.use <- factor(x = group.use)
    names(x = group.use) <- cells
    plot <- SingleRasterMap(
      data = data,
      disp.min = disp.min,
      disp.max = disp.max,
      feature.order = features,
      cell.order = names(x = sort(x = group.use)),
      group.by = group.use
    )
    if (group.bar) {
      # TODO: Change group.bar to annotation.bar
      pbuild <- ggplot_build(plot = plot)
      if (is.null(cols)) {
        cols <- hue_pal()(length(x = levels(x = group.use)))
      }
      names(x = cols) <- levels(x = group.use)
      y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + 0.25
      y.max <- y.pos + 0.5
      plot <- plot + annotation_raster(
        raster = t(x = cols[sort(x = group.use)]),
        xmin = -Inf,
        xmax = Inf,
        ymin = y.pos,
        ymax = y.max
      ) +
        coord_cartesian(ylim = c(0, y.max), clip = 'off')
      if (label) {
        x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
        x.divs <- pbuild$layout$panel_params[[1]]$x.major
        x <- data.frame(group = sort(x = group.use), x = x.divs)
        label.x.pos <- tapply(X = x$x, INDEX = x$group, FUN = median) * x.max
        label.x.pos <- data.frame(group = names(x = label.x.pos), label.x.pos)
        plot <- plot + geom_text(
          stat = "identity",
          data = label.x.pos,
          aes_string(label = 'group', x = 'label.x.pos'),
          y = y.max + y.max * 0.03 * 0.5,
          angle = angle,
          hjust = hjust,
          size = size
        )
        plot <- suppressMessages(plot + coord_cartesian(
          ylim = c(0, y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use))) * size),
          clip = 'off')
        )
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}
# Set a default value if an object is null
#
# @param lhs An object to set if it's null
# @param rhs The value to provide if x is null
#
# @return rhs if lhs is null, else lhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}
# A single heatmap from ggplot2 using geom_raster
#
# @param data A matrix or data frame with data to plot
# @param cell.order ...
# @param feature.order ...
# @param cols A vector of colors to use
# @param disp.min Minimum display value (all values below are clipped)
# @param disp.max Maximum display value (all values above are clipped)
# @param limits A two-length numeric vector with the limits for colors on the plot
# @param group.by A vector to group cells by, should be one grouping identity per cell
#
#' @importFrom ggplot2 ggplot aes_string geom_raster scale_fill_gradient
#' scale_fill_gradientn theme element_blank labs geom_point guides guide_legend
#
SingleRasterMap <- function(
  data,
  cell.order = NULL,
  feature.order = NULL,
  colors = PurpleAndYellow(),
  disp.min = -2.5,
  disp.max = 2.5,
  limits = NULL,
  group.by = NULL
) {
  data <- MinMax(data = data, min = disp.min, max = disp.max)
  data <- Melt(x = t(x = data))
  colnames(x = data) <- c('Feature', 'Cell', 'Expression')
  if (!is.null(x = feature.order)) {
    data$Feature <- factor(x = data$Feature, levels = unique(x = feature.order))
  }
  if (!is.null(x = cell.order)) {
    data$Cell <- factor(x = data$Cell, levels = unique(x = cell.order))
  }
  if (!is.null(x = group.by)) {
    data$Identity <- group.by[data$Cell]
  }
  limits <- limits %||% c(min(data$Expression), max(data$Expression))
  if (length(x = limits) != 2 || !is.numeric(x = limits)) {
    stop("limits' must be a two-length numeric vector")
  }
  plot <- ggplot(data = data) +
    geom_raster(mapping = aes_string(x = 'Cell', y = 'Feature', fill = 'Expression')) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_fill_gradientn(limits = limits, colors = colors) +
    labs(x = NULL, y = NULL, fill = group.by %iff% 'Expression') +
    WhiteBackground() + NoAxes(keep.text = TRUE)
  # if (!is.null(x = group.by)) {
  #   plot <- plot + geom_point(
  #     mapping = aes_string(x = 'Cell', y = 'Feature', color = 'Identity'),
  #     alpha = 0
  #   ) +
  #     guides(color = guide_legend(override.aes = list(alpha = 1)))
  # }
  return(plot)
}

# Melt a data frame
#
# @param x A data frame
#
# @return A molten data frame
#
Melt <- function(x) {
  if (!is.data.frame(x = x)) {
    x <- as.data.frame(x = x)
  }
  return(data.frame(
    rows = rep.int(x = rownames(x = x), times = ncol(x = x)),
    cols = unlist(x = lapply(X = colnames(x = x), FUN = rep.int, times = nrow(x = x))),
    vals = unlist(x = x, use.names = FALSE)
  ))
}
# Set a default value if an object is NOT null
#
# @param lhs An object to set if it's NOT null
# @param rhs The value to provide if x is NOT null
#
# @return lhs if lhs is null, else rhs
#
# @author Hadley Wickham
# @references https://adv-r.hadley.nz/functions.html#missing-arguments
#
`%iff%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(rhs)
  } else {
    return(lhs)
  }
}

#' Returns colors corresponding to each ident class.
#' @export
Palettes <- function(RET, type.use = 1, var.use = NULL) {
  library(scales)
  library(colorspace)
  if (is.null(var.use)) {
    num.levels <- length(levels(RET@seurat))
  } else {
    num.levels <- length(levels(as.factor(RET@seurat[[var.use]][[var.use]])))
  }
  if (type.use == 1) {
    return(hue_pal()(num.levels))
  } else if (type.use == 2) {
    return(rainbow_hcl(num.levels, c = 90, l = 80))
  } else if (type.use == 3) {
    return(rainbow_hcl(num.levels, c = 100, l = 80))
  } else {
    return(rainbow(num.levels))
  }
  # col.palettes <- c(hue_pal()(length(levels(RET@seurat))),
                      # rainbow_hcl(length(levels(RET@seurat)), c = 90, l = 80),
                      # rainbow_hcl(length(levels(RET@seurat)), c = 100, l = 80),
                      # rainbow(length(levels(RET@seurat))))
  # return(col.palettes[as.integer(type.use)])
}
#' Prints table as percentage values
#'
#' @param tbl table
#' 
#' @return none
#'
#' @examples
#' RET <- makeMarkerHeatmaps(RET, marker.lists)
#' printMarkerHeatmaps(RET)
#'
#' @export
percent.table <- function(tbl) {
  props <- tbl/rowSums(tbl)
  return(t(round(props * 100, 1)))
}


#' Summarizes 2 variables
#' 
#' @param var.1 first list of values with length n 
#' @param var.2 second list of values with length m
#' @param transpose whether to flip results
#' @param do.percent whether to use percent or raw counts (default)
#' @param do.format whether to return percent-formatted strings or numeric
#' 
#' @importFrom dplyr mutate_all
#' @importFrom scales percent
#' 
#' @export
breakdownTable <- function(var.1, 
                           var.2, 
                           transpose = F, 
                           do.percent = T, 
                           do.format = T) {
  sum.tbl <- table(var.1, var.2)
  if (do.percent) {
    sum.tbl <- percent.table(sum.tbl) 
  } 
  if (transpose) sum.tbl <- t(sum.tbl)
  MAT <- as.data.frame.matrix(sum.tbl)
  
  if (do.percent & do.format) {
    new.MAT <- dplyr::mutate_all(MAT, funs(scales::percent(., accuracy = 0.1, scale = 1)))
    rownames(new.MAT) <- rownames(MAT)
    return(new.MAT)
  } else if (!do.percent) {
    return(t(MAT))
  }
}

#' Prints markers distinguishing clusters in a Seurat object
#'
#' @param RET list containing Seurat object and plots
#' 
#' @return list containing clustered Seurat object and TSNE plots
#'
#' @examples
#' printSeuratMarkers(RET)
#'
#' @export
printSeuratMarkers <- function(RET) {
  for (ident in levels(Idents(RET@seurat))) {
    print(ident)
    print(prettyPrintMarkers(RET@meta.list$all.markers.full, ident))
  }
}



#' Makes heatmaps for Seurat object
#'
#'
#' @param RET gcdSeurat object
#' @param marker.lists list of markers (may be genes or module scores)
#' 
#' @return list containing Seurat object and marker heatmaps
#'
#' @examples
#' makeMarkerHeatmaps(RET, marker.lists)
#'
#' @export
makeMarkerHeatmaps <- function(RET, marker.lists) {
  require("Seurat")
  require("ggplot2")
  marker.plots <- list(violin = list(), feature = list(), dotplot = list())
  for (marked.group in names(marker.lists)) {
    markers <- marker.lists[[marked.group]][marker.lists[[marked.group]] %in% rownames(RET@seurat)]
    
    # RET@plots$markers[[marked.group]]$heatmap <- DoHeatmap(SubsetData(RET@seurat, max.cells.per.ident = 500), 
    # features = marker.lists[[marked.group]]) + ggtitle(marked.group)
    marker.plots$violin[[marked.group]] <- VlnPlot(RET@seurat, features = markers, pt.size = 0	) #+ ggplot2::ggtitle(marked.group)
    # marker.plots$feature[[marked.group]] <- FeaturePlot(RET@seurat, features = markers, reduction = "umap") #+ ggplot2::ggtitle(marked.group)
    marker.plots$dotplot[[marked.group]] <- DotPlot(RET@seurat, features = markers) + ggplot2::ggtitle(marked.group)
  }
  RET@plots$markers <- marker.plots
  return(RET)
}

#' Prints heatmaps for Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' 
#' @return none
#'
#' @examples
#' RET <- makeMarkerHeatmaps(RET, marker.lists)
#' printMarkerHeatmaps(RET)
#'
#' @export
printMarkerHeatmaps <- function(RET) {
  for(name in names(RET@plots$markers)) {
    print(plot_grid(RET@plots$markers[[name]]$violin, RET@plots$markers[[name]]$feature, ncol=1))
  }
}

# Plot a single expression by identity on a plot
#
# @param type Make either a 'ridge' or 'violin' plot
# @param data Data to plot
# @param idents Idents to use
# @param sort Sort identity classes (on the x-axis) by the average
# expression of the attribute being potted
# @param y.max Maximum Y value to plot
# @param adjust Adjust parameter for geom_violin
# @param cols Colors to use for plotting
# @param log plot Y axis on log scale
#
# @return A ggplot-based Expression-by-Identity plot
#
# @import ggplot2
#' @importFrom stats rnorm
#' @importFrom utils globalVariables
#' @importFrom ggridges geom_density_ridges theme_ridges stat_density_ridges
#' @importFrom ggplot2 ggplot aes_string theme labs geom_violin geom_jitter ylim
#' scale_fill_manual scale_y_log10 scale_x_log10 scale_y_discrete scale_x_continuous waiver
#' @importFrom cowplot theme_cowplot
#'
SingleExIPlot <- function(
  data,
  idents,
  split = NULL,
  type = 'violin',
  sort = FALSE,
  y.max = NULL,
  adjust = 1,
  pt.size = 0,
  cols = NULL,
  log = FALSE,
  show.quantiles = T,
  quantiles = 2
) {
  set.seed(seed = 42)
  if (!is.data.frame(x = data) || ncol(x = data) != 1) {
    stop("'SingleExIPlot requires a data frame with 1 column")
  }
  feature <- colnames(x = data)
  data$ident <- idents
  if ((is.character(x = sort) && nchar(x = sort) > 0) || sort) {
    data$ident <- factor(
      x = data$ident,
      levels = names(x = rev(x = sort(
        x = tapply(
          X = data[, feature],
          INDEX = data$ident,
          FUN = mean
        ),
        decreasing = grepl(pattern = paste0('^', tolower(x = sort)), x = 'decreasing')
      )))
    )
  }
  if (log) {
    noise <- rnorm(n = length(x = data[, feature])) / 200
    data[, feature] <- data[, feature] + 1
  } else {
    noise <- rnorm(n = length(x = data[, feature])) / 100000
  }
  if (all(data[, feature] == data[, feature][1])) {
    warning(paste0("All cells have the same value of ", feature, "."))
  } else{
    data[, feature] <- data[, feature] + noise
  }
  axis.label <- ifelse(test = log, yes = 'Log Expression Level', no = 'Expression Level')
  y.max <- y.max %||% max(data[, feature])
  if (is.null(x = split) || type != 'violin') {
    vln.geom <- geom_violin
    fill <- 'ident'
  } else {
    data$split <- split
    vln.geom <- geom_split_violin
    fill <- 'split'
  }
  switch(
    EXPR = type,
    'violin' = {
      x <- 'ident'
      y <- paste0("`", feature, "`")
      xlab <- 'Identity'
      ylab <- axis.label
      vln.quantiles <- NULL
      if (show.quantiles) {
        vln.quantiles <- seq(0, 1, length.out = quantiles + 1)[2:quantiles]
      }
      geom <- list(
        vln.geom(scale = 'width', adjust = adjust, trim = TRUE,
                 draw_quantiles = vln.quantiles),
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      )
      jitter <- geom_jitter(height = 0, size = pt.size)
      log.scale <- scale_y_log10()
      axis.scale <- ylim
    },
    'ridge' = {
      x <- paste0("`", feature, "`")
      y <- 'ident'
      xlab <- axis.label
      ylab <- 'Identity'
      geom <- list(
        geom_density_ridges(scale = 4),
        theme_ridges(),
        scale_y_discrete(expand = c(0.01, 0)),
        scale_x_continuous(expand = c(0, 0))
      )
      if (show.quantiles) {
        geom <- append(geom, 
                       stat_density_ridges(quantile_lines = T, 
                                           scale = 4, 
                                           quantiles = quantiles))
      }
      jitter <- geom_jitter(width = 0, size = pt.size)
      log.scale <- scale_x_log10()
      axis.scale <- function(...) {
        invisible(x = NULL)
      }
    },
    stop("Unknown plot type: ", type)
  )
  plot <- ggplot(
    data = data,
    mapping = aes_string(x = x, y = y, fill = fill)[c(2, 3, 1)]
  ) +
    labs(x = xlab, y = ylab, title = feature, fill = NULL) +
    theme_cowplot()
  plot <- do.call(what = '+', args = list(plot, geom))
  plot <- plot + if (log) {
    log.scale
  } else {
    axis.scale(min(data[, feature]), y.max)
  }
  if (pt.size > 0) {
    plot <- plot + jitter
  }
  if (!is.null(x = cols)) {
    if (!is.null(x = split)) {
      idents <- unique(x = as.vector(x = data$ident))
      splits <- unique(x = as.vector(x = data$split))
      labels <- if (length(x = splits) == 2) {
        splits
      } else {
        unlist(x = lapply(
          X = idents,
          FUN = function(pattern, x) {
            x.mod <- gsub(
              pattern = paste0(pattern, '.'),
              replacement = paste0(pattern, ': '),
              x = x,
              fixed = TRUE
            )
            x.keep <- grep(pattern = ': ', x = x.mod, fixed = TRUE)
            x.return <- x.mod[x.keep]
            names(x = x.return) <- x[x.keep]
            return(x.return)
          },
          x = unique(x = as.vector(x = data$split))
        ))
      }
      if (is.null(x = names(x = labels))) {
        names(x = labels) <- labels
      }
    } else {
      labels <- unique(x = as.vector(x = data$ident))
    }
    plot <- plot + scale_fill_manual(values = cols, labels = labels)
  }
  return(plot)
}

#' Single cell ridge plot
#'
#' Draws a ridge plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @param object Seurat object
#' @param features Features to plot (gene expression, metrics, PC scores,
#' anything that can be retreived by FetchData)
#' @param cols Colors to use for plotting
#' @param idents Which classes to include in the plot (default is all)
#' @param sort Sort identity classes (on the x-axis) by the average
#' expression of the attribute being potted, can also pass 'increasing' or 'decreasing' to change sort direction
#' @param assay Name of assay to use, defaults to the active assay
#' @param group.by Group (color) cells in different ways (for example, orig.ident)
#' @param y.max Maximum y axis value
#' @param same.y.lims Set all the y-axis limits to the same values
#' @param log plot the feature axis on log scale
#' @param ncol Number of columns if multiple plots are displayed
#' @param combine Combine plots into a single gg object; note that if TRUE; themeing will not work when plotting multiple features
#' @param slot Use non-normalized counts data for plotting
#' @param ... Extra parameters passed on to \code{\link{CombinePlots}}
#'
#' @return A ggplot object
#'
#' @export
#'
#' @examples
#' RidgePlot(object = pbmc_small, features = 'PC1')
#'
GCD.RidgePlot <- function(
  object,
  features,
  cols = NULL,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  group.by = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  log = FALSE,
  ncol = NULL,
  combine = TRUE,
  slot = 'data',
  show.quantiles = TRUE,
  quantiles = 2,
  ...
) {
  return(ExIPlot(
    object = object,
    type = 'ridge',
    features = features,
    idents = idents,
    ncol = ncol,
    sort = sort,
    assay = assay,
    y.max = y.max,
    same.y.lims = same.y.lims,
    cols = cols,
    group.by = group.by,
    log = log,
    combine = combine,
    slot = slot,
    show.quantiles = show.quantiles,
    quantiles = quantiles,
    ...
  ))
}
# Plot feature expression by identity
#
# Basically combines the codebase for VlnPlot and RidgePlot
#
# @param object Seurat object
# @param plot.type Plot type, choose from 'ridge' or 'violin'
# @param features Features to plot (gene expression, metrics, PC scores,
# anything that can be retreived by FetchData)
# @param idents Which classes to include in the plot (default is all)
# @param ncol Number of columns if multiple plots are displayed
# @param sort Sort identity classes (on the x-axis) by the average expression of the attribute being potted
# @param y.max Maximum y axis value
# @param same.y.lims Set all the y-axis limits to the same values
# @param adjust Adjust parameter for geom_violin
# @param pt.size Point size for geom_violin
# @param cols Colors to use for plotting
# @param group.by Group (color) cells in different ways (for example, orig.ident)
# @param split.by A variable to split the plot by
# @param log plot Y axis on log scale
# @param combine Combine plots using cowplot::plot_grid
# @param slot Use non-normalized counts data for plotting
# @param ... Extra parameters passed to \code{\link{CombinePlots}}
#
#' @importFrom scales hue_pal
#
ExIPlot <- function(
  object,
  features,
  type = 'violin',
  idents = NULL,
  ncol = NULL,
  sort = FALSE,
  assay = NULL,
  y.max = NULL,
  same.y.lims = FALSE,
  adjust = 1,
  cols = NULL,
  pt.size = 0,
  group.by = NULL,
  split.by = NULL,
  log = FALSE,
  combine = TRUE,
  slot = 'data',
  show.quantiles = FALSE,
  quantiles = NULL,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  ncol <- ncol %||% ifelse(
    test = length(x = features) > 9,
    yes = 4,
    no = min(length(x = features), 3)
  )
  data <- FetchData(object = object, vars = features, slot = slot)
  features <- colnames(x = data)
  if (is.null(x = idents)) {
    cells <- colnames(x = object)
  } else {
    cells <- names(x = Idents(object = object)[Idents(object = object) %in% idents])
  }
  data <- data[cells, , drop = FALSE]
  idents <- if (is.null(x = group.by)) {
    Idents(object = object)[cells]
  } else {
    object[[group.by, drop = TRUE]][cells]
  }
  if (!is.factor(x = idents)) {
    idents <- factor(x = idents)
  }
  if (is.null(x = split.by)) {
    split <- NULL
  } else {
    split <- object[[split.by, drop = TRUE]][cells]
    if (!is.factor(x = split)) {
      split <- factor(x = split)
    }
    if (is.null(x = cols)) {
      cols <- hue_pal()(length(x = levels(x = idents)))
    } else if (cols == 'interaction') {
      split <- interaction(idents, split)
      cols <- hue_pal()(length(x = levels(x = idents)))
    } else {
      cols <- Col2Hex(cols)
    }
    cols <- Interleave(cols, InvertHex(hexadecimal = cols))
    names(x = cols) <- sort(x = levels(x = split))
  }
  if (same.y.lims && is.null(x = y.max)) {
    y.max <- max(data)
  }
  plots <- lapply(
    X = features,
    FUN = function(x) {
      return(SingleExIPlot(
        type = type,
        data = data[, x, drop = FALSE],
        idents = idents,
        split = split,
        sort = sort,
        y.max = y.max,
        adjust = adjust,
        cols = cols,
        pt.size = pt.size,
        log = log,
        show.quantiles = show.quantiles,
        quantiles = quantiles
      ))
    }
  )
  if (combine) {
    combine.args <- list(
      'plots' = plots,
      'ncol' = ncol
    )
    combine.args <- c(combine.args, list(...))
    if (!'legend' %in% names(x = combine.args)) {
      combine.args$legend <- 'none'
    }
    plots <- do.call(what = 'CombinePlots', args = combine.args)
  }
  return(plots)
}

#' Single cell violin plot
#'
#' Draws a violin plot of single cell data (gene expression, metrics, PC
#' scores, etc.)
#'
#' @inheritParams GCD.RidgePlot
#' @param pt.size Point size for geom_violin
#' @param split.by A variable to split the violin plots by,
#' @param adjust Adjust parameter for geom_violin
#'
#' @return A ggplot object
#'
#' @export
#'
#'
#' @examples
#' VlnPlot(object = pbmc_small, features = 'PC1')
#' VlnPlot(object = pbmc_small, features = 'LYZ', split.by = 'groups')
#'
GCD.VlnPlot <- function(
  object,
  features,
  cols = NULL,
  pt.size = 1,
  idents = NULL,
  sort = FALSE,
  assay = NULL,
  group.by = NULL,
  split.by = NULL,
  adjust = 1,
  y.max = NULL,
  same.y.lims = FALSE,
  log = FALSE,
  ncol = NULL,
  combine = TRUE,
  slot = 'data',
  show.quantiles = T,
  quantiles = 2,
  ...
) {
  return(ExIPlot(
    object = object,
    type = 'violin',
    features = features,
    idents = idents,
    ncol = ncol,
    sort = sort,
    assay = assay,
    y.max = y.max,
    same.y.lims = same.y.lims,
    adjust = adjust,
    pt.size = pt.size,
    cols = cols,
    group.by = group.by,
    split.by = split.by,
    log = log,
    slot = slot,
    combine = combine,
    show.quantiles = show.quantiles,
    quantiles = quantiles,
    ...
  ))
}
