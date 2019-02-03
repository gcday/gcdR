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
Palettes <- function(RET, type.use = 1) {
  library(scales)
  library(colorspace)
  print(type.use)
  print(class(type.use))
  if (type.use == 1) {
    return(hue_pal()(length(levels(RET@seurat))))
  } else if (type.use == 2) {
    return(rainbow_hcl(length(levels(RET@seurat)), c = 90, l = 80))
  } else if (type.use == 3) {
    return(rainbow_hcl(length(levels(RET@seurat)), c = 100, l = 80))
  } else {
    return(rainbow(length(levels(RET@seurat))))
  }
  # col.palettes <- c(hue_pal()(length(levels(RET@seurat))),
                      # rainbow_hcl(length(levels(RET@seurat)), c = 90, l = 80),
                      # rainbow_hcl(length(levels(RET@seurat)), c = 100, l = 80),
                      # rainbow(length(levels(RET@seurat))))
  # return(col.palettes[as.integer(type.use)])
}