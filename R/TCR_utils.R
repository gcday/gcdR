#' Uses expression gating to label T cells.
#'
#'
#' @param object Seurat object
#' @param CD20_cutoff threshold value for CD20+ (default is 0)
#' @param CD3E_cutoff threshold value for CD3+ (default is 0)
#' 
#' @return Seurat object with object$CD4.CD8 filled in 
#'
#' @examples
#' annotateCD4.CD8.TCells(object)
#'
#' @export
#' 
annotateCD3.CD20.cells <- function(object, 
                                   CD20_cutoff = 0, 
                                   CD3E_cutoff = 0) {
  CD3E_expr <- FetchData(object = object, vars = c("CD3E"))
  CD3E.pos <- which(x = CD3E_expr > CD3E_cutoff)
  
  CD20_expr <- FetchData(object = object, vars = c("MS4A1"))
  CD20.pos <- which(x = CD20_expr > CD20_cutoff)
  CD3E.CD20.pos <- intersect(CD3E.pos, CD20.pos)
 
  object$CD3.CD20 <- "CD3E- CD20-"
  object$CD3.CD20[CD3E.pos] <- "CD3E+ CD20-"
  object$CD3.CD20[CD20.pos] <- "CD3E- CD20+"
  object$CD3.CD20[CD3E.CD20.pos] <- "CD3E+ CD20+"
  return(object)
}
#' Uses expression gating to label T cells.
#'
#'
#' @param object Seurat object
#' @param CD4_cutoff threshold value for CD4+
#' @param CD8A_cutoff threshold value for CD8+
#' @param CD3D_cutoff threshold value for CD3+ 
#' 
#' @return Seurat object with object$CD4.CD8 filled in 
#'
#' @examples
#' annotateCD4.CD8.TCells(object)
#'
#' @export
#' 
annotateCD4.CD8.TCells <- function(object, 
                                   CD4_cutoff = 0.5, 
                                   CD8A_cutoff = 0.5, 
                                   CD3D_cutoff = 0.5) {
  CD3D_expr <- FetchData(object = object, vars = c("CD3D"))
  CD3D.pos <- which(x = CD3D_expr > CD3D_cutoff)
  
  CD3D_expr <- FetchData(object = object, vars = c("CD3D"))
  CD3D.pos <- which(x = CD3D_expr > CD3D_cutoff)
  
  CD4_expr <- FetchData(object = object, vars = c("CD4"))
  CD4.pos <- which(x = CD4_expr > CD4_cutoff)
  
  CD8A_expr <- FetchData(object = object, vars = c("CD8A"))
  CD8A.pos <- which(x = CD8A_expr > CD8A_cutoff)
  CD4.pos <- intersect(CD4.pos, CD3D.pos)
  CD8A.pos <- intersect(CD8A.pos, CD3D.pos)
  CD4.CD8.pos <- intersect(CD4.pos, CD8A.pos)
  object$CD4.CD8 <- "CD3-"
  object$CD4.CD8[CD3D.pos] <- "CD3+ DN"
  object$CD4.CD8[CD4.pos] <- "CD4"
  object$CD4.CD8[CD8A.pos] <- "CD8"
  object$CD4.CD8[CD4.CD8.pos] <- "CD4.CD8"
  return(object)
}
#' Applies TCR annotation to Seurat object
#'
#'
#' @param object Seurat object
#' @param barcode_to_clonotype list mapping barcode to clonotypes
#' @param barcode_to_TRAV list mapping barcodes to TRAV genes
#' @param barcode_to_TRBV list mapping barcodes to TRBV genes
#' 
#' @return Seurat object with object$TCR_clonotype, object$TRA_V, 
#' and object$TRB_V filled in 
#'
#' @examples
#' addTCRToSeurat(object, barcode_to_clonotype)
#'
#' @export
#' 
addTCRToSeurat <- function(object, barcode_to_clonotype,
                           barcode_to_TRAV, barcode_to_TRBV) {
  object$TCR_clonotype <- sapply(X = rownames(object[[]]), 
                                     FUN = function(x) {
                                       if (x %in% names(barcode_to_clonotype)) {
                                         return(barcode_to_clonotype[[x]])
                                       } else {
                                         return("none")
                                       }
                                     })
  object$TRA_V <- sapply(X = rownames(object[[]]), 
                             FUN = function(x) {
                               if (x %in% names(barcode_to_TRAV)) {
                                 return(barcode_to_TRAV[[x]])
                               } else {
                                 return("none")
                               }
                             })
  object$TRB_V <- sapply(X = rownames(object[[]]), 
                             FUN = function(x) {
                               if (x %in% names(barcode_to_TRBV)) {
                                 return(barcode_to_TRBV[[x]])
                               } else {
                                 return("none")
                               }
                             })
  return(object)
}

#' Pair of dimensional reduction plots highlighting cells of a given clonotype 
#' @importFrom cowplot plot_grid
#' @importFrom Seurat DimPlot NoLegend
#' @export
plotClonotypeIntVsSCT <- function(object1, prefix1, object2, prefix2, clonotype, clonotype_prefix, clonotype_field = "TCR_clonotype") {
  cells.highlight.1 <- list(
    rownames(object1[[]])[which(object1[[clonotype_field]][[clonotype_field]] == paste0(clonotype_prefix, clonotype))])
  names(cells.highlight.1) <- paste("Clonotype", clonotype)
  cells.highlight.2 <- list(
    rownames(object2[[]])[which(object2[[clonotype_field]][[clonotype_field]] == paste0(clonotype_prefix, clonotype))])
  names(cells.highlight.2) <- paste("Clonotype", clonotype)
  plot.1 <- suppressMessages(DimPlot(object = object1, 
                                     cells.highlight = cells.highlight.1,
                                     sizes.highlight = 0.75) + 
                               labs(title = paste("Clonotype", clonotype, prefix1)) + 
                               NoLegend() + 
                               scale_color_manual(values = c("#bcbcbc", "red")))
  plot.2 <- suppressMessages(DimPlot(object = object2, 
                                     cells.highlight = cells.highlight.2,
                                     sizes.highlight = 0.75) + 
                               labs(title = paste("Clonotype", clonotype, prefix2)) + 
                               NoLegend() + 
                               scale_color_manual(values = c("#bcbcbc", "red")))
  return(plot_grid(plot.1, plot.2))
}

#' Dimensional reduction plot highlighting cells of a given clonotype based on their CD4.CD8 status
#' 
#' @importFrom Seurat DimPlot NoLegend

#' @export

plotClonotypeCD4.CD8 <- function(object, title, clonotype, clonotype_prefix, ...) {
  clonotype.cells <- object[[]][which(object$TCR_clonotype == paste0(clonotype_prefix, 
                                                                     clonotype)),]
  # CD4.CD8.values <- c("CD3-", "CD3+ DN", "CD4", "CD4.CD8", "CD8")
  CD4.CD8.values <- c("CD8", "CD4", "CD4.CD8", "CD3-", "CD3+ DN")
  
  cells.highlight <- list()
  color.pal <- pal_nejm()(length(CD4.CD8.values))
  cols.highlight <- list()
  for (i in 1:length(CD4.CD8.values)) {
    val <- CD4.CD8.values[i]
    matching.cells <- rownames(clonotype.cells)[which(clonotype.cells$CD4.CD8 == val)]
    if (length(matching.cells) > 0) {
      cells.highlight[[paste0(val, " (", length(matching.cells), " cells)" )]] <- matching.cells
      cols.highlight[[paste0(val, " (", length(matching.cells), " cells)" )]] <-  color.pal[i]
    }
  }
  
  return(suppressMessages(DimPlot(object = object, 
                                  cells.highlight = cells.highlight, 
                                  cols.highlight = c(unlist(cols.highlight)),
                                  ...) + 
                            labs(title = paste("Clonotype", clonotype, title)) +
                            scale_color_manual(values = c(Unselected = "#bcbcbc", unlist(cols.highlight)))))
}

#' Dimensional reduction plot highlighting cells of a given clonotype based on their CD4.CD8 status
#' 
#' @importFrom Seurat DimPlot NoLegend

#' @export

plotClonotypeTimePoints <- function(object, title, clonotype, clonotype_prefix, 
                                    timepoint.field = "timepoint", site.field = "site",
                                    site.levels = NULL,
                                    timepoint.levels = NULL,
                                    ...
                                    ) {
  clonotype.cells <- object[[]][which(object$TCR_clonotype == paste0(clonotype_prefix, 
                                                                     clonotype)),]
  timepoint.levels <- timepoint.levels %||% levels(as.factor(object[[timepoint.field]][[timepoint.field]]))
  # site.levels <- site.levels %||% levels(as.factor(object[[site.field]][[site.field]]))
  
  # CD4.CD8.values <- c("CD3-", "CD3+ DN", "CD4", "CD4.CD8", "CD8")
  # CD4.CD8.values <- c("CD8", "CD4", "CD4.CD8", "CD3-", "CD3+ DN")
  
  cells.highlight <- list()
  color.pal <- pal_nejm()(length(timepoint.levels))
  cols.highlight <- list()
  for (i in 1:length(timepoint.levels)) {
    val <- timepoint.levels[i]
    matching.cells <- rownames(clonotype.cells)[which(clonotype.cells[[timepoint.field]] == val)]
    if (length(matching.cells) > 0) {
      cells.highlight[[paste0(val, " (", length(matching.cells), " cells)" )]] <- matching.cells
      cols.highlight[[paste0(val, " (", length(matching.cells), " cells)" )]] <-  color.pal[i]
    }
  }
  
  return(suppressMessages(DimPlot(object = object, 
                                  cells.highlight = cells.highlight, 
                                  cols.highlight = c(unlist(cols.highlight)),
                                  shape.by = site.field,
                                  ...) + 
                            labs(title = paste("Clonotype", clonotype, title)) +
                            scale_color_manual(values = c(Unselected = "#bcbcbc", unlist(cols.highlight)))))
}

#' Splits cell name barcodes into usable dataframe.
#' 
#' @importFrom tidyr separate
#' @importFrom dendextend get_leaves_attr 
#' 
#' @return dataframe with relevant fields: notably $fixed_cellnames,
#' containing properly formatted cell barcodes. 
#' @export
#' 
#' 
separateCellnamesVDJ <- function(tree = NULL, cellnames = NULL) {
  if (is.null(cellnames)) {
    tree <- as.dendrogram(tree)
    labels.tree <- get_leaves_attr(as.dendrogram(tree), attribute = "label")
    cellnames <- grep("GERM", x = labels.tree, value = T, invert = T)
  }
  tumor.df <- data.frame(cellnames)
  tumor.df <- tidyr::separate(tumor.df, col = 1, 
                              into = c("patient", "site", "timepoint", "barcode", "barcode_suffix", "contig", "contig_num"), 
                              remove = F, fill = "right", extra = "drop")
  tumor.df$patient.site <- paste0(tumor.df$patient, "_", tumor.df$site)
  tumor.df$patient.site.timepoint <- paste0(tumor.df$patient, "_", tumor.df$site, "_", tumor.df$timepoint)
  tumor.df$patient.timepoint <- paste0(tumor.df$patient, "_", tumor.df$timepoint)
  tumor.df$site.timepoint <- paste0(tumor.df$site, "_", tumor.df$timepoint)
  tumor.df$fixed_cellnames <- paste0(tumor.df$patient, "_", tumor.df$site, "_", tumor.df$timepoint, "-", tumor.df$barcode) # "-", tumor.df$barcode_suffix)
  tumor.df$fixed_contig_names <- gsub(pattern = " ", replacement = "_", x = tumor.df$cellnames)
    # paste0(tumor.df$patient, "_", tumor.df$site, "_", tumor.df$timepoint, "-", tumor.df$barcode, "_", 
                                        # tumor.df$contig, "_", tumor.df$contig_num) # "-", tumor.df$barcode_suffix)
  
  return(tumor.df)
}

#' @export 
fixIgphymlTreeLeaves <- function(tree, new.labels = NULL) {
  tree <- as.dendrogram(tree)
  labels.tree <- get_leaves_attr(as.dendrogram(tree), attribute = "label")
  
  # tumor.df <- separateCellnamesVDJ(tree)
  if (!is.null(new.labels)) {
    germ.pos <- grep("GERM", x = labels.tree)
    new.labels <- c(new.labels[1:germ.pos],
                    labels.tree[germ.pos],
                    new.labels[germ.pos:length(new.labels)])
  } else {
    new.labels <- gsub(pattern = " ", replacement = "_", x = labels.tree)
  }
  
  # new.labels <- c(new.labels[]labels.tree[germ.pos]
  # new.labels <-  c(dend.cex[1:grep("GERM", x = labels.tree) - 1],
  #                  leaf.cex,
  #                  dend.cex[grep("GERM", x = labels.tree):length(dend.cex)])
  # # tree <- set_labels(tree, tumor.df$fixed_contig_names)
  tree <- set_labels(tree, new.labels)
  return(tree)
}

#' Dendrogram plot
#' 
#' @importFrom dendextend get_leaves_attr get_branches_heights colored_bars
#' @importFrom ggsci pal_jama pal_aaas pal_jco
#' @importFrom RColorBrewer brewer.pal
#' 
#' @export
#' 
igphyml.treePlot <- function(tree, title.text = "", cells.highlight = NULL,
                             cols.highlight = NULL,
                             highlight.name = "Highlight",
                             highlight.legend = T,
                             highlight.legend.pos = "topright",
                             site.colors = NULL,
                             timepoint.colors = NULL,
                             ymin = NULL, ymax = NULL,
                             highlight.leaf = F,
                             leaf.pch = 17,
                             leaf.cex = 2,
                             show.pbmc = T,
                             show.timepoint = F,
                             hide.yaxt = T,
                             use.ylim = T,
                             legend.location = "topleft",
                             bars.yscale = NULL,
                             germ.icon = T) {
  tree <- as.dendrogram(tree)
  
  labels.tree <- get_leaves_attr(as.dendrogram(tree), attribute = "label")
  new.cex <- rep(0, length(labels.tree))
  if (any(grepl("GERM", x = labels.tree))) {
    new.cex[grep("GERM", x = labels.tree)] <- leaf.cex
  }
  dend.pch <- rep(leaf.pch, length(labels.tree))
  tree <- set(tree, "leaves_pch", c(leaf.pch)) %>%  # node point type
    set("leaves_cex", new.cex) %>%  # node point size
    set("leaves_col", c("blue"))
  
  # included.tumor <- grep("GERM", x = labels.tree, value = T, invert = T)
  tumor.df <- separateCellnamesVDJ(tree = tree)
  # brewer.pal(n, name)
  # site.colors <- pal_jama()(3)
  site.colors <- site.colors %||% RColorBrewer::brewer.pal(3, "Set1")
  # patient.colors <- ggsci::pal_npg()(8)
  if (is.null(names(site.colors))) names(site.colors) <- c("lesA", "lesB", "PBMC")
  if (!show.pbmc) site.colors["PBMC"] <- "#D3D3D3"
  timepoint.colors <- timepoint.colors %||% RColorBrewer::brewer.pal(3, "Dark2")
  if(is.null(names(timepoint.colors))) names(timepoint.colors) <- c("pre", "wk2", "wk6")

  
  
  if (is.character(x = cells.highlight)) {
    cells.highlight <- list(cells.highlight)
  }
  else if (is.data.frame(x = cells.highlight) || !is.list(x = cells.highlight)) {
    cells.highlight <- as.list(x = cells.highlight)
  }
  if (length(cells.highlight) != 0 & is.null(x = names(x = cells.highlight))) {
    names(cells.highlight) <- paste0("Group_", 1L:length(x = cells.highlight))
  }
  # names(cells.highlight) <- names.highlight
  
  cols.highlight <- cols.highlight %||% ggsci::pal_igv()(length(cells.highlight))
  dend.colors.site <- list()
  dend.colors.timepoint <- list()
  dend.colors.highlight <- list()
  dend.cex <- list()
  names(cols.highlight) <- names(cells.highlight)
  
  for (i in 1:nrow(tumor.df)) {
    dend.colors.site[i] <- site.colors[[tumor.df[i, "site"]]]
    dend.colors.timepoint[i] <- timepoint.colors[[tumor.df[i, "timepoint"]]]

    if (length(cells.highlight) != 0) {
      highlight.group <- which(sapply(cells.highlight, `%in%`, 
                                    x = tumor.df[i, "fixed_cellnames"]))

      dend.colors.highlight[i] <- ifelse(length(highlight.group) != 0,
                                         cols.highlight[highlight.group[1]],
                                         "#D3D3D3")
      dend.cex[i] <- ifelse(length(highlight.group) != 0,
                            leaf.cex,
                            0)
      # dend.colors.highlight[i] <- ifelse(tumor.df[i, "fixed_cellnames"] %in% cells.highlight,
      #                                    "#FF0000",
      #                                    "#808080")
      
    }
    # dend.colors.patient[i] <- patient.colors[[tumor.df[i, "patient"]]]
  }
  dend.pch <- rep(leaf.pch, length(labels.tree))
  if (any(grepl("GERM", x = labels.tree))) {
    germ.index <- grep("GERM", x = labels.tree)
    # edge case: when germline is the very last/first cell in the object
    if (germ.index == length(labels.tree)) {
      dend.colors.site <- c(dend.colors.site, "#808080")
      dend.colors.timepoint <- c(dend.colors.timepoint, "#808080")
      dend.colors.highlight <- c(dend.colors.highlight, "#808080")
      dend.cex <- c(dend.cex, leaf.cex)
    } else if (germ.index == 1) {
      dend.colors.site <- c("#808080", dend.colors.site)
      dend.colors.timepoint <- c("#808080", dend.colors.timepoint)
      dend.colors.highlight <- c("#808080", dend.colors.highlight)
      dend.cex <- c(leaf.cex, dend.cex)
    } else {
      dend.colors.site <- c(dend.colors.site[1:germ.index-1],
                            "#808080",
                            dend.colors.site[germ.index:length(dend.colors.site)])
      dend.colors.timepoint <- c(dend.colors.timepoint[1:germ.index - 1],
                                 "#808080",
                                 dend.colors.timepoint[germ.index:length(dend.colors.timepoint)])
      
      dend.colors.highlight <- c(dend.colors.highlight[1:germ.index - 1],
                                 "#808080",
                                 dend.colors.highlight[germ.index:length(dend.colors.highlight)])
      dend.cex <- c(dend.cex[1:germ.index - 1], 
                    leaf.cex,
                    dend.cex[germ.index:length(dend.cex)])
    }
    # message("present")
    
  }
  if (any(grepl("GERM", x = labels.tree))) {
    if (germ.icon) {
      dend.pch[grep("GERM", x = labels.tree)] <- 8
    }
  }
  
  if (highlight.leaf) {
    # unique(unlist(cells.highlight))
    # new.cex <- rep(0, length(labels.tree))
    # set("leaves_col", c("blue"))
    tree <- set(tree, "leaves_pch", c(leaf.pch)) %>%  # node point type
      set("leaves_cex", unlist(dend.cex)) %>%  # node point size
      set("leaves_col", unlist(dend.colors.highlight)) %>% 
      set("leaves_pch", dend.pch)
  }
  # tree <- assign_values_to_leaves_nodePar(tree, value = dend.colors.site, nodePar = "leaves_col")
  
  labels(tree) <- rep("", length(labels.tree))
  
  ymin <- ymin %||% min(get_branches_heights(tree))
  ymax <- ymax %||% (max(get_branches_heights(tree)) + 0.05 * max(get_branches_heights(tree)))
  if (use.ylim) {
    plot(tree,
         ylim = c(ymin, ymax),
         yaxt = ifelse(test = hide.yaxt,
                       yes = "n", no = "s"), ann=FALSE)
  } else {
    plot(tree,
         yaxt = ifelse(test = hide.yaxt,
                       yes = "n", no = "s"), ann=FALSE)
  }
  
  if (!show.pbmc) site.colors <- site.colors[1:2]
  bar.colors <- cbind(dend.colors.site)
  # print(dend.colors.site)
  # print(length(dend.colors.site))
  bar.names <- "Site"
  if (show.timepoint) {
    bar.colors <- cbind(bar.colors, dend.colors.timepoint)
    bar.names <- c(bar.names, "Timepoint")
  }
  if (length(cells.highlight) != 0) {
    if (!highlight.leaf) {
      bar.colors <- cbind(bar.colors, dend.colors.highlight)
      bar.names <- c(bar.names, highlight.name)
    }
  }
  
  if (legend.location == "topleft") {
    legend.site <- legend("topleft", legend = names(site.colors), fill = site.colors, title = "Site")
    if (show.timepoint) {
      legend.timepoint <- legend(x = legend.site$rect$left + legend.site$rect$w, 
                                 y = legend.site$rect$top, 
                                 legend = names(timepoint.colors), 
                                 fill = timepoint.colors, title = "Timepoint"
      )
    } else {
      legend.timepoint <- legend.site
    }
    
    if (length(cells.highlight) != 0) {
      
      if (highlight.legend) {
        legend(x = legend.timepoint$rect$left + legend.timepoint$rect$w, 
               y = legend.timepoint$rect$top, 
               legend = names(cells.highlight), 
               col = cols.highlight, title = highlight.name,
               pch = leaf.pch, pt.cex = leaf.cex)
      }
      
      
      
      
      # bars.cols <- cbind(dend.colors.site, dend.colors.highlight)
      
      # if (!highlight.leaf) {
      #   colored_bars(cbind(dend.colors.site, dend.colors.highlight), tree, rowLabels = c(tumor.df[1, "patient"], highlight.name))
      # } else {
      #   colored_bars(dend.colors.site, tree, rowLabels = c(tumor.df[1, "patient"]))
      # }
    }
  }
  # } else {
  #   legend.site <- legend("topright", legend = names(site.colors), fill = site.colors, title = "Site")
  #   if (show.timepoint) {
  #     legend.timepoint <- legend(x = legend.site$rect$left + legend.site$rect$w, 
  #                                y = legend.site$rect$top, 
  #                                legend = names(timepoint.colors), 
  #                                fill = timepoint.colors, title = "Timepoint"
  #     )
  #   } else {
  #     legend.timepoint <- legend.site
  #   }
  #   bar.colors <- cbind(dend.colors.site)
  #   # print(dend.colors.site)
  #   # print(length(dend.colors.site))
  #   bar.names <- "Site"
  #   if (show.timepoint) {
  #     bar.colors <- cbind(bar.colors, dend.colors.timepoint)
  #     bar.names <- c(bar.names, "Timepoint")
  #   }
  #   if (length(cells.highlight) != 0) {
  #     
  #     if (highlight.legend) {
  #       legend(x = legend.timepoint$rect$left + legend.timepoint$rect$w, 
  #              y = legend.timepoint$rect$top, 
  #              legend = names(cells.highlight), 
  #              col = cols.highlight, title = highlight.name,
  #              pch = leaf.pch, pt.cex = leaf.cex)
  #     }
  #     
  #     if (!highlight.leaf) {
  #       bar.colors <- cbind(bar.colors, dend.colors.highlight)
  #       bar.names <- c(bar.names, highlight.name)
  #     }
  #     
  #     
  #     
  #     # bars.cols <- cbind(dend.colors.site, dend.colors.highlight)
  #     
  #     # if (!highlight.leaf) {
  #     #   colored_bars(cbind(dend.colors.site, dend.colors.highlight), tree, rowLabels = c(tumor.df[1, "patient"], highlight.name))
  #     # } else {
  #     #   colored_bars(dend.colors.site, tree, rowLabels = c(tumor.df[1, "patient"]))
  #     # }
  #   }
  # }
  # 
  # } else {
  #     # colored_bars(dend.colors.site, tree, rowLabels = c(tumor.df[1, "patient"]))
  # }
  # print(length(bar.colors))
  # print(length(labels(tree)))
  if (!is.null(bars.yscale)) {
    colored_bars(bar.colors, rowLabels = bar.names, yscale = bars.yscale)
  } else {
    colored_bars(bar.colors, rowLabels = bar.names)
  }
  title(title.text)
  # return(bar.colors)
  
}

#' Adds glycosylation site amino acid position and number to a clones data table.
#' 
#' @importFrom seqinr translate c2s
#' @importFrom stringr str_count str_locate_all 
#' @importFrom future.apply future_sapply 
#' 
#' @export
#' 

countGlycosylationSites <- function(clones.tbl, new.mode = F) {
  if (!new.mode) {
    clones.tbl$NOGAP_SEQUENCE <- future_sapply(clones.tbl$SEQUENCE_IMGT,
                                               FUN = function(igh.seq) {
                                                 return(c2s(translate(gsub(pattern = "\\.",
                                                                           replacement = "", 
                                                                           x = s2c(igh.seq)))))
                                               })
    
    clones.tbl$GLYCOSYLATION_SITES <- future_sapply(clones.tbl$NOGAP_SEQUENCE, 
                                                    FUN = function(nogap.seq) {
                                                      num.glycosylation.sites <- stringr::str_count(nogap.seq, pattern = "N[[:upper:]][ST]")
                                                      if (num.glycosylation.sites != 0) {
                                                        glycosylation.sites <- stringr::str_locate_all(nogap.seq, pattern = "N[[:upper:]][ST]")
                                                        return(paste(unlist(glycosylation.sites[[1]][,1]), sep = ",", collapse = ", "))
                                                      } else {
                                                        return("")
                                                      }})
    clones.tbl$NUM_GLYCOSYLATION_SITES <- future_sapply(clones.tbl$NOGAP_SEQUENCE, 
                                                        FUN = function(nogap.seq) {stringr::str_count((nogap.seq), pattern = "N[[:upper:]][ST]")})
    
  } else {
    out.list <- future_sapply(clones.tbl$SEQUENCE_IMGT, FUN = function(igh.seq) {
      nogap.seq <- c2s(translate(gsub(pattern = "\\.",
                                      replacement = "", 
                                      x = s2c(igh.seq))))
      num.glycosylation.sites <- stringr::str_count(nogap.seq, pattern = "N[[:upper:]][ST]")
      if (num.glycosylation.sites != 0) {
        glycosylation.sites <- stringr::str_locate_all(nogap.seq, pattern = "N[[:upper:]][ST]")
        glycos.str <- paste(unlist(glycosylation.sites[[1]][,1]), sep = ",", collapse = ", ")
      } else {
        glycos.str <- ""
      }
      return(c(nogap.seq, num.glycosylation.sites, glycos.str))
    })
    clones.tbl$NOGAP_SEQUENCE <- unlist(out.list[1:length(clones.tbl$SEQUENCE_IMGT) * 3 - 2])
    clones.tbl$NUM_GLYCOSYLATION_SITES <- unlist(out.list[1:length(clones.tbl$SEQUENCE_IMGT) * 3 - 1])
    clones.tbl$GLYCOSYLATION_SITES <- unlist(out.list[1:length(clones.tbl$SEQUENCE_IMGT) * 3])
  }
    
  
  
  return(clones.tbl)
}

#' Adds glycosylation site amino acid position and number to a clones data table.
#' 
#' @importFrom progressr progressor 
#' 
#' @export
#' 


parseBuildTreesLogfile <- function(file.path, patient = NULL) {
  logfile <- readLines(file.path)
  
  nvec <- length(logfile)
  breaks <- which(! nzchar(logfile))
  nbreaks <- length(breaks)
  if (breaks[nbreaks] < nvec) {
    breaks <- c(breaks, nvec + 1L)
    nbreaks <- nbreaks + 1L
  }
  if (nbreaks > 0L) {
    chunks <- mapply(function(a,b) paste(logfile[a:b], collapse = "\n"),
                     c(1L, 1L + breaks[-nbreaks]),
                     breaks - 1L)
  }
  pass.clones <- list()
  parent_to_cloneid <- list()
  p <- progressr::progressor(steps = length(chunks))
  for (i in 1:length(chunks)) {
    # for (i in 1:2) {
    
    # if ("ID")
    fields <- list(ID = "", 
                   CLONE = "", 
                   PASS = FALSE, 
                   END_MASKED = "", 
                   SEQ_IN = "",
                   SEQ_IMGT = "",
                   SEQ_MASKED = "",
                   IN_FRAME = "",
                   MASKED = "",
                   FRAMESHIFTS = "",
                   FAIL = "",
                   COLLAPSETO = "",
                   COLLAPSEFROM = "",
                   DUPLICATE = FALSE)
    # record.lines <- strsplit(x = chunks[i], split = "\n", fixed = T)
    for (line in strsplit(x = chunks[i], split = "\n", fixed = T)[[1]]) {
      # print(line)
      split.line <- strsplit(x = str_trim(line), split = "> ", fixed = T)[[1]]
      split.line[1] <- gsub(pattern = "-", replacement = "_", x = split.line[1])
      if (!split.line[1] %in% names(fields)) {
        print(split.line[1])
      }
      if(split.line[1] %in% c("PASS", "DUPLICATE")) {
        fields[[split.line[1]]] <- as.logical(split.line[2])
      } else {
        fields[[split.line[1]]] <- split.line[2]
      } 
      
      # print(split.line)
      # split.line[1]
      # print()
      # print(str_trim(line))
    }
    p()
    # if (i %% 1000) message(i / length(chunks))
    if (fields$PASS) {
      if (!fields$ID %in% names(pass.clones)) {
        pass.clones[[fields$ID]] <- c(fields$ID)
        parent_to_cloneid[[fields$ID]] <- fields$CLONE
      }
    } else if (fields$DUPLICATE) {
      # if ()
      parent.clone <- gsub(pattern = "Duplication of ", replacement = "", x = fields$FAIL, fixed = T)
      parent.clone <- gsub(pattern = "Collapsed with ", replacement = "", x = parent.clone, fixed = T)
      
      if (parent.clone %in% names(pass.clones)) {
        pass.clones[[parent.clone]] <- c(pass.clones[[parent.clone]], fields$ID)
      } else {
        pass.clones[[parent.clone]] <- c(parent.clone, fields$ID)
        parent_to_cloneid[[parent.clone]] <- fields$CLONE
        
      }
    }
  }
  patient <- patient %||% strsplit(x = names(pass.clones)[1], split = "_")[[1]][[1]]
  # clone <- fields$CLONE
  new_pass_clones <- list()
  new_pass_df <- NULL
  for (i in 1:length(pass.clones)) {
    clone <- parent_to_cloneid[[names(pass.clones)[i]]]
    subclone <- paste(patient, clone, i, sep = "_")
    # new_pass_clones[[subclone]] <- pass.clones[[i]]
    new_pass_clones[[subclone]] <- separateCellnamesVDJ(cellnames = pass.clones[[i]])
    new_pass_clones[[subclone]]$ig_subclone <- subclone
    
    new_pass_clones[[subclone]]$clone <- paste(patient, clone, sep = "_") 
    new_pass_df <- rbind(new_pass_df, new_pass_clones[[subclone]])
  }
  
  return(new_pass_df)
  # return(new_pass_clones)
}