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


