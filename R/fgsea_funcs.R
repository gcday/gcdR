#' Convers Seurat marker output to DE gene ranks
#'
#'
#' @param de.table DE genes table, returned by \code{\link[Seurat]{FindAllMarkers}}
#' 
#' @return ranks to be used for \code{\link{gseaPerCluster}}
#'
#' @examples
#' seuratDEresToGSEARanks(de.table)
#'
#' @export
seuratDEresToGSEARanks <- function(de.table) {
  # de.markers <- de.table
  de.markers <- mutate(de.table, p_val = ifelse(p_val == 0, 1E-300, p_val))
  de.markers <- filter(de.markers, p_val < 0.01)
  # de.markers[!(is.na(de.markers$p_val)) & de.markers$p_val == 0,]$p_val <- 
  # de.markers <- de.markers[!duplicated(de.markers$gene),]
  de.markers <- mutate(de.markers, sign = ifelse(abs(avg_logFC) > 0.3, sign(avg_logFC), sign(pct.1 - pct.2)))
  de.markers <- mutate(de.markers, rank_metric = -log10(p_val) * sign)
  # de.markers <- mutate(de.markers, rank_metric = -log10(p_val) * sign(avg_logFC))

  de.markers <- de.markers[!duplicated(de.markers$rank_metric),]
  de.markers <- de.markers[!duplicated(de.markers$gene),]
  ranks <- data.frame(de.markers$rank_metric)
  ranks <- setNames(ranks$de.markers.rank_metric,
                    de.markers$gene)
  return(ranks)
}

#' Finds enriched pathways in each cluster using fgsea
#'
#'
#' @param all.markers table of DE genes between conditions, returned by \code{\link[Seurat]{FindAllMarkers}}
#' @param dbs pathway databases 
#' 
#' @return list containing ranks and fgsea output for each ident
#'
#' @examples
#' gseaPerCluster(all.markers, pathways)
#'
#' @export
gseaPerCluster <- function(ranks, pathways) {
  library(stats)
  fgseaRes <- list()
  for (ident in names(ranks)) {
  	message(ident)
    fgsea.out <- fgsea::fgsea(pathways = pathways, stats = ranks[[ident]],
                               minSize=15,
                               maxSize=500,
                               nperm=20000, 
                               nproc = parallel::detectCores())
    fgseaRes[[ident]] <- list()
    fgseaRes[[ident]][["ranks"]] <- ranks[[ident]]
    fgseaRes[[ident]][["fgsea"]] <- fgsea.out
  }
  return(fgseaRes)
}
#' Wrapper running fgsea
#'
#'
#' @param RET list containing Seurat object
#' @param pathways.list list of pathway sets to analyze
#' @param prefixes list of prefixes to trim from name of each pathway
#' @param mouse whether to convert gene names prior to fgsea
#' @return list containing ranks and fgsea output for each ident
#'
#' @examples
#' fgseaWrapper(RET, pathways.list)
#'
#' @export
fgseaWrapper <- function(RET, pathways.list, prefixes = NULL, mouse = FALSE) {
  # RET$fgsea <- list(pathways = list(), results = list(), prefixes=list(), ranks=list(), plots=list())
  # RET$plots$fgsea <- list()
  markers.name = "all.markers.quick"
  if ("all.markers.full" %in% names(RET@meta.list)) {
  	markers.name = "all.markers.full"
  }
  if (!markers.name %in% names(RET@meta.list)) {
  	warning("Missing marker genes")
  }
  markers <- RET@meta.list[[markers.name]]
  # if (mouse) {
  #   markers <- convertMouseGeneList(markers)
  # }


  # for (ident in levels(as.factor(RET[[markers]]$cluster))) { 
  #   RET$fgsea$ranks[[ident]] <- seuratDEresToGSEARanks(filter(RET[[markers]], cluster == ident))
  # }
  RET@meta.list$fgsea <- fgseaFromMarkers(markers, pathways.list, prefixes)
  
  # for (i in 1:length(pathways.list)) {
  #   path.name <- names(pathways.list)[i]
  #   pathways <- pathways.list[[path.name]]
  #   RET$fgsea$results[[path.name]] <- gseaPerCluster(RET$fgsea$ranks, c(pathways))
  #   RET$fgsea$pathways[[path.name]] <- pathways
  #   RET$fgsea$prefixes[[path.name]] <- ifelse(is.null(prefixes), "", prefixes[i])
  #   RET$plots$fgsea[[path.name]] <- gseaTopPlots(fgseaRes = RET$fgsea$results[[path.name]],
  #                                                pathways = pathways, 
  #                                                path.name = path.name, 
  #                                                term.prefix = RET$fgsea$prefixes[[path.name]])
  # }
  return(RET)
}
#' Wrapper running fgsea
#'
#'
#' @param markers DE marker table
#' @param pathways.list list of pathway sets to analyze
#' @param prefixes list of prefixes to trim from name of each pathway
#'
#' @return list containing ranks and fgsea output for each ident
#'
#' @examples
#' fgseaFromMarkers(markers, pathways.list)
#'
#' @export
fgseaFromMarkers <- function(markers, pathways.list, prefixes = NULL) {
  fgsea <- list(pathways = list(), results = list(), prefixes=list(), ranks=list(), plots=list())
  for (ident in levels(as.factor(markers$cluster))) { 
    fgsea$ranks[[ident]] <- seuratDEresToGSEARanks(filter(markers, cluster == ident))
  }

  for (i in 1:length(pathways.list)) {
    path.name <- names(pathways.list)[i]
    pathways <- pathways.list[[path.name]]
    fgsea$results[[path.name]] <- gseaPerCluster(fgsea$ranks, c(pathways))
    fgsea$pathways[[path.name]] <- pathways
    fgsea$prefixes[[path.name]] <- ifelse(is.null(prefixes), "", prefixes[i])
    fgsea$plots[[path.name]] <- gseaTopPlots(fgseaRes = fgsea$results[[path.name]],
                                                 pathways = pathways, 
                                                 path.name = path.name, 
                                                 term.prefix = fgsea$prefixes[[path.name]])
  }
  return(fgsea)
}


#' Plots enriched and depleted pathways in each cluster
#'
#'
#' @param fgseaRes list containing ranks and fgsea output for each ident, output of \link{gseaPerCluster}
#' @param path.name name of pathway database 
#' @param pathways list of pathways
#' @param term.prefix prefix to remove from pathway name for printing
#' @param pval.thresh only include pathways with padj below this value
#' @param num.to.print print at most this many pathways 
#' 
#' @return list containing ranks and fgsea output for each ident
#'
#' @examples
#' gseaTopPlots(all.markers, pathways)
#'
#' @export
gseaTopPlots <- function(fgseaRes, path.name, pathways, term.prefix = "GO_", pval.thresh = 0.01, num.to.print = 10) {
  require("fgsea")
  require("dplyr")
  require("grid")
  require("gridExtra")
  fgsea.plots <- list()
  for (ident in names(fgseaRes)) {
    fgs <- fgseaRes[[ident]]
    sigPathways <- filter(fgs$fgsea, padj < pval.thresh) %>% arrange(pval)
    collapsedSigPathways <- collapsePathways(sigPathways, pathways, fgs$ranks)
    # mainPathways <- sigPathways
    mainPathways <- filter(sigPathways, pathway %in% collapsedSigPathways$mainPathways)
    topPathwaysUp <- filter(mainPathways, NES > 0) %>%
      arrange(desc(NES)) %>% top_n(n = -num.to.print, wt = pval)
    topPathwaysDown <- filter(mainPathways, NES < 0)%>%
      arrange(NES) %>% top_n(n = -num.to.print, wt = pval)
    fgsea.plots[[ident]] <- list()
    fgsea.plots[[ident]]$up <- arrangeGrob(top=textGrob(paste("Up pathways from", path.name), 
                                                        gp = gpar(fontvsize=15, fontface="bold")),
                                           GCD.plotGseaTable(pathways[topPathwaysUp$pathway], 
                                                             fgs$ranks, fgs$fgsea,
                                                             term.prefix = term.prefix, 
                                                             num.to.print = num.to.print,
                                                             gseaParam = 1))
    
    fgsea.plots[[ident]]$down <- arrangeGrob(top=textGrob(paste("Down pathways from", 
                                                                path.name), 
                                                          gp = gpar(fontvsize=15, fontface="bold")),
                                             GCD.plotGseaTable(pathways[topPathwaysDown$pathway], 
                                                               fgs$ranks, fgs$fgsea, 
                                                               term.prefix=term.prefix, 
                                                               num.to.print = num.to.print,
                                                               gseaParam = 1))
  }
  return(fgsea.plots)
}
#' Plots enriched and depleted pathways for a cluster 
#'
#'
#' @param fgseaRes fgsea results
#' @param stats fgsea ranks
#' @param pathways list of pathways
#' @param gseaParam param for fgsea algorithm
#' @param term.prefix prefix to remove from pathway name for printing
#' @param pval.thresh only include pathways with padj below this value
#' @param num.to.print print at most this many pathways
#' @param colwidths ex. c(10, 3, 0.8, 1.2)
#' 
#' @return list containing ranks and fgsea output for each ident
#'
#' @examples
#' GCD.plotGseaTable(all.markers, pathways)
#'
#' @export
GCD.plotGseaTable <- function (pathways, stats, fgseaRes, gseaParam = 1, term.prefix, num.to.print,
                               colwidths = c(10, 4, 0.8, 1.6))
{
  require("grid")
  require("gridExtra")
  pathways <- head(pathways, num.to.print)
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    list(textGrob(substr(gsub(term.prefix,"",pn),1,55), just = "right", x = unit(0.95, "npc")), 
         ggplot() + geom_segment(aes(x = p, xend = p, y = 0, 
                                     yend = statsAdj[p]), size = 0.2) + 
           scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) + 
           scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
           xlab(NULL) + ylab(NULL) + 
           theme(panel.background = element_blank(), axis.line = element_blank(), 
                 axis.text = element_blank(), axis.ticks = element_blank(), 
                 panel.grid = element_blank(), axis.title = element_blank(), 
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0, "null"), 4)), 
         textGrob(sprintf("%.2f", annotation$NES)), 
         textGrob(sprintf("%.2e", annotation$pval)))
  })
  pads <- list()
  for (i in 0:(num.to.print-length(ps))) {
    pads[[i + 1]] <- list(nullGrob())
  }
  rankPlot <- ggplot() + geom_blank() + 
    scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) + 
    scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
    xlab(NULL) + ylab(NULL) + 
    theme(panel.background = element_blank(), axis.line = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), 
          panel.grid = element_blank(), 
          axis.title = element_blank(), 
          plot.margin = unit(c(0, 0, 0.5, 0), "npc"), panel.spacing = unit(c(0, 0, 0, 0), "npc"))
  
  arrangeGrob(grobs = c(lapply(c("Pathway", "Gene ranks", "NES", "pval"), textGrob, 
                               gp = gpar(fontface="bold")),
                        unlist(ps, recursive = FALSE), 
                        list(nullGrob(), rankPlot, nullGrob(), nullGrob()), 
                        unlist(pads, recursive=FALSE)), ncol = 4, widths = colwidths)
}


#' Prints fgsea plots for Seurat object
#'
#'
#' @param RET list containing Seurat object
#' @param redraw whether plots should be redrawn before plotting
#' @param path.dbs.to.check if given, limit plotting to only these pathway sets
#'
#' @return list containing ranks and fgsea output for each ident
#'
#' @examples
#' printGseaPlots(RET, pathways.list)
#'
#' @export
printGseaPlots <- function(RET, redraw = F, path.dbs.to.check = NULL, ...) {
  require("grid")
  require("gridExtra")
  if (is.null(path.dbs.to.check)) {
    path.dbs.to.check <- names(RET@meta.list$fgsea$plots)
  }
  for (pathway in path.dbs.to.check) {
    if (redraw) {
      RET@meta.list$plots$fgsea[[pathway]] <- gseaTopPlots(fgseaRes = RET@meta.list$fgsea$results[[pathway]],
                                                 pathways = RET@meta.list$fgsea$pathways[[pathway]],
                                                 path.name = pathway,
                                                 term.prefix = RET@meta.list$fgsea$prefixes[[pathway]], ...)
    }
    for (ident in names(RET@meta.list$fgsea$plots[[pathway]])) {
      grid.arrange(top=textGrob(ident, gp = gpar(fontvsize=20, fontface="bold")),
                   RET@meta.list$fgsea$plots[[pathway]][[ident]]$up,
                   RET@meta.list$fgsea$plots[[pathway]][[ident]]$down)
    }
  }
}

#' Prints fgsea plots for Seurat object
#'
#'
#' @param RET list containing Seurat object
#' @param redraw whether plots should be redrawn before plotting
#' @param path.dbs.to.check if given, limit plotting to only these pathway sets
#'
#' @return list containing ranks and fgsea output for each ident
#'
#' @examples
#' printGseaPlots(RET, pathways.list)
#'
#' @export
printGseaPlotsFGSEA <- function(fgsea, redraw = F, path.dbs.to.check = NULL, ...) {
  require("grid")
  require("gridExtra")
  if (is.null(path.dbs.to.check)) {
    path.dbs.to.check <- names(fgsea$plots)
  }
  for (pathway in path.dbs.to.check) {
    if (redraw) {
      fgsea$plots[[pathway]] <- gseaTopPlots(fgseaRes = fgsea$results[[pathway]],
                                                 pathways = fgsea$pathways[[pathway]],
                                                 path.name = pathway,
                                                 term.prefix = fgsea$prefixes[[pathway]], ...)
    }
    for (ident in names(fgsea$plots[[pathway]])) {
      grid.arrange(top=textGrob(ident, gp = gpar(fontvsize=20, fontface="bold")),
                   fgsea$plots[[pathway]][[ident]]$up,
                   fgsea$plots[[pathway]][[ident]]$down)
    }
  }
}

#' Correlates genes in leading edges of enriched genes
#'
#'
#' @param RET list containing Seurat object
#' @param pval.thresh significance threshold for plotting
#' @param path.dbs.to.check if given, limit plotting to only these pathway sets
#'
#' @return list containing ranks and fgsea output for each ident
#'
#' @examples
#' correlateFgseaLeadingEdges(RET, pval.thresh)
#'
#' @export
correlateFgseaLeadingEdges <- function(RET, pval.thresh = 0.001, path.dbs.to.check = NULL, mouse = FALSE) {
  require("dplyr")
  leading.edges <- list()
  if(is.null(path.dbs.to.check)) {
    path.dbs.to.check <- names(RET@meta.list$fgsea$plots)
  }
  for (pathway in path.dbs.to.check) {
    fgseaRes <- RET@meta.list$fgsea$results[[pathway]]
    for (ident in names(fgseaRes)) {
      # print(ident)
      fgs <- fgseaRes[[ident]]
      sigPathways <- dplyr::filter(fgs$fgsea, padj < pval.thresh)
      collapsedSigPathways <- collapsePathways(sigPathways, RET@meta.list$fgsea$pathways[[pathway]], fgs$ranks)
      mainPathways <- filter(sigPathways, 
        pathway %in% collapsedSigPathways$mainPathways)
      # print(sigPathways)
      topPathwaysUp <- filter(mainPathways, NES > 0) %>%
        arrange(desc(NES)) %>% top_n(n = -10, wt = pval)
      topPathwaysDown <- filter(mainPathways, NES < 0) %>%
        arrange(NES) %>% top_n(n = -10, wt = pval)
      topUpLists <- topPathwaysUp$leadingEdge
      # names(topUpLists) <- topPathwaysUp$pathway
      topDownLists <- topPathwaysDown$leadingEdge

      leadingEdgeLists <- c(topUpLists, topDownLists) #mainPathways$leadingEdge
      names(leadingEdgeLists) <- c(topPathwaysUp$pathway, topPathwaysDown$pathway)
      # names(leadingEdgeLists) <- mainPathways$pathway
      # print(leadingEdgeLists)
      for (path.name in names(leadingEdgeLists)) {
        # print(path.name)
        if (path.name %in% names(leading.edges)) {
          leading.edges[[path.name]] <- unique(c(unlist(leading.edges[[path.name]]), 
                                                 leadingEdgeLists[[path.name]]))
        } else {
          leading.edges[[path.name]] <- unique(c(unlist(leadingEdgeLists[[path.name]])))
        }
        # if (mouse) {
        #   leading.edges[[path.name]] <- convertHumanGeneList(leading.edges[[path.name]])
        # }
      }
    }
  }

  RET <- scoreModules(RET, leading.edges)
  # RET <- correlateVars(RET, c(names(leading.edges)))
  return(RET)
}
