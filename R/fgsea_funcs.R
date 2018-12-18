#' Finds enriched pathways using fgsea
#'
#'
#' @param RET list containing Seurat object
#' @param dbs pathway databases 
#' 
#' @return list containing Seurat object and marker heatmaps
#'
#' @examples
#' makeMarkerHeatmaps(RET, marker.lists)
#'
#' @export
findEnrichedPathways <- function(RET, dbs) {
  RET[["enriched"]] <- getSigUpPathways(RET$all.markers, dbs)
  return(RET)
}
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
  # de.markers[!(is.na(de.markers$p_val)) & de.markers$p_val == 0,]$p_val <- 
  # de.markers <- de.markers[!duplicated(de.markers$gene),]
  de.markers <- mutate(de.markers, rank_metric = -log10(p_val) * sign(avg_logFC))
  # de.markers <- de.markers[!duplicated(de.markers$rank_metric),]
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
gseaPerCluster <- function(all.markers, pathways) {
  require(fgsea)
  require(dplyr)
  fgseaRes <- list()
  for (ident in levels(as.factor(all.markers$cluster))) { 
    ranks <- seuratDEresToGSEARanks(filter(all.markers, cluster == ident))
    fgsea.out <- fgsea(pathways = pathways, stats = ranks,
                               minSize=15,
                               maxSize=500,
                               nperm=100000)
    fgseaRes[[ident]] <- list()
    fgseaRes[[ident]][["ranks"]] <- ranks
    fgseaRes[[ident]][["fgsea"]] <- fgsea.out
  }
  return(fgseaRes)
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
    sigPathways <- filter(fgs$fgsea, padj < pval.thresh) %>% arrange(padj)
    collapsedSigPathways <- collapsePathways(sigPathways, pathways, fgs$ranks)
    mainPathways <- filter(sigPathways, pathway %in% collapsedSigPathways$mainPathways)
    # topPathwaysUp <- filter(fgs$fgsea, padj < pval.thresh, NES > 0) %>%
    #   arrange(desc(NES)) %>% top_n(n = 10, wt = NES)
    # topPathwaysDown <- filter(fgs$fgsea, padj < pval.thresh, NES < 0) %>%
    #   arrange(NES) %>% top_n(n = -10, wt = NES)
    # topPathwaysUp <- filter(mainPathways, NES > 0) %>%
    #   arrange(desc(NES)) %>% top_n(n = 10, wt = NES)
    # topPathwaysDown <- filter(mainPathways, NES < 0)%>%
    #   arrange(NES) %>% top_n(n = -10, wt = NES)
    topPathwaysUp <- filter(mainPathways, NES > 0) %>%
      arrange(desc(NES)) %>% top_n(n = -num.to.print, wt = padj)
    topPathwaysDown <- filter(mainPathways, NES < 0)%>%
      arrange(NES) %>% top_n(n = -num.to.print, wt = padj)
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
                               colwidths = c(10, 3, 0.8, 1.2))
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
                                     yend = statsAdj[p]), size = 0.2 ) + 
           scale_x_continuous(limits = c(0, length(statsAdj)), expand = c(0, 0)) + 
           scale_y_continuous(limits = c(-1, 1), expand = c(0, 0)) + 
           xlab(NULL) + ylab(NULL) + 
           theme(panel.background = element_blank(), axis.line = element_blank(), 
                 axis.text = element_blank(), axis.ticks = element_blank(), 
                 panel.grid = element_blank(), axis.title = element_blank(), 
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0, "null"), 4)), 
         textGrob(sprintf("%.2f", annotation$NES)), 
         textGrob(sprintf("%.1e", annotation$padj)))
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
  
  arrangeGrob(grobs = c(lapply(c("Pathway", "Gene ranks", "NES", "padj"), textGrob, 
                               gp = gpar(fontface="bold")),
                        unlist(ps, recursive = FALSE), 
                        list(nullGrob(), rankPlot, nullGrob(), nullGrob()), 
                        unlist(pads, recursive=FALSE)), ncol = 4, widths = colwidths)
}

#' Wrapper running fgsea
#'
#'
#' @param RET list containing Seurat object
#' @param pathways.list list of pathway sets to analyze
#' @param prefixes list of prefixes to trim from name of each pathway
#'
#' @return list containing ranks and fgsea output for each ident
#'
#' @examples
#' fgseaWrapper(RET, pathways.list)
#'
#' @export
fgseaWrapper <- function(RET, pathways.list, prefixes = NULL) {
  RET$fgsea <- list(pathways = list(), results = list(), prefixes=list())
  RET$plots$fgsea <- list()
  for (i in 1:length(pathways.list)) {
    path.name <- names(pathways.list)[i]
    pathways <- pathways.list[[path.name]]
    RET$fgsea$results[[path.name]] <- gseaPerCluster(RET$all.markers.full, c(pathways))
    RET$fgsea$pathways[[path.name]] <- pathways
    RET$fgsea$prefixes[[path.name]] <- ifelse(is.null(prefixes), "", prefixes[i])
    RET$plots$fgsea[[path.name]] <- gseaTopPlots(fgseaRes = RET$fgsea$results[[path.name]],
                                                 pathways = pathways, 
                                                 path.name = path.name, 
                                                 term.prefix = RET$fgsea$prefixes[[path.name]])
  }
  return(RET)
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
    path.dbs.to.check <- names(RET$plots$fgsea)
  }
  for (pathway in path.dbs.to.check) {
    if (redraw) {
      RET$plots$fgsea[[pathway]] <- gseaTopPlots(fgseaRes = RET$fgsea$results[[pathway]],
                                                 pathways = RET$fgsea$pathways[[pathway]],
                                                 path.name = pathway,
                                                 term.prefix = RET$fgsea$prefixes[[pathway]], ...)
    }
    for (ident in names(RET$plots$fgsea[[pathway]])) {
      grid.arrange(top=textGrob(ident, gp = gpar(fontvsize=20, fontface="bold")),
                   RET$plots$fgsea[[pathway]][[ident]]$up,
                   RET$plots$fgsea[[pathway]][[ident]]$down)
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
correlateFgseaLeadingEdges <- function(RET, pval.thresh = 0.01, path.dbs.to.check = NULL) {
  require("dplyr")
  leading.edges <- list()
  if(is.null(path.dbs.to.check)) {
    path.dbs.to.check <- names(RET$plots$fgsea)
  }
  for (pathway in path.dbs.to.check) {
    for (ident in names(RET$fgsea$results[[pathway]])) {
      fgs <- RET$fgsea$results[[pathway]][[ident]]
      sigPathways <- filter(fgs$fgsea, padj < pval.thresh)
      leadingEdgeLists <- sigPathways$leadingEdge
      names(leadingEdgeLists) <- sigPathways$pathway
      for (path.name in names(leadingEdgeLists)) {
        if (path.name %in% names(leading.edges)) {
          leading.edges[[path.name]] <- list(unique(c(unlist(leading.edges[[path.name]]), 
                                                 leadingEdgeLists[[path.name]])))
        } else {
          leading.edges[[path.name]] <- list(leadingEdgeLists[[path.name]])
        }
      }
    }
  }
  # print(leading.edges)
  RET <- scoreModules(RET, leading.edges)
  # RET <- correlateVars(RET, c(names(leading.edges)))
  return(RET)
}

#' Identifies light chain identity of sample
#'
#'
#' @param RET list containing Seurat object
#' 
#' @return updated \code{RET}
#'
#' @examples
#' RET <- findLikelyLightChainIdent(RET)
#'
#' @export
findLikelyLightChainIdent <- function(RET) {
  require("dplyr")
  require("tibble")
  avg.exprs <- as.data.frame(rowMeans(RET$seurat@data))
  colnames(avg.exprs) <-c("avg_expr")
  avg.exprs <- rownames_to_column(avg.exprs, var = "gene")
  IG.C.genes <- grep(pattern = "^IG[LK]C", x = avg.exprs$gene, value = TRUE)
  IG.V.genes <- grep(pattern = "^IG[LK]V", x = avg.exprs$gene, value = TRUE)
  IG.C.exprs <- filter(avg.exprs, gene %in% IG.C.genes)
  IG.V.exprs <- filter(avg.exprs, gene %in% IG.V.genes)
  top.C <- top_n(IG.C.exprs, n=1, wt = avg_expr)
  top.V <- top_n(IG.V.exprs, n=1, wt = avg_expr)
  RET$IG.light.V <- top.V$gene
  RET$IG.light.C <- top.C$gene
  return(RET)
}