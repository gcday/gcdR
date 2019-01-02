#' Performs canonical correction analysis (CCA)
#'
#'
#' @param all.treatments Seurat object containing cells in both conditions 
#' @param condition.1 first condition of interest
#' @param condition.2 second condition of interest
#' @return list containing Seurat objects and plots
#'
#' @examples
#' seuratCCAWrapper(all.treatments, condition.1, condition.2)
#'
#' @export
seuratCCAWrapper <- function(all.treatments, condition.1, condition.2) {
  require("Seurat")
  c1.counts <- GetAssayData(subset(all.treatments, idents = condition.1), slot = "counts")
  # c1.counts <- SubsetData(all.treatments, ident.use = c(condition.1), do.clean=TRUE)@raw.data
  COND1.RET <- seuratStartFromCounts(c1.counts, project = condition.1, mito.prefix = "^mt-")

  # COND1.RET <- seuratFilterWrapper(c1.counts, project.name = condition.1, mito.prefix = "^mt-")
  COND1.RET$seurat@meta.data$treatment <- condition.1
  
  c2.counts <- GetAssayData(subset(all.treatments, idents = condition.2), slot = "counts")

  # c2.counts <- SubsetData(all.treatments, ident.use = c(condition.2), do.clean=TRUE)@raw.data
  # COND2.RET <- seuratFilterWrapper(c2.counts, project.name = condition.2, mito.prefix = "^mt-")
  COND2.RET <- seuratStartFromCounts(c2.counts, project = condition.2, mito.prefix = "^mt-")

  COND2.RET$seurat@meta.data$treatment <- condition.2
  
  COND1.genes <- head(VariableFeatures(COND1.RET$seurat), 3000)
  COND2.genes <- head(VariableFeatures(COND2.RET$seurat), 3000)
  genes.use <- unique(c(COND1.genes, COND2.genes))
  genes.use <- intersect(genes.use, rownames(COND1.RET$seurat@scale.data))
  genes.use <- intersect(genes.use, rownames(COND2.RET$seurat@scale.data))
  
  C1.C2.combined <- RunCCA(COND1.RET$seurat, COND2.RET$seurat, features = genes.use, num.cc = 30)
  RET <- list()
  RET$cond.1 <- condition.1
  RET$cond.2 <- condition.2
  PLOTS <- list()
  
  PLOTS[["p1"]] <- DimPlot(object = C1.C2.combined, reduction = "cca", group.by = "treatment", 
                pt.size = 0.5, do.return = TRUE)
  
  PLOTS[["p2"]] <- VlnPlot(object = C1.C2.combined, features = c("CC1"), group.by = "treatment",
                do.return = TRUE)
  
  PLOTS[["p3"]] <- MetageneBicorPlot(C1.C2.combined, grouping.var = "treatment", dims.eval = 1:30,
                          display.progress = TRUE)
  
  RET[["plots"]] <- PLOTS
  RET[["CCA"]] <- C1.C2.combined
  RET <- seuratAlignClusterWrapper(RET)
  # saveRDS(RET, paste(condition.1, condition.2, "CCA", "RDS", sep = "."))
  return(RET)
}
#' Aligns clusters between two conditions
#'
#'
#' @param RET list containing Seurat objects for cond1 and cond2 and plots
#' @param dims.align dimensions to be used for alignment and clustering
#' @param resolution resolution to be used for clustering
#' @return list containing aligned Seurat object and TSNE plot of the CCA
#'
#' @examples
#' seuratAlignClusterWrapper(RET)
#'
#' @export
seuratAlignClusterWrapper <- function(RET, dims.align = 1:20, resolution = 0.6) {
  require("Seurat")
  RET$CCA <- AlignSubspace(RET$CCA, reduction.type = "cca", grouping.var = "treatment", 
                                  dims.align = dims.align)
  
  RET$plots[["p4"]] <- VlnPlot(object = RET$CCA,
                               features.plot = "ACC1",
                               group.by = "treatment",
                               do.return = TRUE)
  RET$plots[["p5"]] <- VlnPlot(object = RET$CCA,
                               features.plot = "ACC2", 
                               group.by = "treatment",
                           do.return = TRUE)
  RET$CCA <- RunTSNE(RET$CCA, reduction.use = "cca.aligned", 
                     dims.use = dims.align, do.fast = T)
  RET$CCA <- FindClusters(RET$CCA, reduction.type = "cca.aligned", 
                              resolution = resolution, dims.use = dims.align)
  return(RET)
}
#' Finds conserved markers for all clusters
#'
#'
#' @param RET list containing aligned Seurat object
#' 
#' @return list containing aligned Seurat object and TSNE plot of the CCA
#'
#' @examples
#' seuratAlignClusterWrapper(RET)
#'
#' @export
seuratConservedMarkerWrapper <- function(RET) {
  RET[["conserved.markers"]] <- list()
  for (i in levels(RET$CCA@ident)){
    print(i)
    RET$conserved.markers[[i]] <- FindConservedMarkers(RET$CCA, ident.1 = i, grouping.var = "treatment")
  }
  return(RET)
}

# #' Finds DE genes in each cluster between conditions 
# #'
# #'
# #' @param RET list containing aligned Seurat object
# #' 
# #' @return list containing aligned Seurat object and TSNE plot of the CCA
# #'
# #' @examples
# #' seuratAlignClusterWrapper(RET)
# #'
# #' @export
# seuratMarkersBetweenConditions <- function(RET) {
#   require("Seurat")
#   require("tibble")
#   require("dplyr")
#   RET$CCA@meta.data$celltype.treatment <- paste0(RET$CCA@ident,
#                                                  "_",
#                                                  RET$CCA@meta.data$treatment)
#   RET$CCA <- StashIdent(RET$CCA, save.name = "celltype")
#   RET$CCA <- SetAllIdent(RET$CCA, id = "celltype.treatment")
  
#   library(tibble)
#   RET$response_markers <- NULL
#   for (ident in levels(as.factor(RET$CCA@meta.data$celltype))) {
#     print(ident)
#     RET$response_markers <- rbind(RET$response_markers,
#                          rownames_to_column(FindMarkers(RET$CCA,
#                                                         ident.1 = paste0(ident, "_", RET$cond.2),
#                                                         ident.2 = paste0(ident, "_", RET$cond.1)),
#                                             var="gene") %>% mutate(cluster = ident))
#   }
#   RET$CCA <- SetAllIdent(RET$CCA, id = "celltype")
#   return(RET)
# }

#' Makes TSNE plot after CCA 
#'
#'
#' @param RET list containing aligned Seurat object
#' 
#' @return list containing Seurat object and TSNE plots (tsne.treatment and tsne.ident) of the CCA
#'
#' @examples
#' seuratTSNEPlotCCA(RET)
#'
#' @export
seuratTSNEPlotCCA <- function(RET) {
  RET@meta.list$plots$tsne.treatment <- TSNEPlot(RET$CCA, do.return = T, pt.size = 0.5, group.by = "treatment")
  RET@meta.list$plots$tsne.ident <- TSNEPlot(RET$CCA, do.label = T, do.return = T, pt.size = 0.5)
  return(RET)
  # plot_grid(p1, p2)
}







