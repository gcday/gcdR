#' Gene expression markers for all identity classes
#'
#' Modified to add progress indicator for Shiny apps.
#' Finds markers (differentially expressed genes) for each of the identity classes in a dataset
#'
#' @param node A node to find markers for and all its children; requires
#' \code{\link{BuildClusterTree}} to have been run previously; replaces \code{FindAllMarkersNode}
#' @param return.thresh Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)
#'
#' @return Matrix containing a ranked list of putative markers, and associated
#' statistics (p-values, ROC score, etc.)
#'
#' @importFrom ape drop.tip
#' @importFrom stats setNames
#'
#' @export
#'
#' @aliases FindAllMarkersNodeShiny
#'
#' @examples
#' # Find markers for all clusters
#' all.markers <- FindAllMarkers(object = pbmc_small)
#' head(x = all.markers)
#'
#' # Pass a value to node as a replacement for FindAllMarkersNode
#' pbmc_small <- BuildClusterTree(object = pbmc_small)
#' all.markers <- FindAllMarkers(object = pbmc_small, node = 4)
#' head(x = all.markers)
#'
FindAllMarkersShiny <- function(
  object,
  assay = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  return.thresh = 1e-2,
  ...
) {
  MapVals <- function(vec, from, to) {
    vec2 <- setNames(object = to, nm = from)[as.character(x = vec)]
    vec2[is.na(x = vec2)] <- vec[is.na(x = vec2)]
    return(unname(obj = vec2))
  }
  if ((test.use == "roc") && (return.thresh == 1e-2)) {
    return.thresh <- 0.7
  }
  if (is.null(x = node)) {
    idents.all <- sort(x = unique(x = Idents(object = object)))
  } else {
    tree <- Tool(object = object, slot = 'BuildClusterTree')
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' before finding markers on nodes")
    }
    descendants <- DFT(tree = tree, node = node, include.children = TRUE)
    all.children <- sort(x = tree$edge[, 2][!tree$edge[, 2] %in% tree$edge[, 1]])
    descendants <- MapVals(
      vec = descendants,
      from = all.children,
      to = tree$tip.label
    )
    drop.children <- setdiff(x = tree$tip.label, y = descendants)
    keep.children <- setdiff(x = tree$tip.label, y = drop.children)
    orig.nodes <- c(
      node,
      as.numeric(x = setdiff(x = descendants, y = keep.children))
    )
    tree <- drop.tip(phy = tree, tip = drop.children)
    new.nodes <- unique(x = tree$edge[, 1, drop = TRUE])
    idents.all <- (tree$Nnode + 2):max(tree$edge)
  }
  genes.de <- list()
  for (i in 1:length(x = idents.all)) {
    if (verbose) {
      incProgress(1/length(x = idents.all), detail = paste0("Calculating cluster ", idents.all[i]))
      message("Calculating cluster ", idents.all[i])
    }
    genes.de[[i]] <- tryCatch(
      expr = {
        FindMarkers(
          object = object,
          assay = assay,
          ident.1 = if (is.null(x = node)) {
            idents.all[i]
          } else {
            tree
          },
          ident.2 = if (is.null(x = node)) {
            NULL
          } else {
            idents.all[i]
          },
          features = features,
          logfc.threshold = logfc.threshold,
          test.use = test.use,
          min.pct = min.pct,
          min.diff.pct = min.diff.pct,
          verbose = verbose,
          only.pos = only.pos,
          max.cells.per.ident = max.cells.per.ident,
          random.seed = random.seed,
          latent.vars = latent.vars,
          min.cells.feature = min.cells.feature,
          min.cells.group = min.cells.group,
          pseudocount.use = pseudocount.use,
          ...
        )
      },
      error = function(cond) {
        return(NULL)
      }
    )
  }
  gde.all <- data.frame()
  for (i in 1:length(x = idents.all)) {
    if (is.null(x = unlist(x = genes.de[i]))) {
      next
    }
    gde <- genes.de[[i]]
    if (nrow(x = gde) > 0) {
      if (test.use == "roc") {
        gde <- subset(
          x = gde,
          subset = (myAUC > return.thresh | myAUC < (1 - return.thresh))
        )
      } else if (is.null(x = node) || test.use %in% c('bimod', 't')) {
        gde <- gde[order(gde$p_val, -gde$avg_logFC), ]
        gde <- subset(x = gde, subset = p_val < return.thresh)
      }
      if (nrow(x = gde) > 0) {
        gde$cluster <- idents.all[i]
        gde$gene <- rownames(x = gde)
      }
      if (nrow(x = gde) > 0) {
        gde.all <- rbind(gde.all, gde)
      }
    }
  }
  if ((only.pos) && nrow(x = gde.all) > 0) {
    return(subset(x = gde.all, subset = avg_logFC > 0))
  }
  rownames(x = gde.all) <- make.unique(names = as.character(x = gde.all$gene))
  if (nrow(x = gde.all) == 0) {
    warning("No DE genes identified", call. = FALSE)
  }
  if (!is.null(x = node)) {
    gde.all$cluster <- MapVals(
      vec = gde.all$cluster,
      from = new.nodes,
      to = orig.nodes
    )
  }
  return(gde.all)
}

#' Finds DE genes in each cluster between conditions 
#'
#'
#' @param RET list containing aligned Seurat object
#' @param cond.var name of variable storing condition values in Seurat object's metadata.
#' @param cond.1 first condition
#' @param cond.2 second condition 
#'
#' @return dataframe containing DE markers between conditions
#'
#' @examples
#' markersBetweenConditions(RET)
#'
#' @export
markersBetweenConditions <- function(RET, cond.var, cond.1, cond.2) {
  cond.markers <- seuratMarkersBetweenConditions(RET@seurat, cond.var, cond.1, cond.2)
  if (!"markers" %in% names(RET@meta.list)) {
    RET@meta.list$markers <- list()
  }
  RET@meta.list$markers[[paste0(cond.1, "_vs_", cond.2)]] <- cond.markers   
  return(RET)
}


#' Finds DE genes in each cluster between conditions 
#'
#'
#' @param SRT list containing aligned Seurat object
#' @param cond.var name of variable storing condition values in Seurat object's metadata.
#' @param cond.1 first condition
#' @param cond.2 second condition 
#'
#' @return dataframe containing DE markers between conditions
#'
#' @examples
#' seuratMarkersBetweenConditions(SRT)
#'
#' @export
seuratMarkersBetweenConditions <- function(SRT, cond.var, cond.1, cond.2) {
  SRT$celltype.cond <- paste0(Idents(SRT), "_", SRT[[cond.var]][[cond.var]])
  cond.levels <- levels(SRT[[cond.var]])
  old.idents <- levels(as.factor(Idents(SRT)))
  Idents(SRT) <- "celltype.cond"
  new.idents <- levels(as.factor(Idents(SRT)))
  num.idents <- length(levels(old.idents))
  cond.markers <- NULL
  for (old.id in old.idents) {
    id.1 <- paste0(old.id, "_", cond.1)
    id.2 <- paste0(old.id, "_", cond.2)
    message(id.1)
    if (id.1 %in% new.idents && id.2 %in% new.idents) {
      SRT.subset <- subset(SRT, idents = c(id.1, id.2))
      if (ncol(subset(SRT.subset, idents = c(id.1))) > 3 && ncol(subset(SRT.subset, idents = c(id.2))) > 3) {
        markers <- FindMarkers(SRT.subset, ident.1 = id.1, ident.2 = id.2, 
                          min.pct = 0, logfc.threshold = 0.05, verbose = T)
        markers <- tibble::rownames_to_column(markers, var="gene")
        cond.markers <- rbind(cond.markers, dplyr::mutate(markers, cluster = old.id))
      }
    }
  }
  return(cond.markers)
}

#' Wrapper function to find markers between clusters for a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param do.fast whether fast markers should be computed (default: True)
#' @param do.full whether full markers should be computed (default: True)
#' @return list containing clustered Seurat object and TSNE plots
#'
#' @examples
#' seuratAllMarkers(RET)
#'
#' @export
seuratAllMarkers <- function(RET, do.fast = TRUE, do.full = FALSE) {
  require(Seurat)
  if (do.fast) {
    message("Calculating fast markers...")
    RET@markers$all.markers.quick <- FindAllMarkers(RET@seurat)
  }
  if (do.full) {
    message("Calculating full markers...")
    RET@markers$all.markers.full <- FindAllMarkers(RET@seurat, logfc.threshold = 0.05)
  }
  return(RET)
}