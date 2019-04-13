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
#' @export
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
  idents.use = NULL,
  shiny = FALSE,
  metadata.add = NULL,
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
    if (!is.null(idents.use)) {
      idents.all <- idents.all[idents.all %in% idents.use]
    }
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
      message("Calculating cluster ", idents.all[i])
    } 
    if (shiny) {
      incProgress(1/length(x = idents.all), detail = paste0("Calculating cluster ", idents.all[i]))
    }
    genes.de[[i]] <- tryCatch(
      expr = {
        FindMarkers.Seurat_GCD(
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
          metadata.add = metadata.add,
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
#' @param fast.threshold whether full markers should be computed (default: True)
#' @return list containing clustered Seurat object and TSNE plots
#'
#' @examples
#' seuratAllMarkers(RET)
#'
#' @export
seuratAllMarkers <- function(RET, 
                             do.fast = TRUE, 
                             do.full = FALSE, 
                             fast.threshold = 0.25, 
                             shiny = F, 
                             idents.use = NULL) {
  require(Seurat)
  if (do.fast) {
    message("Calculating fast markers...")
    # RET@markers[[RET@meta.list$active.ident]]
    RET <- setMarkersByCluster(RET, "quick", 
      FindAllMarkersShiny(RET@seurat, logfc.threshold = fast.threshold, shiny = shiny, idents.use = idents.use), idents.use = idents.use)
  }
  if (do.full) {
    message("Calculating full markers...")
    RET <- setMarkersByCluster(RET, "full", 
      FindAllMarkersShiny(RET@seurat, logfc.threshold = 0.05, shiny = shiny, idents.use = idents.use), idents.use = idents.use)
  }
  return(RET)
}
#' Wrapper function to find module markers between clusters for a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param do.fast whether fast markers should be computed (default: True)
#' @param do.full whether full markers should be computed (default: True)
#' @param fast.threshold whether full markers should be computed (default: True)
#' @return gcdSeurat
#'
#' @examples
#' seuratModuleMarkers(RET)
#'
#' @export
seuratModuleMarkers <- function(RET, modules.use = NULL, 
                                markers.name = "modules",
                                idents.use = NULL, shiny = F,
                                logfc.threshold = 0.05) {
  modules.use <- modules.use %||% names(RET@meta.list$modules)
  
  RET <- setMarkersByCluster(RET, markers.name, 
                             FindAllMarkersShiny(RET@seurat, 
                                                 logfc.threshold = lofgc.threshold, 
                                                 metadata.add = modules.use,
                                                 features = modules.use,
                                                 min.pct = 0,
                                                 shiny = shiny, 
                                                 idents.use = idents.use), idents.use = idents.use)
}
                                

#' Adds markers to cluster list 
#'
#' @param RET gcdSeurat
#' @param markers.type name for markers (e.g. "quick" or "full")
#' @param markers returned data frame from Seurat FindAllMarkers
#' @return gcdSeurat
#'
#' @examples
#' setMarkersByCluster(RET, markers.type, markers)
#'
#' @export
setMarkersByCluster <- function(RET, markers.type, markers, idents.use) {
  ident.name <- ActiveIdent(RET)
  if (is.null(idents.use)) {
    idents.use <- levels(RET@seurat)
  }
  if (! markers.type %in% names(RET@markers)) {
    RET@markers[[markers.type]] <- list()
  } 
  if (ident.name %in% names(RET@markers[[markers.type]])) {
    if (length(intersect(names(RET@markers[[markers.type]][[ident.name]]), 
        levels(RET@seurat))) < 1) {
      RET@markers[[markers.type]][[ident.name]] <- list()
    }
  } else {
      RET@markers[[markers.type]][[ident.name]] <- list()
  } 
  # message(names(RET@markers[[markers.type]][[ident.name]]))
  for (ident in idents.use) {
    # message("Ident ", ident)
    ident.markers <- dplyr::filter(markers, cluster == ident)
    # if (nrow(ident.markers) != 0) {
    RET@markers[[markers.type]][[ident.name]][[ident]] <- ident.markers
    # } else {
      # RET@markers[[markers.type]][[ident.name]][[ident]] <- list()
    
  }
  return(RET)
}

#' FindMarkers
#' @param ident.1 Identity class to define markers for; pass an object of class
#' \code{phylo} or 'clustertree' to find markers for a node in a cluster tree;
#' passing 'clustertree' requires \code{\link{BuildClusterTree}} to have been run
#' @param ident.2 A second identity class for comparison; if \code{NULL},
#' use all other cells for comparison; if an object of class \code{phylo} or
#' 'clustertree' is passed to \code{ident.1}, must pass a node to find markers for
#' @param assay Assay to use in differential expression testing
#' @param reduction Reduction to use in differential expression testing - will test for DE on cell embeddings
#' @param group.by Regroup cells into a different identity class prior to performing differential expression (see example)
#' @param subset.ident Subset a particular identity class prior to regrouping. Only relevant if group.by is set (see example)
#'
#' @importFrom methods is
#'
#' @export
#'
#'
# `@method FindMarkers_GCD Seurat

FindMarkers.Seurat_GCD <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  reduction = NULL, 
  features = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  metadata.add = NULL,
  ...
) {
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) & !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  data.slot <- ifelse(
    test = test.use %in% c("negbinom", "poisson", "DESeq2"),
    yes = 'counts',
    no = 'data'
  )
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <-  GetAssayData(object = object[[assay]], slot = data.slot)
  } else {
    if (data.slot == "counts") {
      stop("The following tests cannot be used when specifying a reduction as they assume a count model: negbinom, poisson, DESeq2")
    }
    data.use <- t(x = Embeddings(object = object, reduction = reduction))
  }
  if (!is.null(metadata.add)) {
    message("Adding vars")
    metadata.vals <- FetchData(object, vars = metadata.add)
    data.use <- rbind(data.use, t(metadata.vals))
  }
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  } else if (ident.1 == 'clustertree' || is(object = ident.1, class2 = 'phylo')) {
    if (is.null(x = ident.2)) {
      stop("Please pass a node to 'ident.2' to run FindMarkers on a tree")
    }
    tree <- if (is(object = ident.1, class2 = 'phylo')) {
      ident.1
    } else {
      Tool(object = object, slot = 'BuildClusterTree')
    }
    if (is.null(x = tree)) {
      stop("Please run 'BuildClusterTree' or pass an object of class 'phylo' as 'ident.1'")
    }
    ident.1 <- tree$tip.label[GetLeftDescendants(tree = tree, node = ident.2)]
    ident.2 <- tree$tip.label[GetRightDescendants(tree = tree, node = ident.2)]
  }
  if (length(x = as.vector(x = ident.1)) > 1 &&
      any(as.character(x = ident.1) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(x = !as.character(x = ident.1) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- WhichCells(object = object, idents = ident.1)
  }
  # if NULL for ident.2, use all other cells
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.2) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      ident.2 <- setdiff(x = colnames(x = data.use), y = ident.1)
    } else {
      ident.2 <- WhichCells(object = object, idents = ident.2)
    }
  }
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(ident.1, ident.2)
    )
  }
  de.results <- FindMarkers(
    object = data.use,
    cells.1 = ident.1,
    cells.2 = ident.2,
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
  return(de.results)
}