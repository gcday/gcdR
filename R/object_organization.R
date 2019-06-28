
#' gcdSeurat
#'
#' Replacing the simple list I had before with this, hopefully improving memory usage.
#'
#' @slot seurat Unnormalized data such as raw counts or TPMs
#' @slot meta.list Normalized expression data
#' 
#' @import Seurat
#' @name gcdSeurat
#' @rdname gcdSeurat
#' @exportClass gcdSeurat
#'
setClass("gcdSeurat", 
         slots = c(
           seurat = "Seurat", 
           meta.list = "list",
           plots = "list",
           info = "list",
           markers = "list",
           fgsea = "list",
           enrichr = "list"
         )
)




#' Reads Seurat/gcdSeurat object and converts it to updated gcdSeurat
#'
#'
#' @param path Filepath to .rds object storing a gcdSeurat or Seurat object
#' @param object In-memory gcdSeurat/Seurat object.
#' @param shiny Print shiny progress updates?
#' 
#' @return gcdSeurat object
#' 
#' @export
#' 
readSeuratRDStoGCD <- function(path = "", object = NULL, shiny = FALSE) {
  if (!is.null(object)) {
    IN.OBJ <- object
  } else {
    IN.OBJ <- readRDS(path)
  }
  RET <- new(Class = "gcdSeurat")
  if (class(IN.OBJ) == "gcdSeurat") {
    RET <- gcdUpdateSeurat(RET = IN.OBJ, shiny = shiny)
  } else {
    RET <- gcdUpdateSeurat(RET = NULL, SRT = IN.OBJ, shiny = shiny)
  }
  # if (!"active.ident" %in% names(RET@meta.list)) {
  #   RET <- findActiveIdent(RET)
  # }
  return(RET)
}

#' Reads Seurat/gcdSeurat object and converts it to updated gcdSeurat
#'
#'
#' @param RET gcdSeurat object
#' @return RET with meta.list modified to inclue active.ident
#'
#' @export
#'
gcdUpdateSeurat <- function(RET = NULL, SRT = NULL, shiny = TRUE) {
  if (is.null(RET)) {
     RET <- new(Class = "gcdSeurat")
  } else {
    SRT <- RET@seurat
  }
  dims <- NULL
  if (shiny) incProgress(amount = 0, message = "Checking Seurat object version...")
  if (utils::compareVersion(as.character(SRT@version), "3.0") != 1) {
    # pre Seurat 3.0, must update
    message("Updating Seurat object to 3.0...")
    if (shiny) {
      incProgress(amount = 0, message = "Updating Seurat object to 3.0...")
    }
    dims <- c(SRT@calc.params$RunTSNE$dims.use)
    SRT <- UpdateSeuratObject(SRT)
  }
  else {
    commands.use <- grep(pattern = "FindNeighbors", x = names(SRT@commands), value = T)
    if (length(commands.use) == 0) {
      commands.use <- names(SRT@commands)
    }
    for (name in commands.use) {
      message("command name: ", name)
      if ("dims" %in% names(SRT@commands[[name]]@params)) {
        dims <- SRT@commands[[name]]@params$dims
        message("Using dims (", dims, ") from ", name, " command.")
        break
      }
    }
    if (is.null(dims)) {
      dims <- c(1:10)
      message("Couldn't find dims in command history, using ", dims, " instead.")
    }
  }
  RET@seurat <- SRT
  RET@meta.list$dims <- dims
  return(RET)
}


#' Reads Seurat/gcdSeurat object and converts it to updated gcdSeurat
#'
#'
#' @param RET gcdSeurat object
#' @return RET with meta.list modified to inclue active.ident
#'
#' @export
#'
gcdRunUMAP <- function(RET, shiny = TRUE, check.tree = T) {
  
  UMAP.present <- T
  if (shiny) incProgress(amount = 0.3, message = "Checking for UMAP...")
  if (!"umap" %in% names(RET@seurat@reductions)) {
    if (reticulate::py_module_available(module = 'umap')) {
      if (shiny) incProgress(amount = 0, message = "Running UMAP...")

      RET@seurat <- RunUMAP(RET@seurat, dims = c(RET@meta.list$dims))
    } else if (shiny) {
      showModal(modalDialog(title = "Warning!",
                            "Missing umap package, only t-SNE will be available. Please install through pip (e.g. pip install umap-learn).",
                            easyClose = T))
      UMAP.present <- F
    }
  }
  if (UMAP.present) {
    RET@meta.list$reductions <- list(UMAP = "umap", `t-SNE` = "tsne")
  } else {
    RET$meta.list$reductions <- list(`t-SNE` = "tsne")
  }
  if (check.tree) {
    if (is.null(Tool(RET@seurat, slot = "BuildClusterTree"))) {
      if (shiny) incProgress(amount = 0.2, message = "Building cluster tree...")
      RET@seurat <- BuildClusterTree(RET@seurat, dims = RET@meta.list$dims)
    }
  }
 
  # RET@seurat <- BuildClusterTree(RET@seurat, dims = dims)
  return(RET)
}

#' Renames idents of cells in a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param new.idents list containing new idents (in same order as current)
#' @return list containing renamed Seurat object and re-plotted TSNE 
#'
#' @examples
#' renameIdents(RET)
#'
#' @export
renameAllIdents <- function(RET, new.idents) {
  for (i in 1:length(new.idents)) {
    if (levels(Idents(RET@seurat))[i] != new.idents[i]) {
      RET <- gcdRenameIdent(RET, levels(Idents(RET@seurat))[i], new.idents[i])
    }
  }
  return(RET)
  # levels(Idents(RET@seurat)) <- new.idents
  # if ("fgsea" %in% names(RET@meta.list)) {
  #   names(RET@meta.list$fgsea$ranks) <- new.idents
  #   for (pathway in names(RET@meta.list$fgsea$results)) {
  #     names(RET@meta.list$fgsea$results[[pathway]]) <- new.idents
  #   }
  # }
  # if (RET@meta.list$active.ident %in% names(RET@markers$quick)) {
  #   levels(RET@markers$quick[[RET@meta.list$active.ident]]) <- new.idents
  # }
  # if (RET@meta.list$active.ident %in% names(RET@markers$full)) {
  #   levels(RET@markers$full[[RET@meta.list$active.ident]]) <- new.idents
  # }
  # RET@plots$UMAP <- DimPlot(RET@seurat, label = T, reduction = "umap")
}

#' Renames selected groups of cells in a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param ident.to.rename SINGLE ident to rename
#' @param new.label new name for ident
#' @return list containing renamed Seurat object and re-plotted TSNE 
#'
#' @examples
#' renameIdents(RET)
#'
#' @export
gcdRenameIdent <- function(RET, ident.to.rename, new.label) {
  active.ident <- suppressWarnings(ActiveIdent(RET))
  levels(RET@seurat)[levels(RET@seurat) == ident.to.rename] <- new.label
  
  RET@seurat[[active.ident]][[active.ident]] <- Idents(RET@seurat)
  # levels(RET@seurat[[ActiveIdent(RET)]][[ActiveIdent(RET)]])[levels(RET@seurat[[ActiveIdent(RET)]][[ActiveIdent(RET)]]) == ident.to.rename] <- new.label
  
  # cluster.tree <- Tool(RET@seurat, slot = "BuildClusterTree")
  if (!is.null(Tool(RET@seurat, slot = "BuildClusterTree"))) {
    RET@seurat@tools$BuildClusterTree$tip.label[RET@seurat@tools$BuildClusterTree$tip.label == ident.to.rename] <- new.label
    # Tool(RET@seurat, slot = "BuildClusterTree") <- cluster.tree
  }

  for (marker.names in names(RET@markers)) {
    # message("renaming ", marker.names)
    if (active.ident %in% names(RET@markers[[marker.names]])) {
      # RET@markers[[marker.names]][[active.ident]]
      names(RET@markers[[marker.names]][[active.ident]])[names(RET@markers[[marker.names]][[active.ident]]) == ident.to.rename] <- new.label
    }
  }
  if ("results" %in% names(RET@fgsea)) {
    if (active.ident %in% names(RET@fgsea$results)) {
      names(RET@fgsea$results[[active.ident]]$ranks)[names(RET@fgsea$results[[active.ident]]$ranks) == ident.to.rename] <- new.label
      names(RET@fgsea$results[[active.ident]]$output)[names(RET@fgsea$results[[active.ident]]$output) == ident.to.rename] <- new.label
    }
  }
  if (active.ident %in% names(RET@enrichr)) {
    if (ident.to.rename %in% names(RET@enrichr[[active.ident]])) {
      names(RET@enrichr[[active.ident]])[names(RET@enrichr[[active.ident]]) == ident.to.rename] <- new.label
    }
  }
  ident.expr.summaries <- names(RET@meta.list$precalc.ident.expr)
  for (expr.names in ident.expr.summaries[tidyselect::starts_with(match = active.ident, vars = ident.expr.summaries)]) {
    colnames(RET@meta.list$precalc.ident.expr[[expr.names]]$expr_counts)[colnames(RET@meta.list$precalc.ident.expr[[expr.names]]$expr_counts) == ident.to.rename] <- new.label
    colnames(RET@meta.list$precalc.ident.expr[[expr.names]]$log_expr)[colnames(RET@meta.list$precalc.ident.expr[[expr.names]]$log_expr) == ident.to.rename] <- new.label
  }
  return(RET)
  # levels(DATA$orig.RET@seurat)[levels(DATA$orig.RET@seurat) %in% c(input$GRPS.TO.RENAME)] <- input$NEW.GRPS.NAME
}

#' Renames selected groups of cells in a Seurat object
#' @rdname Idents
#' @export
#' @method CellGroups gcdSeurat
#'
CellGroups.gcdSeurat <- function(object, ...) {
    meta.data <- slot(slot(object = object, name = "seurat"), name = "meta.data")
        
    # metadata.counts <- apply(meta.data, 2, function(x) length(table(x)))
    # ident.choices <- colnames(meta.data)[metadata.counts<=100]
    return(colnames(meta.data)[!sapply(meta.data, is.numeric)])
    # for (choice in SRTmeta.data)
      # RET@meta.list$ident.choices) {
    # if(isTRUE(all.equal(as.character(Idents(RET@seurat)),
                        # as.character(SRT@seurat@meta.data[[choice]])))) {
      # active.ident <- choice
      # break
    # }
}
#' @export
#'
CellGroups <- function(object, ...) {
  UseMethod(generic = 'CellGroups', object = object)
}

#' Finds and sets current active ident.
#'
#'
#' @param RET gcdSeurat object
#' @return RET with meta.list modified to inclue active.ident
#' 
#' @method ActiveIdent gcdSeurat
#' @export
#' 
ActiveIdent.gcdSeurat <- function(object, ...) {
  active.ident <- object@meta.list$active.ident
  # active.ident <- NULL
  cell.groups <- CellGroups(object)
  exact.matches <- cell.groups[sapply(X = cell.groups,
                                      FUN = function(choice) {
                                        isTRUE(all.equal(as.character(Idents(object@seurat)),
                                                         as.character(object@seurat@meta.data[[choice]])))
                                      })]
  if ("seurat_clusters" %in% exact.matches & length(exact.matches) > 1) {
    exact.matches <- exact.matches[exact.matches != "seurat_clusters"]
  }
  if (length(exact.matches) == 0) {
    warning("Couldn't find active ident in metadata!")
    return(NULL)
  }
  return(exact.matches[1])
  # for (choice in CellGroups(object)) {
  #   if() {
  #     active.ident <- choice
  #     object@meta.list$active.ident <- active.ident
  #     break
  #   }
  # }
  # if (is.null(active.ident)) {
  #   
  # }
  # return(object@meta.list$active.ident)
}

#' @export
#'
ActiveIdent <- function(object, ...) {
  UseMethod(generic = 'ActiveIdent', object = object)
}


#' @export
#' 
IdentsSaved <- function(object, ...) {
  UseMethod(generic = 'IdentsSaved', object = object)
}

#' Checks whether ident has been saved
#'
#'
#' @param object gcdSeurat
#' @return whether idents must be saved
#' 
#' @method IdentsSaved gcdSeurat
#' @export
#' 
IdentsSaved.gcdSeurat <- function(object, ...) {
  # active.ident <- object@meta.list$active.ident
  active.ident <- NULL
  for (choice in CellGroups(object)) {
    if(isTRUE(all.equal(as.character(Idents(object@seurat)),
                        as.character(object@seurat@meta.data[[choice]])))) {
      active.ident <- choice
      # object@meta.list$active.ident <- active.ident
      break
    }
  }
  return(!is.null(active.ident)) 

  # (is.null(active.ident)) {
    # return(FALSE)
    # warning("Couldn't find active ident in metadata!")
  # }
  # return(object@meta.list$active.ident)
}

#' Checks whether ident has been saved
#'
#'
#' @param object gcdSeurat
#' @return whether idents must be saved
#' 
#' @method IdentsSaved gcdSeurat
#' @export
#' 
ClosestIdent <- function(object, ...) {
  # active.ident <- object@meta.list$active.ident
  closest.ident <- NULL
  closest.diff <- Inf
  ident.levels <- levels(Idents(object@seurat))
  for (choice in CellGroups(object)) {
    choice.levels <- levels(object@seurat@meta.data[[choice]])
    if (length(ident.levels) != length(choice.levels)) {
      next
    }
    sym.diff <- union(setdiff(ident.levels, choice.levels), setdiff(choice.levels, ident.levels))
    if (length(sym.diff) < closest.diff) {
      closest.diff <- length(sym.diff)
      closest.ident <- choice
    } else if (length(sym.diff) == closest.diff) {
      closest.ident <- c(closest.ident, choice)
    }
  }
  message(closest.ident, " ", closest.diff)
  return(closest.ident) 
}

#' Finds and sets current active ident.
#'
#'
#' @param RET gcdSeurat object
#' @return RET with meta.list modified to inclue active.ident
#' 
#' @export
#' 
setActiveIdent <- function(RET) {
  
  active.ident <- NULL
  for (choice in CellGroups(RET)) {
    if(isTRUE(all.equal(as.character(Idents(RET@seurat)),
                        as.character(RET@seurat@meta.data[[choice]])))) {
      active.ident <- choice
      RET@meta.list$active.ident <- active.ident
      break
    }
  }
  if (is.null(active.ident)) {
    warning("Couldn't find active ident in metadata!")
  }
  else {
    RET@meta.list$active.ident <- active.ident
  }
  return(RET)
}
#' Finds and sets current active ident.
#'
#'
#' @param RET gcdSeurat object
#' @return RET with meta.list modified to inclue active.ident
#' 
#' @export
#' 
tryOverwriteIdent <- function(RET, OVERWRITE.IDENT) {
  if (is.null(ClosestIdent(RET))) {
    return(RET)
  }
  RET@seurat[[OVERWRITE.IDENT]] <- Idents(RET@seurat)
  for (marker.type in names(RET@markers)) {
    if (ActiveIdent(RET) %in% names(RET@markers[[marker.type]])) {
      RET@markers[[marker.type]][[OVERWRITE.IDENT]] <- RET@markers[[marker.type]][[ActiveIdent(RET)]]
      # names(RET@markers[[marker.type]][[OVERWRITE.IDENT]]) <- levels(Idents(RET@seurat))
      }
  }
  return(RET)
}


#' Finds best markers available for current idents.
#'
#' @param RET gcdSeurat object
#' @return markers slots containing entries matching ActiveIdent(RET)
#' 
#' @method ActiveMarkers gcdSeurat
#' @export
#' 
ActiveMarkers.gcdSeurat <- function(object, verbose = F, ...) {
  active.ident <- ActiveIdent(object)
  # Taking the set of markers with the most total genes observed 
  valid.choices <- NULL
  for (choice in names(object@markers)) {
    if (active.ident %in% names(object@markers[[choice]])) {
      valid.choices <- c(choice, valid.choices)
    }
  }
  return(valid.choices)
}

#' Finds best markers available for current idents.
#'
#' @param RET gcdSeurat object
#' @return RET with meta.list modified to inclue active.ident
#' 
#' @method BestMarkers gcdSeurat
#' @export
#' 
BestMarkers.gcdSeurat <- function(object, verbose = F, slot = NULL, ...) {
  max.markers <- 0
  best.markers <- NULL
  active.ident <- ActiveIdent(object)
  # Taking the set of markers with the most total genes observed 
  if (!is.null(slot)) {
    return(object@markers[[slot]][[active.ident]])
  }
  for (choice in names(object@markers)) {
    if (active.ident %in% names(object@markers[[choice]])) {
      num.markers <- sum(unlist(lapply(object@markers[[choice]][[active.ident]], nrow)))
      if (num.markers >= max.markers) best.markers <- choice
    }
  }
  if (is.null(best.markers)) {
    warning("Couldn't find active markers in metadata!")
  }
  return(object@markers[[best.markers]][[active.ident]])
}

#' @export
#'
BestMarkers <- function(object, ...) {
  UseMethod(generic = 'BestMarkers', object = object)
}

#' @export
#'
ActiveMarkers <- function(object, ...) {
  UseMethod(generic = 'ActiveMarkers', object = object)
}


#' Returns fgsea results for current idents. 
#'
#' @param RET gcdSeurat object
#' @return RET with meta.list modified to inclue active.ident
#' 
#' @importFrom dplyr %>%
#' @method GseaRes gcdSeurat
#' @export
#' 
GseaRes.gcdSeurat  <- function(object, ...) {
  max.markers <- 0
  best.markers <- NULL
  if ("results" %in% names(object@fgsea)) {
    if (ActiveIdent(object) %in% names(object@fgsea$results)) {
      return(object@fgsea$results[[ActiveIdent(object)]])
    } 
  }
  warning("Couldn't find fgsea in metadata!")
  return(NULL)
}

#' @export
#'
GseaRes <- function(object, ...) {
  UseMethod(generic = 'GseaRes', object = object)
}