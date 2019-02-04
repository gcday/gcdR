
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
           fgsea = "list"
         )
)




#' Reads Seurat/gcdSeurat object and converts it to updated gcdSeurat
#'
#'
#' @param RET gcdSeurat object
#' @return RET with meta.list modified to inclue active.ident
#' 
#' @export
#' 
readSeuratRDStoGCD <- function(path, shiny = FALSE) {
  IN.OBJ <- readRDS(path)
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
    for (name in names(SRT@commands)) {
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
gcdRunUMAP <- function(RET, shiny = TRUE) {
  
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
  if (shiny) incProgress(amount = 0.2, message = "Building cluster tree...")
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
  levels(Idents(RET@seurat)) <- new.idents
  if ("fgsea" %in% names(RET@meta.list)) {
    names(RET@meta.list$fgsea$ranks) <- new.idents
    for (pathway in names(RET@meta.list$fgsea$results)) {
      names(RET@meta.list$fgsea$results[[pathway]]) <- new.idents
    }
  }
  if (RET@meta.list$active.ident %in% names(RET@markers$quick)) {
    levels(RET@markers$quick[[RET@meta.list$active.ident]]) <- new.idents
  }
  if (RET@meta.list$active.ident %in% names(RET@markers$full)) {
    levels(RET@markers$full[[RET@meta.list$active.ident]]) <- new.idents
  }
  # RET@plots$UMAP <- DimPlot(RET@seurat, label = T, reduction = "umap")
  return(RET)
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
  levels(RET@seurat)[levels(RET@seurat) == ident.to.rename] <- new.label
  active.ident <- suppressWarnings(ActiveIdent(RET))
  # cluster.tree <- Tool(RET@seurat, slot = "BuildClusterTree")
  if (!is.null(Tool(RET@seurat, slot = "BuildClusterTree"))) {
    RET@seurat@tools$BuildClusterTree$tip.label[RET@seurat@tools$BuildClusterTree$tip.label == ident.to.rename] <- new.label
    # Tool(RET@seurat, slot = "BuildClusterTree") <- cluster.tree
  }

  for (marker.names in names(RET@markers)) {
    # message("renaming ", marker.names)
    if (active.ident %in% names(RET@markers[[marker.names]])) {
      RET@markers[[marker.names]][[active.ident]]
      names(RET@markers[[marker.names]][[active.ident]])[names(RET@markers[[marker.names]][[active.ident]]) == ident.to.rename] <- new.label
    }
    # for (ident.name in names(RET@markers[[marker.names]])) {
      # levels(RET@markers[[marker.names]][[ident.name]]$cluster)[levels(RET@markers[[marker.names]][[ident.name]]$cluster) %in% c(idents.to.rename)] <- new.label
    # }
    # levels(RET@markers[[marker.names]]$cluster)[levels(RET@markers[[marker.names]]$cluster) %in% c(idents.to.rename)] <- new.label
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
        
    metadata.counts <- apply(meta.data, 2, function(x) length(table(x)))
    ident.choices <- colnames(meta.data)[metadata.counts<=100]
    return(ident.choices)
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
  for (choice in CellGroups(object)) {
    if(isTRUE(all.equal(as.character(Idents(object@seurat)),
                        as.character(object@seurat@meta.data[[choice]])))) {
      active.ident <- choice
      object@meta.list$active.ident <- active.ident
      break
    }
  }
  if (is.null(active.ident)) {
    warning("Couldn't find active ident in metadata!")
  }
  return(object@meta.list$active.ident)
}

#' @export
#'
ActiveIdent <- function(object, ...) {
  UseMethod(generic = 'ActiveIdent', object = object)
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
          # print(names(DATA$orig.RET@markers[[marker.type]][[input$OVERWRITE.IDENT]]))
      names(RET@markers[[marker.type]][[OVERWRITE.IDENT]]) <- levels(Idents(RET@seurat))
      }
  }
  return(RET)
}
