
#' Filters modules
#'
#' @param modules named list of modules (containing gene names) to be added
#'
#' @export
#'
filterModules <- function(modules, control.pool) {
  filtered.modules <- list()
  for (module.name in names(modules)) {
    if (class(modules[[module.name]]) == "character") {
      modules[[module.name]] <- list(modules[[module.name]])
    }
    modules[[module.name]] <- c(unlist(modules[[module.name]]))
    features = modules[[module.name]]
    features = features[features %in% control.pool]
    if (length(features) <= 1 | length(features) > length(control.pool) / 24) {
      next
    }
    filtered.modules[[module.name]] <- features
  }
  return(filtered.modules)
}

#' Adds module scores to a Seurat object
#'
#'
#' @param RET list containing Seurat object
#' @param modules named list of modules (containing gene names) to be added
#'
#' @return list containing Seurat object and TSNE plots (tsne.treatment and tsne.ident) of the CCA
#' @examples
#' scoreModules(RET, modules)
#' 
#' @importFrom Seurat GetAssayData 
#'
#' @export

scoreModules <- function(RET, modules, verbose = T, min.cells = 5) {
  if (!'modules' %in% names(RET@meta.list)) {
    RET@meta.list$modules <- list()
  }
  # i <- 0
  assay.data <- GetAssayData(RET@seurat)
  if (verbose) message("Filtering control pool")
  cell.counts <- apply(assay.data, 1, function(x) {sum(x != 0)})
  # data.avg <- Matrix::rowMeans(x = assay.data[rownames(RET@seurat), ])
  # control.pool <- names(data.avg[data.avg != 0])
  control.pool <- names(cell.counts[cell.counts >= min.cells])
  filtered.modules <- filterModules(modules, control.pool)
  
  RET@seurat <- AddModuleScoreGCD(RET@seurat, 
                                  features = filtered.modules,
                                  pool = control.pool, 
                                  verbose = verbose,
                                  feature.names = names(filtered.modules))
  for (i in 1:length(filtered.modules)) {
    module.name <- names(filtered.modules)[i]
    # message(i, module.name)
    # RET@seurat[[module.name]] <- RET@seurat[[paste0("FilteredModule", i)]]
    # RET@seurat@meta.data <- dplyr::select(RET@seurat@meta.data, -matches(paste0("FilteredModule", i)))
    RET@meta.list$modules[[module.name]] <- filtered.modules[[module.name]]
  }
  
  # RET@meta.list$modules <- append(RET@meta.list$modules, filtered.modules)
  return(RET)
}
# scoreModules <- function(RET, modules, verbose = T) {
#   if (!'modules' %in% names(RET@meta.list)) {
#     RET@meta.list$modules <- list()
#   }
#   i <- 0
#   assay.data <- GetAssayData(RET@seurat)
#   data.avg <- Matrix::rowMeans(x = assay.data[rownames(RET@seurat), ])
#   control.pool <- names(data.avg[data.avg != 0])
#   filtered.modules <- filterModules(modules, control.pool)
#   for (module.name in names(filtered.modules)) {
#     i <- i + 1
#     features <- filtered.modules[[module.name]]
#     RET@meta.list$modules[[module.name]] <- features
#     if (verbose) message("Scoring ", module.name, " (", i, " of ", length(modules), ")")
#     RET@seurat <- AddModuleScore(RET@seurat, features = features,
#                                  name = module.name, pool = control.pool,
#                                  ctrl = min(vapply(X = features, FUN = length,
#                                                    FUN.VALUE = numeric(length = 1))))
#     RET@seurat@meta.data[[module.name]] <- RET@seurat@meta.data[[paste0(module.name, "1")]]
#     RET@seurat@meta.data <- dplyr::select(RET@seurat@meta.data, -matches(paste0(module.name, "1")))
#   }
#   return(RET)
# }

#' Calculate module scores for featre expression programs in single cells
#'
#' Calculate the average expression levels of each program (cluster) on single cell level,
#' subtracted by the aggregated expression of control featre sets.
#' All analyzed featres are binned based on averaged expression, and the control featres are
#' randomly selected from each bin.
#'
#' @param object Seurat object
#' @param features Featre expression programs in list
#' @param pool List of features to check expression levels agains, defaults to \code{rownames(x = object)}
#' @param nbin Number of bins of aggregate expression levels for all analyzed features
#' @param ctrl Number of control features selected from the same bin per analyzed feature
#' @param k Use feature clusters returned from DoKMeans
#' @param assay Name of assay to use
#' @param name Name for the expression programs
#' @param seed Set a random seed
#'
#' @return Returns a Seurat object with module scores added to object meta data
#'
# @importFrom Hmisc cut2
#' @importFrom ggplot2 cut_number
#' @importFrom Matrix rowMeans colMeans
#' @importFrom future.apply future_lapply
#'
#' @references Tirosh et al, Science (2016)
#'
#' @export
#'
#' @examples
#' \dontrun{
#' cd_features <- list(c(
#'   'CD79B',
#'   'CD79A',
#'   'CD19',
#'   'CD180',
#'   'CD200',
#'   'CD3D',
#'   'CD2',
#'   'CD3E',
#'   'CD7',
#'   'CD8A',
#'   'CD14',
#'   'CD1C',
#'   'CD68',
#'   'CD9',
#'   'CD247'
#' ))
#' pbmc_small <- AddModuleScore(
#'   object = pbmc_small,
#'   features = cd_features,
#'   ctrl = 5,
#'   name = 'CD_Features'
#' )
#' head(x = pbmc_small[])
#' }
#'
AddModuleScoreGCD <- function(
  object,
  features,
  pool = NULL,
  nbin = 24,
  ctrl = 100,
  k = FALSE,
  assay = NULL,
  name = 'Cluster',
  seed = 1,
  verbose = F,
  do.parallel = T,
  check.parallel = F,
  feature.names = NULL,
  check.features = F,
  bin_size = 256
) {
  set.seed(seed = seed)
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = 'k')
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster == i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    if (check.features) {
    features <- lapply(
      X = features,
      FUN = function(x) {
        return(intersect(x = x, y = rownames(x = object)))
      }
    )
    }
    cluster.length <- length(x = features)
  }
  if (!all(LengthCheck(values = features))) {
    warning(paste(
      'Could not find enough features in the object from the following feature lists:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'Attempting to match case...'
    ))
    features <- lapply(
      X = features.old,
      FUN = CaseMatch,
      match = rownames(x = object)
    )
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste(
      'The following feature lists do not have enough features present in the object:',
      paste(names(x = which(x = !LengthCheck(values = features)))),
      'exiting...'
    ))
  }
  
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg, n = nbin, labels = FALSE, right = FALSE)
  # data.cut <- as.numeric(x = Hmisc::cut2(x = data.avg, m = round(x = length(x = data.avg) / (nbin + 1))))
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)
  for (i in 1:cluster.length) {
    features.use <- features[[i]]
    for (j in 1:length(x = features.use)) {
      ctrl.use[[i]] <- c(
        ctrl.use[[i]],
        names(x = sample(
          x = data.cut[which(x = data.cut == data.cut[features.use[j]])],
          size = ctrl,
          replace = FALSE
        ))
      )
    }
  }
  ctrl.use <- lapply(X = ctrl.use, FUN = unique)

  
  if (do.parallel) {
    if (verbose) message("Starting control gene calculations...")
    # ctrl.scores.par <- future.apply::future_lapply(ctrl.use, function(features.use, assay.data) {
    #   Matrix::colMeans(x = assay.data[features.use, ])
    # }, assay.data)
    ctrl.scores.par <- future.apply::future_lapply(ctrl.use, function(features.use) {
      Matrix::colMeans(x = assay.data[features.use, ])
    })
    mat.ctrl <- matrix(unlist(ctrl.scores.par), nrow = length(ctrl.use), byrow = T)
    if (verbose) message("Starting feature gene calculations...")
    # features.scores.par <- future.apply::future_lapply(features, function(features.use, assay.data) {
    #   Matrix::colMeans(x = assay.data[features.use, , drop = FALSE])
    # }, assay.data)
    features.scores.par <- future.apply::future_lapply(features, function(features.use) {
      Matrix::colMeans(x = assay.data[features.use, , drop = FALSE])
    })
    mat.features <- matrix(unlist(features.scores.par), nrow = cluster.length, byrow = T)
  }
  if (!do.parallel | check.parallel) {
    ctrl.scores <- matrix(
      data = numeric(length = 1L),
      nrow = length(x = ctrl.use),
      ncol = ncol(x = object)
    )
    if (verbose) message("Before ctrl colMeans")
    for (i in 1:length(ctrl.use)) {
      features.use <- ctrl.use[[i]]
      ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, ])
    }
    features.scores <- matrix(
      data = numeric(length = 1L),
      nrow = cluster.length,
      ncol = ncol(x = object)
    )
    if (verbose) message("Before feature colMeans")
    for (i in 1:cluster.length) {
        features.use <- features[[i]]
        data.use <- assay.data[features.use, , drop = FALSE]
        features.scores[i, ] <- Matrix::colMeans(x = data.use)
    }
  }
  
  if (do.parallel) {
    if (check.parallel) {
      if (!all(ctrl.scores == mat.ctrl) & !all(features.scores == mat.features)) {
        stop(paste(
          'The parallel and sequential objects are not equal!',
          "Control matrices match:",
          all(ctrl.scores == mat.ctrl),
          "Feature matrices match:",
          all(features.scores == mat.features),
          'exiting...'
        ))
      }
    }
    ctrl.scores <- mat.ctrl
    features.scores <- mat.features
  }
  message("cluster length: ", cluster.length)
  
  feature.names <- feature.names %||% paste0(name, 1:cluster.length)
  
  features.scores.use <- features.scores - ctrl.scores
  rownames(x = features.scores.use) <- feature.names
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  gc(verbose = FALSE)
  DefaultAssay(object = object) <- assay.old
  return(object)
}