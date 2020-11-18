# Get the number of threads provided by the current plan
#
# @return The number of threads (workers) for the current future plan, or 1 if no workers detected
#
#' @importFrom future plan
#' @export
PlanThreads <- function() {
  nthreads <- eval(expr = formals(fun = plan())$workers)
  return(nthreads %||% 1)
}

# Check the length of components of a list
#
# @param values A list whose components should be checked
# @param cutoff A minimum value to check for
#
# @return a vector of logicals
#
LengthCheck <- function(values, cutoff = 0) {
  return(vapply(
    X = values,
    FUN = function(x) {
      return(length(x = x) > cutoff)
    },
    FUN.VALUE = logical(1)
  ))
}

#' Reads in a list of markers by cell type
#'
#'
#' @param RET list containing Seurat object and plots
#' 
#' @return none
#'
#' @examples
#' readMarkersYAML(yaml.path)
#' 
#' 
#' @importFrom yaml read_yaml
#'
#' @export
readMarkersYAML <- function(yaml.path) {
  raw.list <- read_yaml(yaml.path)
  markers.list <- list()
  alt.names <- list()
  for (group.name in names(raw.list)) {
    markers <- raw.list[[group.name]]
    markers.list[[group.name]] <- NULL
    for (marker.line in markers) {
      splits <- unlist(strsplit(marker.line[[1]], split = "/", fixed = T))
      if (length(splits) == 1) {
        markers.list[[group.name]] <- c(markers.list[[group.name]], trimws(splits[[1]]))
      } else {
        markers.list[[group.name]] <- c(markers.list[[group.name]], trimws(splits[[1]]))
        alt.names[[trimws(splits[[1]])]] <-  trimws(splits[[2]])
      }
      
    }
    
  }
  return(list(markers.list = markers.list, alt.names = alt.names))
}

#' Calculates average expression and counts cells with 
#' nonzero expression for each ident. 
#' 
#' @param object A Seurat object
#' @param assay Assay to use in differential expression testing
#' @param slot Slot to pull data from
#' @param verbose Print message when starting calculation for each ident
#'
#' @importFrom future.apply future_apply
#' 
#' 
#' @return List of two matrices: log_expr, containing mean of exmp1-transformed 
#' expression values for each feature across the cells in each ident; 
#' and expr_counts, containing the count of cells in each ident expressing 
#' each feature
#' 
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ident.summary <- SummarizeExprAcrossIdents(pbmc_small)
#' }
#'
SummarizeExprAcrossIdentsOld <- function(object, 
                                      assay = NULL, 
                                      features = NULL, 
                                      slot = 'data', 
                                      verbose = F) {
  # trying this as recommended in https://github.com/HenrikBengtsson/future.apply/issues/39
  # goal: disable future package's high-overhead checking for global variables
  opts <- options()
  on.exit(options(opts))
  options(future.globals.maxSize = +Inf)
  
  assay <- assay %||% DefaultAssay(object = object)
  data.use <-  GetAssayData(object = object[[assay]], slot = slot)
  features <- features %||% rownames(x = data.use)
  data.use <- data.use[features,]
  # holding averaged expression of all cells in each ident
  avg_expr_df <- data.frame(row.names = features)
  # number of cells with >= thresh.min for each feature per ident
  expr_count_df <- data.frame(row.names = features)
  thresh.min <- 0
  for (ident in levels(Idents(object))) {
    if (verbose) message("Averaging expression for cluster: ", ident)
    cells.ident <- WhichCells(object = object, idents = ident)
    avg_expr_df[[ident]] <- future.apply::future_apply(
      X = data.use[, cells.ident, drop = FALSE],
      MARGIN = 1,
      FUN = function(x) {
        return(mean(x = expm1(x = x)))
        # return(log(x = mean(x = expm1(x = x)) + pseudocount.use))
      }
    )
    expr_count_df[[ident]] <- 
      rowSums(x = data.use[, cells.ident, drop = FALSE] > thresh.min)
  }
  return(list(log_expr = as.matrix(x = avg_expr_df), expr_counts = as.matrix(x = expr_count_df)))
}