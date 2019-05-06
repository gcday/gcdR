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