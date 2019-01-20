#' Launches Seurat viewer app
#' @examples
#' runCellViewer()
#' @export
runCellViewer <- function() {
  appDir <- system.file("shiny_app", package = "gcdR")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `gcdR`.", call. = FALSE)
  }
  shiny::runApp(appDir, display.mode = "normal")
}