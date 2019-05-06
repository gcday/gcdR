#' Launches Seurat viewer app
#' 
#' @param ... Arguments passed to shiny::runApp
#' @examples
#' runCellViewer()
#' 
#' @importFrom shiny reactiveValues runApp
#' @export
runCellViewer <- function(object = NULL, ...) {
  appDir <- system.file("shiny_app", package = "gcdR")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `gcdR`.", call. = FALSE)
  }
  DATA <- reactiveValues() 
  
  # shiny::runApp(appDir, display.mode = "normal", ...)
  
  
  source(file.path(appDir, "server.R"), local = T)
  source(file.path(appDir, "ui.R"), local = T)
  if (!is.null(object)) {
    DATA$INPUT.OBJECT <- object
  }
  shiny::shinyApp(ui = ui, server = server,  options = list("display.mode" = "normal"), ...)
}