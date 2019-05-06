library( shiny )
library(shinydashboard)
library(shinyBS)

library( Seurat )
library(scales)
library(knitr)
library(kableExtra)
library(future)
# library(gcdR)
library(ggtree)
library(yaml)
library(dplyr)
library(DT)
library(colorspace)

library(shinyFiles)
library(fs)

options(shiny.maxRequestSize = 120000*1024^2)


mydebug <- function(msg="[DEBUG]") {
  DEBUG <- FALSE
  if (DEBUG) {
    print(sprintf("%s - %s - %s", msg1, as.character(Sys.time()), as.character(deparse(sys.calls()[[sys.nframe()-1]]))))
  }
}

appDir <- system.file("shiny_app", package = "gcdR")
if (appDir == "") {
  stop("Could not find example directory. Try re-installing `gcdR`.", call. = FALSE)
}

server <- shinyServer(function(input, output, session){
  source(file.path(appDir, "server", "de_markers.R"), local = TRUE)$value
  source(file.path(appDir, "server", "fgsea_panel.R"), local = TRUE)$value
  source(file.path(appDir, "server", "initial_load_and_download.R"), local = TRUE)$value
  source(file.path(appDir, "server", "overview_panel.R"), local = TRUE)$value
  source(file.path(appDir, "server", "recluster.R"), local = TRUE)$value
  source(file.path(appDir, "server", "renaming.R"), local = TRUE)$value
  source(file.path(appDir, "server", "gene_expression_panel.R"), local = TRUE)$value
  source(file.path(appDir, "server", "marker_lists_panel.R"), local = TRUE)$value
  
  

  
  output$NAMES <- renderUI({
    req(input$FILTER)
    GROUP.NAMES <- levels(DATA$orig.RET@seurat@meta.data[,input$FILTER])
    selectInput("SELECTION", "Filter by Sample/Cluster/Type", GROUP.NAMES, multiple = T, selectize = T)
  })
  
  output$DIM.REDUC.CHOICE <- renderUI({
    req(DATA$orig.RET)
    radioGroupButtons("DIM.REDUC", "Reduction",
                 choices = DATA$orig.RET@meta.list$reductions, selected = DATA$orig.RET@meta.list$reductions[1])
  })
})
