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


server <- shinyServer(function(input, output, session){
  DATA <- reactiveValues() 
  source(file.path("server", "de_markers.R"), local = TRUE)$value
  source(file.path("server", "fgsea_panel.R"), local = TRUE)$value
  source(file.path("server", "initial_load_and_download.R"), local = TRUE)$value
  source(file.path("server", "overview_panel.R"), local = TRUE)$value
  source(file.path("server", "recluster.R"), local = TRUE)$value
  source(file.path("server", "renaming.R"), local = TRUE)$value
  source(file.path("server", "gene_expression_panel.R"), local = TRUE)$value
  source(file.path("server", "marker_lists_panel.R"), local = TRUE)$value
  
  

  
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
