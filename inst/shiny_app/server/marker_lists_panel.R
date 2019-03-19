
observeEvent(input$MARKERS.LIST.PATH$datapath, {
  DATA$markers.list <- read_yaml(input$MARKERS.LIST.PATH$datapath)
})

multiPanelSeuratFigs <- function(features, fig.type = "violin", n.row = 2, n.col = 2, height = 600) {
  plots.list <- NULL
  if (fig.type == "violin") {
    plots.list <- VlnPlot(DATA$RET@seurat,
                          features = c(features),
                          cols = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE)), pt.size = 0,
                          combine = F)
    for (i in 1:length(plots.list)) {
      plots.list[[i]] <- plots.list[[i]] + NoLegend()
    }
  }
  else if (fig.type == "feature") {
    plots.list <- FeaturePlot(DATA$RET@seurat,
                              features = c(features),
                              reduction = input$DIM.REDUC,
                              combine = F)
  }
  return(renderPlot({
    cowplot::plot_grid(plotlist = plots.list,
                       nrow = n.row, ncol = n.col)
  }, width="auto",height=height))
}

plotMultiFeatures <- function(markers.name, fig.type = "violin", n.row = 2, n.col = 2) {
  markers <- DATA$markers.list[[markers.name]]
  markers <- markers[markers %in% rownames(DATA$RET@seurat)]
  if (length(markers) == 0) {
    return(NULL)
  }
  split.markers <- split(markers, 
                         ceiling(seq_along(markers)/(n.row * n.col)))
  if (fig.type == "violin") {
    height <- 600
  } else if (fig.type == "feature") {
    height <- 900
  }
  return(do.call(tabPanel,
                 args = c(title = markers.name,
                          purrr::map(names(split.markers),
                                     .f = function(chunk.name){
                                       multiPanelSeuratFigs(split.markers[[chunk.name]],
                                                            fig.type = fig.type, 
                                                            n.row = n.row, n.col = n.col,
                                                            height = height)
                                     }
                          )
                 )
  )
  )
}

output$MARKER.SETS <- renderUI({
  # fluidPage(
  do.call(tabBox,
          c(purrr::map(names(DATA$markers.list), plotMultiFeatures,
                       fig.type = input$MARKERS.TYPE),
            width = 12, height = "auto"))
  # tabBox(width = 12,
  #        tabPanel("Violin", uiOutput("MARKER.PLOTS")),
  #        tabPanel("Feature", uiOutput("MARKER.FEATURE"))
  #       )
  
})

output$MARKER.DOTPLOTS <- renderUI({
  req(DATA$RET, DATA$markers.list)
  do.call(what = shiny::tabsetPanel,
          args = purrr::map(names(DATA$markers.list), 
                            .f = function(markers.name) {
                              tabPanel(title=markers.name,
                                       shiny::renderPlot({
                                         DotPlot(DATA$orig.RET@seurat, 
                                                 features = c(DATA$markers.list[[markers.name]]))
                                       }, width=600,height=400))
                            }
          )
  )
})