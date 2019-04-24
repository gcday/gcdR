
observeEvent(input$MARKERS.LIST.PATH$datapath, {
  DATA$markers.list <- readMarkersYAML(input$MARKERS.LIST.PATH$datapath)
})

multiPanelSeuratFigs <- function(features, fig.type = "violin", n.row = 2, 
                                 n.col = 2, height = 600, alt.names = NULL, vln.show.dots = F) {
  plots.list <- NULL
  if (fig.type == "feature") {
    plots.list <- FeaturePlot(DATA$RET@seurat,
                              features = c(features),
                              reduction = input$DIM.REDUC,
                              combine = F)
  } else {
    if (fig.type == "violin") {
      plots.list <- GCD.VlnPlot(DATA$RET@seurat,
                            features = c(features),
                            cols = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE)), 
                            pt.size = ifelse(vln.show.dots, 1, 0),
                            combine = F, slot = input$DATA.SLOT)
    } else if (fig.type == "ridge") {
      plots.list <- GCD.RidgePlot(DATA$RET@seurat,
                            features = c(features),
                            cols = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE)),
                            combine = F, slot = input$DATA.SLOT)
    }
    else {
      warning("Failed to create multi-panel plot, invalid fig.type ", fig.type)
      return(NULL)
    }
    for (i in 1:length(plots.list)) {
      plots.list[[i]] <- plots.list[[i]] + NoLegend()
    }
  }
  if (!is.null(alt.names)) {
    for (i in 1:length(features)) {
      if (features[[i]] %in% names(alt.names)) {
        plots.list[[i]] <- plots.list[[i]] + 
          ggtitle(paste0(features[[i]], " (", alt.names[[features[[i]]]], ")"))
      }
    }
  }
  return(renderPlot({
    cowplot::plot_grid(plotlist = plots.list,
                       nrow = n.row, ncol = n.col)
  }, width="auto",height=height))
}

plotMultiFeatures <- function(markers.name, fig.type = "violin", n.row = 2, n.col = 2) {
  markers <- DATA$markers.list$markers.list[[markers.name]]
  markers <- markers[markers %in% rownames(DATA$RET@seurat)]
  if (length(markers) == 0) {
    return(NULL)
  }
  if (fig.type == "heatmap") {
    return(tabPanel(title = markers.name,
           renderPlot({
      gcdDoHeatmap(DATA$orig.RET@seurat, 
                   cells = WhichCells(DATA$orig.RET@seurat, 
                                     downsample = input$DIFFEXP.HEATMAP.CELLS),
                   cols = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE)),
                   features = markers, slot = input$DATA.SLOT)},
      width="auto",height=900))
    )
  } else if (fig.type == "dotplot") {
    return(tabPanel(title = markers.name,
                    renderPlot({
                      DotPlot(DATA$orig.RET@seurat, 
                              features = markers)},
                      width="auto",height=100 + 30 * length(levels(DATA$orig.RET@seurat))))
    )
  }
  split.markers <- split(markers, 
                         ceiling(seq_along(markers)/(n.row * n.col)))
  if (fig.type == "violin") {
    height <- 600
  } else if (fig.type == "feature") {
    height <- 900
  } else if (fig.type == "ridge") {
    height <- 900
  }
  return(do.call(tabPanel,
                 args = c(title = markers.name,
                          purrr::map(names(split.markers),
                                     .f = function(chunk.name){
                                       multiPanelSeuratFigs(split.markers[[chunk.name]],
                                                            fig.type = fig.type, 
                                                            n.row = n.row, n.col = n.col,
                                                            height = height, 
                                                            alt.names = DATA$markers.list$alt.names,
                                                            vln.show.dots = input$VIOLIN.SHOW.DOTS)
                                     }
                          )
                 )
  )
  )
}

output$MARKER.SETS <- renderUI({
  # fluidPage(
  do.call(tabBox,
          c(purrr::map(names(DATA$markers.list$markers.list), plotMultiFeatures,
                       fig.type = input$MARKERS.TYPE),
            width = 12, height = "auto"))
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