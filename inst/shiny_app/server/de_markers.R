library(htmlwidgets)
observeEvent(input$DO.MARKERS, {
  if (!IdentsSaved(DATA$orig.RET)) {
    if (!is.null(ClosestIdent(DATA$orig.RET))) {
      DATA$orig.RET <- tryOverwriteIdent(DATA$orig.RET, ClosestIdent(DATA$orig.RET))
      updateGroupOptions()
    } else {
      showModal(newNameModal(warn = T))
    }
  }
  showModal({
    modalDialog(
      "Select groups for calculating DE genes (if none selected, all will be used--which may be very slow! (>30 min)).",
      tipify(selectInput("IDENTS.TO.USE", label = "Groups for DE gene calculation",
                         choices = levels(DATA$orig.RET@seurat), selectize = T, multiple = T),
             "Compares cells in each group to cells in ALL other groups, regardless of what you choose here!"),
      popify(sliderInput("DE.LOGFC", "Log-FC threshold",
                         min = 0, max = 2,
                         value = 0.5, step = 0.05), 
             "Limits testing to genes whose avg. expression in in-group cells differs from avg. expression in outgroup cells by at least this value, speeding up the testing (but potentially missing important signals!). 0.5 works well for finding markers between very distinct populations (e.g. B vs. T cells), but lower values are needed when comparing expression in e.g. tumor subsets."),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("DO.MARKERS.OK", "Ok"))
    )
  })
})
observeEvent(input$DO.MARKERS.OK, {
  req(DATA$orig.RET, DATA$orig.RET@markers)
  withProgress(message = 'DE markers', value = 0, {
    DATA$orig.RET <- seuratAllMarkers(DATA$orig.RET, shiny = T, 
                                      fast.threshold = input$DE.LOGFC, idents.use = input$IDENTS.TO.USE)
  })
  DATA$RET <- DATA$orig.RET 
  removeModal()
})


output$DE.PANEL <- renderUI({
  fluidPage(
    fluidRow(box(width = "100%",
                 height = "auto",
                 uiOutput("DE.MARKERS"))
    ))
})
makeDEMarkerPanel <- function(ident, only.pos = F) {
  markers.list <- BestMarkers(DATA$orig.RET, slot = input$SELECTED.DE.SLOT)
  top_markers <- markers.list[[ident]]
  if (only.pos) {
    top_markers <- dplyr::filter(top_markers, 
                                 avg_logFC > 0)
  }
  tabPanel(title=ident,
           renderDT({
             datatable(top_markers %>% 
               dplyr::arrange(p_val_adj, -avg_logFC) %>%
               dplyr::select(-p_val, -cluster) %>%
               mutate(avg_logFC = round(avg_logFC, 2)),
               options = list(autoWidth = T, scrollX = T)) %>% 
               formatSignif(columns = c("p_val_adj", "avg_logFC"), digits = 3)
           }),
           #, height = "550px"), 
           # options = list(pageLength = 10, autoWidth = T)),
           renderPlot({
             top_markers <- dplyr::filter(top_markers, 
                                          avg_logFC > 0)
             if (nrow(top_markers) > 25) {
               markers <- top_markers$gene[1:25]
             } else {
               markers <- top_markers$gene
             }
             gcdDoHeatmap(DATA$orig.RET@seurat, 
                          cells = WhichCells(DATA$orig.RET@seurat, downsample = input$DIFFEXP.HEATMAP.CELLS),
                          cols = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE)),
                          features = markers, slot = input$DATA.SLOT)
           }, height = 650, width = "auto")
  )
}


output$DE.MARKERS.SLOT <- renderUI({
  req(DATA$orig.RET, ActiveMarkers(DATA$orig.RET))
  
  marker.slots <- ActiveMarkers(DATA$orig.RET)
  selectInput("SELECTED.DE.SLOT", 
              "Marker slots", marker.slots, multiple = F,
              selectize = T, selected = BestMarkers(DATA$orig.RET))
})

output$DE.MARKERS <- renderUI({
  
  req(DATA$orig.RET,  isTruthy(BestMarkers(DATA$orig.RET)), input$SELECTED.DE.SLOT)
  
  markers.list <- BestMarkers(DATA$orig.RET, slot = input$SELECTED.DE.SLOT)
  # req(length(intersect(names(markers.list), levels(DATA$orig.RET@seurat))) >= 1)
  
  do.call(tabBox,
          c(width = 12,
            height = "auto",
            purrr::map(names(markers.list),
                       makeDEMarkerPanel,
                       only.pos = input$DE.ONLY.POS
                       )
          )
  )
})