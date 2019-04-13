reclusterModal <- function(failed = FALSE) {
  modalDialog(
    tipify(sliderInput("RECLUSTER.RES", "Resolution",
                       min = 0.1, max = 3,
                       value = 0.8, step = 0.05), "Higher resolution values result in more clusters"),
    if (failed)
      div(tags$b("Error reassigning ident name.", style = "color: red;")),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("OK.RECLUSTER", "OK")
    )
  )
}

observeEvent(input$OK.RECLUSTER, {
  if (input$RECLUSTER.RES) {
    withProgress(message = 'Clustering', value = 0, {
      DATA$orig.RET <- seuratClusterWrapper(DATA$orig.RET, resolution = input$RECLUSTER.RES, dims = DATA$orig.RET@meta.list$dims)
    })
    
    DATA$orig.RET <- setActiveIdent(DATA$orig.RET)
    updateGroupOptions()
    DATA$RET <- DATA$orig.RET
    removeModal()
  } else {
    showModal(reclusterModal(failed = TRUE))
  }
})

observeEvent(input$RECLUSTER, {
  if (!IdentsSaved(DATA$orig.RET)) {
    if (!is.null(ClosestIdent(DATA$orig.RET))) {
      DATA$orig.RET <- tryOverwriteIdent(DATA$orig.RET, ClosestIdent(DATA$orig.RET))
    } else {
      showModal(newNameModal(warn = T))
    }
  } else {
    showModal(reclusterModal())
  }
})


