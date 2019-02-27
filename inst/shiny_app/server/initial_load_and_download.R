addCellCycleScoring <- function(RET) {
  cc.genes.species <- cc.genes
  if (DATA$mouse) {
    cc.genes.species <- mouse.cc.genes
  }
  if (!"Phase" %in% colnames( RET@seurat@meta.data )) {
    incProgress(amount = 0.2, message = "Adding cell cycle scores...")
    RET <- gcdCellCycleScoring(RET, cc.genes.species$s.genes, 
                               cc.genes.species$g2m.genes, set.ident = F)
  }
  return(RET)
}
updateGroupOptions <- function() {
  choices <- CellGroups(DATA$orig.RET)
  updateSelectizeInput(session, 'CLUSTER', 'Cluster by',
                       choices = choices,
                       selected = ActiveIdent(DATA$orig.RET), server = T)
  updateSelectizeInput(session, 'FILTER', 'Filter by',
                       choices = c(None="", choices),
                       selected = input$FILTER, server = T)
  split.choices <- choices[choices != ActiveIdent(DATA$orig.RET)]
  
  updateSelectizeInput(session, 'SPLIT.BY', 'Optional split value',
                       choices = c(None="", split.choices),
                       selected = ifelse(isTruthy(input$SPLIT.BY), input$SPLIT.BY, split.choices[1]))
}

openInitialRDS <- function() {
  RET <- readSeuratRDStoGCD(input$IMAGE$datapath, shiny = T)
  RET <- gcdRunUMAP(RET)
  DATA$mouse <- "Sdc1" %in% rownames(RET@seurat)
  
  RET <- addCellCycleScoring(RET)
  
  DATA$orig.RET <- RET
  
  DATA$RET <- RET
  DATA$ACTIVE.FILTER <- F
  
  numeric.meta.data <- colnames(DATA$orig.RET@seurat@meta.data)[sapply(DATA$orig.RET@seurat@meta.data, function(x) is.numeric(x))]
  
  updateSelectizeInput(session, 'GENE','Gene', 
                       choices = list("Metadata" = numeric.meta.data,
                                      "Genes" = sort(rownames(DATA$orig.RET@seurat))),
                       selected = numeric.meta.data[1], server = T)
  if (!IdentsSaved(DATA$orig.RET)) {
    if (!is.null(ClosestIdent(DATA$orig.RET))) {
      DATA$orig.RET <- tryOverwriteIdent(DATA$orig.RET, ClosestIdent(DATA$orig.RET))
      updateGroupOptions()
    } else {
      showModal(newNameModal())
    }
  } else {
    DATA$orig.RET <- setActiveIdent(DATA$orig.RET)
    updateGroupOptions()
  }
  message(ActiveIdent(DATA$orig.RET))
}

observeEvent(input$IMAGE$datapath, {
  withProgress(message = 'Loading', value = 0, {
    openInitialRDS()
  })
})

output$DOWNLOAD.DATA <- downloadHandler(filename = function() {return(input$IMAGE$name)},
                                        content = function(file) {
                                          saveRDS(DATA$orig.RET, file)
                                        }
)