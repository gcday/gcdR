


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

openInitialRDS <- function(rds.path) {
  RET <- readSeuratRDStoGCD(rds.path, shiny = T)
  RET <- gcdRunUMAP(RET)
  # DATA$mouse <- "Sdc1" %in% rownames(RET@seurat)
  DATA$species <- guessGCDSeuratSpecies(RET)
  
  RET <- guessSpeciesAddCellCycleScoring(RET)
  
  DATA$orig.RET <- RET
  
  DATA$RET <- RET
  DATA$ACTIVE.FILTER <- F
  
  numeric.meta.data <- colnames(DATA$orig.RET@seurat@meta.data)[sapply(DATA$orig.RET@seurat@meta.data, function(x) is.numeric(x))]
  
  updateSelectizeInput(session, 'GENE', 'Gene', 
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
    openInitialRDS(input$IMAGE$datapath )
  })
})
volumes <- c(Home = fs::path_home(), 
             OneDrive = "/mnt/c/Users/grady/OneDrive - Leland Stanford Junior University",
             Mayo = "/mnt/c/Users/grady/OneDrive - Leland Stanford Junior University/Mayo/Mouse sc-RNAseq/",
             getVolumes()())
shinyFileChoose(input,'file', session=session, roots=volumes, 
                defaultRoot = "OneDrive",
                defaultPath = "Levy Lab/scRNA-seq/",
                filetypes = c("rds")
)
shinyFileSave(input, "serverside_save", session = session, 
              roots = volumes, 
              defaultRoot = "OneDrive",
              defaultPath = "Levy Lab/scRNA-seq/")


# shinyFileChoose(input, 'LOCAL.FILE', roots = volumes, session = session)
observeEvent(input$file, {
  # req(isTruthy(input$file))
  inFile <- parseFilePaths(roots=volumes, input$file)
  # print(input$file)
  # print(inFile)
  req(isTruthy(inFile$datapath))
  withProgress(message = 'Loading', value = 0, {
    openInitialRDS(inFile$datapath)
  })
  
}, ignoreInit = T, ignoreNULL = T)




observeEvent(input$serverside_save, {
  # req(isTruthy(input$file))
  outFile <- parseSavePath(roots=volumes, input$serverside_save)
  # print(input$file)
  # print(inFile)
  req(isTruthy(outFile$datapath))
  withProgress(message = 'Saving', value = 0, {
    saveRDS(DATA$orig.RET, outFile$datapath)
    # openInitialRDS(inFile$datapath)
  })
  
}, ignoreInit = T, ignoreNULL = T)
#


# observeEvent(parseFilePaths(volumes, input$LOCAL.FILE), {
#   # (input$LOCAL.FILE, {
#   # 
#   # print(input$LOCAL.FILE)
#  
#   local.file <- parseFilePaths(volumes, input$LOCAL.FILE)
#   print(local.file)
#   req(local.file$datapath)
#   withProgress(message = 'Loading', value = 0, {
#     openInitialRDS(local.file$datapath)
#   })
# })

output$DOWNLOAD.DATA <- downloadHandler(filename = function() {return(input$IMAGE$name)},
                                        content = function(file) {
                                          saveRDS(DATA$orig.RET, file)
                                        }
)