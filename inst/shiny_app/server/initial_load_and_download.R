


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

updateSplitVals <- function() {
  req(input$SPLIT.BY, input$CLUSTER, DATA$RET)
  
  DATA$do.percent <- input$SPLIT.TYPE == "percent"
  split.vals <- as.factor(DATA$RET@seurat[[input$SPLIT.BY]][[input$SPLIT.BY]])
  DATA$split.dt.1 <- breakdownTable(split.vals, Idents(DATA$RET@seurat), transpose = T)
  DATA$split.dt.2 <- breakdownTable(Idents(DATA$RET@seurat), 
                                    split.vals)
  
  
  summary.table <- table(Idents(DATA$RET@seurat), split.vals)
  # if (input$SPLIT.TYPE == "percent") summary.table <- gcdR::percent.table(summary.table)
  
  
  
  revised.dt <- NULL
  for (i in 1:length(levels(DATA$RET@seurat))) {
    for (j in 1:length(levels(split.vals))) {
      val <- summary.table[i + (j - 1) * length(levels(DATA$RET@seurat))]
      row <- levels(DATA$RET@seurat)[i]
      col <- levels(split.vals)[j]
      new.dt <- data.frame(row, col, val)
      revised.dt <- rbind(revised.dt, new.dt)
    }
  }
  colnames(revised.dt) <- c(input$CLUSTER, input$SPLIT.BY, "count")
  # DATA$split.datatable <-
  
  # var.1 <- 

  DATA$split.datatable <- revised.dt
}

openInitialRDS <- function(object) {

  RET <- readSeuratRDStoGCD(object = object, shiny = T)
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
    # DATA$orig.RET <- setActiveIdent(DATA$orig.RET)
    # choices <- CellGroups(DATA$orig.RET)
    # split.choices <- CellGroups(DATA$orig.RET)
    # [choices != ActiveIdent(DATA$orig.RET)]
    # split.selected <- ifelse(isTruthy(input$SPLIT.BY), input$SPLIT.BY, split.choices[1])
    updateGroupOptions()
  }
  updateSplitVals()
  message("finished!", ActiveIdent(DATA$orig.RET))
}

observeEvent(input$IMAGE$datapath, {
  withProgress(message = 'Loading', value = 0, {
    openInitialRDS(object = readRDS(file = input$IMAGE$datapath ))
  })
})
volumes <- getOption(x = "gcdR.cellViewer.volumes", 
                     default = c(Home = fs::path_home(),
                                 shinyFiles::getVolumes()()))
defaultRoot <- getOption(x = "gcdR.cellViewer.defaultRoot", 
                        default = NULL)
defaultPath <- getOption(x = "gcdR.cellViewer.defaultPath", 
                         default = "")
                        

             
             # OneDrive = "/mnt/c/Users/grady/OneDrive - Leland Stanford Junior University",
             # Mayo = "/mnt/c/Users/grady/OneDrive - Leland Stanford Junior University/Mayo/Mouse sc-RNAseq/",
             # getVolumes()())
shinyFileChoose(input,'file', session=session, roots=volumes, 
                defaultRoot = defaultRoot,
                defaultPath = defaultPath,
                filetypes = c("rds")
)
shinyFileSave(input, "serverside_save", session = session, 
              roots = volumes, 
              defaultRoot = defaultRoot,
              defaultPath = defaultPath)
              
              # defaultRoot = "OneDrive",
              # defaultPath = "Levy Lab/scRNA-seq/")


# shinyFileChoose(input, 'LOCAL.FILE', roots = volumes, session = session)
observeEvent(input$file, {
  # req(isTruthy(input$file))
  inFile <- parseFilePaths(roots=volumes, input$file)
  # print(input$file)
  # print(inFile)
  req(isTruthy(inFile$datapath))
  withProgress(message = 'Loading', value = 0, {
    openInitialRDS(object = readRDS(file = inFile$datapath))
  })
  
}, ignoreInit = T, ignoreNULL = T)




observeEvent(input$serverside_save, {
  # req(isTruthy(input$file))
  outFile <- parseSavePath(roots=volumes, input$serverside_save)
  # print(input$file)
  # print(inFile)
  req(isTruthy(outFile$datapath))
  withProgress(message = 'Saving', value = 0, {
    saveRDS(DATA$orig.RET, outFile$datapath, compress = F)
    # openInitialRDS(inFile$datapath)
  })
  
}, ignoreInit = T, ignoreNULL = T)
#
observeEvent(DATA$INPUT.OBJECT, ignoreNULL = T,
  {
  req(DATA$INPUT.OBJECT)
  withProgress(message = 'Loading', value = 0, {
    openInitialRDS(object = DATA$INPUT.OBJECT)
    DATA$INPUT.OBJECT <- NULL
  })
})


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