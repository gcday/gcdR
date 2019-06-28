observeEvent(input$CLUSTER, ignoreInit = T, {
  req(input$CLUSTER, DATA$orig.RET)
  must.update <- F
  # means current group label assignments must be saved
  if (!IdentsSaved(DATA$orig.RET)) {
    if (!is.null(ClosestIdent(DATA$orig.RET))) {
      DATA$orig.RET <- tryOverwriteIdent(DATA$orig.RET, ClosestIdent(DATA$orig.RET))
    } else {
      must.update <- T
      DATA$orig.RET@seurat[[paste0("new_", ClosestIdent(DATA$orig.RET))]] <- Idents(object = DATA$orig.RET@seurat)
    }
  }
  
  if (input$CLUSTER != ActiveIdent(DATA$orig.RET)) {
    if(class(DATA$orig.RET@seurat[[input$CLUSTER]][[input$CLUSTER]]) == "character") {
      # Getting the sort order correct (e.g. 0, 1, 2,...) for numeric idents stored as chars, hopefully
      Idents(DATA$orig.RET@seurat) <- as.factor(type.convert(DATA$orig.RET@seurat@meta.data[[input$CLUSTER]]))
      # message(levels(as.factor(type.convert(DATA$orig.RET@seurat@meta.data[[input$CLUSTER]]))))
    } else if (class(DATA$orig.RET@seurat[[input$CLUSTER]][[input$CLUSTER]]) == "factor") {
      Idents(DATA$orig.RET@seurat) <- input$CLUSTER
    } else {
      Idents(DATA$orig.RET@seurat) <- type.convert(DATA$orig.RET@seurat@meta.data[[input$CLUSTER]], as.is = T)
    }
    message("Updating Seurat idents... ", input$CLUSTER)
    
    DATA$orig.RET@seurat <- BuildClusterTree(DATA$orig.RET@seurat, dims = DATA$orig.RET@meta.list$dims)
    DATA$RET <- DATA$orig.RET
    if (must.update) updateGroupOptions()
    
  } else {
    message("Unnecessary update attempt.")
  }
})

observeEvent(input$SELECTION, ignoreInit = T, {
  if (!is.null(input$SELECTION)) {
    DATA$ACTIVE.FILTER <- T
    # print(length(input$SELECTION))
    pass.cells <- which(DATA$orig.RET@seurat@meta.data[,input$FILTER] %in% input$SELECTION)
    CELLS.USE <- colnames(DATA$orig.RET@seurat)[pass.cells]
    DATA$RET@seurat <- subset(DATA$orig.RET@seurat, cells = CELLS.USE)
    message("Updating Seurat filtering...")
  } else {
    DATA$RET <- DATA$orig.RET
    DATA$ACTIVE.FILTER <- NULL
    message("Unnecessary filter update attempt.")
  }
}, ignoreNULL = F)

observeEvent(input$FILTER, {
  if (is.null(input$FILTER)) {
    if (!is.null(DATA$ACTIVE.FILTER)) {
      DATA$RET <- DATA$orig.RET
    }
    DATA$ACTIVE.FILTER <- NULL
  } else if (input$FILTER == "") {
    if (!is.null(DATA$ACTIVE.FILTER)) {
      DATA$RET <- DATA$orig.RET
    }
    DATA$ACTIVE.FILTER <- NULL
  } else {
    print("Unnecessary filter update attempt.")
  }
}, ignoreNULL = F)
output$DIM.REDUC <- renderPlot({
  req(DATA$RET)
  dim.plt <- DimPlot(DATA$RET@seurat,
                     label = input$LABELS,
                     label.size = 6,
                     repel = T,
                     cols = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE)),
                     reduction = input$DIM.REDUC)
  return(dim.plt)
})

observeEvent(input$SPLIT.BY, ignoreNULL = T, ignoreInit = T, {
  
  updateSplitVals()
})
output$DIM.REDUC.SPLIT <- renderPlot({
  req(DATA$RET, input$SPLIT.BY)
  SPLIT.RET <- DATA$RET
  if(class(SPLIT.RET@seurat[[input$SPLIT.BY]][[input$SPLIT.BY]]) == "character") {
    # Getting the sort order correct (e.g. 0, 1, 2,...) for numeric idents stored as chars, hopefully
    Idents(SPLIT.RET@seurat) <- as.factor(type.convert(SPLIT.RET@seurat@meta.data[[input$SPLIT.BY]]))
    message(levels(as.factor(type.convert(SPLIT.RET@seurat@meta.data[[input$SPLIT.BY]]))))
  } else if (class(SPLIT.RET@seurat[[input$SPLIT.BY]][[input$SPLIT.BY]]) == "factor") {
    Idents(SPLIT.RET@seurat) <- input$SPLIT.BY
  } else {
    Idents(SPLIT.RET@seurat) <- type.convert(SPLIT.RET@seurat@meta.data[[input$SPLIT.BY]], as.is = T)
  }
  dim.plt <- DimPlot(SPLIT.RET@seurat,
                     group.by = input$SPLIT.BY,
                     reduction = input$DIM.REDUC,
                     cols = Palettes(SPLIT.RET, type.use = as.integer(input$SPLIT.PALETTE)))
  return(dim.plt)
})


output$SPLIT.SUMMARY.1 <- renderUI({
  req(DATA$RET, input$SPLIT.BY, input$CLUSTER, DATA$split.datatable, DATA$split.dt.1, DATA$split.dt.2)
  fluidPage(
    fluidRow(
      box(title = paste0(input$SPLIT.BY, " makeup across ", input$CLUSTER),
          width = NULL,
          height = "auto",
          div(style = 'overflow-x: scroll; overflow-y: scroll', 
              renderDT(DATA$split.dt.1, 
                       height = "auto",
                       options = list(paging = F, autoWidth = T, searching = F, 
                                      info = F, ordering = F)
              )
          )
      )
    ),
    fluidRow(
      renderPlot({
        plot <- ggplot(DATA$split.datatable,
                       mapping=aes_string(x=input$CLUSTER, y="count", fill=input$SPLIT.BY)[c(2, 3, 1)])
        return(plot + geom_bar(stat="identity", position = ifelse(DATA$do.percent, "fill", "stack")) +
          scale_fill_manual(values = Palettes(DATA$orig.RET, 
                                              as.integer(input$SPLIT.PALETTE), 
                                              var.use = input$SPLIT.BY)) + SeuratAxes() + NoGrid() + RotatedAxis())
        })
    ),
    fluidRow(
      box(title = paste0(input$CLUSTER, " makeup across ", input$SPLIT.BY),
          width = NULL,
          height = "auto",
          div(style = 'overflow-x: scroll; overflow-y: scroll',
              renderDT(DATA$split.dt.2, 
                       height = "auto",
                       options = list(paging = F, autoWidth = T, searching = F, 
                                      info = F, ordering = F)
              )
          )
      )
    ),
    fluidRow(
      renderPlot({
        plot <- ggplot(DATA$split.datatable,
                       aes_string(x=input$SPLIT.BY, 
                                  y="count", 
                                  fill = input$CLUSTER)[c(2, 3, 1)])
        return(plot + 
          geom_bar(stat="identity", position = ifelse(DATA$do.percent, "fill", "stack")) +
          scale_fill_manual(values = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE))) + 
          SeuratAxes() + NoGrid() + RotatedAxis())
        })
    )
  )
})

output$CLUSTER.TREE <- renderPlot({
  req(DATA$orig.RET)
  if (is.null(Tool(DATA$orig.RET@seurat, slot = "BuildClusterTree"))) {
    DATA$orig.RET@seurat <- BuildClusterTree(DATA$orig.RET@seurat, dims = DATA$orig.RET@meta.list$dims)
  }
  data.tree <- Tool(object = DATA$orig.RET@seurat, slot = "BuildClusterTree")
  ape::plot.phylo(x = data.tree, direction = "rightwards")
  ape::nodelabels()
})