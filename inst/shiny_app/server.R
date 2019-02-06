library(gcdR)
mydebug <- function(msg="[DEBUG]") {
  DEBUG <- FALSE
  if (DEBUG) {
    print(sprintf("%s - %s - %s", msg1, as.character(Sys.time()), as.character(deparse(sys.calls()[[sys.nframe()-1]]))))
  }
}


server <- function( input, output, session){
  DATA <- reactiveValues() 
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

  newNameModal <- function(failed = FALSE, warn = FALSE) {
    req(DATA$orig.RET)
    print(ClosestIdent(DATA$orig.RET))
    print(CellGroups(DATA$orig.RET))
    modalDialog(
      
      selectInput("OVERWRITE.IDENT", "Select an ident column to overwrite with current group assignments. ", 
                  choices = c(None = '', CellGroups(DATA$orig.RET), selected = ClosestIdent(DATA$orig.RET)),
                   multiple = F),
      textInput("NEW.IDENTS", "Or, enter new name for current group assignments.",
                placeholder = paste0('new_', input$CLUSTER)),
      if (failed)
        div(tags$b("Group name already exists!", style = "color: red;")),
      if (warn)
        div(tags$b("Current idents have been modified! Please save before continuing or they will be lost.", style = "color: orange;")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("OK.SAVE.GROUPS", "OK")
      )
    )
  }
  observeEvent(input$OK.SAVE.GROUPS, {
    req(isTruthy(input$NEW.IDENTS) | isTruthy(input$OVERWRITE.IDENT))
    
    if (input$OVERWRITE.IDENT != "") {
      DATA$orig.RET <- tryOverwriteIdent(DATA$orig.RET, input$OVERWRITE.IDENT)
      DATA$RET <- DATA$orig.RET
      updateGroupOptions()
      removeModal()
    } else if (input$NEW.IDENTS != "" & !input$NEW.IDENTS %in% colnames(DATA$orig.RET@seurat@meta.data)) {
      DATA$orig.RET@seurat[[input$NEW.IDENTS]] <- Idents(DATA$orig.RET@seurat)
      for (marker.type in names(DATA$orig.RET@markers)) {
        DATA$orig.RET@markers[[marker.type]][[input$NEW.IDENTS]] <- DATA$orig.RET@markers[[marker.type]][[ActiveIdent(DATA$orig.RET)]]
      }
      DATA$RET <- DATA$orig.RET
      updateGroupOptions()
      removeModal()
    } else {
      showModal(newNameModal(failed = TRUE))
    }
  })
  
  observeEvent(input$RENAME.IDENT, {
    showModal(newNameModal())
  })
  
  
  newGroupNameModal <- function(failed = FALSE) {
    modalDialog(
      selectInput("GRPS.TO.RENAME", label = "Group to rename",
                  choices = levels(DATA$orig.RET@seurat), selectize = T, multiple = F
      ),
      textInput("NEW.GRPS.NAME", label = "New label for selected group"),
      if (failed)
        div(tags$b("Error reassigning ident name.", style = "color: red;")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ok.grp.rename", "OK")
      )
    )
  }
  observeEvent(input$CHANGE.GRP.NAME, {
    showModal(newGroupNameModal())
  })
  
  observeEvent(input$ok.grp.rename, {
    req(input$NEW.GRPS.NAME)
    if (length(input$GRPS.TO.RENAME) >= 1) {
      print("doing it")
      DATA$orig.RET <- gcdRenameIdent(DATA$orig.RET, input$GRPS.TO.RENAME, input$NEW.GRPS.NAME)
      print(levels(DATA$orig.RET@seurat))
      DATA$RET <- DATA$orig.RET
      removeModal()
    } else {
      showModal(newGroupNameModal(failed = TRUE))
    }
  })
  
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
                  choices = levels(DATA$orig.RET@seurat), selectize = T, multiple = T),  "Compares cells in each group to cells in ALL other groups, regardless of what you choose here!"),
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
  
  output$DOWNLOAD.DATA <- downloadHandler(filename = function() {return(input$IMAGE$name)},
    content = function(file) {saveRDS(DATA$orig.RET, file)})
  

  output$NAMES <- renderUI({
    req(input$FILTER)
    GROUP.NAMES <- levels(DATA$orig.RET@seurat@meta.data[,input$FILTER])
    selectInput("SELECTION", "Filter by Sample/Cluster/Type", GROUP.NAMES, multiple = T, selectize = T)
  })
  
  output$DIM.REDUC.CHOICE <- renderUI({
    req(DATA$orig.RET)
    radioButtons("DIM.REDUC", "Reduction",
                 choices = DATA$orig.RET@meta.list$reductions, selected = DATA$orig.RET@meta.list$reductions[1])
  })

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
        message(levels(as.factor(type.convert(DATA$orig.RET@seurat@meta.data[[input$CLUSTER]]))))
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
      print(length(input$SELECTION))
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
                       reduction = input$DIM.REDUC)
    return(dim.plt)
  })
  

  doPlotsSingleGene <- function(gene, vln.width = "auto", vln.height = 400,
                                feat.width = "auto", feat.height = 500,
                                table.height = "200px") {
    tabPanel(title = gene,
      fluidRow(
        box(
          width = 5, status = "primary",
          renderPlot({
            VlnPlot(DATA$RET@seurat, features = c(gene),
                    cols = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE)), pt.size = 0) + NoLegend()
          }, width = vln.width, height = vln.height),
          style=paste0("padding-top:", (feat.height - vln.height ) / 2 ,"px;")
        ),
        box(
          status = "primary", width = 7,
          renderPlot({
            FeaturePlot(DATA$RET@seurat, features = c(gene),
                        reduction = input$DIM.REDUC)
          }, width=feat.width, height=feat.height)
        )
      ),
      fluidRow(
        box(width = 12, status = "primary",
          div(style = 'overflow-x: scroll; overflow-y: scroll',
            renderDT({geneSummaryTable(gene)},
                   height = table.height,
                   options = list(paging = F, autoWidth = T,
                                  searching = F, info = F,
                                  ordering = F))
          )
        )
      )
    )
  }

  output$GENE.PLOTS <- renderUI({
    req(input$GENE, DATA$RET)
    do.call(tabBox,
              c(purrr::map(input$GENE, doPlotsSingleGene),
                width = 12, height = "auto"))
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
  
  breakdownTable <- function(var.1, var.2, transpose = F) {
    if (transpose) {
      MAT <- as.data.frame.matrix(t(percent.table(table(var.1, var.2))))
    } else {
      MAT <- as.data.frame.matrix(percent.table(table(var.1, var.2)))
    }
    new.MAT <- dplyr::mutate_all(MAT, funs(scales::percent(., accuracy = 0.1, scale = 1)))
    rownames(new.MAT) <- rownames(MAT)
    return(new.MAT)
  }
    
  output$SPLIT.SUMMARY.1 <- renderUI({
    req(DATA$RET, input$SPLIT.BY)
    split.vals <- DATA$orig.RET@seurat[[input$SPLIT.BY]][[input$SPLIT.BY]]
    dt.1 <- breakdownTable(split.vals, Idents(DATA$RET@seurat), T)
    dt.2 <- breakdownTable(Idents(DATA$RET@seurat), split.vals)
    fluidPage(
      fluidRow(
        box(title = paste0(input$SPLIT.BY, " makeup across ", input$CLUSTER),
            width = NULL,
            height = "auto",
            div(style = 'overflow-x: scroll; overflow-y: scroll', 
              renderDT(dt.1, 
                     height = "auto",
                     options = list(paging = F, autoWidth = T, searching = F, 
                                    info = F, ordering = F)
              )
            )
        )
      ),
      fluidRow(
        box(title = paste0(input$CLUSTER, " makeup across ", input$SPLIT.BY),
            width = NULL,
            height = "auto",
            div(style = 'overflow-x: scroll; overflow-y: scroll',
                renderDT(dt.2, 
                         height = "auto",
                         options = list(paging = F, autoWidth = T, searching = F, 
                                        info = F, ordering = F)
                )
            )
        )
      )
    )
  })
  
  geneSummaryTable <- function(gene) {
    EXPR.VALS <- FetchData(DATA$RET@seurat, vars = c(gene))[[gene]]
    GROUP.NAMES <- levels(Idents(DATA$RET@seurat))
    IDENTS <- Idents(DATA$RET@seurat)
    CELL.COUNTS <- table( IDENTS )
    CELL.COUNTS <- CELL.COUNTS[match( GROUP.NAMES, names( CELL.COUNTS ))]
    EXPR.CELLS <- tapply( EXPR.VALS , IDENTS , function(x) sum( x > 0 ))
    EXPR.CELLS <- EXPR.CELLS[match( GROUP.NAMES, names( EXPR.CELLS ))]
    EXPR.FRAC <- scales::percent(c(EXPR.CELLS / CELL.COUNTS), accuracy = 0.1)
    CELLS.FRAC <- scales::percent(c(EXPR.CELLS / sum(EXPR.CELLS)), accuracy = 0.1)
    MEAN.EXPR <- format( round( tapply( EXPR.VALS , IDENTS , mean ) ,2), nsmall=2)
    MEAN.EXPR <- MEAN.EXPR[match( GROUP.NAMES, names( MEAN.EXPR ))]
    XXX <- rbind(CELL.COUNTS, CELLS.FRAC, EXPR.CELLS, EXPR.FRAC, MEAN.EXPR)
    rownames( XXX ) <- c('Cell Count', 'Pct of Total Cells', 'Expressing Cell Count', 'Pct of Expressing Cell',
                         'Mean Expression All Cells')
    return(XXX)
  }
  output$DE.PANEL <- renderUI({
    fluidPage(
      fluidRow(box(width = "100%",
          height = "auto",
          uiOutput("DE.MARKERS"))
      ))
  })

  output$DE.MARKERS <- renderUI({
    req(DATA$orig.RET, DATA$orig.RET@markers$quick[[input$CLUSTER]])

    markers.list <- DATA$orig.RET@markers$quick[[ActiveIdent(DATA$orig.RET)]]
    req(length(intersect(names(markers.list), levels(DATA$orig.RET@seurat))) >= 1)
    
    do.call(tabBox,
            c(width = 12,
              height = "auto",
            purrr::map(names(markers.list),
                              function(ident){
                  tabPanel(title=ident,
                       renderDT({
                         dplyr::filter(markers.list[[ident]], 
                                       avg_logFC > 0) %>% 
                           dplyr::arrange(p_val_adj, -avg_logFC) %>%
                           dplyr::select(-p_val, -cluster) %>%
                           mutate(avg_logFC = round(avg_logFC, 2),
                                  p_val_adj = formatC(p_val_adj, format = "e", digits = 2))
                       }, height = "550px", options = list(pageLength = 10, autoWidth = T)),
                       renderPlot({
                         top_markers <- dplyr::filter(markers.list[[ident]], 
                                       avg_logFC > 0) %>% 
                           dplyr::arrange(p_val_adj, -avg_logFC)
                         if (nrow(top_markers) > 25) {
                           markers <- top_markers$gene[1:25]
                         } else {
                           markers <- top_markers$gene
                         }
                         gcdDoHeatmap(DATA$orig.RET@seurat, 
                                   cells = WhichCells(DATA$orig.RET@seurat, downsample = 250),
                                   cols = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE)),
                                   features = markers, slot = "data")
                       }, height = 650, width = "auto")
                       )
            }))
    )
  })
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
                              .f = function(markers.name){
                                tabPanel(title=markers.name,
                                         shiny::renderPlot({
                                           DotPlot(DATA$orig.RET@seurat, 
                                                   features = c(DATA$markers.list[[markers.name]]))
                                         }, width=600,height=400))
                              }
            )
    )
  })
}
