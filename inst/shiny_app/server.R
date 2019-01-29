library(gcdR)
mydebug <- function(msg="[DEBUG]") {
  DEBUG <- FALSE
  if (DEBUG) {
    print(sprintf("%s - %s - %s", msg1, as.character(Sys.time()), as.character(deparse(sys.calls()[[sys.nframe()-1]]))))
  }
}


server <- function( input, output, session){
  DATA <- reactiveValues() 
  
  updateSeuratAndRunUMAP <- function(SRT) {
    # dims = NULL
    dims <- NULL
    if (utils::compareVersion(as.character(SRT@version), "3.0") != 1) {
      # pre Seurat 3.0, must update
      message("Updating Seurat object to 3.0...")
      dims <- c(SRT@calc.params$RunTSNE$dims.use)
      SRT <- UpdateSeuratObject(SRT)
    } else {
      for (name in names(SRT@commands)) {
        message("command name: ", name)
        if ("dims" %in% names(SRT@commands[[name]]@params)) {
          dims <- SRT@commands[[name]]@params$dims
          message("Using dims (", dims, ") from ", name, " command.")
          break
        }
      }
      if (is.null(dims)) {
        dims <- c(1:10)
        message("Couldn't find dims in command history, using ", dims, " instead.")
      }
    }
    UMAP.present <- T
    if (!"umap" %in% names(SRT@reductions)) {
      if (reticulate::py_module_available(module = 'umap')) {
        SRT <- RunUMAP(SRT, dims = c(dims))
      } else {
        showModal(modalDialog(title = "Warning!",
        "Missing umap package, only t-SNE will be available. Please install through pip (e.g. pip install umap-learn).",
        easyClose = T))
        UMAP.present <- F
      }
    }
    DATA$reductions <- ifelse(UMAP.present, c(UMAP = "umap", `t-SNE` = "tsne"), c(`t-SNE` = "tsne"))
    DATA$dims <- dims
    SRT <- BuildClusterTree(SRT, dims = dims)
    return(SRT)
  }
  
  addCellCycleScoring <- function(RET) {
    cc.genes.species <- cc.genes
    if (DATA$mouse) {
      cc.genes.species <- mouse.cc.genes
    }
    if (!"Phase" %in% colnames( RET@seurat@meta.data )) {
      RET <- gcdCellCycleScoring(RET, cc.genes.species$s.genes, 
                                 cc.genes.species$g2m.genes, set.ident = F)
    }
    return(RET)
  }
  
  openInitialRDS <- function() {
    IN.OBJ <- readRDS(input$IMAGE$datapath)
    RET <- new(Class = "gcdSeurat")
    if (class(IN.OBJ) == "gcdSeurat") {
      RET <- IN.OBJ
    } else {
      RET@seurat <- IN.OBJ
    }
    RET@seurat <- updateSeuratAndRunUMAP(RET@seurat)
    DATA$mouse <- "Sdc1" %in% rownames(RET@seurat)
    
    RET <- addCellCycleScoring(RET)

    COUNTS <- apply( RET@seurat@meta.data, 2, function(x) length(table(x)))
    DATA$CHOICES <- colnames( RET@seurat@meta.data )[COUNTS<=100]
    selected.ident <- NULL
    for (choice in DATA$CHOICES) {
      if(isTRUE(all.equal(as.character(Idents(RET@seurat)),
                          as.character(RET@seurat@meta.data[[choice]])))) {
        selected.ident <- choice
        break
      }
    }
    DATA$selected.ident <- selected.ident
    # if ("Sdc1" %in% rownames(RET@seurat))
    
    numeric.meta.data <- colnames(RET@seurat@meta.data)[sapply(RET@seurat@meta.data, function(x) is.numeric(x))]
    
    DATA$GENES.OPT <- c(numeric.meta.data, sort(rownames(RET@seurat)))
    # radioButtons("DIM.REDUC", "Reduction",  
    #              choices = c(`t-SNE` = "tsne"), selected = "tsne"),
    # if ("umap" %in% names(RET@seurat@reductions)) {
    #   updateRadioButtons(session, "DIM.REDUC", "Reduction",  
    #                      c(UMAP = "umap", `t-SNE` = "tsne"), selected = "umap")
    # }
    
    updateSelectizeInput(session, 'GENE','Gene', 
                         choices = list("Metadata" = numeric.meta.data,
                                        "Genes" = sort(rownames(RET@seurat))),
                                        selected = numeric.meta.data[1], server = T)
    updateSelectizeInput(session, 'CLUSTER','Cluster by', 
                         choices = DATA$CHOICES, 
                         selected = DATA$selected.ident, server = T)
    updateSelectizeInput(session, 'FILTER','Filter by', 
                         choices = c(None="", DATA$CHOICES),
                         selected = NULL, server = T)
    split.choices <- DATA$CHOICES[DATA$CHOICES != input$CLUSTER]
    
    selected.split <- ifelse(length(split.choices)>= 1, split.choices[1], NULL)
    
    updateSelectizeInput(session, 'SPLIT.BY', 'Optional split value', 
                         choices = c(None="", split.choices), 
                         selected = selected.split)
    
      
    DATA$orig.RET <- RET
    DATA$RET <- RET
    DATA$ACTIVE.FILTER <- F
    DATA$COLS <- rainbow_hcl(length(levels(Idents(DATA$orig.RET@seurat))), c = 100, l = 80) 
    
    #DATA$COLS <- list(hue_pal()(length(levels(Idents(DATA$orig.RET@seurat)))),
                      # rainbow_hcl(length(levels(Idents(DATA$orig.RET@seurat))), c = 90, l = 80),
                      # rainbow_hcl(length(levels(Idents(DATA$orig.RET@seurat))), c = 100, l = 80))
    message(DATA$selected.ident)
  }
  
  observeEvent(input$IMAGE$datapath, {openInitialRDS()})

  newNameModal <- function(failed = FALSE) {
    modalDialog(
      textInput("NEW.IDENTS", "Enter new name for current group assignments.",
                placeholder = paste0('new_', input$CLUSTER)),
      if (failed)
        div(tags$b("Group name already exists!", style = "color: red;")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("OK.SAVE.GROUPS", "OK")
      )
    )
  }
  observeEvent(input$OK.SAVE.GROUPS, {
    if (input$NEW.IDENTS != "" & !input$NEW.IDENTS %in% colnames(DATA$orig.RET@seurat@meta.data)) {
      DATA$orig.RET@seurat[[input$NEW.IDENTS]] <- Idents(DATA$orig.RET@seurat)
      DATA$CHOICES <- c(input$NEW.IDENTS, DATA$CHOICES)
      DATA$selected.ident <- input$NEW.IDENTS
      updateSelectInput(session, 'CLUSTER', 'Cluster by', 
                        choices = DATA$CHOICES, selected = DATA$selected.ident, server = T)
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
      selectInput("GRPS.TO.RENAME", label = "Group(s) to rename",
                  choices = levels(DATA$orig.RET@seurat), selectize = T, multiple = T
      ),
      textInput("NEW.GRPS.NAME", label = "New label for selected group(s)"),
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
      DATA$orig.RET <- gcdRenameIdents(DATA$orig.RET, input$GRPS.TO.RENAME, input$NEW.GRPS.NAME)
      print(levels(DATA$orig.RET@seurat))
      DATA$orig.RET@seurat <- BuildClusterTree(DATA$orig.RET@seurat, dims = DATA$dims)
      DATA$RET <- DATA$orig.RET
      removeModal()
    } else {
      showModal(newGroupNameModal(failed = TRUE))
    }
  })
  observeEvent(input$DO.MARKERS, {
    req(DATA$orig.RET, DATA$orig.RET@markers)
    withProgress(message = 'DE markers', value = 0, {
      DATA$orig.RET@markers$all.markers.quick <- FindAllMarkersShiny(DATA$orig.RET@seurat, 
                                                                     logfc.threshold = 0.5)
    })
    # DATA$orig.RET <- seuratAllMarkers(DATA$orig.RET)
  })
  
  output$DOWNLOAD.DATA <- downloadHandler(filename = function() {return(input$IMAGE$name)},
    content = function(file) {saveRDS(DATA$orig.RET, file)})
  
  output$SPLIT.VAR <- renderUI({
    if (is.null( input$IMAGE ) | is.null(input$CLUSTER) ) return(NULL)
    message("updating SPLIT.BY")
    choices <- DATA$CHOICES[DATA$CHOICES != input$CLUSTER]
    selected.split <- ifelse(length(choices) >= 1, choices[1], NULL)
    selectInput('SPLIT.BY', 'Optional split value', choices = c(None="", choices), selected = selected.split, selectize = T)
  })

  output$NAMES <- renderUI({
    req(input$FILTER)
    GROUP.NAMES <- levels(DATA$orig.RET@seurat@meta.data[,input$FILTER])
    selectInput("SELECTION", "Filter by Sample/Cluster/Type", GROUP.NAMES, multiple = T, selectize = T)
  })
  # output$MAXX <- renderUI({
  #   req(input$EXPR.LIM)
  #   if (input$EXPR.LIM) {
  #     sliderInput("MAXX", "Log-Scale Maximum",
  #                 min = 0, max = 10,
  #                 value = 5, step = 0.05)
  #   }
  # })
  output$DIM.REDUC.CHOICE <- renderUI({
    req(DATA$reductions)
    radioButtons("DIM.REDUC", "Reduction",  
                 choices = c(DATA$reductions), selected = DATA$reductions[1])
  })

  observeEvent(input$CLUSTER, ignoreInit = T, {
    if(!isTRUE(all.equal(as.character(Idents(DATA$orig.RET@seurat)),
                        as.character(DATA$orig.RET@seurat@meta.data[[input$CLUSTER]])))) {
      if(class(DATA$orig.RET@seurat[[input$CLUSTER]][[input$CLUSTER]]) == "character") {
        # Getting the sort order correct (e.g. 0, 1, 2,...) for numeric idents, hopefully
        Idents(DATA$orig.RET@seurat) <- as.factor(type.convert(DATA$orig.RET@seurat@meta.data[[input$CLUSTER]]))
      } else {
        Idents(DATA$orig.RET@seurat) <- type.convert(DATA$orig.RET@seurat@meta.data[[input$CLUSTER]], as.is = T)
      }
      DATA$COLS <- rainbow_hcl(length(levels(Idents(DATA$orig.RET@seurat))), c = 100, l = 80) 
      message("Updating Seurat idents...")
      DATA$RET <- DATA$orig.RET
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
                       cols = DATA$COLS, 
                       reduction = input$DIM.REDUC)
    return(dim.plt)
  })
  # output$MULTIFEATURE.PLOT <- renderUI({
  #   req(input$GENE, DATA$RET)
  #   do.call(what = shiny::fluidPage,
  #           args = purrr::map(input$GENE, .f = function(gene){
  #             fluidRow(title=gene,
  #                      renderPlot({
  #                        FeaturePlot(DATA$RET@seurat, features = c(gene),
  #                                    reduction = input$DIM.REDUC)
  #                        }, width=800,height=600))
  #           })
  #   )
  # })
  # output$MULTIGENE.PLOT <- renderUI({
  #   req(input$GENE, DATA$RET)
  #   do.call(what = shiny::fluidPage,
  #           args = purrr::map(input$GENE, .f = function(gene){
  #             shiny::fluidRow(title=gene,
  #                             shiny::renderPlot({
  #                               VlnPlot(DATA$RET@seurat, features = c(gene),
  #                                       cols = DATA$COLS, pt.size = 0) + NoLegend()
  #                               }, width=600,height=400))
  #             })
  #   )
  # })
  doPlotsSingleGene <- function(gene, vln.width = "auto", vln.height = 400,
                                feat.width = "auto", feat.height = 500, 
                                table.height = "200px") {
    tabPanel(title = gene,
      fluidRow(
        box(
          width = 5, status = "primary",
          renderPlot({
            VlnPlot(DATA$RET@seurat, features = c(gene),
                    cols = DATA$COLS, pt.size = 0) + NoLegend()
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
                width = 12, height = "900px"))
  })
  
  # output$MULTIGENE.PLOT <- renderUI({
  #   req(input$GENE, DATA$RET)
  #   do.call(shiny::fluidPage,
  #           c(tags$head(
  #             tags$style(HTML('overflow-x: scroll; overflow-y: scroll'))),
  #             purrr::map(input$GENE, doPlotsSingleGene)))
  #                           # purrr::map(input$GENE, .f = doPlotsSingleGene))))
  #           # args = purrr::map(input$GENE, .f = doPlotsSingleGene))
  #   # c(rbind(,b))
  # })
  
  output$CLUSTER.TREE <- renderPlot({
    req(DATA$orig.RET)
    cowplot::plot_grid(ggtree(Tool(DATA$orig.RET@seurat, slot = "BuildClusterTree")) + geom_tree() + theme_tree() + geom_tiplab())
  })
  
  breakdownTable <- function(var.1, var.2, transpose = F) {
    if (transpose) {
      MAT <- as.data.frame.matrix(t(percent.table(table(var.1, var.2))))
    } else {
      MAT <- as.data.frame.matrix(percent.table(table(var.1, var.2)))
    }
    MAT <- format(round(MAT, 2), nsmall = 1, width = 5)
    return(MAT)
  }
    
  output$SPLIT.SUMMARY.1 <- renderUI({
    req(DATA$RET, input$SPLIT.BY)
    split.vals <- DATA$orig.RET@seurat[[input$SPLIT.BY]][[input$SPLIT.BY]]
    dt.1 <- breakdownTable(split.vals, Idents(DATA$RET@seurat), T)
    dt.2 <- breakdownTable(Idents(DATA$RET@seurat), split.vals)
    fluidPage(
      paste0("Makeup of ", input$SPLIT.BY, " groups split between ", input$CLUSTER, " group"),
      fluidRow(
        renderDT(dt.1, 
                 height = "150px", 
                 options = list(paging = F, autoWidth = T, searching = F, 
                                info = F, ordering = F)
                 )),
      paste0("Makeup of ", input$CLUSTER, " groups split between ", input$SPLIT.BY, " group"),
      fluidRow(
        renderDT(dt.2, 
                  height = "150px", 
                  options = list(paging = F, autoWidth = T, searching = F,
                                info = F, ordering = F)))
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

  output$DE.MARKERS <- renderUI({
    req(DATA$orig.RET, DATA$orig.RET@markers)
    markers.name <- NULL
    if ("all.markers.full" %in% names(DATA$orig.RET@markers)) {
      markers.name <- "all.markers.full"
    } else if ("all.markers.quick" %in% names(DATA$orig.RET@markers)) {
      markers.name <- "all.markers.quick"
    }
    req(markers.name)
    
    do.call(tabsetPanel,
            purrr::map(levels(DATA$orig.RET@markers[[markers.name]]$cluster), 
                              function(ident){
                  tabPanel(title=ident,
                       renderDT({
                         dplyr::filter(DATA$orig.RET@markers[[markers.name]], avg_logFC > 0, cluster == ident) %>% 
                           dplyr::arrange(p_val_adj) %>%
                           dplyr::select(-p_val) %>%
                           mutate(avg_logFC = round(avg_logFC, 2), 
                                  p_val_adj = formatC(p_val_adj, format = "e", digits = 2))
                       }, height = "550px", options = list(pageLength = 10, autoWidth = T)))
            })
    )
  })
  observeEvent(input$MARKERS.LIST.PATH$datapath, {
    DATA$markers.list <- read_yaml(input$MARKERS.LIST.PATH$datapath)
  })
  # newNameModal <- function(failed = FALSE) {
  #   modalDialog(
  #     textInput("NEW.IDENTS", "Enter new name for current group assignments.",
  #               placeholder = paste0('new_', input$CLUSTER)),
  #     if (failed)
  #       div(tags$b("Group name already exists!", style = "color: red;")),
  #     footer = tagList(
  #       modalButton("Cancel"),
  #       actionButton("OK.SAVE.GROUPS", "OK")
  #     )
  #   )
  # }
  # observeEvent(input$EDIT.MARKERS, {
  #   
  # })
  
  output$MARKER.SETS <- renderUI({
    fluidPage(
      fluidRow(
        column(width = 3, 
             fileInput('MARKERS.LIST.PATH', 'Choose marker list file',
                       multiple=FALSE,
                       accept=c('.yaml'))),
        column(width = 3, 
             downloadButton("DOWNLOAD.MARKERS.LIST", "Save marker list")),
        column(width = 3,
               wellPanel(actionButton("EDIT.MARKERS", "Edit markers")))),
      uiOutput("MARKER.DOTPLOTS"),
      uiOutput("MARKER.FIGS")
      )
  })
  
  plotMultiFeatures <- function(markers.name, n.row = 2, n.col = 2, width = 800, height = 600) {
    markers <- DATA$markers.list[[markers.name]]
    split.markers <- split(markers, 
                           ceiling(seq_along(markers)/(n.row * n.col)))
    message("chunk1", split.markers)
    return(do.call(fluidPage,
                   purrr::map(names(split.markers),
                             .f = function(markers.chunk){
                               plt.width <- width * min(length(split.markers[[markers.chunk]]), n.col) / n.col
                               plt.height <- height * min(ceiling(length(split.markers[[markers.chunk]]) / n.col), n.row) / n.row
                               fluidRow(renderPlot({
                                 VlnPlot(DATA$RET@seurat,
                                         features = c(split.markers[[markers.chunk]]),
                                         ncol = n.col, cols = DATA$COLS,
                                         pt.size = 0) + NoLegend()
                                  # FeaturePlot(DATA$RET@seurat, 
                                  #                     features = c(split.markers[[markers.chunk]]),
                                  #                     reduction = input$DIM.REDUC,
                                  #                     ncol = n.col)
                                        }, width=plt.width,height=plt.height))
                                 

                             })))
  }
  output$MARKER.FIGS <- renderUI({
    req(DATA$RET, DATA$markers.list)
    do.call(tabsetPanel,
            purrr::map(names(DATA$markers.list), 
                       .f = function(markers.name){
                         tabPanel(title=markers.name,
                                  renderUI(plotMultiFeatures(markers.name)))
                       }
            )
    )
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
  output$DOWNLOAD.MARKERS.LIST <- downloadHandler(filename = function() {return(input$MARKERS.LIST.PATH$name)},
                                          content = function(file) {write_yaml(DATA$markers.list, file)})
  # navbarPage(
}
