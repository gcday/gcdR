mydebug <- function(msg="[DEBUG]") {
  DEBUG <- FALSE
  if (DEBUG) {
    print(sprintf("%s - %s - %s", msg1, as.character(Sys.time()), as.character(deparse(sys.calls()[[sys.nframe()-1]]))))
  }
}




server <- function( input, output , session){
  DATA <- reactiveValues() 
  
  updateSeuratAndRunUMAP <- function(SRT) {
    # dims = NULL
    DATA$dims <- NULL
    if (utils::compareVersion(as.character(SRT@version), "3.0") != 1) {
      # pre Seurat 3.0, must update
      message("Updating Seurat object to 3.0...")
      DATA$dims <- c(SRT@calc.params$RunTSNE$dims.use)
      SRT <- UpdateSeuratObject(SRT)
    } else {
      DATA$dims <- SRT@commands$RunTSNE.pca$dims
    }
    if (!"umap" %in% names(SRT@reductions)) {
      if (reticulate::py_module_available(module = 'umap')) {
        SRT <- RunUMAP(SRT, dims = c(DATA$dims))
      } else {
        showModal(modalDialog(title = "Warning!",
        "Missing umap package, only t-SNE will be available. Please install through pip (e.g. pip install umap-learn).",
        easyClose = T))
      }
    } 
    return(SRT)
  }
  
  openInitialRDS <- function() {
    SRT <- readRDS( input$IMAGE$datapath)
    SRT <- updateSeuratAndRunUMAP(SRT)
    COUNTS <- apply( SRT@meta.data, 2, function(x) length(table( x )))
    DATA$CHOICES <- colnames( SRT@meta.data )[COUNTS<=100]
    selected.ident <- NULL
    for (choice in DATA$CHOICES) {
      if(isTRUE(all.equal(as.character(Idents(SRT)),
                          as.character(SRT@meta.data[[choice]])))) {
        selected.ident <- choice
        break
      }
    }
    DATA$selected.ident <- selected.ident
    
    numeric.meta.data <- colnames(SRT@meta.data)[sapply(SRT@meta.data, function(x) is.numeric(x))]
    
    DATA$GENES.OPT <- c(numeric.meta.data, sort(rownames(SRT)))
    if ("umap" %in% names(SRT@reductions)) {
      updateRadioButtons(session, "DIM.REDUC", "Reduction",  
                         c(UMAP = "umap", `t-SNE` = "tsne"), selected = "umap")
    }
    
    updateSelectizeInput(session, 'GENE','Gene', 
                         choices = DATA$GENES.OPT,
                         selected = NULL, server = T)
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
    
    DATA$COLS <- hue_pal()(length(levels(Idents(SRT))))
    # DATA$COLS <- rainbow_hcl(length(levels(Idents(SRT)))) 
    # DATA$COLS <-  rainbow(length(levels(Idents(SRT))))
    DATA$orig.seurat <- SRT
    DATA$seurat <- SRT
    DATA$ACTIVE.FILTER <- F
    message(DATA$selected.ident)
  }
  
  observeEvent(input$IMAGE$datapath, {openInitialRDS()})
  observeEvent(input$do.cell.cycle, {
    cc.genes.species <- cc.genes
    if ("Sdc1" %in% rownames(DATA$orig.seurat)) cc.genes.species <- mouse.cc.genes
    DATA$orig.seurat <- gcdSeuratCellCycleScoring(DATA$orig.seurat, 
                                                  cc.genes.species$s.genes, 
                                                  cc.genes.species$g2m.genes, set.ident = F)
    
    COUNTS <- apply(DATA$orig.seurat@meta.data, 2, function(x) length(table(x)))
    DATA$CHOICES <- colnames(DATA$orig.seurat@meta.data)[COUNTS<=100]
    
    numeric.meta.data <- colnames(DATA$orig.seurat@meta.data)[sapply(DATA$orig.seurat@meta.data, function(x) is.numeric(x))]
    DATA$GENES.OPT <- c(numeric.meta.data, sort(rownames(DATA$orig.seurat)))
    DATA$seurat <- DATA$orig.seurat
  })
  
  newNameModal <- function(failed = FALSE) {
    modalDialog(
      textInput("NEW.IDENTS", "Enter new name for current group assignments.",
                placeholder = paste0('new_', input$CLUSTER)
      ),
      if (failed)
        div(tags$b("Group name already exists!", style = "color: red;")),

      footer = tagList(
        modalButton("Cancel"),
        actionButton("ok", "OK")
      )
    )
  }
  observeEvent(input$RENAME.IDENT, {
    showModal(newNameModal())
  })
  
  observeEvent(input$CHANGE.GRP.NAME, {
    showModal(newGroupNameModal())
  })
  newGroupNameModal <- function(failed = FALSE) {
    modalDialog(
      selectInput("GRPS.TO.RENAME", label = "Group(s) to rename",
                  choices = levels(DATA$orig.seurat), selectize = T, multiple = T
      ),
      textInput("NEW.GRPS.NAME", label = "New label for selected group(s)",
                placeholder = paste(input$GRPS.TO.RENAME, sep = ",")),
      if (failed)
        div(tags$b("Error reassigning ident name.", style = "color: red;")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("ok.grp.rename", "OK")
      )
    )
  }
  observeEvent(input$ok.grp.rename, {
    req(input$NEW.GRPS.NAME)
    if (length(input$GRPS.TO.RENAME) >= 1) {
      levels(DATA$orig.seurat)[levels(DATA$orig.seurat) %in% c(input$GRPS.TO.RENAME)] <- input$NEW.GRPS.NAME
      DATA$seurat <- DATA$orig.seurat
      removeModal()
    } else {
      showModal(newGroupNameModal(failed = TRUE))
    }
  })
  
  
  observeEvent(input$ok, {
    if (input$NEW.IDENTS != "" & !input$NEW.IDENTS %in% colnames(DATA$orig.seurat@meta.data)) {
      DATA$orig.seurat[[input$NEW.IDENTS]] <- Idents(DATA$orig.seurat)
      DATA$CHOICES <- c(input$NEW.IDENTS, DATA$CHOICES)
      DATA$selected.ident <- input$NEW.IDENTS
      updateSelectInput(session, 'CLUSTER', 'Cluster by', 
                        choices = DATA$CHOICES, selected = DATA$selected.ident, server = T)
      removeModal()
    } else {
      showModal(newNameModal(failed = TRUE))
    }rst
  })
  
  observeEvent(input$DO.MARKERS, {
    req(DATA$orig.seurat)
    DATA$all.markers <- FindAllMarkers(DATA$orig.seurat, only.pos = T)
    output$DE.MARKERS <- renderDataTable({DATA$all.markers})
  })

  output$DOWNLOAD <- downloadHandler(
    filename = function() {
      return(input$IMAGE$name)
    },
    content = function(file) {
      saveRDS(DATA$orig.seurat, file)
  })
  output$SPLIT.VAR <- renderUI({
    if (is.null( input$IMAGE ) | is.null(input$CLUSTER) ) return(NULL)
    message("updating SPLIT.BY")
    choices <- DATA$CHOICES[DATA$CHOICES != input$CLUSTER]

    selected.split <- ifelse(length(choices)>= 1, choices[1], NULL)
    selectInput('SPLIT.BY','Optional split value', choices = c(None="", choices), selected = selected.split, selectize = T)
  })

  output$NAMES <- renderUI({
    req(input$FILTER)
    GROUP.NAMES <- levels(DATA$orig.seurat@meta.data[,input$FILTER])
    selectInput("SELECTION", "Filter by Sample/Cluster/Type", GROUP.NAMES, multiple = T, selectize = T)
  })
  output$MAXX <- renderUI({
    req(input$EXPR.LIM)
    if (input$EXPR.LIM) {
      sliderInput("MAXX", "Log-Scale Maximum",
                  min = 0, max = 10,
                  value = 5, step = 0.05)
    }
  })

  observeEvent(input$CLUSTER, ignoreInit = T, {
    if(!isTRUE(all.equal(as.character(Idents(DATA$orig.seurat)),
                        as.character(DATA$orig.seurat@meta.data[[input$CLUSTER]])))) {
      if(class(DATA$orig.seurat[[input$CLUSTER]][[input$CLUSTER]]) == "character") {
        # Getting the sort order correct (e.g. 0, 1, 2,...) for numeric idents, hopefully
        Idents(DATA$orig.seurat) <- as.factor(type.convert(DATA$orig.seurat@meta.data[[input$CLUSTER]]))
      } else {
        Idents(DATA$orig.seurat) <- type.convert(DATA$orig.seurat@meta.data[[input$CLUSTER]], as.is = T)
      }
      DATA$COLS <-  hue_pal()(length(levels(Idents(DATA$orig.seurat))))
      # DATA$COLS <- rainbow_hcl(length(levels(Idents(SRT))))
      # DATA$COLS <-  rainbow(length(levels(Idents(DATA$orig.seurat))))
      
      message("Updating Seurat idents...")
      DATA$seurat <- DATA$orig.seurat
    } else {
      message("Unnecessary update attempt.")
    }
  })
  
  observeEvent(input$SELECTION, {
    if (!is.null(input$SELECTION)) {
      DATA$ACTIVE.FILTER <- T
      print(length(input$SELECTION))
      pass.cells <- which(DATA$orig.seurat@meta.data[,input$FILTER] %in% input$SELECTION)
      CELLS.USE <- colnames(DATA$orig.seurat)[pass.cells]
      DATA$seurat <- subset(DATA$orig.seurat, cells = CELLS.USE)
      message("Updating Seurat filtering...")
    } else {
      DATA$seurat <- DATA$orig.seurat
      DATA$ACTIVE.FILTER <- NULL
      message("Unnecessary filter update attempt.")
    }
  }, ignoreNULL = F)
  
  observeEvent(input$FILTER, {
    if (is.null(input$FILTER)) {
      if (!is.null(DATA$ACTIVE.FILTER)) {
        DATA$seurat <- DATA$orig.seurat
      }
      DATA$ACTIVE.FILTER <- NULL
    } else if (input$FILTER == "") {
      if (!is.null(DATA$ACTIVE.FILTER)) {
        DATA$seurat <- DATA$orig.seurat
      }
      DATA$ACTIVE.FILTER <- NULL
    } else {
      print("Unnecessary filter update attempt.")
    }
  }, ignoreNULL = F)
  
  output$DIM.REDUC <- renderPlot({
    req(DATA$seurat)
    dim.plt <- DimPlot(DATA$seurat,
                       label = input$LABELS,
                       label.size = 6,
                       repel = T,
                       cols = DATA$COLS, 
                       reduction = input$DIM.REDUC)
    cowplot::plot_grid(dim.plt)
  })
  
  output$MULTIFEATURE.PLOT <- renderUI({
    req(input$GENE, DATA$seurat)
    do.call(what = shiny::fluidPage,
            args = purrr::map(input$GENE, .f = function(gene){
              fluidRow(title=gene,
                       renderPlot({
                         FeaturePlot(DATA$seurat, features = c(gene),
                                     reduction = input$DIM.REDUC)
                         }, width=800,height=600))
            })
    )
  })
  
  output$MULTIGENE.PLOT <- renderUI({
    req(input$GENE, DATA$seurat)
    do.call(what = shiny::fluidPage,
            args = purrr::map(input$GENE, .f = function(gene){
              shiny::fluidRow(title=gene,
                              shiny::renderPlot({
                                VlnPlot(DATA$seurat, features = c(gene),
                                        cols = DATA$COLS, pt.size = 0) + NoLegend()
                                }, width=600,height=400))
              })
    )
  })
  output$CLUSTER.TREE <- renderPlot({
    req(DATA$orig.seurat)
    print(DATA$dims)
    DATA$orig.seurat <- BuildClusterTree(DATA$orig.seurat, dims = DATA$dims)
    cowplot::plot_grid(ggtree(Tool(DATA$orig.seurat, slot = "BuildClusterTree")) + geom_tree() + theme_tree() + geom_tiplab())
  })

  
  output$SPLIT.SUMMARY.1 <- renderUI({
    req(DATA$seurat, input$SPLIT.BY)
    fluidPage(
      h4(paste0("Fraction of cells in ", input$SPLIT.BY, " group belonging to ", input$CLUSTER, " group")),
      fluidRow({
        renderTable({
            MAT <- as.data.frame.matrix(t(percent.table(table(DATA$orig.seurat[[input$SPLIT.BY]][[input$SPLIT.BY]], 
                                                              Idents(DATA$seurat)))))
            MAT <- format( round( MAT, 2), nsmall=1, width = 5)
            MAT <- cbind( rownames( MAT ), MAT )
            colnames( MAT )[1] <- '.'
            return(MAT)
          }, digits = -2)}),
      h4(paste0("Fraction of cells in ", input$CLUSTER, " group belonging to ", input$SPLIT.BY, " group")),
      fluidRow({renderTable({
            MAT <- as.data.frame.matrix((percent.table(table(Idents(DATA$seurat), 
                                                             DATA$orig.seurat[[input$SPLIT.BY]][[input$SPLIT.BY]]))))
            MAT <- format(round(MAT, 2), nsmall = 1, width = 5)
            MAT <- cbind(rownames(MAT), MAT)
            colnames( MAT )[1] <- '.'
            return(MAT)
        }, digits = -2)})
  )
  })
  geneSummaryTable <- function(gene) {
    XX <- FetchData(DATA$seurat, vars = c(gene))[[gene]]
    GROUP.NAMES <- levels(Idents(DATA$seurat))
    XY <- DATA$seurat@meta.data[,input$CLUSTER]
    INDEX.X <- 1:nrow(DATA$seurat@meta.data)
    COUNT2 <- table( XY[INDEX.X] )
    COUNT2 <- COUNT2[match( GROUP.NAMES, names( COUNT2 ))]
    COUNT1 <- tapply( XX[INDEX.X] , XY[INDEX.X] , function(x) sum( x > 0 ))
    COUNT1 <- COUNT1[match( GROUP.NAMES, names( COUNT1 ))]
    COUNT3 <- format( round( 100*COUNT1 / COUNT2, 1), nsmall=1)
    COUNT4 <- format( round( 100*COUNT2 / sum(COUNT2), 1), nsmall=1)
    COUNT6 <- format( round( tapply( XX[INDEX.X] , XY[INDEX.X] , mean ) ,2), nsmall=2 )
    COUNT6 <- COUNT6[match( GROUP.NAMES, names( COUNT6 ))]
    XXX <- rbind( COUNT2, COUNT4, COUNT1, COUNT3, COUNT6 )
    rownames( XXX ) <- c('Cell Count','Pct of Total Cells','Expressing Cell Count','Pct of Expressing Cell',
                         'Mean Expression All Cells')
    XXX <- cbind( rownames( XXX ), XXX )
    colnames( XXX )[1] <- '.'
    return( XXX )
  }

  output$GENE.SUMMARY <- renderUI({
    req(input$GENE, DATA$seurat)
    do.call(what = shiny::fluidPage,
            args= purrr::map(input$GENE, .f = function(gene){
              shiny::fluidRow(h2(gene),
                              renderTable({geneSummaryTable(gene)}))
            })
    )
  })
}