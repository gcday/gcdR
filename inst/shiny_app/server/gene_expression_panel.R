doPlotsSingleGene <- function(gene, vln.width = "auto", vln.height = 400,
                              feat.width = "auto", feat.height = 500,
                              table.height = "200px") {
  if (input$GENE.VIOLIN.SHOW) {
    pt.size <- 1
  } else {
    pt.size <- 0
  }
  tabPanel(title = gene,
           fluidRow(
             box(
               width = 5, status = "primary",
               renderPlot({
                 VlnPlot(DATA$RET@seurat, features = c(gene),
                         cols = Palettes(DATA$orig.RET, as.integer(input$COLOR.PALETTE)), 
                         pt.size = pt.size) + NoLegend()
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