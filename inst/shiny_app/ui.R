# command + A  to highlight all
# command + enter to run code
# will bring up external window with Shiny App

library( shiny )

library( Seurat )
library(scales)
library(knitr)
library(kableExtra)
library(gcdR)
library(future)
library(ggtree)
plan(multicore, workers = parallel::detectCores() - 1)



ui <- shinyUI(fluidPage(
  
  #headerPanel('scRNA-Seq Gene and Cluster Viewer'),
  titlePanel('scRNA-Seq Gene and Cluster Viewer'),
  # sidebarPanel(
  fluidRow(  
    column(3,
           fileInput('IMAGE','Choose RDS File',
                       multiple=FALSE,
                       accept=c('.RDS')),
           conditionalPanel("input$IMAGE$name", {
             downloadButton("DOWNLOAD", "Save RDS file")
           }),
           # checkboxInput("EXPR.LIM", 
           #               "Force log-scale max for gene plots", FALSE),
           # uiOutput("MAXX"),
           checkboxInput("LABELS", "Include Labels", FALSE),
           radioButtons("DIM.REDUC", "Reduction",  
                          choices = c(`t-SNE` = "tsne"), selected = "tsne")),
    column(3,
           selectizeInput('GENE','Gene', choices = NULL,
                          selected=NULL, multiple=T),
           selectizeInput('CLUSTER','Group cells by', choices = NULL, selected = NULL),
           selectizeInput('SPLIT.BY', 'Alt. group by (optional)', 
                          choices = NULL, selected = NULL),
           selectizeInput('FILTER','Filter by', choices = c(None=""), selected = NULL),
           uiOutput("NAMES")),
    column(3, 
           actionButton("do.cell.cycle", "CC scoring"),
           actionButton("DO.MARKERS", "Find DE markers"),
           actionButton("CHANGE.GRP.NAME", "Rename group label(s)"),
           actionButton("RENAME.IDENT", "Save group labels")
           )),
  mainPanel(
    tabsetPanel(tabPanel("Dim.plot", plotOutput("DIM.REDUC", width='800px', height='600px')),
                tabPanel("Feat.plot", uiOutput("MULTIFEATURE.PLOT")),
                tabPanel("MultiGene", uiOutput("MULTIGENE.PLOT")),
                tabPanel("Group breakdown", uiOutput('SPLIT.SUMMARY.1')),
                tabPanel("Gene stats", uiOutput("GENE.SUMMARY")),
                tabPanel("Cluster tree", plotOutput("CLUSTER.TREE")),
                tabPanel("DE markers", dataTableOutput("DE.MARKERS")))
  )
))



options(shiny.maxRequestSize = 1000*1024^2)

# shinyApp(ui = ui, server = server)



