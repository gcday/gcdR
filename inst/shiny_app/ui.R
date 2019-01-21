# command + A  to highlight all
# command + enter to run code
# will bring up external window with Shiny App

library( shiny )

library( Seurat )
library(scales)
library(knitr)
library(kableExtra)
library(gcdR)
# library(future)
# plan(multiprocess, workers = 15)



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
           # c(UMAP = "umap", TSNE = "tsne"), 
                          # selected = "umap")),
             
      # column(3,
             # ,
             # selectInput("DIM.REDUC", "Dimensionality reduction", 
                         # choices = c("tsne", "umap"), selected = "umap", selectize = T)),
             # radioButtons('LABELS', 'Include Labels',
                          # c(yes='YES',no='NO'),'YES')),
      # column(3,
      #        # sliderInput("MAXX", "Log-Scale Maximum",
      #        #               min = 0, max = 10,
      #        #               value = 5, step = 0.05),
      #        ,
    column(3,
           selectizeInput('GENE','Gene', choices = NULL,
                          selected=NULL, multiple=T),
           selectizeInput('CLUSTER','Group cells by', choices = NULL, selected = NULL),
           # uiOutput("CLUSTERS"),
           
           selectizeInput('SPLIT.BY', 'Alt. group by (optional)', 
                          choices = NULL, selected = NULL),
           selectizeInput('FILTER','Filter by', choices = c(None=""), selected = NULL),
           # conditionalPanel("input$FILTER != null",  selectizeInput("")
           # uiOutput("SPLIT.VAR",),
           # uiOutput("FILTERBY"),
           uiOutput("NAMES")),
    column(3, 
           actionButton("do.cell.cycle", "CC scoring"),
           actionButton("CHANGE.GRP.NAME", "Rename group label(s)"),
           actionButton("RENAME.IDENT", "Save group labels")
           )),
  mainPanel(
    tabsetPanel(tabPanel("Dim.plot", plotOutput("DIM.REDUC", width='800px',height='600px')),
                tabPanel("Feat.plot", uiOutput("MULTIFEATURE.PLOT")),
                         # plotOutput("FEAT.PLOT", width='800px',height='600px')),
                tabPanel("MultiGene", uiOutput("MULTIGENE.PLOT")),
                tabPanel("Group breakdown", uiOutput('SPLIT.SUMMARY.1')),
                tabPanel("Gene stats", uiOutput("GENE.SUMMARY")))
                # tabPanel("Summary", tableOutput('SUMMARY')))
  )
))



options(shiny.maxRequestSize = 1000*1024^2)

# shinyApp(ui = ui, server = server)



