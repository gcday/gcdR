library( shiny )
library(shinydashboard)
library(shinyBS)

library( Seurat )
library(scales)
library(knitr)
library(kableExtra)
library(future)
# library(gcdR)
library(ggtree)
library(yaml)
library(dplyr)
library(DT)
library(colorspace)
plan(multicore, workers = parallel::detectCores() - 1)


sidebar <- dashboardSidebar(
  tags$style(".skin-blue .sidebar a { color: #444; }"),
  
  fileInput('IMAGE','Choose RDS File',
            multiple=FALSE,
            accept=c('.RDS')),
  div(style = "margin-left: 10px; margin-bot: 5px;",
  downloadButton("DOWNLOAD.DATA", "Save RDS file")),
  sidebarMenu(id = "sidebarTabs",
    menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
    conditionalPanel("input.sidebarTabs == 'overview'",
                      div(style = "margin-left: 10px; margin-top: 5px;",
                       checkboxInput("LABELS", "Include Labels", FALSE),
                       actionButton("RECLUSTER", "Re-cluster", icon = icon("calculator")),
                       tipify(
                         selectizeInput('SPLIT.BY', 'Secondary split value', 
                                      choices = NULL, selected = NULL),
                          "Used for group breakdown table")
                     )),
    menuItem("Gene expression", icon = icon("dna"), tabName = "genePanel"),
    conditionalPanel("input.sidebarTabs == 'genePanel'",
                     div(style = "margin-left: 10px;",
                         tipify(selectizeInput("GENE", "Gene", choices = NULL,
                                        selected=NULL, multiple=T), 
                                "Enter one or more genes", placement = "bottom")
                     )),
    menuItem("DE genes", icon = icon("chart-bar"), tabName = "DEPanel"),
    conditionalPanel("input.sidebarTabs == 'DEPanel'",
                     div(style = "margin-left: 10px;",
                         actionButton("DO.MARKERS", "Find DE markers", icon = icon("calculator"))
                     )),
    menuItem("Marker sets", icon = icon("th"), tabName = "markerPanel"),
    conditionalPanel("input.sidebarTabs == 'markerPanel'",
                     div(style = "margin-left: 10px;",
                         fileInput('MARKERS.LIST.PATH', 'Choose marker list file',
                                   multiple=FALSE,
                                   accept=c('.yaml')),
                         radioButtons("MARKERS.TYPE", "Plot type", 
                                      choices = c(Violin = "violin", Feature = "feature"),
                                      selected = "violin")
                     ))
  ),
  uiOutput("DIM.REDUC.CHOICE"),
  
  tipify(
    selectizeInput('CLUSTER', 'Group cells by', choices = NULL, selected = NULL),
    "Reccomended: group cells by clusters here, then use \"Secondary split value\" under Overview to analyze breakdown by sample, cell cycle phase, etc. "
  )  ,
  
  selectizeInput('FILTER','Filter by', choices = c(None=""), selected = NULL),
  uiOutput("NAMES"),
  
  actionButton("CHANGE.GRP.NAME", "Rename group label"),
  actionButton("RENAME.IDENT", "Save current groups"),
  radioButtons("COLOR.PALETTE", "Color palette",  
                                       choices = c(`1` = 1, `2` = 2, `3` = 3, `4` = 4), selected = 1)
  
)
body <- dashboardBody(
  tabItems(
    tabItem(tabName = "overview",
            fluidRow(
              div(class = "col-sm-12 col-md-10 col-lg-8",
                  box(title = "Dim Plot",
                      width='100%',
                      height = "auto",
                      plotOutput("DIM.REDUC",
                                 width='100%', height='600px')
                  )
              )
            ),
            fluidRow(
              div(class = "col-sm-12 col-md-10 col-lg-8",
                  box(title = "Dim Plot Split",
                      width='100%',
                      height = "auto",
                      plotOutput("DIM.REDUC.SPLIT",
                                 width='100%', height='600px')
                  )
              )
            ),
            fluidRow(
              div(class = "col-sm-12 col-md-9 col-lg-7",
              box(title = "Cluster Tree",
                  width='100%',
                  height = "auto",
                  plotOutput("CLUSTER.TREE", width='100%', height='350px')
              ))
            ),
            fluidRow(
              box(title = "Group breakdown",
                  width = 12,
                  height = "auto",
                  div(style = 'overflow-x: scroll;', 
                      uiOutput('SPLIT.SUMMARY.1'))
                  )
            )
    ),
    tabItem(tabName = "genePanel",
            fluidRow(
                box(title = "Gene analysis",
                    width = 12,
                    height = "1200px",
                    uiOutput("GENE.PLOTS")
                )
              )
    ),
    tabItem(tabName = "markerPanel",
            fluidRow(
              box(title = "Markers",
                  width = 12,
                  height = "auto",
                  uiOutput("MARKER.SETS"))
            )
    ),
    tabItem(tabName = "DEPanel",
            fluidRow(
              box(title = "DE genes",
                  width = 12,
                  height = "auto",
                  uiOutput("DE.MARKERS")
              )
          )
    )
    )
)

  
  


    



options(shiny.maxRequestSize = 1000*1024^2)

ui <- dashboardPage(
  dashboardHeader(title = 'scRNA-Seq Gene and Cluster Viewer'),
  sidebar,
  body
)
# ui <- shinyUI(fluidPage(
#   
#   #headerPanel('scRNA-Seq Gene and Cluster Viewer'),
#   titlePanel('scRNA-Seq Gene and Cluster Viewer'),
#   # sidebarPanel(
#   fluidRow(  
#     column(3,
#            fileInput('IMAGE','Choose RDS File',
#                        multiple=FALSE,
#                        accept=c('.RDS')),
#            conditionalPanel("input.IMAGE.name", {
#              downloadButton("DOWNLOAD.DATA", "Save RDS file")
#            }),
#            # checkboxInput("EXPR.LIM", 
#            #               "Force log-scale max for gene plots", FALSE),
#            # uiOutput("MAXX"),
#            checkboxInput("LABELS", "Include Labels", FALSE),
#            radioButtons("DIM.REDUC", "Reduction",  
#                           choices = c(`t-SNE` = "tsne"), selected = "tsne")),
#            
#     column(3,
#            # radioButtons("MULTIGENE.MODE", "Viewing mode",  
#            #              choices = c(`Single gene` = "single", `Multi-gene` = "multi"), selected = "single"),
#            selectizeInput("GENE", "Gene", choices = NULL,
#                           selected=NULL, multiple=T)),
#            # uiOutput("GENE.LIST")),
#     column(3,
#            selectizeInput('CLUSTER', 'Group cells by', choices = NULL, selected = NULL),
#            selectizeInput('SPLIT.BY', 'Alt. group by (optional)', 
#                           choices = NULL, selected = NULL),
#            selectizeInput('FILTER','Filter by', choices = c(None=""), selected = NULL),
#            uiOutput("NAMES")),
#            # radioButtons("COLOR.PALETTE", "Palette",  
#                         # choices = c(`1` = 1, `2` = 2, `3` = 3), selected = 1)),
#     column(3, 
#            actionButton("do.cell.cycle", "CC scoring"),
#            actionButton("DO.MARKERS", "Find DE markers"),
#            actionButton("CHANGE.GRP.NAME", "Rename group label(s)"),
#            actionButton("RENAME.IDENT", "Save group labels")
#            )),
#   # mainPanel(
#   tabsetPanel(tabPanel("Dim.plot", plotOutput("DIM.REDUC", width='800px', height='600px')),
#               # tabPanel("Feat.plot", uiOutput("MULTIFEATURE.PLOT")),
#               tabPanel("MultiGene", uiOutput("MULTIGENE.PLOT")),
#               tabPanel("Group breakdown", uiOutput('SPLIT.SUMMARY.1')),
#               tabPanel("Gene stats", uiOutput("GENE.SUMMARY")),
#               tabPanel("Cluster tree", plotOutput("CLUSTER.TREE", width='600px', height='450px')),
#               tabPanel("DE markers", uiOutput("DE.MARKERS")),
#               tabPanel("Marker sets", uiOutput("MARKER.SETS"))
#                        # renderUI({
#                        
#                         # uiOutput("MARKER.SETS"))
#               )
#                        
#                        
#                        # uiOutput("MARKER.SETS")))
#   # )
# ))




# shinyApp(ui = ui, server = server)



