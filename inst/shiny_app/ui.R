# command + A  to highlight all
# command + enter to run code
# will bring up external window with Shiny App

library( shiny )
library(shinydashboard)

library( Seurat )
library(scales)
library(knitr)
library(kableExtra)
library(future)
library(gcdR)
library(ggtree)
library(yaml)
library(dplyr)
library(DT)
library(colorspace)
plan(multicore, workers = parallel::detectCores() - 1)




sidebar <- dashboardSidebar(
  fileInput('IMAGE','Choose RDS File',
            multiple=FALSE,
            accept=c('.RDS')),
  conditionalPanel("input.IMAGE.name", {
    downloadButton("DOWNLOAD.DATA", "Save RDS file")
  }),
  checkboxInput("LABELS", "Include Labels", FALSE),
  uiOutput("DIM.REDUC.CHOICE"),
  # radioButtons("DIM.REDUC", "Reduction",  
  #              choices = c(`t-SNE` = "tsne"), selected = "tsne"),
  selectizeInput("GENE", "Gene", choices = NULL,
  selected=NULL, multiple=T),
  selectizeInput('CLUSTER', 'Group cells by', choices = NULL, selected = NULL),
  selectizeInput('SPLIT.BY', 'Alt. group by (optional)', 
                 choices = NULL, selected = NULL),
  selectizeInput('FILTER','Filter by', choices = c(None=""), selected = NULL),
  uiOutput("NAMES"),
  actionButton("DO.MARKERS", "Find DE markers"),
  actionButton("CHANGE.GRP.NAME", "Rename group label(s)"),
  actionButton("RENAME.IDENT", "Save group labels")
  # sidebarMenu(
  #   menuItem("Gene Plots", tabName = "GENE.PLOTS")
  # )
)
body <- dashboardBody(
  fluidRow(
    box(title = "Dim Plot",
        width = 10,
        height = "650px",
        collapsible = T,
        plotOutput("DIM.REDUC", 
                   width='auto', height='600px')
    )
  ),
  fluidRow(
    box(title = "Cluster Tree",
        width = 5,
        height = "400px",
        collapsible = T,
        plotOutput("CLUSTER.TREE", width='auto', height='300px')
    )
  ),
  fluidRow(
    box(title = "Gene analysis",
        width = 12,
        height = "1000px",
        collapsible = T,
        uiOutput("GENE.PLOTS")
    )
  ),
  
    # tabBox(
    #   title = "Overview", width = 10,
    #   # The id lets us use input$tabset1 on the server to find the current tab
    #   id = "tabset1", height = "650px",
    #   selected = "Dim.plot",
    #   tabPanel("Dim.plot", "First tab content", 
    #            plotOutput("DIM.REDUC", 
    #                       width='auto', height='600px')),
    #   tabPanel("Cluster.tree", "Cluster tree",
    #            plotOutput("CLUSTER.TREE", width='auto', height='450px'))
    #   # tabPanel("Tab2", "Tab content 2")
    # )
  # ),
  br(),
  fluidRow(
    tabBox(
      title = "Gene figures",
      id = "GENE.PLOTS",
      height = "600px",
      width = 12, 
      selected = "Genes",
      tabPanel("Genes", div(style = 'overflow-x: scroll; overflow-y: scroll',
                            uiOutput("MULTIGENE.PLOT"))),
      # tabPanel("Genes", div(style = 'overflow-x: scroll; overflow-y: scroll', 
      #                       uiOutput("MULTIGENE.PLOT"))),
      tabPanel("Marker sets", div(style = 'overflow-x: scroll; overflow-y: scroll', 
                                  uiOutput("MARKER.SETS")))
      # tabPanel("Tab1", "Tab content 1"),
      # tabPanel("Tab2", "Tab content 2"),
    )),
  fluidRow(
    tabBox(
      title = "Group/gene stats", width = 12,
      height = "450px",
      selected = "Group breakdown",
      tabPanel("Group breakdown",
               div(style = 'overflow-x: scroll; overflow-y: scroll', uiOutput('SPLIT.SUMMARY.1'))),
      tabPanel("Gene stats",
               div(style = 'overflow-x: scroll; overflow-y: scroll', uiOutput("GENE.SUMMARY"))),
      tabPanel("DE markers",
               div(style = 'overflow-x: scroll; overflow-y: scroll', uiOutput("DE.MARKERS"))))
  )
  # fluidRow(
  #   uiOutput()
  # )
  # fluidRow(
  #   tabItems(
  #     tabItem(tabName = "GENE.PLOTS",
  #             fluidRow(
  #               tabBox(width = 12,height="500",
  #                      tabPanel("Genes",
  #                               uiOutput("nlp_sentences_tree")))))),
  #   tabBox(
  #     title = "Gene figures",
  #     id = "GENE.PLOTS",
  #     height = "600px",
  #     width = 12, 
  #     selected = "Genes",
  #     tabPanel("Genes", div(style = 'overflow-x: scroll; overflow-y: scroll',
  #                           uiOutput("MULTIGENE.PLOT"))),
  #     # tabPanel("Genes", div(style = 'overflow-x: scroll; overflow-y: scroll', 
  #     #                       uiOutput("MULTIGENE.PLOT"))),
  #     tabPanel("Marker sets", div(style = 'overflow-x: scroll; overflow-y: scroll', 
  #                                 uiOutput("MARKER.SETS")))
  #     # tabPanel("Tab1", "Tab content 1"),
  #     # tabPanel("Tab2", "Tab content 2"),
  #   )),
  # fluidRow(
  #   tabBox(
  #     title = "Group/gene stats", width = 12,
  #     height = "450px",
  #     selected = "Group breakdown",
  #     tabPanel("Group breakdown", 
  #              div(style = 'overflow-x: scroll; overflow-y: scroll', uiOutput('SPLIT.SUMMARY.1'))),
  #     tabPanel("Gene stats", 
  #              div(style = 'overflow-x: scroll; overflow-y: scroll', uiOutput("GENE.SUMMARY"))),
  #     tabPanel("DE markers", 
  #              div(style = 'overflow-x: scroll; overflow-y: scroll', uiOutput("DE.MARKERS"))))
  # )
    
)

# body <- dashboardBody(
#   fluidRow(
#     tabBox(
#       title = "Overview",
#       # The id lets us use input$tabset1 on the server to find the current tab
#       id = "tabset1", height = "650px",
#       selected = "Dim.plot",
#       tabPanel("Dim.plot", "First tab content", 
#                plotOutput("DIM.REDUC", 
#                           width='800px', height='600px')),
#       tabPanel("Cluster.tree", "Cluster tree",
#                plotOutput("CLUSTER.TREE", width='600px', height='450px'))
#       # tabPanel("Tab2", "Tab content 2")
#     ),
#     tabBox(
#       title = "Group/gene stats",
#       height = "250px",
#       selected = "Group breakdown",
#       tabPanel("Group breakdown", uiOutput('SPLIT.SUMMARY.1')),
#       tabPanel("Gene stats", uiOutput("GENE.SUMMARY")),
#       tabPanel("DE markers", uiOutput("DE.MARKERS"))
#     )
#   ),
#   fluidRow(
#     tabBox(
#       title = "Gene figures",
#       height = "600px",
#       width = 12, 
#       selected = "Genes",
#       tabPanel("Genes", uiOutput("MULTIGENE.PLOT")),
#       tabPanel("Marker sets", uiOutput("MARKER.SETS"))
#       # tabPanel("Tab1", "Tab content 1"),
#       # tabPanel("Tab2", "Tab content 2"),
#     )
#   )
# )
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



