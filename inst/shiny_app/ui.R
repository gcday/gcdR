library(shiny)
library(shinydashboard)
library(shinyBS)
library(Seurat)
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
library(shinyFiles)
library(shinyWidgets)
library(fs)

options(shiny.maxRequestSize = 120000*1024^2)



plan(multicore, workers = parallel::detectCores() - 1)


sidebar <- dashboardSidebar(
  tags$style(".skin-blue .sidebar a { color: #444; }"),
  tags$style(".pretty {
    white-space: inherit;
             width: 14rem;
             }

             .pretty .state label{
             text-indent: 0;
             padding-left: 2rem;
             }

             .pretty .state label:after,
             .pretty .state label:before{
                top: calc(50% - ((1em + 2px)/2));
             }
             .pretty .state label:after,
             .pretty .state label:before,
             .pretty.p-icon .state .icon {
                top: calc(50% - ((1em + 2px)/2));
             }"),
  # tags$style(".pretty {
  #   white-space: normal;
  #            }
  #            .pretty .state label::after, .pretty .state label::before
  #            {
  #            top:-2px;
  #            }"),
  
  fileInput('IMAGE','Choose RDS File',
            multiple=FALSE,
            accept=c('.RDS')),
  shinyFilesButton('file', 'Load Dataset', 'Please select a dataset', FALSE),
  # shinySaveButton('save.file', "Save dataset", "Please select a name", "")
  # conditionalPanel("input$file",
  shinySaveButton("serverside_save", "Save dataset", "Save file as...", filetype = list(RDS = "rds")),
  # , filetype = list(RDS = "rds")),
  # ),
  
  div(style = "margin-left: 10px; margin-bot: 5px;",
      # shinyFilesButton("LOCAL.FILE", "Local file", "Please select a file", multiple = FALSE),
    downloadButton("DOWNLOAD.DATA", "Save RDS file")
    ),
  sidebarMenu(id = "sidebarTabs",
    menuItem("Overview", tabName = "overview", icon = icon("dashboard")),
    conditionalPanel("input.sidebarTabs == 'overview'",
                      div(style = "margin-left: 10px; margin-top: 5px;",
                       checkboxInput("LABELS", "Include Labels", FALSE),
                       actionButton("RECLUSTER", "Re-cluster", icon = icon("calculator")),
                       tipify(selectizeInput('SPLIT.BY', 'Secondary split value', 
                                      choices = NULL, selected = NULL),
                          "Used for group breakdown table; reccomended choices include \"Phase\" (cell cycle) and sample", placement = "top"),
                       radioGroupButtons("SPLIT.TYPE", "Split type", choices = c(Counts = "counts", Percent = "percent"), selected = "percent"), 
                       radioGroupButtons("SPLIT.PALETTE", "Color palette for split", 
                                         direction = "horizontal", 
                                         choices = c(`1` = 1, `2` = 2, `3` = 3, `4` = 4), 
                                         selected = 2)
                     )),
    menuItem("Gene expression", icon = icon("dna"), tabName = "genePanel"),
    conditionalPanel("input.sidebarTabs == 'genePanel'",
                     div(style = "margin-left: 10px;",
                         tipify(selectizeInput("GENE", "Gene", choices = NULL,
                                        selected=NULL, multiple=T), 
                                "Enter one or more genes", placement = "top"),
                         checkboxInput("GENE.VIOLIN.SHOW",
                                       "Show cells as dots in violin plot", 
                                       value = F)
                     )),
    menuItem("DE genes", icon = icon("balance-scale"), tabName = "DEPanel"),
    conditionalPanel("input.sidebarTabs == 'DEPanel'",
                     div(style = "margin-left: 10px;",
                         actionButton("DO.MARKERS", "Find DE markers", icon = icon("calculator"))
                     )),
    menuItem("GSEA", icon = icon("chart-bar"), tabName = "GSEAPanel"),
    conditionalPanel("input.sidebarTabs == 'GSEAPanel'",
                     div(style = "margin-left: 10px; ",
                         actionButton("DO.FGSEA", "Run FGSEA", icon = icon("calculator")),
                          prettyCheckbox("FGSEA.ONLY.POS",
                                      "Show only enriched pathways (positive ES)",
                                      value = T, status = "primary")
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
      div(style = "margin-left: 10px;",
          uiOutput("DIM.REDUC.CHOICE"),
          tipify(
            selectizeInput('CLUSTER', 'Group cells by', choices = NULL, selected = NULL),
            "Reccomended: group cells by clusters here, then use \"Secondary split value\" under Overview to analyze breakdown by sample, cell cycle phase, etc.", placement = "top"
          ),
          selectizeInput('FILTER','Filter by', choices = c(None=""), selected = NULL),
          uiOutput("NAMES"),
          
          actionButton("CHANGE.GRP.NAME", "Rename group"),
          actionButton("RENAME.IDENT", "Save current groups"),
          radioGroupButtons("COLOR.PALETTE", "Color palette", 
                               direction = "horizontal", 
                               choices = c(`1` = 1, `2` = 2, `3` = 3, `4` = 4), 
                            selected = 1)
          # radioButtons("COLOR.PALETTE", "Color palette",  
          #              choices = c(`1` = 1, `2` = 2, `3` = 3, `4` = 4), selected = 1)
      )
  


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
    ),
    tabItem(tabName = "GSEAPanel",
            fluidRow(
              box(title = "GSEA",
                  width = 12,
                  height = "auto",
                  uiOutput("FGSEA.PANEL")
              )
            )
    )
    )
)

  
  


    



ui <- shinyUI(dashboardPage(
  dashboardHeader(title = 'scRNA-Seq Gene and Cluster Viewer'),
  sidebar,
  body
))