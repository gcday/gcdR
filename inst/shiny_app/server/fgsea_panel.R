output$FGSEA.PANEL <- renderUI({
  fluidPage(
    fluidRow(box(width = "100%",
                 height = "auto",
                 uiOutput("FGSEA"))
    ))
})


makeFGSEAPanel <- function(pathway.list.name, ident, only.pos) {
  req(GseaRes(DATA$orig.RET)$output[[ident]][[pathway.list.name]])
  tabPanel(title=pathway.list.name,
           renderDT({
             results.dt <- GseaRes(DATA$orig.RET)$output[[ident]][[pathway.list.name]]
             if (only.pos) {
               results.dt <- filter(results.dt, ES > 0)
             }
             results.dt <- mutate(results.dt, 
                                  leadingEdge = lapply(leadingEdge, paste, collapse = ","))
   
             return(datatable(
               #dplyr::select(results.dt, -pval) %>%
                                dplyr::arrange(results.dt, padj),
                              options = list(autoWidth = T, scrollX = T,
                                             columnDefs = list(
                                               list(targets = 1,
                                                render = DT::JS(
                                                  "function(data, type, row, meta) {",
                                                  "return type === 'display' && data.length > 25 ?",
                                                  "'<span title=\"' + data + '\">' + data.substr(0, 25) + '...</span>' : data;",
                                                  "}")),
                                               list(targets = 7,
                                                    render = DT::JS(
                                                      "function(data, type, row, meta) {",
                                                      "return type === 'display' && data.length > 50 ?",
                                                      "'<span title=\"' + data + '\">' + data.substr(0, 50) + '...</span>' : data;",
                                                      "}"))
                                               
                                               
                              ))) %>%
                      formatSignif(columns = c("padj", "ES", "NES"),
                                   digits = 3))
             
           }))
}

output$FGSEA <- renderUI({
  req(DATA$orig.RET, GseaRes(DATA$orig.RET))
  
  markers.list <- GseaRes(DATA$orig.RET)$output
  # req(length(intersect(names(markers.list), levels(DATA$orig.RET@seurat))) >= 1)
  
  do.call(tabBox,
          c(width = 12,
            height = "auto",
            purrr::map(names(markers.list),
                       function(ident) {
                         tabPanel(title=ident,
                                  do.call(tabBox,
                                          c(width = 12,
                                            height = "auto",
                                            purrr::map(names(GseaRes(DATA$orig.RET)$output[[ident]]),
                                                       makeFGSEAPanel, ident, input$FGSEA.ONLY.POS)
                                          )))}
                       
            )))
  
})

output$ENRICHR.PANEL <- renderUI({
  fluidPage(
    fluidRow(box(width = "100%",
                 height = "auto",
                 uiOutput("ENRICHR"))
    ))
})

output$ENRICHR <- renderUI({
  req(DATA$orig.RET, DATA$orig.RET@enrichr[[ActiveIdent(DATA$orig.RET)]])
  enrichr.list <- DATA$orig.RET@enrichr[[ActiveIdent(DATA$orig.RET)]]
  # markers.list <- GseaRes(DATA$orig.RET)$output
  # req(length(intersect(names(markers.list), levels(DATA$orig.RET@seurat))) >= 1)
  
  do.call(tabBox,
          c(width = 12,
            height = "auto",
            purrr::map(names(enrichr.list),
                       function(ident) {
                         tabPanel(title=ident,
                                  do.call(tabBox,
                                          c(width = 12,
                                            height = "auto",
                                            purrr::map(names(enrichr.list[[ident]]), 
                                                       function(db.name) {
                                                        tabPanel(title = db.name, 
                                                                 renderDT({
                                                                   datatable(enrichr.list[[ident]][[db.name]],
                                                                             options = list(autoWidth = T, scrollX = T))
                                                                 })
                                                                 )
                                                       }
                                                      )
                                          )
                                  )
                         )
                                            
                          }))
  )

})

observeEvent(input$DO.FGSEA, {
  gene.sets <- list("H: hallmark gene sets" = list("H"),
                    "C1: positional gene sets" = list("C1"),
                    "C2: curated gene sets" = list(`CGP: chemical and genetic perturbations` = "CGP", 
                                                   `CP: Canonical pathways` = "CP", 
                                                   `CP:BIOCARTA: BioCarta gene sets` = "CP:BIOCARTA",
                                                   `CP:KEGG: KEGG gene sets` = "CP:KEGG",
                                                   `CP:REACTOME: Reactome gene sets` = "CP:REACTOME"),
                    "C3: motif gene sets" = list(`MIR: microRNA targets` = "MIR",
                                                 `TFT: transcription factor targets` = "TFT"),
                    "C4: computational gene sets" = list(`CGN: cancer gene neighborhoods` = "CGN",
                                                         `CM: cancer modules` = "CM"),
                    "C5: GO gene sets" = list(`BP: GO biological process` = "BP",
                                              `CC: GO cellular component` = "CC",
                                              `MF: GO molecular function` = "MF"),
                    "C6: oncogenic signatures" = list("C6"),
                    "C7: immunologic signatures" = list("C7"))
  
  
  showModal({
    modalDialog(
      "Select groups for calculating DE genes (if none selected, all will be used--which may be very slow! (>30 min)).",
      selectInput("MSIGDBR.SPECIES", label = "Species", 
                  choices = c("Homo sapiens", "Mus musculus"), selected = DATA$species),
      selectInput("MSIGDBR.SETS", label = "Gene sets", 
                  choices = gene.sets, selectize = T, multiple = T),
      textInput("FGSEA.NAME", label = "Label for selected pathways", value = "Selected pathways"),
      
      
      #                                     ""),
      #                      selected = numeric.meta.data[1], server = T)
      # 
      # tipify(selectInput("IDENTS.TO.USE", label = "Groups for DE gene calculation",
      #                    choices = levels(DATA$orig.RET@seurat), selectize = T, multiple = T),
      #        "Compares cells in each group to cells in ALL other groups, regardless of what you choose here!"),
      # popify(sliderInput("DE.LOGFC", "Log-FC threshold",
      #                    min = 0, max = 2,
      #                    value = 0.5, step = 0.05), 
      #        "Limits testing to genes whose avg. expression in in-group cells differs from avg. expression in outgroup cells by at least this value, speeding up the testing (but potentially missing important signals!). 0.5 works well for finding markers between very distinct populations (e.g. B vs. T cells), but lower values are needed when comparing expression in e.g. tumor subsets."),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("DO.FGSEA.OK", "Ok"))
    )
  })
})

observeEvent(input$DO.FGSEA.OK, {
  req(DATA$orig.RET, BestMarkers(DATA$orig.RET))
  withProgress(message = 'FGSEA', value = 0, {
    library(msigdbr)
    pathways.df <- msigdbr(species = input$MSIGDBR.SPECIES)
    pathways.list <- list()
    pathways.list[[input$FGSEA.NAME]] <- filter(pathways.df, gs_cat %in% input$MSIGDBR.SETS | gs_subcat %in% input$MSIGDBR.SETS) %>% split(x = .$gene_symbol, f = .$gs_name)
    DATA$orig.RET <- fgseaWrapper(DATA$orig.RET, pathways.list, shiny = T, nproc = 4)
  })
  DATA$RET <- DATA$orig.RET 
  removeModal()
})
