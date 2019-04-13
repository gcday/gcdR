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
    print("renaming")
    DATA$orig.RET <- gcdRenameIdent(DATA$orig.RET, input$GRPS.TO.RENAME, input$NEW.GRPS.NAME)
    print(levels(DATA$orig.RET@seurat))
    DATA$RET <- DATA$orig.RET
    removeModal()
  } else {
    showModal(newGroupNameModal(failed = TRUE))
  }
})

output$rename.panel <- renderUI({
  req(DATA$orig.RET)
  return(div(
             do.call(div, c(style = 'overflow-y: scroll; height: 500px;',
    purrr::map(1:length(levels(DATA$orig.RET@seurat)),
                            function (i) {
                              ident <- levels(DATA$orig.RET@seurat)[i]
                              return(textInput(paste0("dynamicInputRename_", i), width = "60%", label = ident, value = ident))
                            }
          ))),
  actionButton("ok.multi.rename", "Save")
  ))
  # actionButton(ok.multi.rename
})


observeEvent(input$ok.multi.rename, {
  new.namelist <- list()
  for (i in 1:length(levels(DATA$orig.RET@seurat))) {
    new.namelist[i] <- input[[paste0("dynamicInputRename_", i)]]
  }
  if (length(new.namelist) == length(unique(new.namelist))) {
    print("Okay new group names")
    print(c(unlist(new.namelist)))
    DATA$orig.RET <- renameAllIdents(DATA$orig.RET, c(unlist(new.namelist)))
    DATA$RET <- DATA$orig.RET
  }
})
