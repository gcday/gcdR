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
    print("doing it")
    DATA$orig.RET <- gcdRenameIdent(DATA$orig.RET, input$GRPS.TO.RENAME, input$NEW.GRPS.NAME)
    print(levels(DATA$orig.RET@seurat))
    DATA$RET <- DATA$orig.RET
    removeModal()
  } else {
    showModal(newGroupNameModal(failed = TRUE))
  }
})
