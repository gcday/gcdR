#' Uses expression gating to label T cells.
#'
#'
#' @param object Seurat object
#' @param CD20_cutoff threshold value for CD20+ (default is 0)
#' @param CD3E_cutoff threshold value for CD3+ (default is 0)
#' 
#' @return Seurat object with object$CD4.CD8 filled in 
#'
#' @examples
#' annotateCD4.CD8.TCells(object)
#'
#' @export
#' 
annotateCD3.CD20.cells <- function(object, 
                                   CD20_cutoff = 0, 
                                   CD3E_cutoff = 0) {
  CD3E_expr <- FetchData(object = object, vars = c("CD3E"))
  CD3E.pos <- which(x = CD3E_expr > CD3E_cutoff)
  
  CD20_expr <- FetchData(object = object, vars = c("MS4A1"))
  CD20.pos <- which(x = CD20_expr > CD20_cutoff)
  CD3E.CD20.pos <- intersect(CD3E.pos, CD20.pos)
 
  object$CD3.CD20 <- "CD3E- CD20-"
  object$CD3.CD20[CD3E.pos] <- "CD3E+ CD20-"
  object$CD3.CD20[CD20.pos] <- "CD3E- CD20+"
  object$CD3.CD20[CD3E.CD20.pos] <- "CD3E+ CD20+"
  return(object)
}
#' Uses expression gating to label T cells.
#'
#'
#' @param object Seurat object
#' @param CD4_cutoff threshold value for CD4+
#' @param CD8A_cutoff threshold value for CD8+
#' @param CD3D_cutoff threshold value for CD3+ 
#' 
#' @return Seurat object with object$CD4.CD8 filled in 
#'
#' @examples
#' annotateCD4.CD8.TCells(object)
#'
#' @export
#' 
annotateCD4.CD8.TCells <- function(object, 
                                   CD4_cutoff = 0.5, 
                                   CD8A_cutoff = 0.5, 
                                   CD3D_cutoff = 0.5) {
  CD3D_expr <- FetchData(object = object, vars = c("CD3D"))
  CD3D.pos <- which(x = CD3D_expr > CD3D_cutoff)
  
  CD3D_expr <- FetchData(object = object, vars = c("CD3D"))
  CD3D.pos <- which(x = CD3D_expr > CD3D_cutoff)
  
  CD4_expr <- FetchData(object = object, vars = c("CD4"))
  CD4.pos <- which(x = CD4_expr > CD4_cutoff)
  
  CD8A_expr <- FetchData(object = object, vars = c("CD8A"))
  CD8A.pos <- which(x = CD8A_expr > CD8A_cutoff)
  CD4.pos <- intersect(CD4.pos, CD3D.pos)
  CD8A.pos <- intersect(CD8A.pos, CD3D.pos)
  CD4.CD8.pos <- intersect(CD4.pos, CD8A.pos)
  object$CD4.CD8 <- "CD3-"
  object$CD4.CD8[CD3D.pos] <- "CD3+ DN"
  object$CD4.CD8[CD4.pos] <- "CD4"
  object$CD4.CD8[CD8A.pos] <- "CD8"
  object$CD4.CD8[CD4.CD8.pos] <- "CD4.CD8"
  return(object)
}
#' Applies TCR annotation to Seurat object
#'
#'
#' @param object Seurat object
#' @param barcode_to_clonotype list mapping barcode to clonotypes
#' @param barcode_to_TRAV list mapping barcodes to TRAV genes
#' @param barcode_to_TRBV list mapping barcodes to TRBV genes
#' 
#' @return Seurat object with object$TCR_clonotype, object$TRA_V, 
#' and object$TRB_V filled in 
#'
#' @examples
#' addTCRToSeurat(object, barcode_to_clonotype)
#'
#' @export
#' 
addTCRToSeurat <- function(object, barcode_to_clonotype,
                           barcode_to_TRAV, barcode_to_TRBV) {
  object$TCR_clonotype <- sapply(X = rownames(object[[]]), 
                                     FUN = function(x) {
                                       if (x %in% names(barcode_to_clonotype)) {
                                         return(barcode_to_clonotype[[x]])
                                       } else {
                                         return("none")
                                       }
                                     })
  object$TRA_V <- sapply(X = rownames(object[[]]), 
                             FUN = function(x) {
                               if (x %in% names(barcode_to_TRAV)) {
                                 return(barcode_to_TRAV[[x]])
                               } else {
                                 return("none")
                               }
                             })
  object$TRB_V <- sapply(X = rownames(object[[]]), 
                             FUN = function(x) {
                               if (x %in% names(barcode_to_TRBV)) {
                                 return(barcode_to_TRBV[[x]])
                               } else {
                                 return("none")
                               }
                             })
  return(object)
}

#' Pair of dimensional reduction plots highlighting cells of a given clonotype 
#' @importFrom cowplot plot_grid
#' @importFrom Seurat DimPlot NoLegend
#' @export
plotClonotypeIntVsSCT <- function(object1, prefix1, object2, prefix2, clonotype, clonotype_prefix, clonotype_field = "TCR_clonotype") {
  cells.highlight.1 <- list(
    rownames(object1[[]])[which(object1[[clonotype_field]][[clonotype_field]] == paste0(clonotype_prefix, clonotype))])
  names(cells.highlight.1) <- paste("Clonotype", clonotype)
  cells.highlight.2 <- list(
    rownames(object2[[]])[which(object2[[clonotype_field]][[clonotype_field]] == paste0(clonotype_prefix, clonotype))])
  names(cells.highlight.2) <- paste("Clonotype", clonotype)
  plot.1 <- suppressMessages(DimPlot(object = object1, 
                                     cells.highlight = cells.highlight.1,
                                     sizes.highlight = 0.75) + 
                               labs(title = paste("Clonotype", clonotype, prefix1)) + 
                               NoLegend() + 
                               scale_color_manual(values = c("#bcbcbc", "red")))
  plot.2 <- suppressMessages(DimPlot(object = object2, 
                                     cells.highlight = cells.highlight.2,
                                     sizes.highlight = 0.75) + 
                               labs(title = paste("Clonotype", clonotype, prefix2)) + 
                               NoLegend() + 
                               scale_color_manual(values = c("#bcbcbc", "red")))
  return(plot_grid(plot.1, plot.2))
}

#' Dimensional reduction plot highlighting cells of a given clonotype based on their CD4.CD8 status
#' 
#' @importFrom Seurat DimPlot NoLegend

#' @export

plotClonotypeCD4.CD8 <- function(object, title, clonotype, clonotype_prefix, ...) {
  clonotype.cells <- object[[]][which(object$TCR_clonotype == paste0(clonotype_prefix, 
                                                                     clonotype)),]
  # CD4.CD8.values <- c("CD3-", "CD3+ DN", "CD4", "CD4.CD8", "CD8")
  CD4.CD8.values <- c("CD8", "CD4", "CD4.CD8", "CD3-", "CD3+ DN")
  
  cells.highlight <- list()
  color.pal <- pal_nejm()(length(CD4.CD8.values))
  cols.highlight <- list()
  for (i in 1:length(CD4.CD8.values)) {
    val <- CD4.CD8.values[i]
    matching.cells <- rownames(clonotype.cells)[which(clonotype.cells$CD4.CD8 == val)]
    if (length(matching.cells) > 0) {
      cells.highlight[[paste0(val, " (", length(matching.cells), " cells)" )]] <- matching.cells
      cols.highlight[[paste0(val, " (", length(matching.cells), " cells)" )]] <-  color.pal[i]
    }
  }
  
  return(suppressMessages(DimPlot(object = object, 
                                  cells.highlight = cells.highlight, 
                                  cols.highlight = c(unlist(cols.highlight)),
                                  ...) + 
                            labs(title = paste("Clonotype", clonotype, title)) +
                            scale_color_manual(values = c(Unselected = "#bcbcbc", unlist(cols.highlight)))))
}

#' Dimensional reduction plot highlighting cells of a given clonotype based on their CD4.CD8 status
#' 
#' @importFrom Seurat DimPlot NoLegend

#' @export

plotClonotypeTimePoints <- function(object, title, clonotype, clonotype_prefix, 
                                    timepoint.field = "timepoint", site.field = "site",
                                    site.levels = NULL,
                                    timepoint.levels = NULL,
                                    ...
                                    ) {
  clonotype.cells <- object[[]][which(object$TCR_clonotype == paste0(clonotype_prefix, 
                                                                     clonotype)),]
  timepoint.levels <- timepoint.levels %||% levels(as.factor(object[[timepoint.field]][[timepoint.field]]))
  # site.levels <- site.levels %||% levels(as.factor(object[[site.field]][[site.field]]))
  
  # CD4.CD8.values <- c("CD3-", "CD3+ DN", "CD4", "CD4.CD8", "CD8")
  # CD4.CD8.values <- c("CD8", "CD4", "CD4.CD8", "CD3-", "CD3+ DN")
  
  cells.highlight <- list()
  color.pal <- pal_nejm()(length(timepoint.levels))
  cols.highlight <- list()
  for (i in 1:length(timepoint.levels)) {
    val <- timepoint.levels[i]
    matching.cells <- rownames(clonotype.cells)[which(clonotype.cells[[timepoint.field]] == val)]
    if (length(matching.cells) > 0) {
      cells.highlight[[paste0(val, " (", length(matching.cells), " cells)" )]] <- matching.cells
      cols.highlight[[paste0(val, " (", length(matching.cells), " cells)" )]] <-  color.pal[i]
    }
  }
  
  return(suppressMessages(DimPlot(object = object, 
                                  cells.highlight = cells.highlight, 
                                  cols.highlight = c(unlist(cols.highlight)),
                                  shape.by = site.field,
                                  ...) + 
                            labs(title = paste("Clonotype", clonotype, title)) +
                            scale_color_manual(values = c(Unselected = "#bcbcbc", unlist(cols.highlight)))))
}





#' Adds glycosylation site amino acid position and number to a clones data table.
#' 
#' @importFrom seqinr translate c2s
#' @importFrom stringr str_count str_locate_all 
#' @importFrom future.apply future_sapply 
#' 
#' @export
#' 

countGlycosylationSites <- function(clones.tbl, new.mode = F) {
  if (!new.mode) {
    clones.tbl$NOGAP_SEQUENCE <- future_sapply(clones.tbl$SEQUENCE_IMGT,
                                               FUN = function(igh.seq) {
                                                 return(c2s(translate(gsub(pattern = "\\.",
                                                                           replacement = "", 
                                                                           x = s2c(igh.seq)))))
                                               })
    
    clones.tbl$GLYCOSYLATION_SITES <- future_sapply(clones.tbl$NOGAP_SEQUENCE, 
                                                    FUN = function(nogap.seq) {
                                                      num.glycosylation.sites <- stringr::str_count(nogap.seq, pattern = "N[[:upper:]][ST]")
                                                      if (num.glycosylation.sites != 0) {
                                                        glycosylation.sites <- stringr::str_locate_all(nogap.seq, pattern = "N[[:upper:]][ST]")
                                                        return(paste(unlist(glycosylation.sites[[1]][,1]), sep = ",", collapse = ", "))
                                                      } else {
                                                        return("")
                                                      }})
    clones.tbl$NUM_GLYCOSYLATION_SITES <- future_sapply(clones.tbl$NOGAP_SEQUENCE, 
                                                        FUN = function(nogap.seq) {stringr::str_count((nogap.seq), pattern = "N[[:upper:]][ST]")})
    
  } else {
    out.list <- future_sapply(clones.tbl$SEQUENCE_IMGT, FUN = function(igh.seq) {
      nogap.seq <- c2s(translate(gsub(pattern = "\\.",
                                      replacement = "", 
                                      x = s2c(igh.seq))))
      num.glycosylation.sites <- stringr::str_count(nogap.seq, pattern = "N[[:upper:]][ST]")
      if (num.glycosylation.sites != 0) {
        glycosylation.sites <- stringr::str_locate_all(nogap.seq, pattern = "N[[:upper:]][ST]")
        glycos.str <- paste(unlist(glycosylation.sites[[1]][,1]), sep = ",", collapse = ", ")
      } else {
        glycos.str <- ""
      }
      return(c(nogap.seq, num.glycosylation.sites, glycos.str))
    })
    clones.tbl$NOGAP_SEQUENCE <- unlist(out.list[1:length(clones.tbl$SEQUENCE_IMGT) * 3 - 2])
    clones.tbl$NUM_GLYCOSYLATION_SITES <- unlist(out.list[1:length(clones.tbl$SEQUENCE_IMGT) * 3 - 1])
    clones.tbl$GLYCOSYLATION_SITES <- unlist(out.list[1:length(clones.tbl$SEQUENCE_IMGT) * 3])
  }
    
  
  
  return(clones.tbl)
}

#' Adds glycosylation site amino acid position and number to a clones data table.
#' 
#' @importFrom progressr progressor 
#' 
#' @export
#' 


parseBuildTreesLogfile <- function(file.path, patient = NULL) {
  logfile <- readLines(file.path)
  
  nvec <- length(logfile)
  breaks <- which(! nzchar(logfile))
  nbreaks <- length(breaks)
  if (breaks[nbreaks] < nvec) {
    breaks <- c(breaks, nvec + 1L)
    nbreaks <- nbreaks + 1L
  }
  if (nbreaks > 0L) {
    chunks <- mapply(function(a,b) paste(logfile[a:b], collapse = "\n"),
                     c(1L, 1L + breaks[-nbreaks]),
                     breaks - 1L)
  }
  pass.clones <- list()
  parent_to_cloneid <- list()
  p <- progressr::progressor(steps = length(chunks))
  for (i in 1:length(chunks)) {
    # for (i in 1:2) {
    
    # if ("ID")
    fields <- list(ID = "", 
                   CLONE = "", 
                   PASS = FALSE, 
                   END_MASKED = "", 
                   SEQ_IN = "",
                   SEQ_IMGT = "",
                   SEQ_MASKED = "",
                   IN_FRAME = "",
                   MASKED = "",
                   FRAMESHIFTS = "",
                   FAIL = "",
                   COLLAPSETO = "",
                   COLLAPSEFROM = "",
                   DUPLICATE = FALSE)
    # record.lines <- strsplit(x = chunks[i], split = "\n", fixed = T)
    for (line in strsplit(x = chunks[i], split = "\n", fixed = T)[[1]]) {
      # print(line)
      split.line <- strsplit(x = str_trim(line), split = "> ", fixed = T)[[1]]
      split.line[1] <- gsub(pattern = "-", replacement = "_", x = split.line[1])
      if (!split.line[1] %in% names(fields)) {
        print(split.line[1])
      }
      if(split.line[1] %in% c("PASS", "DUPLICATE")) {
        fields[[split.line[1]]] <- as.logical(split.line[2])
      } else {
        fields[[split.line[1]]] <- split.line[2]
      } 
      
      # print(split.line)
      # split.line[1]
      # print()
      # print(str_trim(line))
    }
    p()
    # if (i %% 1000) message(i / length(chunks))
    if (fields$PASS) {
      if (!fields$ID %in% names(pass.clones)) {
        pass.clones[[fields$ID]] <- c(fields$ID)
        parent_to_cloneid[[fields$ID]] <- fields$CLONE
      }
    } else if (fields$DUPLICATE) {
      # if ()
      parent.clone <- gsub(pattern = "Duplication of ", replacement = "", x = fields$FAIL, fixed = T)
      parent.clone <- gsub(pattern = "Collapsed with ", replacement = "", x = parent.clone, fixed = T)
      
      if (parent.clone %in% names(pass.clones)) {
        pass.clones[[parent.clone]] <- c(pass.clones[[parent.clone]], fields$ID)
      } else {
        pass.clones[[parent.clone]] <- c(parent.clone, fields$ID)
        parent_to_cloneid[[parent.clone]] <- fields$CLONE
        
      }
    }
  }
  patient <- patient %||% strsplit(x = names(pass.clones)[1], split = "_")[[1]][[1]]
  # clone <- fields$CLONE
  new_pass_clones <- list()
  new_pass_df <- NULL
  for (i in 1:length(pass.clones)) {
    clone <- parent_to_cloneid[[names(pass.clones)[i]]]
    subclone <- paste(patient, clone, i, sep = "_")
    # new_pass_clones[[subclone]] <- pass.clones[[i]]
    new_pass_clones[[subclone]] <- separateCellnamesVDJ(cellnames = pass.clones[[i]])
    new_pass_clones[[subclone]]$ig_subclone <- subclone
    
    new_pass_clones[[subclone]]$clone <- paste(patient, clone, sep = "_") 
    new_pass_df <- rbind(new_pass_df, new_pass_clones[[subclone]])
  }
  
  return(new_pass_df)
  # return(new_pass_clones)
}