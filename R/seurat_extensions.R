#' Adds module scores to a Seurat object
#'
#'
#' @param RET list containing Seurat object
#' @param modules named list of modules (containing gene names) to be added
#' 
#' @return list containing Seurat object and TSNE plots (tsne.treatment and tsne.ident) of the CCA
#'
#' @examples
#' scoreModules(RET, modules)
#'
#' @export
scoreModules <- function(RET, modules) {
  require("Seurat")
  if (!'modules' %in% names(RET)) {
    RET$modules <- list()
  }
  i <- 1
  assay.data <- GetAssayData(RET$seurat)
  data.avg <- Matrix::rowMeans(x = assay.data[rownames(RET$seurat), ])
  control.pool <- names(data.avg[data.avg != 0])

  for (module.name in names(modules)) {
    RET$modules[[module.name]] <- modules[[module.name]]
    print(paste("Scoring ", module.name))
    RET$seurat <- AddModuleScore(RET$seurat, features = modules[[module.name]],
                                 ctrl = 100, name = module.name, pool = control.pool)
    RET$seurat@meta.data[[module.name]] <- RET$seurat@meta.data[[paste0(module.name, "1")]]
    print(paste("Completed", i, "of", length(modules)))
    i <- i + 1
  }
  return(RET)
}

#' Correlates modules/genes 
#'
#'
#' @param RET list containing Seurat object
#' @param vars.to.corr variables to correlate ()
#' @param only.var if True, only consider named variables and not all modules
#' 
#' @return list containing Seurat object and TSNE plots (tsne.treatment and tsne.ident) of the CCA
#'
#' @examples
#' correlationAcrossIdents(RET, vars.to.corr, True)
#'
#' @export
correlationAcrossIdents <- function(RET, vars.to.corr, only.var = TRUE) {
  require("dplyr")
  require("tibble")
  all.correlations <- list()
  exp_matrix <- as.matrix(GetAssayData(RET$seurat))
  if (only.var) {
    exp.sd <- apply(exp_matrix, 1, function(x){sd(as.numeric(x))})
    exp_matrix <- exp_matrix[exp.sd != 0,]
  }
  # genes.to.use <- rownames(exp_matrix)
  exp_matrix <- rbind(exp_matrix,
                      as.matrix(t(FetchData(RET$seurat,
                                            vars = names(RET$modules)))))
  exp_df <- data.frame(exp_matrix)
  
  # combining gene expression and module scores into a single matrix
  # vars.all = c(genes.to.use, names(RET$modules))
  idents.exprs = list()
  for (ident in levels(as.factor(Idents(RET$seurat)))) {
    print(paste("Subsetting for ident:", ident))
    idents.exprs[[ident]] <- as.matrix(dplyr::select(exp_df, 
                                              one_of(names(Idents(RET$seurat)[Idents(RET$seurat)== ident]))))
    ident.sd <- apply(idents.exprs[[ident]], 1, function(x){sd(as.numeric(x))})
    idents.exprs[[ident]] <- idents.exprs[[ident]][ident.sd != 0,]
  }
  all.var.corrs <- list()
  
  if (!'module.corrs' %in% names(RET)) {
    print("erasing old corrs")
    RET$module.corrs <- list()
  }
  i <- 1
  for (var.name in vars.to.corr) {
  	print(paste0("Calculating correlations for ", var.name, " [", i, " of ", length(vars.to.corr), "]"))
    # print(exp_matrix[var.name,])
    var.corrs <- apply(exp_matrix, 1, function(x){cor(as.numeric(exp_matrix[var.name,]), x)})
    var.df <- data.frame(var.corrs, row.names = names(var.corrs))
    colnames(var.df) <- c("All cells")
    var.df <- rownames_to_column(var.df, var = "gene")
    i <- i + 1
    for (ident in names(idents.exprs)) {
      # print(paste("Ident:", ident))
      ident.expr <- idents.exprs[[ident]]
      ident.corrs <- apply(ident.expr, 1, function(x){cor(as.numeric(ident.expr[var.name,]), x)})
      var.ident.df <- data.frame(ident.corrs, row.names = names(ident.corrs))
      colnames(var.ident.df) <- c(ident)
      var.ident.df <- rownames_to_column(var.ident.df, var = "gene")
      var.df <- left_join(var.df, var.ident.df, by = "gene")
    }
    RET$module.corrs[[var.name]] <- var.df
  }
  return(RET)
}

# #' returns list of pathways
# #'
# #'
# #' @param RET list containing Seurat object
# #' @param vars.to.corr variables to correlate ()
# #' @param only.var if True, only consider named variables and not all modules
# #' 
# #' @return list containing Seurat object and TSNE plots (tsne.treatment and tsne.ident) of the CCA
# #'
# #' @examples
# #' correlationAcrossIdents(RET, vars.to.corr, True)
# #'
# #' @export
# getPathwayLists <- function() {
#   require("fgsea")
#   hallmark.pathways <- gmtPathways("scRNASEQ/MSIGDB/h.all.v6.2.symbols.gmt")
#   pos.pathways <- gmtPathways("scRNASEQ/MSIGDB/c1.all.v6.2.symbols.gmt")
#   CP.pathways <- gmtPathways("scRNASEQ/MSIGDB/c2.cp.v6.2.symbols.gmt")
#   motif.pathways <- gmtPathways("scRNASEQ/MSIGDB/c3.all.v6.2.symbols.gmt")
#   GO.BP.pathways <- gmtPathways("scRNASEQ/MSIGDB/c5.bp.v6.2.symbols.gmt")
#   return(list(Hallmark = hallmark.pathways,
#               Positional = pos.pathways,
#               `Canonical Pathways` = CP.pathways, 
#               `microRNA & TF motifs` = motif.pathways,
#               `GO biological process` = GO.BP.pathways))
# }