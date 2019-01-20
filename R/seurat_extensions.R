#' gcdSeurat
#'
#' Replacing the simple list I had before with this, hopefully improving memory usage.
#'
#' @slot seurat Unnormalized data such as raw counts or TPMs
#' @slot meta.list Normalized expression data
#' 
#' @import Seurat
#' @name gcdSeurat
#' @rdname gcdSeurat
#' @exportClass gcdSeurat
#'
setClass("gcdSeurat", 
  slots = c(
    seurat = "Seurat", 
    meta.list = "list",
    info = "list"
  )
)



#' Adds module scores to a Seurat object
#'
#'
#' @param RET list containing Seurat object
#' @param modules named list of modules (containing gene names) to be added
#' 
#' @return list containing Seurat object and TSNE plots (tsne.treatment and tsne.ident) of the CCA
#' @examples
#' scoreModules(RET, modules)
#'
#' @export
scoreModules <- function(RET, modules) {
  require("Seurat")
  if (!'modules' %in% names(RET@meta.list)) {
    RET@meta.list$modules <- list()
  }
  i <- 0
  assay.data <- GetAssayData(RET@seurat)
  data.avg <- Matrix::rowMeans(x = assay.data[rownames(RET@seurat), ])
  control.pool <- names(data.avg[data.avg != 0])
  for (module.name in names(modules)) {
    
    i <- i + 1

    # message(class(modules[[module.name]]))
    if (class(modules[[module.name]]) == "character") {
    	modules[[module.name]] = list(modules[[module.name]])
    }
    features = modules[[module.name]]
    features[[1]] = features[[1]][features[[1]] %in% control.pool]
    if (length(features[[1]]) <= 1) {
    	next
    }
    RET@meta.list$modules[[module.name]] <- features
    message("Scoring ", module.name, " (", i, " of ", length(modules), ")")
    RET@seurat <- AddModuleScore(RET@seurat, features = features,
                                 name = module.name, pool = control.pool, 
                                 ctrl = min(vapply(X = features, FUN = length, 
            FUN.VALUE = numeric(length = 1))))
    RET@seurat@meta.data[[module.name]] <- RET@seurat@meta.data[[paste0(module.name, "1")]]
    # message(paste("Completed", i, "of", length(modules)))
    
  }
  return(RET)
}

#' Performs cell cycle scoring
#'
#'
#' @param object Seurat object
#' @param s.features list of features (genes) associated with S phase 
#' @param g2m.features list of features (genes) associated with G2/M phase 
#' @param set.ident whether to set idents of Seurat object 
#' 
#' @return list containing Seurat object and TSNE plots (tsne.treatment and tsne.ident) of the CCA
#'
#' @examples
#' gcdCellCycleScoring(RET, modules)
#'
#' @export
gcdCellCycleScoring <- function (RET, s.features, g2m.features, set.ident = FALSE) 
{
    RET@seurat <- gcdSeuratCellCycleScoring(RET@seurat, s.features, g2m.features, set.ident)
    if (!'modules' %in% names(RET@meta.list)) {
   		RET@meta.list$modules <- list()
  	}
    RET@meta.list$modules["S.Score"] = list(s.features)
    RET@meta.list$modules["G2M.Score"] = list(g2m.features)
    return(RET)
}

#' Performs cell cycle scoring
#'
#'
#' @param object Seurat object
#' @param s.features list of features (genes) associated with S phase 
#' @param g2m.features list of features (genes) associated with G2/M phase 
#' @param set.ident whether to set idents of Seurat object 
#' 
#' @return list containing Seurat object and TSNE plots (tsne.treatment and tsne.ident) of the CCA
#'
#' @examples
#' gcdCellCycleScoring(RET, modules)
#'
#' @export
gcdSeuratCellCycleScoring <- function (object, s.features, g2m.features, set.ident = FALSE) 
{
    name <- "Cell Cycle"
    features <- list(S.Score = s.features, G2M.Score = g2m.features)
    assay.data <- GetAssayData(object)
  	data.avg <- Matrix::rowMeans(x = assay.data[rownames(object), ])
    control.pool <- names(data.avg[data.avg != 0])

    object.cc <- AddModuleScore(object = object, features = features, 
        name = name, ctrl = min(vapply(X = features, FUN = length, 
            FUN.VALUE = numeric(length = 1))), pool = control.pool)
    cc.columns <- grep(pattern = name, x = colnames(x = object.cc[[]]), 
        value = TRUE)
    cc.scores <- object.cc[[cc.columns]]
    rm(object.cc)
    gc(verbose = FALSE)
    assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, 
        first = "S", second = "G2M", null = "G1") {
        if (all(scores < 0)) {
            return(null)
        }
        else {
            if (length(which(x = scores == max(scores))) > 1) {
                return("Undecided")
            }
            else {
                return(c(first, second)[which(x = scores == max(scores))])
            }
        }
    })
    cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
        by = 0)
    colnames(x = cc.scores) <- c("rownames", "S.Score", "G2M.Score", 
        "Phase")
    rownames(x = cc.scores) <- cc.scores$rownames
    cc.scores <- cc.scores[, c("S.Score", "G2M.Score", "Phase")]
    object[[colnames(x = cc.scores)]] <- cc.scores
    if (set.ident) {
        object[["old.ident"]] <- Idents(object = object)
        Idents(object = object) <- "Phase"
    }
    object$CC.Difference <- object$S.Score - object$G2M.Score
    return(object)
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
  require("future")
  require("future.apply")
  all.correlations <- list()
  exp_matrix <- as.matrix(GetAssayData(RET@seurat))
  if (only.var) {
  	message("Subsetting to variable genes only")
    exp.sd <- apply(exp_matrix, 1, function(x){sd(as.numeric(x))})
    exp_matrix <- exp_matrix[exp.sd != 0,]
  }
  # genes.to.use <- rownames(exp_matrix)
  exp_matrix <- rbind(exp_matrix,
                      as.matrix(t(FetchData(RET@seurat,
                                            vars = names(RET@meta.list$modules)))))
  exp_df <- data.frame(exp_matrix)
  
  # combining gene expression and module scores into a single matrix
  # vars.all = c(genes.to.use, names(RET@meta.list$modules))
  idents.exprs = list()
  for (ident in levels(as.factor(Idents(RET@seurat)))) {
    message("Subsetting for ident: ", ident)
    idents.exprs[[ident]] <- as.matrix(dplyr::select(exp_df, 
                                              one_of(names(Idents(RET@seurat)[Idents(RET@seurat)== ident]))))
    ident.sd <- apply(idents.exprs[[ident]], 1, function(x){sd(as.numeric(x))})
    idents.exprs[[ident]] <- idents.exprs[[ident]][ident.sd != 0,]
  }
  all.var.corrs <- list()
  
  if (!'module.corrs' %in% names(RET@meta.list)) {
    message("erasing old corrs")
    RET@meta.list$module.corrs <- list()
  }
  i <- 1
  for (var.name in vars.to.corr) {
  	message("Calculating correlations for ", var.name, " [", i, " of ", length(vars.to.corr), "]")
    # correlations across all idents
    var.corrs <- apply(exp_matrix, 1, function(x){cor(as.numeric(exp_matrix[var.name,]), x)})
    var.df <- data.frame(var.corrs, row.names = names(var.corrs))
    colnames(var.df) <- c("All cells")
    var.df <- rownames_to_column(var.df, var = "gene")
    i <- i + 1
    var.ident.corrs <- list()

    for (ident in names(idents.exprs)) {
      ident.expr <- idents.exprs[[ident]]
      if (var.name %in% row.names(ident.expr)) {
      	ident.corrs <- apply(ident.expr, 1, function(x){cor(as.numeric(ident.expr[var.name,]), x)})
	      var.ident.df <- data.frame(ident.corrs, row.names = names(ident.corrs))
	      colnames(var.ident.df) <- c(ident)
	      var.ident.df <- rownames_to_column(var.ident.df, var = "gene")
	      var.df <- left_join(var.df, var.ident.df, by = "gene")
      }
    }
    RET@meta.list$module.corrs[[var.name]] <- var.df
  }
  return(RET)
}

getCorrLists <- function(vars.to.corr, exp_matrix, ident.name = "All cells") {
	calcCorrs <- function(var.name) {
  	var.expr <- as.numeric(exp_matrix[var.name,])
  	var.corrs <- apply(exp_matrix, 1, function(x){cor(var.expr, x)})
    var.df <- data.frame(var.corrs, row.names = names(var.corrs))
    colnames(var.df) <- ident.name
    var.df <- rownames_to_column(var.df, var = "gene")
    return(var.df)
  }
	module.corrs <- future_lapply(vars.to.corr, calcCorrs)
	names(module.corrs) <- vars.to.corr
	return(module.corrs)
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
#' correlateVariables(RET, vars.to.corr, True)
#'
#' @export
correlateVariables <- function(RET, vars.to.corr, only.var = TRUE, split.by = NULL) {
  require("dplyr")
  require("tibble")
  require("future.apply")
  all.correlations <- list()
  exp_matrix <- as.matrix(GetAssayData(RET@seurat))
  if (only.var) {
  	message("Subsetting to variable genes only")
  	SRT <- FindVariableFeatures(RET@seurat,	selection.method = "vst", nfeatures = 5000)
  	vars <- VariableFeatures(SRT)
  	vars <- unique(c(vars, vars.to.corr))
  	exp_matrix <- as.matrix(t(FetchData(RET@seurat, vars = vars)))
    # exp.sd <- apply(exp_matrix, 1, function(x){sd(as.numeric(x))})
    # exp_matrix <- exp_matrix[exp.sd != 0,]
  }
  if ("modules" %in% names(RET@meta.list)) {
  	exp_matrix <- rbind(exp_matrix,
                      as.matrix(t(FetchData(RET@seurat,
                                            vars = names(RET@meta.list$modules)))))
  }

	module.corrs <- getCorrLists(vars.to.corr, exp_matrix)
	if (!is.null(split.by)) {
		RET@seurat$old.idents <- Idents(RET@seurat)
		Idents(RET@seurat) <- split.by
	}
	for (ident in levels(RET@seurat)) {
		message("Calculating corrs for ident:",ident)
  	cell.names <- WhichCells(RET@seurat, idents = c(ident))
  	ident.exprs <- exp_matrix[, c(cell.names)]
  	ident.exprs.sd <- apply(ident.exprs, 1, function(x){sd(as.numeric(x))})
    ident.exprs <- ident.exprs[ident.exprs.sd != 0,]
    expressed.vars <- vars.to.corr[vars.to.corr %in% rownames(ident.exprs)]
    if (length(expressed.vars) >= 1) {
    	ident.corrs <- getCorrLists(expressed.vars, ident.exprs, ident.name = ident)
  		for (var in expressed.vars) {
  			module.corrs[[var]] <- dplyr::left_join(module.corrs[[var]], ident.corrs[[var]], by = "gene")
  		}
  	}
  }
  if (!is.null(split.by)) {
		Idents(RET@seurat) <- RET@seurat$old.idents
	}
  RET@meta.list$module.corrs <- module.corrs
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


#' Identifies light chain identity of sample
#'
#'
#' @param RET list containing Seurat object
#' @param sample.variable name distinguishing sample idents
#'
#' @return updated \code{RET}
#'
#' @examples
#' RET <- findLikelyLightChainIdent(RET)
#'
#' @export
findLikelyLightChainIdent <- function(RET, sample.variable = NULL, mouse = FALSE, tumor.idents = NULL) {
	ig.list <- list()
	# if (!is.null(tumor.idents)) {
	# 	print(tumor.idents)
	# 	SRT <- subset(RET@seurat, idents = tumor.idents)
	# }
	ig.genes <- list(heavy.V = list(), heavy.C = list(), light.V = list(), light.C = list())
	if (!is.null(sample.variable)) {
		RET@seurat$old.idents <- Idents(RET@seurat)
		Idents(RET@seurat) <- sample.variable
		for (sample in levels(Idents(RET@seurat))) {
			sub.SRT <- subset(RET@seurat, idents = c(sample))
			Idents(sub.SRT) <- "old.idents"
			igs <- seuratFindLikelyLightChainIdent(sub.SRT, mouse, tumor.idents)
			ig.list[[sample]] <- igs
			ig.genes$heavy.V <- unique(c(unlist(ig.genes$heavy.V), igs$heavy.V))
			ig.genes$heavy.C <- unique(c(unlist(ig.genes$heavy.C), igs$heavy.C))
			ig.genes$light.V <- unique(c(unlist(ig.genes$light.V), igs$light.V))
			ig.genes$light.C <- unique(c(unlist(ig.genes$light.C), igs$light.C))

		}	
		Idents(RET@seurat) <- RET@seurat$old.idents
	} else {
		ig.list[["all"]] <- seuratFindLikelyLightChainIdent(RET@seurat, mouse, tumor.idents)
	}
	RET@meta.list$ig.genes <- ig.list
	RET@meta.list$all.igs <- ig.genes
	return(RET)
}

#' Identifies light chain identity of sample
#'
#'
#' @param SRT Seurat object
#' @param mouse whether mouse (e.g. Igh...) or human (IGH...) gene names should be assumed
#' @param tumor.idents identities of clusters containing tumor cells
#' 
#' @return IG.light.V, IG.light.C
#'
#' @examples
#' RET <- findLikelyLightChainIdent(RET)
#'
#' @export
seuratFindLikelyLightChainIdent <- function(SRT, mouse = FALSE, tumor.idents = NULL) {
	if (mouse) {
		IGH.genes <- c("Igha", "Ighd", "Ighe", "Ighg1", "Ighg2b", "Ighg2c", "Ighg3", "Ighm")
		light.pattern.V <- "^Ig[lk]v"
		light.pattern.C <- "^Ig[lk]c"
		IGHV.pattern <- "^Ighv"
	} else {
		IGH.genes <- c("IGHA1", "IGHA2", "IGHD", "IGHE", 
										"IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM")
		light.pattern.V <- "^IG[LK]V"
		light.pattern.C <- "^IG[LK]C"
		IGHV.pattern <- "^IGHV"
	}
	if (!is.null(tumor.idents)) {
		SRT <- subset(SRT, idents = tumor.idents)
	}
  avg.exprs <- as.data.frame(Matrix::rowMeans(GetAssayData(SRT, slot = "counts")))
  colnames(avg.exprs) <- c("avg_expr")
  avg.exprs <- tibble::rownames_to_column(avg.exprs, var = "gene")
  IG.C.genes <- grep(pattern = light.pattern.C, x = avg.exprs$gene, value = TRUE)
  IGH.C.exprs <- dplyr::filter(avg.exprs, gene %in% IGH.genes)

  IGH.V.genes <- grep(pattern = IGHV.pattern, x = avg.exprs$gene, value = TRUE) 
  IGH.V.exprs <- dplyr::filter(avg.exprs, gene %in% IGH.V.genes)

  IG.V.genes <- grep(pattern = light.pattern.V, x = avg.exprs$gene, value = TRUE)
  IG.C.exprs <- dplyr::filter(avg.exprs, gene %in% IG.C.genes)
  IG.V.exprs <- dplyr::filter(avg.exprs, gene %in% IG.V.genes)
  # print(IG.V.exprs)
  top.C <- dplyr::top_n(IG.C.exprs, n=3, wt = avg_expr) %>% dplyr::arrange(desc(avg_expr))
  top.V <- dplyr::top_n(IG.V.exprs, n=3, wt = avg_expr) %>% dplyr::arrange(desc(avg_expr))
  heavy.C <- dplyr::top_n(IGH.C.exprs, n=3, wt = avg_expr) %>% dplyr::arrange(desc(avg_expr))
  heavy.V <- dplyr::top_n(IGH.V.exprs, n=3, wt = avg_expr) %>% dplyr::arrange(desc(avg_expr))
  message(top.C)
  message(top.V)
  message(heavy.C)
  message(heavy.V)
  top.C <- dplyr::top_n(IG.C.exprs, n=1, wt = avg_expr)
  top.V <- dplyr::top_n(IG.V.exprs, n=1, wt = avg_expr)
  heavy.C <- dplyr::top_n(IGH.C.exprs, n=1, wt = avg_expr)
  heavy.V <- dplyr::top_n(IGH.V.exprs, n=1, wt = avg_expr)
  return(list(light.V = top.V$gene, light.C = top.C$gene, heavy.C = heavy.C$gene, heavy.V = heavy.V$gene))
}

#' Finds DE genes in each cluster between conditions 
#'
#'
#' @param RET list containing aligned Seurat object
#' @param cond.var name of variable storing condition values in Seurat object's metadata.
#' @param cond.1 first condition
#' @param cond.2 second condition 
#'
#' @return dataframe containing DE markers between conditions
#'
#' @examples
#' markersBetweenConditions(RET)
#'
#' @export
markersBetweenConditions <- function(RET, cond.var, cond.1, cond.2) {
  cond.markers <- seuratMarkersBetweenConditions(RET@seurat, cond.var, cond.1, cond.2)
  if (!"markers" %in% names(RET@meta.list)) {
    RET@meta.list$markers <- list()
  }
  RET@meta.list$markers[[paste0(cond.1, "_vs_", cond.2)]] <- cond.markers   
  return(RET)
}


#' Finds DE genes in each cluster between conditions 
#'
#'
#' @param SRT list containing aligned Seurat object
#' @param cond.var name of variable storing condition values in Seurat object's metadata.
#' @param cond.1 first condition
#' @param cond.2 second condition 
#'
#' @return dataframe containing DE markers between conditions
#'
#' @examples
#' seuratMarkersBetweenConditions(SRT)
#'
#' @export
seuratMarkersBetweenConditions <- function(SRT, cond.var, cond.1, cond.2) {
  SRT$celltype.cond <- paste0(Idents(SRT), "_", SRT[[cond.var]][[cond.var]])
  cond.levels <- levels(SRT[[cond.var]])
  old.idents <- levels(as.factor(Idents(SRT)))
  Idents(SRT) <- "celltype.cond"
  new.idents <- levels(as.factor(Idents(SRT)))
  num.idents <- length(levels(old.idents))
  cond.markers <- NULL
  for (old.id in old.idents) {
    id.1 <- paste0(old.id, "_", cond.1)
    id.2 <- paste0(old.id, "_", cond.2)
    message(id.1)
    # cond.markers[[old.id]] <- FindMarkers()
    if (id.1 %in% new.idents && id.2 %in% new.idents) {
      SRT.subset <- subset(SRT, idents = c(id.1, id.2))
      if (ncol(subset(SRT.subset, idents = c(id.1))) > 3 && ncol(subset(SRT.subset, idents = c(id.2))) > 3) {
        markers <- FindMarkers(SRT.subset, ident.1 = id.1, ident.2 = id.2, 
                          min.pct = 0, logfc.threshold = 0.05, verbose = T)
        markers <- tibble::rownames_to_column(markers, var="gene")
        cond.markers <- rbind(cond.markers, dplyr::mutate(markers, cluster = old.id))
      }
    }
  }
  return(cond.markers)
}