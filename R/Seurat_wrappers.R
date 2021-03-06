#' Creates table of marker genes differentially expressed within cells belonging to ident
#'
#'
#' @param markers data frame returned by \code{\link[Seurat]{FindAllMarkers}}
#' @param ident cluster/group of interest
#'
#' @return Well-formatted data frame
#'
#' @examples
#' prettyPrintMarkers()
#'
#' @export
prettyPrintMarkers <- function(markers, ident) {
  # require("knitr")
  return(filter(markers, cluster == ident, p_val_adj < 1E-25) %>%
    top_n(n = -10, wt = p_val_adj) %>% arrange(p_val_adj) %>%
    dplyr::select(-p_val, -cluster) %>%
    mutate(avg_logFC = round(avg_logFC, 2), p_val_adj = formatC(p_val_adj, format = "e", digits = 2)))
}

#' Creates table of pathways enriched within cells belonging to ident
#'
#' 
#'
#' @param pathways data frame returned by \code{\link[fgsea]{fgsea}}
#' @param ident cluster/group of interest
#'
#' @return Well-formatted data frame
#'
#' @examples
#' prettyPrintPathways()
#'
#' @export
prettyPrintPathways <- function(pathways, ident) {
  return(filter(pathways, cluster == ident) %>% 
           dplyr::select(-Genes, -set, -cluster) %>%
           top_n(-10, Adjusted.P.value) %>% arrange(Adjusted.P.value) %>%
           mutate_if(is.numeric, funs(formatC(., format = "e", digits = 2))))
}

#' Wrapper function to initiate Seurat analysis from count matrices
#'
#'
#' @param read.counts sparse matrix of read counts (returned by \code{\link[Seurat]{Read10X}})
#' @param project name for Seurat object (required)
#' @param min.genes param for filtering
#' @param max.genes param for filtering
#' @param max.UMI param for filtering
#' @param max.mito param for filtering
#' @param mito.prefix param for filtering
#' @return Output of seuratFilterWrapper
#'
#' @examples
#' prettyPrintMarkers()
#'
#' @export
seuratStartFromCounts <- function(read.counts, project, min.genes = 200, 
                                  max.genes = 5000, max.UMI = 30000, 
                                  max.mito = 0.25, mito.prefix = "^MT-") {
  # removes prefixes added when using mixed references (if present)
  row.names(read.counts) <- gsub("hg19_","",row.names(read.counts))
  # row.names(read.counts) <- gsub("mm10_","",row.names(read.counts))
  SRT <- CreateSeuratObject(counts = read.counts, project = project)
  return(seuratFilterWrapper(SRT, min.genes = min.genes, max.genes = max.genes, max.UMI = max.UMI,
                             max.mito = max.mito, mito.prefix = mito.prefix))
}

#' Wrapper function to perform initial filtering and normalization of a Seurat object
#'
#'
#' @param SRT Seurat object 
#' @param min.genes param for filtering
#' @param max.genes param for filtering
#' @param max.UMI param for filtering
#' @param max.mito param for filtering
#' @param mito.prefix prefix for mitochondrial genes ("^MT-" for hg19)
#' @return list containing plots and filtered Seurat object
#'
#' @examples
#' seuratFilterWrapper(SRT, min.genes = min.genes, max.genes = max.genes, max.UMI = max.UMI,
#'    max.mito = max.mito, mito.prefix = mito.prefix)
#'
#' @export
seuratFilterWrapper <- function(SRT, min.genes = 200, max.genes = 5000, min.UMI = 0, max.UMI = 30000, 
                                max.mito = 0.25, mito.prefix = "^MT-", erythro.genes = NULL, max.erythro = 0.01) {
  require("Seurat")
  RET <- list(params = list(), plots = list(), qc = list())
  INFO <- list(params = list(), qc = list())
  INFO$params$min.genes <- min.genes
  INFO$params$max.genes <- max.genes
  INFO$params$max.UMI <- max.UMI
  INFO$params$min.UMI <- min.UMI
  INFO$params$max.mito <- max.mito
  INFO$params$erythro.genes <- erythro.genes
  INFO$params$max.erythro <- max.erythro

  INFO$qc$raw.cell.count <- ncol(GetAssayData(SRT))
  if (!"nUMI" %in% names(SRT@meta.data)) {
  	SRT <- AddMetaData(SRT, SRT$nCount_RNA, col.name = "nUMI")
  }
  if (!"nGene" %in% names(SRT@meta.data)) {
  	SRT <- AddMetaData(SRT, SRT$nFeature_RNA, col.name = "nGene")
  }
  
  mito.genes <- grep(pattern = mito.prefix, x = rownames(GetAssayData(SRT)), value = TRUE)
  percent.mito <- Matrix::colSums(SRT@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(SRT@assays$RNA@counts)
  SRT <- AddMetaData(object = SRT, metadata = percent.mito, col.name = "percent.mito")

  if (!is.null(erythro.genes)) {
  	erythro.fil <- erythro.genes[erythro.genes %in% rownames(GetAssayData(SRT))]
  	percent.erythro <- Matrix::colSums(SRT@assays$RNA@counts[erythro.fil, ]) / Matrix::colSums(SRT@assays$RNA@counts)
 		SRT <- AddMetaData(object = SRT, metadata = percent.erythro, col.name = "percent.erythro")
  	INFO$qc$too.much.erythro <- sum(SRT$percent.erythro > max.erythro)
  	# RET$plots$erythro <- VlnPlot(object = SRT, 
  	# 	pt.size = 0, features = c("percent.erythro")) + ggplot2::ggtitle("Before mito/UMI filtering")
  }

  INFO$qc$too.few.genes <- sum(SRT$nGene < min.genes) 
  INFO$qc$too.many.genes <- sum(SRT$nGene > max.genes) 
  INFO$qc$too.few.UMIs <- sum(SRT$nUMI < min.UMI) 
  INFO$qc$too.many.UMIs <- sum(SRT$nUMI > max.UMI) 
  INFO$qc$too.much.mito <- sum(SRT$percent.mito > max.mito) 
  
  # RET$plots$pre.mito.UMI.filtering <-  VlnPlot(object = SRT, pt.size = 0, 
  # 	features = c("nGene", "nUMI", "percent.mito"), ncol = 3) + ggplot2::ggtitle("Before mito/UMI filtering")

  # p <-  FeatureScatter(object = SRT, feature1 = "nUMI", feature2 = "percent.mito")
  # p <- p + ggplot2::labs(title = "Before mito/UMI filtering", subtitle = p$labels$title) 
  # RET$plots$pre.filter.UMI.mito <- p
  # 
  # p <-  FeatureScatter(object = SRT, feature1 = "nUMI", feature2 = "nGene")
  # p <- p + ggplot2::labs(title = "Before mito/UMI filtering", subtitle = p$labels$title)
  # RET$plots$pre.filter.UMI.nGene <- p

  if (is.null(erythro.genes)) {
  	cells.keep = colnames(SRT)[SRT$nGene > min.genes & 
  		SRT$nGene < max.genes & 
  		SRT$nUMI < max.UMI & 
  		SRT$nUMI > min.UMI & 
  		SRT$percent.mito < max.mito]
  } else {
  	cells.keep = colnames(SRT)[SRT$nGene > min.genes & 
  		SRT$nGene < max.genes & 
  		SRT$nUMI < max.UMI & 
  		SRT$nUMI > min.UMI & 
  		SRT$percent.mito < max.mito &
  		SRT$percent.erythro < max.erythro]
  }  
  SRT <- subset(SRT, cells = cells.keep)
	INFO$qc$passing.filters <- ncol(SRT)
  # RET$plots$post.mito.UMI.filtering <-  VlnPlot(object = SRT, pt.size = 0,
  #                                          features = c("nGene", "nUMI", "percent.mito"),
  #                                          ncol = 3) + ggplot2::ggtitle("After mito/UMI filtering")
  # p <-  FeatureScatter(object = SRT, feature1 = "nUMI", feature2 = "percent.mito")
  # p <- p + ggplot2::labs(title = "After mito/UMI filtering", subtitle = p$labels$title) 
  # RET$plots$post.filter.UMI.mito <- p
# 
#   p <-  FeatureScatter(object = SRT, feature1 = "nUMI", feature2 = "nGene")
#   p <- p + ggplot2::labs(title = "After mito/UMI filtering", subtitle = p$labels$title)
#   RET$plots$post.filter.UMI.nGene <- p
  SRT <- NormalizeData(object = SRT, normalization.method = "LogNormalize", scale.factor = 10000)
  RET.gcdSeurat <- new(Class = "gcdSeurat", 
  	seurat = SRT, meta.list = RET, info = INFO)
  # RET.gcdSeurat <- seuratVariableWrapper(RET.gcdSeurat)
  # RET.gcdSeurat <- seuratPCAWrapper(RET.gcdSeurat)
  return(RET.gcdSeurat)
}

#' Prints filtering info for a gcdSeurat object, describing how many cells pass each filter. 
#'
#'
#' @param RET Seurat object 
#'
#' @examples
#' seuratFilterWrapper(SRT, min.genes = min.genes, max.genes = max.genes, max.UMI = max.UMI,
#'    max.mito = max.mito, mito.prefix = mito.prefix)
#'
#' @export
gcdPrintSeuratQC <- function(RET) {
  message(sprintf("Out of %d initial cells, %d (%.2f%%) pass all filters.", 
  	RET@info$qc$raw.cell.count, RET@info$qc$passing.filters, 
  	100 * RET@info$qc$passing.filters / RET@info$qc$raw.cell.count))
  message(sprintf("%d (%.2f%%) of all cells have under %d detected UMIs.", 
  	RET@info$qc$too.few.UMIs, 
  	100 * RET@info$qc$too.few.UMIs / RET@info$qc$raw.cell.count, 
  	RET@info$params$min.UMI))
  message(sprintf("%d (%.2f%%) of all cells have over %d detected UMIs.", 
  	RET@info$qc$too.many.UMIs, 
  	100 * RET@info$qc$too.many.UMIs / RET@info$qc$raw.cell.count, 
  	RET@info$params$max.UMI))
  message(sprintf("%d (%.2f%%) of all cells have under %d detected genes.", 
  	RET@info$qc$too.few.genes, 
  	100 * RET@info$qc$too.few.genes / RET@info$qc$raw.cell.count, 
  	RET@info$params$min.genes))
  message(sprintf("%d (%.2f%%) of all cells have over %d detected genes.", 
  	RET@info$qc$too.many.genes, 
  	100 * RET@info$qc$too.many.genes / RET@info$qc$raw.cell.count, 
  	RET@info$params$max.genes))
  message(sprintf("%d (%.2f%%) of all cells have mitochondrial read fractions greater than %.2f.", 
  	RET@info$qc$too.much.mito, 
  	100 * RET@info$qc$too.much.mito / RET@info$qc$raw.cell.count, 
  	RET@info$params$max.mito))
  if (!is.null(RET@info$params$erythro.genes)) {
  	message(sprintf("%d (%.2f%%) of all cells have erythrocytic read fractions greater than %.2f.", 
  	RET@info$qc$too.much.erythro, 
  	100 * RET@info$qc$too.much.erythro / RET@info$qc$raw.cell.count, 
  	RET@info$params$max.erythro))
  }
}

#' Wrapper function to find variable genes in Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @return list containing scaled Seurat object and plot of variable genes
#'
#' @examples
#' seuratVariableWrapper(RET)
#'
#' @importFrom Seurat ScaleData FindVariableFeatures NormalizeData
#'
#' 
#' @export
seuratVariableWrapper <- function(RET, nfeatures = 2500, vars.to.regress = NULL, do.normalize = F,
                                  scale.all.features = TRUE) {
  if (do.normalize) {
    RET@seurat <- NormalizeData(RET@seurat, normalization.method = "LogNormalize", scale.factor = 10000)
  }
	RET@seurat <- FindVariableFeatures(object = RET@seurat, selection.method = "vst", nfeatures = nfeatures)
  # RET@plots[["variable.genes"]] <- VariableFeaturePlot(object = RET@seurat)
	if (scale.all.features) {
	  features.scale <- rownames(RET@seurat)
	} else {
	  features.scale <- NULL
	}
	message(features.scale)
	
	message(length(features.scale))
  RET@seurat <- ScaleData(RET@seurat, vars.to.regress = vars.to.regress, features = features.scale)
  return(RET)
}

#' Wrapper function to perform PCA and generate PCA plots for a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param do.jackstraw whether JackStraw analysis is performed
#'
#' @return list containing post-PCA Seurat object and PCA plots
#'
#' @examples
#' seuratPCAWrapper(RET)
#'
#' @export
seuratPCAWrapper <- function(RET, do.jackstraw = FALSE, pcs.compute = 30) {
  RET@seurat <- RunPCA(object =  RET@seurat, pc.genes =  RET@seurat@var.genes, do.print = FALSE,
                  pcs.compute = pcs.compute)
  if (do.jackstraw) {
  	RET@seurat <- JackStraw(RET@seurat, dims = pcs.compute)
  	RET@seurat <- ScoreJackStraw(RET@seurat, dims = pcs.compute)
  }
  # RET@plots$pc.elbow <- ElbowPlot(RET@seurat, ndims = 50)
  # RET@plots$pca <- DimPlot(RET@seurat, dim.1 = 1, dim.2 = 2, reduction = "pca") 
  return(RET)
}

#' Wrapper function to perform clustering and run TSNE on a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param dims dimensions to be used for clustering
#' @param resolution resolution to be used for clustering
#'
#' @return list containing clustered Seurat object and TSNE plots
#'
#' @examples
#' seuratClusterWrapper(RET)
#' 
#' @importFrom Seurat FindNeighbors FindClusters RunTSNE RunUMAP BuildClusterTree
#'
#' @export
seuratClusterWrapper <- function(RET, dims = NULL, resolution = 0.8, do.TSNE = T, do.UMAP = T, do.tree = T) {
  if (is.null(dims)) {
    if (!is.null(RET@meta.list$dims)) {
      dims <- RET@meta.list$dims
    } else {
      dims <- 1:10
    }
  } else {
    RET@meta.list$dims <- dims
  }
  RET@seurat <- FindNeighbors(object = RET@seurat, dims = dims)

RET@seurat <- FindClusters(RET@seurat, resolution = resolution)
  if (do.TSNE) {
   	RET@seurat <- RunTSNE(RET@seurat, dims = dims)
  	# RET@plots[["TSNE"]] <- DimPlot(RET@seurat, label = T, reduction = "tsne", repel = T)
  }
  if (do.UMAP) {
    RET@seurat <- RunUMAP(RET@seurat, dims = dims)
    # RET@plots[["UMAP"]] <- DimPlot(RET@seurat, label = T, reduction = "umap", repel = T)
  } 
  if (do.tree) {
    RET@seurat <- BuildClusterTree(RET@seurat, dims = dims)
  }
  return(RET)
}
 
#' Subsets a gcdSeurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param cells cells to subset
#'
#' @return subsetted gcdSeurat
#'
#' @examples
#' gcdSubsetSeurat(RET)
#'
#' @export
gcdSubsetSeurat <- function(RET, ...) {
  RET@seurat <- subset(RET@seurat, ...)
  RET@markers <- list()
  RET@fgsea <- list()
  if (!is.null(Tool(RET@seurat, slot = "BuildClusterTree"))) {
    Tool(RET@seurat, slot = "BuildClusterTree") <- NULL
  }
  RET@seurat <- NormalizeData(RET@seurat)
  return(RET)
}





