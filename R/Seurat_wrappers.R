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

convertHumanGeneList <- function(genes){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = genes, 
                   mart = human, 
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, 
                   uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  # print(head(humanx))
  return(humanx)
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
seuratFilterWrapper <- function(SRT, min.genes = 200, max.genes = 5000, max.UMI = 30000, 
                                max.mito = 0.25, mito.prefix = "^MT-") {
  require("Seurat")
  RET <- list()
  PLOTS <- list()
  RET[["raw.cell.count"]] <- ncol(SRT@assays[[SRT@active.assay]]@data)
  RET[["min.genes"]] <- min.genes
  RET[["max.genes"]] <- max.genes
  RET[["max.UMI"]] <- max.UMI
  SRT <- AddMetaData(SRT, SRT@meta.data$nCount_RNA, col.name = "nUMI")
  SRT <- AddMetaData(SRT, SRT@meta.data$nFeature_RNA, col.name = "nGene")


  oldCellCount <- ncol(SRT@assays[[SRT@active.assay]]@data )
  SRT <- SubsetData(SRT, subset.name = "nGene", low.threshold = min.genes)
  RET[["under.min.genes"]] <-  oldCellCount - ncol(SRT@assays[[SRT@active.assay]]@data )
  oldCellCount <- ncol(SRT@assays[[SRT@active.assay]]@data)
  SRT <- SubsetData(SRT, subset.name = "nGene", high.threshold = max.genes)
  RET[["over.max.genes"]] <-  oldCellCount - ncol(SRT@assays[[SRT@active.assay]]@data )
  
  mito.genes <- grep(pattern = mito.prefix, x = rownames(x = SRT@assays[[SRT@active.assay]]@data), value = TRUE)
  percent.mito <- Matrix::colSums(SRT@assays$RNA@counts[mito.genes, ]) / Matrix::colSums(SRT@assays$RNA@counts)
  SRT <- AddMetaData(object = SRT, metadata = percent.mito, col.name = "percent.mito")
  
  PLOTS$pre.mito.UMI.filtering <-  VlnPlot(object = SRT,
                                           features = c("nGene", "nUMI", "percent.mito"),
                                           ncol = 3) + ggtitle("Before mito/UMI filtering")
  PLOTS$pre.filter.UMI.mito <- FeatureScatter(object = SRT, feature1 = "nUMI", feature2 = "percent.mito") + ggtitle("Before mito/UMI filtering")
  PLOTS$pre.filter.UMI.nGene <- FeatureScatter(object = SRT, feature1 = "nUMI", feature2 = "nGene") + ggtitle("Before mito/UMI filtering")

  oldCellCount <- ncol(SRT@assays[[SRT@active.assay]]@data )
  SRT <- SubsetData(SRT, subset.name = "nUMI", high.threshold = max.UMI)
  RET[["over.max.UMI"]] <- oldCellCount - ncol(SRT@assays[[SRT@active.assay]]@data )
  
  oldCellCount <- ncol(SRT@assays[[SRT@active.assay]]@data )
  SRT <- SubsetData(SRT, subset.name = "percent.mito", high.threshold = max.mito)

  RET[["over.max.mito"]] <- oldCellCount - ncol(SRT@assays[[SRT@active.assay]]@data )
  PLOTS$post.mito.UMI.filtering <-  VlnPlot(object = SRT,
                                           features = c("nGene", "nUMI", "percent.mito"),
                                           ncol = 3) + ggtitle("After mito/UMI filtering")
  PLOTS$post.filter.UMI.mito <- FeatureScatter(object = SRT, feature1 = "nUMI", feature2 = "percent.mito") + ggtitle("After mito/UMI filtering")
  PLOTS$post.filter.UMI.nGene <- FeatureScatter(object = SRT, feature1 = "nUMI", feature2 = "nGene") + ggtitle("After mito/UMI filtering")


  SRT <- NormalizeData(object = SRT, normalization.method = "LogNormalize", scale.factor = 10000)
  RET[["seurat"]] <- SRT
  RET[["plots"]] <- PLOTS
  RET <- seuratVariableWrapper(RET)
  RET <- seuratPCAWrapper(RET)
  return(RET)
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
#' @export
seuratVariableWrapper <- function(RET) {
  # RET[["seurat"]] <- FindVariableGenes(object = RET[["seurat"]], 
  #                                      do.plot = FALSE,
  #                                      mean.function = ExpMean,
  #                                      dispersion.function = LogVMR,
  #                                      x.low.cutoff = 0.1, x.high.cutoff = 5,
  # #                                      y.cutoff = 0.5)
  # RET[["seurat"]] <- FindVariableFeatures(object = RET[["seurat"]], 
  # 	selection.method = 'mean.var.plot', 
  # 	mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))

	RET[["seurat"]] <- FindVariableFeatures(object = RET[["seurat"]], 
        selection.method = "vst", nfeatures = 5000)

  # FindVariableFeatures(object = RET[["seurat"]],
  #                                      dispersion.cutoff = c(0.5, Inf),
  #                                      mean.cutoff = c(0.1, 8))
  #  
                                       # mean.function = FastExpMean,
                                       # dispersion.function = FastLogVMR,
  RET[["plots"]][["variable.genes"]] <- VariableFeaturePlot(object = RET[["seurat"]])
  # RET[["plots"]][["variable.genes"]] <- VariableGenePlot(object = RET[["seurat"]],
  #                                      x.low.cutoff = 0.1, x.high.cutoff = 5,
  #                                      y.cutoff = 0.5)
  RET[["seurat"]] <- ScaleData(RET[["seurat"]], vars.to.regress = c("nCount_RNA", "percent.mito"))
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
seuratPCAWrapper <- function(RET, do.jackstraw = FALSE) {
  RET[["seurat"]] <- RunPCA(object =  RET[["seurat"]], pc.genes =  RET[["seurat"]]@var.genes, do.print = FALSE,
                  pcs.compute = 50)
  if (do.jackstraw) {
  	RET[["seurat"]] <- JackStraw(RET[["seurat"]], dims = 50)
  	RET$seurat <- ScoreJackStraw(RET$seurat, dims = 50)
  }
  RET[["plots"]][["pc.elbow"]] <- ElbowPlot(RET[["seurat"]], ndims = 50)
  RET[["plots"]][["pca"]] <- DimPlot(RET[["seurat"]], dim.1 = 1, dim.2 = 2, reduction = "pca") 
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
#' @export
seuratClusterWrapper <- function(RET, dims = 1:10, resolution = 0.50) {
  require(Seurat)
  RET$seurat <- FindNeighbors(object = RET$seurat, dims = dims)

  RET$seurat <- FindClusters(RET$seurat, resolution = resolution)
  RET$seurat <- RunTSNE(RET$seurat, dims = dims)
  RET$seurat <- RunUMAP(RET$seurat, dims = dims, reduction = "pca")

  RET$plots[["TSNE"]] <- DimPlot(RET$seurat, label = T, reduction = "tsne")
  RET$plots[["UMAP"]] <- DimPlot(RET$seurat, label = T, reduction = "umap")

  return(RET)
}

#' Wrapper function to find markers between clusters for a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param do.fast whether fast markers should be computed (default: True)
#' @param do.full whether full markers should be computed (default: True)
#' @return list containing clustered Seurat object and TSNE plots
#'
#' @examples
#' seuratAllMarkers(RET)
#'
#' @export
seuratAllMarkers <- function(RET, do.fast = TRUE, do.full = FALSE) {
  require(Seurat)
  if (do.fast) {
    print("Calculating fast markers...")
    RET[["all.markers.quick"]] <- FindAllMarkers(RET$seurat)
  }
  if (do.full) {
    print("Calculating full markers...")
    RET[["all.markers.full"]] <- FindAllMarkers(RET$seurat, logfc.threshold = 0, min.pct = 0)
  }
  return(RET)
}

#' Prints markers distinguishing clusters in a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' 
#' @return list containing clustered Seurat object and TSNE plots
#'
#' @examples
#' printSeuratMarkers(RET)
#'
#' @export
printSeuratMarkers <- function(RET) {
  for (ident in levels(Idents(RET$seurat))) {
    print(ident)
    print(kable(prettyPrintMarkers(RET$all.markers.quick, ident)))
  }
}

#' Renames idents of cells in a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param new.idents list containing new idents (in same order as current)
#' @return list containing renamed Seurat object and re-plotted TSNE 
#'
#' @examples
#' renameIdents(RET)
#'
#' @export
renameIdents <- function(RET, new.idents) {
  require(Seurat)
  # ident.lists = list()
  # names(new.idents) <- levels(Idents(RET$seurat))
  # unlist()
  # RET$seurat <- RenameIdents()
  ident.levels <- levels(RET$seurat)
  for (i in 1:length(ident.levels)) {
  		cells.use <- WhichCells(object = RET$seurat, idents = ident.levels[1])
     	RET$seurat <- SetIdent(object = RET$seurat, cells = cells.use, value = new.idents[i])
  }
  # RET$plots[["TSNE"]] <- TSNEPlot(RET$seurat, do.return = TRUE, do.label = TRUE)
  return(RET)
}
#' Renames idents of cells in a Seurat object
#'
#'
#' @param srt.obj Seurat object
#' @param new.idents list containing new idents (in same order as current)
#' @return list containing renamed Seurat object and re-plotted TSNE 
#'
#' @examples
#' renameSeuratIdents(srt.obj, new.idents)
#'
#' @export
renameSeuratIdents <- function(srt.obj, new.idents) {
  require(Seurat)
  for (i in 1:length(levels(srt.obj@ident))) { 
    srt.obj <- RenameIdent(srt.obj, old.ident.name = i - 1,
                              new.ident.name = new.idents[i])
  }
  # RET$plots[["TSNE"]] <- TSNEPlot(RET$seurat, do.return = TRUE, do.label = TRUE)
  return(srt.obj)
}

#' Makes heatmaps for Seurat object
#'
#'
#' @param RET list containing Seurat object
#' @param marker.lists list of markers (may be genes or module scores)
#' 
#' @return list containing Seurat object and marker heatmaps
#'
#' @examples
#' makeMarkerHeatmaps(RET, marker.lists)
#'
#' @export
makeMarkerHeatmaps <- function(RET, marker.lists) {
  require("Seurat")
  require("ggplot2")
  marker.plots <- list(violin = list(), feature = list(), dotplot = list())
  # RET$plots$markers <- list()
  for (marked.group in names(marker.lists)) {
    # marker.plots$[[marked.group]] <- list()
    markers <- marker.lists[[marked.group]][marker.lists[[marked.group]] %in% rownames(RET$seurat)]

    # RET$plots$markers[[marked.group]]$heatmap <- DoHeatmap(SubsetData(RET$seurat, max.cells.per.ident = 500), 
                                              # features = marker.lists[[marked.group]]) + ggtitle(marked.group)
    marker.plots$violin[[marked.group]] <- VlnPlot(RET$seurat, features = markers) #+ ggplot2::ggtitle(marked.group)
    marker.plots$feature[[marked.group]] <- FeaturePlot(RET$seurat, features = markers, reduction = "umap") #+ ggplot2::ggtitle(marked.group)
    marker.plots$dotplot[[marked.group]] <- DotPlot(RET$seurat, features = markers) + ggplot2::ggtitle(marked.group)
  }
  RET$plots$markers <- marker.plots
  return(RET)
}
#' Prints heatmaps for Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' 
#' @return none
#'
#' @examples
#' RET <- makeMarkerHeatmaps(RET, marker.lists)
#' printMarkerHeatmaps(RET)
#'
#' @export
printMarkerHeatmaps <- function(RET) {
  for(name in names(RET$plots$markers)) {
    print(plot_grid(RET$plots$markers[[name]]$violin, RET$plots$markers[[name]]$feature, ncol=1))
  }
}

