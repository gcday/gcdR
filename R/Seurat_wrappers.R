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
  SRT <- AddMetaData(SRT, SRT@meta.data$nCount_RNA, col.name = "nGene")
  SRT <- AddMetaData(SRT, SRT@meta.data$nFeature_RNA, col.name = "nUMI")

  # SRT@meta.data$nUMI <- SRT@meta.data$nCount_RNA
  # SRT@meta.data$nGene <- SRT@meta.data$nFeature_RNA

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
  
  oldCellCount <- ncol(SRT@assays[[SRT@active.assay]]@data )
  SRT <- SubsetData(SRT, subset.name = "nUMI", high.threshold = max.UMI)
  RET[["over.max.UMI"]] <- oldCellCount - ncol(SRT@assays[[SRT@active.assay]]@data )
  
  oldCellCount <- ncol(SRT@assays[[SRT@active.assay]]@data )
  SRT <- SubsetData(SRT, subset.name = "percent.mito", high.threshold = max.mito)

  RET[["over.max.mito"]] <- oldCellCount - ncol(SRT@assays[[SRT@active.assay]]@data )
  PLOTS$post.mito.UMI.filtering <-  VlnPlot(object = SRT,
                                           features = c("nGene", "nUMI", "percent.mito"),
                                           ncol = 3) + ggtitle("After mito/UMI filtering")
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
  #                                      y.cutoff = 0.5)
  RET[["seurat"]] <- FindVariableFeatures(object = RET[["seurat"]],
                                       dispersion.cutoff = c(0.5, Inf),
                                       mean.cutoff = c(0.1, 8))
  #  
                                       # mean.function = FastExpMean,
                                       # dispersion.function = FastLogVMR,
  RET[["plots"]][["variable.genes"]] <- VariableFeaturePlot(object = RET[["seurat"]])
  # RET[["plots"]][["variable.genes"]] <- VariableGenePlot(object = RET[["seurat"]],
  #                                      x.low.cutoff = 0.1, x.high.cutoff = 5,
  #                                      y.cutoff = 0.5)
  RET[["seurat"]] <- ScaleData(RET[["seurat"]])
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
seuratPCAWrapper <- function(RET, do.jackstraw = TRUE) {
  RET[["seurat"]] <- RunPCA(object =  RET[["seurat"]], pc.genes =  RET[["seurat"]]@var.genes, do.print = FALSE,
                  pcs.compute = 50)
  if (do.jackstraw) {
  	RET[["seurat"]] <- JackStraw(RET[["seurat"]], dims = 50)
  }
  RET[["plots"]][["pc.elbow"]] <- PCElbowPlot(RET[["seurat"]], num.pc = 50)
  RET[["plots"]][["pca"]] <- PCAPlot(RET[["seurat"]], dim.1 = 1, dim.2 = 2, do.return = TRUE)
  return(RET)
}

#' Wrapper function to perform clustering and run TSNE on a Seurat object
#'
#'
#' @param RET list containing Seurat object and plots
#' @param dims.use dimensions to be used for clustering
#' @param resolution resolution to be used for clustering
#'
#' @return list containing clustered Seurat object and TSNE plots
#'
#' @examples
#' seuratClusterWrapper(RET)
#'
#' @export
seuratClusterWrapper <- function(RET, dims.use = 1:10, resolution = 0.50) {
  require(Seurat)
  RET$seurat <- FindClusters(RET$seurat, dims.use = dims.use, resolution = resolution)
  RET$seurat <- RunTSNE(RET$seurat, dims.use = dims.use)
  RET$plots[["TSNE"]] <- TSNEPlot(RET$seurat, do.return = TRUE)
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
    RET[["all.markers.full"]] <- FindAllMarkers(RET$seurat, logfc.threshold = 0.05)
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
  for (ident in levels(RET$seurat@ident)) {
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
  for (i in 1:length(levels(RET$seurat@ident))) { 
      RET$seurat <- RenameIdent(object = RET$seurat, old.ident.name = i - 1,
            new.ident.name = new.idents[i])
  }
  RET$plots[["TSNE"]] <- TSNEPlot(RET$seurat, do.return = TRUE, do.label = TRUE)
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
#' Performs canonical correction analysis (CCA)
#'
#'
#' @param all.treatments Seurat object containing cells in both conditions 
#' @param condition.1 first condition of interest
#' @param condition.2 second condition of interest
#' @return list containing Seurat objects and plots
#'
#' @examples
#' seuratCCAWrapper(all.treatments, condition.1, condition.2)
#'
#' @export
seuratCCAWrapper <- function(all.treatments, condition.1, condition.2) {
  require("Seurat")
  c1.counts <- SubsetData(all.treatments, ident.use = c(condition.1), do.clean=TRUE)@raw.data
  COND1.RET <- seuratStartFromCounts(c1.counts, project = condition.1, mito.prefix = "^mt-")

  # COND1.RET <- seuratFilterWrapper(c1.counts, project.name = condition.1, mito.prefix = "^mt-")
  COND1.RET$seurat@meta.data$treatment <- condition.1
  
  c2.counts <- SubsetData(all.treatments, ident.use = c(condition.2), do.clean=TRUE)@raw.data
  # COND2.RET <- seuratFilterWrapper(c2.counts, project.name = condition.2, mito.prefix = "^mt-")
  COND2.RET <- seuratStartFromCounts(c2.counts, project = condition.2, mito.prefix = "^mt-")

  COND2.RET$seurat@meta.data$treatment <- condition.2
  
  COND1.genes <- head(rownames(COND1.RET$seurat@hvg.info), 2000)
  COND2.genes <- head(rownames(COND2.RET$seurat@hvg.info), 2000)
  genes.use <- unique(c(COND1.genes, COND2.genes))
  genes.use <- intersect(genes.use, rownames(COND1.RET$seurat@scale.data))
  genes.use <- intersect(genes.use, rownames(COND2.RET$seurat@scale.data))
  
  C1.C2.combined <- RunCCA(COND1.RET$seurat, COND2.RET$seurat, genes.use = genes.use, num.cc = 30)
  RET <- list()
  RET$cond.1 <- condition.1
  RET$cond.2 <- condition.2
  PLOTS <- list()
  
  PLOTS[["p1"]] <- DimPlot(object = C1.C2.combined, reduction.use = "cca", group.by = "treatment", 
                pt.size = 0.5, do.return = TRUE)
  
  PLOTS[["p2"]] <- VlnPlot(object = C1.C2.combined, features.plot = "CC1", group.by = "treatment",
                do.return = TRUE)
  
  PLOTS[["p3"]] <- MetageneBicorPlot(C1.C2.combined, grouping.var = "treatment", dims.eval = 1:30,
                          display.progress = TRUE)
  
  RET[["plots"]] <- PLOTS
  RET[["CCA"]] <- C1.C2.combined
  RET <- seuratAlignClusterWrapper(RET)
  # saveRDS(RET, paste(condition.1, condition.2, "CCA", "RDS", sep = "."))
  return(RET)
}
#' Aligns clusters between two conditions
#'
#'
#' @param RET list containing Seurat objects for cond1 and cond2 and plots
#' @param dims.align dimensions to be used for aligment and clustering
#' @param resolution resolution to be used for clustering
#' @return list containing aligned Seurat object and TSNE plot of the CCA
#'
#' @examples
#' seuratAlignClusterWrapper(RET)
#'
#' @export
seuratAlignClusterWrapper <- function(RET, dims.align = 1:20, resolution = 0.6) {
  require("Seurat")
  RET$CCA <- AlignSubspace(RET$CCA, reduction.type = "cca", grouping.var = "treatment", 
                                  dims.align = dims.align)
  
  RET$plots[["p4"]] <- VlnPlot(object = RET$CCA,
                               features.plot = "ACC1",
                               group.by = "treatment",
                               do.return = TRUE)
  RET$plots[["p5"]] <- VlnPlot(object = RET$CCA,
                               features.plot = "ACC2", 
                               group.by = "treatment",
                           do.return = TRUE)
  RET$CCA <- RunTSNE(RET$CCA, reduction.use = "cca.aligned", 
                     dims.use = dims.align, do.fast = T)
  RET$CCA <- FindClusters(RET$CCA, reduction.type = "cca.aligned", 
                              resolution = resolution, dims.use = dims.align)
  return(RET)
}
#' Finds conserved markers for all clusters
#'
#'
#' @param RET list containing aligned Seurat object
#' 
#' @return list containing aligned Seurat object and TSNE plot of the CCA
#'
#' @examples
#' seuratAlignClusterWrapper(RET)
#'
#' @export
seuratConservedMarkerWrapper <- function(RET) {
  RET[["conserved.markers"]] <- list()
  for (i in levels(RET$CCA@ident)){
    print(i)
    RET$conserved.markers[[i]] <- FindConservedMarkers(RET$CCA, ident.1 = i, grouping.var = "treatment")
  }
  return(RET)
}

#' Finds DE genes in each cluster between conditions 
#'
#'
#' @param RET list containing aligned Seurat object
#' 
#' @return list containing aligned Seurat object and TSNE plot of the CCA
#'
#' @examples
#' seuratAlignClusterWrapper(RET)
#'
#' @export
seuratMarkersBetweenConditions <- function(RET) {
  require("Seurat")
  require("tibble")
  require("dplyr")
  RET$CCA@meta.data$celltype.treatment <- paste0(RET$CCA@ident,
                                                 "_",
                                                 RET$CCA@meta.data$treatment)
  RET$CCA <- StashIdent(RET$CCA, save.name = "celltype")
  RET$CCA <- SetAllIdent(RET$CCA, id = "celltype.treatment")
  
  library(tibble)
  RET$response_markers <- NULL
  for (ident in levels(as.factor(RET$CCA@meta.data$celltype))) {
    print(ident)
    RET$response_markers <- rbind(RET$response_markers,
                         rownames_to_column(FindMarkers(RET$CCA,
                                                        ident.1 = paste0(ident, "_", RET$cond.2),
                                                        ident.2 = paste0(ident, "_", RET$cond.1)),
                                            var="gene") %>% mutate(cluster = ident))
  }
  RET$CCA <- SetAllIdent(RET$CCA, id = "celltype")
  return(RET)
}

#' Makes TSNE plot after CCA 
#'
#'
#' @param RET list containing aligned Seurat object
#' 
#' @return list containing Seurat object and TSNE plots (tsne.treatment and tsne.ident) of the CCA
#'
#' @examples
#' seuratTSNEPlotCCA(RET)
#'
#' @export
seuratTSNEPlotCCA <- function(RET) {
  RET$plots$tsne.treatment <- TSNEPlot(RET$CCA, do.return = T, pt.size = 0.5, group.by = "treatment")
  RET$plots$tsne.ident <- TSNEPlot(RET$CCA, do.label = T, do.return = T, pt.size = 0.5)
  return(RET)
  # plot_grid(p1, p2)
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
  RET$plots$markers <- list()
  for (marked.group in names(marker.lists)) {
    RET$plots$markers[[marked.group]] <- list()
    RET$plots$markers[[marked.group]]$heatmap <- DoHeatmap(SubsetData(RET$seurat, max.cells.per.ident = 500), 
                                              genes.use = marker.lists[[marked.group]], 
                                              slim.col.label = TRUE, group.label.rot = TRUE, 
                                              do.plot=FALSE) + ggtitle(marked.group)
    RET$plots$markers[[marked.group]]$dotplot <- DotPlot(RET$seurat, genes.plot = marker.lists[[marked.group]],
                                                    do.return=TRUE) + ggtitle(marked.group)
  }
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
    print(plot_grid(RET$plots$markers[[name]]$heatmap, RET$plots$markers[[name]]$dotplot, ncol=1))
  }
}

