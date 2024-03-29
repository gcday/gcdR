% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_wrappers.R
\name{seuratClusterWrapper}
\alias{seuratClusterWrapper}
\title{Wrapper function to perform clustering and run TSNE on a Seurat object}
\usage{
seuratClusterWrapper(
  RET,
  dims = NULL,
  resolution = 0.8,
  do.TSNE = T,
  do.UMAP = T,
  do.tree = T
)
}
\arguments{
\item{RET}{list containing Seurat object and plots}

\item{dims}{dimensions to be used for clustering}

\item{resolution}{resolution to be used for clustering}
}
\value{
list containing clustered Seurat object and TSNE plots
}
\description{
Wrapper function to perform clustering and run TSNE on a Seurat object
}
\examples{
seuratClusterWrapper(RET)

}
