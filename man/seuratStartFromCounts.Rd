% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_wrappers.R
\name{seuratStartFromCounts}
\alias{seuratStartFromCounts}
\title{Wrapper function to initiate Seurat analysis from count matrices}
\usage{
seuratStartFromCounts(
  read.counts,
  project,
  min.genes = 200,
  max.genes = 5000,
  max.UMI = 30000,
  max.mito = 0.25,
  mito.prefix = "^MT-"
)
}
\arguments{
\item{read.counts}{sparse matrix of read counts (returned by \code{\link[Seurat]{Read10X}})}

\item{project}{name for Seurat object (required)}

\item{min.genes}{param for filtering}

\item{max.genes}{param for filtering}

\item{max.UMI}{param for filtering}

\item{max.mito}{param for filtering}

\item{mito.prefix}{param for filtering}
}
\value{
Output of seuratFilterWrapper
}
\description{
Wrapper function to initiate Seurat analysis from count matrices
}
\examples{
prettyPrintMarkers()

}
