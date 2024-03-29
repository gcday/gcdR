% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_wrappers.R
\name{seuratFilterWrapper}
\alias{seuratFilterWrapper}
\title{Wrapper function to perform initial filtering and normalization of a Seurat object}
\usage{
seuratFilterWrapper(
  SRT,
  min.genes = 200,
  max.genes = 5000,
  min.UMI = 0,
  max.UMI = 30000,
  max.mito = 0.25,
  mito.prefix = "^MT-",
  erythro.genes = NULL,
  max.erythro = 0.01
)
}
\arguments{
\item{SRT}{Seurat object}

\item{min.genes}{param for filtering}

\item{max.genes}{param for filtering}

\item{max.UMI}{param for filtering}

\item{max.mito}{param for filtering}

\item{mito.prefix}{prefix for mitochondrial genes ("^MT-" for hg19)}
}
\value{
list containing plots and filtered Seurat object
}
\description{
Wrapper function to perform initial filtering and normalization of a Seurat object
}
\examples{
seuratFilterWrapper(SRT, min.genes = min.genes, max.genes = max.genes, max.UMI = max.UMI,
   max.mito = max.mito, mito.prefix = mito.prefix)

}
