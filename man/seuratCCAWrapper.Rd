% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_CCA_wrappers.R
\name{seuratCCAWrapper}
\alias{seuratCCAWrapper}
\title{Performs canonical correction analysis (CCA)}
\usage{
seuratCCAWrapper(all.treatments, condition.1, condition.2)
}
\arguments{
\item{all.treatments}{Seurat object containing cells in both conditions}

\item{condition.1}{first condition of interest}

\item{condition.2}{second condition of interest}
}
\value{
list containing Seurat objects and plots
}
\description{
Performs canonical correction analysis (CCA)
}
\examples{
seuratCCAWrapper(all.treatments, condition.1, condition.2)

}
