% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/object_organization.R
\name{renameAllIdents}
\alias{renameAllIdents}
\title{Renames idents of cells in a Seurat object}
\usage{
renameAllIdents(RET, new.idents)
}
\arguments{
\item{RET}{list containing Seurat object and plots}

\item{new.idents}{list containing new idents (in same order as current)}
}
\value{
list containing renamed Seurat object and re-plotted TSNE
}
\description{
Renames idents of cells in a Seurat object
}
\examples{
renameIdents(RET)

}
