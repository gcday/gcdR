% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat_wrappers.R
\name{prettyPrintPathways}
\alias{prettyPrintPathways}
\title{Creates table of pathways enriched within cells belonging to ident}
\usage{
prettyPrintPathways(pathways, ident)
}
\arguments{
\item{pathways}{data frame returned by \code{\link[fgsea]{fgsea}}}

\item{ident}{cluster/group of interest}
}
\value{
Well-formatted data frame
}
\description{
Creates table of pathways enriched within cells belonging to ident
}
\examples{
prettyPrintPathways()

}
