% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{makeMarkerHeatmaps}
\alias{makeMarkerHeatmaps}
\title{Makes heatmaps for Seurat object}
\usage{
makeMarkerHeatmaps(RET, marker.lists)
}
\arguments{
\item{RET}{gcdSeurat object}

\item{marker.lists}{list of markers (may be genes or module scores)}
}
\value{
list containing Seurat object and marker heatmaps
}
\description{
Makes heatmaps for Seurat object
}
\examples{
makeMarkerHeatmaps(RET, marker.lists)

}
