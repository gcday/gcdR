% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_exp.R
\name{FindAllMarkersShiny}
\alias{FindAllMarkersShiny}
\alias{FindAllMarkersNodeShiny}
\title{Gene expression markers for all identity classes}
\usage{
FindAllMarkersShiny(
  object,
  assay = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  node = NULL,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  return.thresh = 0.01,
  idents.use = NULL,
  shiny = FALSE,
  metadata.add = NULL,
  ...
)
}
\arguments{
\item{node}{A node to find markers for and all its children; requires
\code{\link{BuildClusterTree}} to have been run previously; replaces \code{FindAllMarkersNode}}

\item{return.thresh}{Only return markers that have a p-value < return.thresh, or a power > return.thresh (if the test is ROC)}
}
\value{
Matrix containing a ranked list of putative markers, and associated
statistics (p-values, ROC score, etc.)
}
\description{
Modified to add progress indicator for Shiny apps.
Finds markers (differentially expressed genes) for each of the identity classes in a dataset
}
\examples{
# Find markers for all clusters
all.markers <- FindAllMarkers(object = pbmc_small)
head(x = all.markers)

# Pass a value to node as a replacement for FindAllMarkersNode
pbmc_small <- BuildClusterTree(object = pbmc_small)
all.markers <- FindAllMarkers(object = pbmc_small, node = 4)
head(x = all.markers)

}
