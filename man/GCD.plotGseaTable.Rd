% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fgsea_funcs.R
\name{GCD.plotGseaTable}
\alias{GCD.plotGseaTable}
\title{Plots enriched and depleted pathways for a cluster}
\usage{
GCD.plotGseaTable(
  pathways,
  stats,
  fgseaRes,
  gseaParam = 1,
  term.prefix,
  num.to.print,
  colwidths = c(10, 4, 0.8, 1.6)
)
}
\arguments{
\item{pathways}{list of pathways}

\item{stats}{fgsea ranks}

\item{fgseaRes}{fgsea results}

\item{gseaParam}{param for fgsea algorithm}

\item{term.prefix}{prefix to remove from pathway name for printing}

\item{num.to.print}{print at most this many pathways}

\item{colwidths}{ex. c(10, 3, 0.8, 1.2)}

\item{pval.thresh}{only include pathways with padj below this value}
}
\value{
list containing ranks and fgsea output for each ident
}
\description{
Plots enriched and depleted pathways for a cluster
}
\examples{
GCD.plotGseaTable(all.markers, pathways)

}
