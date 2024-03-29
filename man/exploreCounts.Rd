% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exploreCounts.R
\name{exploreCounts}
\alias{exploreCounts}
\title{Explore counts structure}
\usage{
exploreCounts(
  object,
  group,
  typeTrans = "VST",
  gene.selection = "pairwise",
  col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen"),
  ggplot_theme = theme_light()
)
}
\arguments{
\item{object}{a \code{DESeqDataSet} from DESeq2 or \code{DGEList} object from edgeR}

\item{group}{factor vector of the condition from which each sample belongs}

\item{typeTrans}{transformation method for PCA/clustering with DESeq2: \code{"VST"} or \code{"rlog"}}

\item{gene.selection}{selection of the features in MDSPlot (\code{"pairwise"} by default)}

\item{col}{colors used for the PCA/MDS (one per biological condition)}

\item{ggplot_theme}{ggplot2 theme function (\code{theme_light()} by default)}
}
\value{
A list containing the dds object and the results object
}
\description{
Explore counts structure: PCA (DESeq2) or MDS (edgeR) and clustering
}
\author{
Hugo Varet
}
