% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarizeResults.edgeR.r
\name{summarizeResults.edgeR}
\alias{summarizeResults.edgeR}
\title{Summarize edgeR analysis}
\usage{
summarizeResults.edgeR(
  out.edgeR,
  group,
  counts,
  alpha = 0.05,
  col = c("lightblue", "orange", "MediumVioletRed", "SpringGreen"),
  log2FClim = NULL,
  padjlim = NULL,
  ggplot_theme = theme_light()
)
}
\arguments{
\item{out.edgeR}{the result of \code{run.edgeR()}}

\item{group}{factor vector of the condition from which each sample belongs}

\item{counts}{matrix of raw counts}

\item{alpha}{significance threshold to apply to the adjusted p-values}

\item{col}{colors for the plots}

\item{log2FClim}{numeric vector containing both upper and lower y-axis limits for all the MA-plots produced (NULL by default to set them automatically)}

\item{padjlim}{numeric value between 0 and 1 for the adjusted p-value upper limits for all the volcano plots produced (NULL by default to set them automatically)}

\item{ggplot_theme}{ggplot2 theme function (\code{theme_light()} by default)}
}
\value{
A list containing: (i) a list of \code{data.frames} from \code{exportResults.edgeR()} and (ii) a table summarizing the number of differentially expressed features
}
\description{
Summarize edgeR analysis: diagnotic plots, dispersions plot, summary of the independent filtering, export results...
}
\author{
Hugo Varet
}
