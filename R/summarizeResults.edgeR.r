#' Summarize edgeR analysis
#'
#' Summarize edgeR analysis: diagnotic plots, dispersions plot, summary of the independent filtering, export results...
#'
#' @param out.edgeR the result of \code{run.edgeR()}
#' @param group factor vector of the condition from which each sample belongs
#' @param counts matrix of raw counts
#' @param alpha significance threshold to apply to the adjusted p-values
#' @param col colors for the plots
#' @param log2FClim numeric vector containing both upper and lower y-axis limits for all the MA-plots produced (NULL by default to set them automatically)
#' @param padjlim numeric value between 0 and 1 for the adjusted p-value upper limits for all the volcano plots produced (NULL by default to set them automatically)
#' @param ggplot_theme ggplot2 theme function (\code{theme_light()} by default)
#' @return A list containing: (i) a list of \code{data.frames} from \code{exportResults.edgeR()} and (ii) a table summarizing the number of differentially expressed features
#' @author Hugo Varet

summarizeResults.edgeR <- function(out.edgeR, group, counts, alpha=0.05,
                                   col=c("lightblue","orange","MediumVioletRed","SpringGreen"),
                                   log2FClim=NULL, padjlim=NULL, ggplot_theme=theme_light()){  
  # create the figures/tables directory if does not exist
  if (!I("figures" %in% dir())) dir.create("figures", showWarnings=FALSE)
  if (!I("tables" %in% dir())) dir.create("tables", showWarnings=FALSE)
  
  # boxplots before and after normalisation
  countsBoxplots(out.edgeR$dge, group=group, col=col, ggplot_theme=ggplot_theme)
  
  # dispersions
  BCVPlot(dge=out.edgeR$dge, ggplot_theme=ggplot_theme)
  
  # exporting results of the differential analysis
  complete <- exportResults.edgeR(out.edgeR=out.edgeR, group=group, counts=counts, alpha=alpha)
  
  # small table with number of differentially expressed features
  nDiffTotal <- nDiffTotal(complete=complete, alpha=alpha)
  cat("Number of features down/up and total:\n")
  print(nDiffTotal, quote=FALSE)
  
  # histograms of raw p-values
  rawpHist(complete=complete, ggplot_theme=ggplot_theme)
  
  # MA-plots
  MAPlot(complete=complete, alpha=alpha, log2FClim=log2FClim, ggplot_theme=ggplot_theme)
  
  # Volcano plots
  volcanoPlot(complete=complete, alpha=alpha, padjlim=padjlim, ggplot_theme=ggplot_theme)
  
  return(list(complete=complete, nDiffTotal=nDiffTotal))
}
