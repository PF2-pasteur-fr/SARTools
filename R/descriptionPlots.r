#' Description plots of the counts
#'
#' Description plots of the counts according to the biological condition
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors for the plots (one per biological condition)
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return PNG files in the "figures" directory and the matrix of the most expressed sequences
#' @author Hugo Varet

descriptionPlots <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), ggplot_theme=theme_gray()){
  # create the figures directory if does not exist
  if (!I("figures" %in% dir())) dir.create("figures", showWarnings=FALSE)
  
  # total number of reads per sample
  barplotTotal(counts=counts, group=group, col=col, ggplot_theme=ggplot_theme)
  
  # percentage of null counts per sample
  barplotNull(counts=counts, group=group, col=col, ggplot_theme=ggplot_theme)
  
  # distribution of counts per sample
  densityPlot(counts=counts, group=group, col=col, ggplot_theme=ggplot_theme)
  
  # features which catch the most important number of reads
  majSequences <- majSequences(counts=counts, group=group, col=col, ggplot_theme=ggplot_theme)
  
  # SERE and pairwise scatter plots
  cat("Matrix of SERE statistics:\n")
  print(tabSERE(counts))
  pairwiseScatterPlots(counts=counts, ggplot_theme=ggplot_theme)

  return(majSequences)
}
