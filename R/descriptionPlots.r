#' Description plots of the counts
#'
#' Description plots of the counts according to the biological condition
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors for the plots (one per biological condition)
#' @return PNG files in the "figures" directory and the matrix of the most expressed sequences
#' @author Hugo Varet

descriptionPlots <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen")){
  # create the figures directory if does not exist
  if (!I("figures" %in% dir())) dir.create("figures", showWarnings=FALSE)
  
  # total number of reads per sample
  barplotTotal(counts=counts, group=group, col=col)
  
  # percentage of null counts per sample
  barplotNull(counts=counts, group=group, col=col)
  
  # distribution of counts per sample
  densityPlot(counts=counts, group=group, col=col)
  
  # features which catch the most important number of reads
  majSequences <- majSequences(counts=counts, group=group, col=col)
  
  # SERE and pairwise scatter plots
  cat("Matrix of SERE statistics:\n")
  print(tabSERE(counts))
  pairwiseScatterPlots(counts=counts, group=group)

  return(majSequences)
}
