#' Pairwise SERE for two samples
#'
#' Compute the SERE coefficient for two samples
#'
#' @param observed \code{matrix} with two columns containing observed counts of two samples
#' @return The SERE coefficient for the two samples
#' @references Schulze, Kanwar, Golzenleuchter et al, SERE: Single-parameter quality control and sample comparison for RNA-Seq, BMC Genomics, 2012
#' @author See paper published

SERE <- function(observed){
  #calculate lambda and expected values
  laneTotals <- colSums(observed)
  total <- sum(laneTotals)
  fullObserved <- observed[rowSums(observed)>0,]
  fullLambda <- rowSums(fullObserved)/total
  fullLhat <- fullLambda > 0
  fullExpected<- outer(fullLambda, laneTotals)

  #keep values
  fullKeep <- which(fullExpected > 0)

  #calculate degrees of freedom (nrow*(ncol -1) >> number of parameters - calculated (just lamda is calculated >> thats why minus 1)
  #calculate pearson and deviance for all values
  oeFull <- (fullObserved[fullKeep] - fullExpected[fullKeep])^2/ fullExpected[fullKeep] # pearson chisq test
  dfFull <- length(fullKeep) - sum(fullLhat!=0)

  sqrt(sum(oeFull)/dfFull)
}