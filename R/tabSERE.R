#' SERE statistics for several samples
#'
#' Compute the SERE statistic for each pair of samples
#'
#' @param counts \code{matrix} of raw counts
#' @return The \code{matrix} of SERE values
#' @author Marie-Agnes Dillies and Hugo Varet

tabSERE <- function(counts){
  sere <- matrix(NA, ncol=ncol(counts), nrow=ncol(counts))
  for (i in 1:ncol(counts)){
    for (j in 1:ncol(counts)){
      sere[i,j] <- SERE(counts[,c(i,j)])
    }
  }
  colnames(sere) <- rownames(sere) <- colnames(counts)
  return(invisible(round(sere, digits=3)))
}
