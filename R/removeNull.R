#' Remove features with null counts in all samples
#'
#' Remove features with null counts in all samples. These features do not contain any information and will not be used for the statistical analysis.
#'
#' @param counts \code{matrix} of raw counts
#' @return The \code{matrix} of counts without features with only null counts
#' @author Marie-Agnes Dillies and Hugo Varet

removeNull <- function(counts){
  return(counts[rowSums(counts) > 0,])
}
