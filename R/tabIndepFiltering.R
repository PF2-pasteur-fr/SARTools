#' Table of the number of features discarded by the independent filtering (if use of DESeq2)
#'
#' Compute the number of features discarded by the independent filtering for each comparison (if use of DESeq2)
#'
#' @param results list of results of \code{results(dds,...)} with chosen parameters
#' @return A \code{matrix} with the threshold and the number of features discarded for each comparison
#' @author Marie-Agnes Dillies and Hugo Varet

tabIndepFiltering <- function(results){
  out <- matrix(NA,ncol=3,nrow=length(names(results)),dimnames=list(names(results),c("Test vs Ref","BaseMean Threshold","# discarded")))
  for (name in names(results)){
	  threshold <- metadata(results[[name]])$filterThreshold
	  out[name,2] <- round(threshold,2)
    use <- results[[name]]$baseMean > threshold
	  out[name,3] <- ifelse(is.na(table(use)["FALSE"]),0,table(use)["FALSE"])
  }
  out[,1] <- gsub("_"," ",rownames(out))
  rownames(out) <- NULL
  return(out)
}
