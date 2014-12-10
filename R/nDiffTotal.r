#' Number of differentially expressed features per comparison
#'
#' Number of down- and up-regulated features per comparison
#'
#' @param complete list of \code{data.frame} containing features results (from \code{exportResults.DESeq2()} or \code{exportResults.edgeR()})
#' @param alpha threshold to apply to the FDR
#' @return A matrix with the number of up, down and total of features per comparison
#' @author Marie-Agnes Dillies and Hugo Varet

nDiffTotal <- function(complete, alpha=0.05){
  nDiffTotal <- matrix(NA,ncol=4,nrow=length(complete),dimnames=list(names(complete),c("Test vs Ref", "# down","# up","# total")))
  for (name in names(complete)){
    complete.name <- complete[[name]]
    if (!is.null(complete.name$betaConv)){
	  nDiffTotal[name,2:3]=c(nrow(complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange<=0),]),
                             nrow(complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange>=0),]))
	} else{
	  nDiffTotal[name,2:3]=c(nrow(complete.name[which(complete.name$padj <= alpha & complete.name$log2FoldChange<=0),]),
                             nrow(complete.name[which(complete.name$padj <= alpha & complete.name$log2FoldChange>=0),]))	
	}
  }
  nDiffTotal[,4] <- nDiffTotal[,2] + nDiffTotal[,3]
  nDiffTotal[,1] <- gsub("_"," ",rownames(nDiffTotal))
  rownames(nDiffTotal) <- NULL
  return(nDiffTotal)
}
