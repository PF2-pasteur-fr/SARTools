#' Box-plots of (normalized) counts distribution per sample
#'
#' Box-plots of raw and normalized counts distributions per sample to assess the effect of the normalization
#'
#' @param object a \code{DESeqDataSet} object from DESeq2 or a \code{DGEList} object from edgeR
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the boxplots (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named countsBoxplots.png in the figures directory containing boxplots of the raw and normalized counts
#' @author Marie-Agnes Dillies and Hugo Varet

countsBoxplots <- function(object, group, col = c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  if (class(object)=="DESeqDataSet"){
    counts <- counts(object)
    counts <- removeNull(counts)
    norm.counts <- counts(object, normalized=TRUE)
    norm.counts <- removeNull(norm.counts)  
  } else{
    counts <- object$counts
    counts <- removeNull(counts)
	tmm <- object$samples$norm.factors
    N <- colSums(object$counts)
    f <- tmm * N/mean(tmm * N)
    norm.counts <- scale(object$counts, center=FALSE, scale=f)
    norm.counts <- removeNull(norm.counts)    
  }

  if (outfile) png(filename="figures/countsBoxplots.png",width=2*min(2200,1800+800*ncol(norm.counts)/10),height=1800,res=300)
    par(mfrow=c(1,2))
	# raw counts
    boxplot(log2(counts+1), col = col[as.integer(group)], las = 2,
	        main = "Raw counts distribution", ylab = expression(log[2] ~ (raw ~ count + 1)))
    legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
	# norm counts
    boxplot(log2(norm.counts+1), col = col[as.integer(group)], las = 2,
	        main = "Normalized counts distribution", ylab = expression(log[2] ~ (norm ~ count + 1)))
    legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
  if (outfile) dev.off()
}
