#' Assess the estimations of the size factors
#'
#' Plots to assess the estimations of the size factors
#'
#' @param dds a \code{DESeqDataSet} object
#' @return Two files in the figures directory: diagSizeFactorsHist.png containing one histogram per sample and diagSizeFactorsTC.png for a plot of the size factors vs the total number of reads
#' @author Marie-Agnes Dillies and Hugo Varet

diagSizeFactorsPlots <- function(dds){
  # histograms
  nrow <- ceiling(sqrt(ncol(counts(dds))))
  ncol <- ceiling(ncol(counts(dds))/nrow)
  png(filename="figures/diagSizeFactorsHist.png", width=300*max(ncol,nrow), height=300*min(ncol,nrow))
    par(mfrow=sort(c(nrow,ncol)))
	geomeans <- exp(rowMeans(log(counts(dds))))
    samples <- colnames(counts(dds))
    counts.trans <- log2(counts(dds)/geomeans)
    xmin <- min(counts.trans[is.finite(counts.trans)],na.rm=TRUE)
    xmax <- max(counts.trans[is.finite(counts.trans)],na.rm=TRUE)
    for (j in 1:ncol(dds)){
  	  hist(log2(counts(dds)[,j]/geomeans), nclass=100, xlab=expression(log[2] ~ (counts/geometric~mean)), las=1, xlim=c(xmin,xmax),
           main=paste0("Size factors diagnostic - Sample ",samples[j]),col="skyblue")
      abline(v = log2(sizeFactors(dds)[j]), col="red", lwd=1.5)
    }
  dev.off()
  # total read counts vs size factors
  png(filename="figures/diagSizeFactorsTC.png",width=400,height=400)  
	plot(sizeFactors(dds), colSums(counts(dds)), pch=19, las=1, xlab="Size factors",
	     ylab="Total number of reads",main="Diagnostic: size factors vs total number of reads")
	abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0), lty=2, col="grey")
  dev.off()
}


