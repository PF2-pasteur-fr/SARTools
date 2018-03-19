#' Assess the estimations of the size factors
#'
#' Plots to assess the estimations of the size factors
#'
#' @param dds a \code{DESeqDataSet} object
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors for the plots
#' @param outfile TRUE to export the figure in a png file
#' @param plots vector of plots to generate
#' @return Two files in the figures directory: diagSizeFactorsHist.png containing one histogram per sample and diagSizeFactorsTC.png for a plot of the size factors vs the total number of reads
#' @author Marie-Agnes Dillies and Hugo Varet

diagSizeFactorsPlots <- function(dds, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), 
                                 outfile=TRUE, plots=c("diag","sf_libsize")){
  # histograms
  if ("diag" %in% plots){
    ncol <- ifelse(ncol(counts(dds))<=4, ceiling(sqrt(ncol(counts(dds)))), 3)
    nrow <- ceiling(ncol(counts(dds))/ncol)
    if (outfile) png(filename="figures/diagSizeFactorsHist.png", width=cairoSizeWrapper(1400*ncol), height=cairoSizeWrapper(1400*nrow), res=300)
    par(mfrow=c(nrow,ncol))
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
    if (outfile) dev.off()
  }
  
  # total read counts vs size factors
  if ("sf_libsize" %in% plots){
    if (outfile) png(filename="figures/diagSizeFactorsTC.png", width=1800, height=1800, res=300)  
    sf <- sizeFactors(dds)
    libsize <- colSums(counts(dds))/1e6
    plot(sf, libsize, pch=16, las=1,
         col = col[as.integer(group)],
         xlab="Size factors", ylab="Total number of reads (millions)",
         main="Diagnostic: size factors vs total number of reads")
    abs <- range(sf); meanAbs <- mean(abs); abs <- abs(abs[2]-abs[1])/25;
    ord <- range(libsize); meanOrd <- mean(ord); ord <- abs(ord[2]-ord[1])/25;
    text(sf - ifelse(sf > meanAbs, abs, -abs), 
         libsize - ifelse(libsize > meanOrd, ord, -ord),
         colnames(dds), col=col[as.integer(group)])
    abline(lm(libsize ~ sf + 0), lty=2, col="grey")
    if (outfile) dev.off()
  }
}


