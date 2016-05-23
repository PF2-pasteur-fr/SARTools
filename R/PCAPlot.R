#' PCA of samples (if use of DESeq2)
#'
#' Principal Component Analysis of samples based on the 500 most variant features on VST- or rlog-counts (if use of DESeq2)
#'
#' @param counts.trans a matrix a transformed counts (VST- or rlog-counts)
#' @param group factor vector of the condition from which each sample belongs
#' @param n number of features to keep among the most variant
#' @param col colors to use (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named PCA.png in the figures directory with a pairwise plot of the three first principal components
#' @author Marie-Agnes Dillies and Hugo Varet

PCAPlot <- function(counts.trans, group, n=min(500,nrow(counts.trans)), 
                    col=c("lightblue","orange","MediumVioletRed","SpringGreen"),
                    outfile=TRUE){
  # PCA on the 500 most variables features
  rv = apply(counts.trans, 1, var, na.rm=TRUE)
  pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE), ][1:n,]))
  prp <- pca$sdev^2 * 100 / sum(pca$sdev^2)
  prp <- round(prp[1:3],2)

  # create figure
  if (outfile) png(filename="figures/PCA.png",width=1800*2,height=1800,res=300)
    par(mfrow=c(1,2))
	# axes 1 et 2
	abs=range(pca$x[,1]); abs=abs(abs[2]-abs[1])/25;
    ord=range(pca$x[,2]); ord=abs(ord[2]-ord[1])/25;
    plot(pca$x[,1], pca$x[,2],
         las = 1, cex = 2, pch = 16, col = col[as.integer(group)],
	     xlab = paste0("PC1 (",prp[1],"%)"), 
	     ylab = paste0("PC2 (",prp[2],"%)"), 
	     main = "Principal Component Analysis - Axes 1 and 2")
    abline(h=0,v=0,lty=2,col="lightgray")
    text(pca$x[,1] - ifelse(pca$x[,1]>0,abs,-abs), pca$x[,2] - ifelse(pca$x[,2]>0,ord,-ord),
         colnames(counts.trans), col=col[as.integer(group)])

	# axes 1 et 3
	abs=range(pca$x[,1]); abs=abs(abs[2]-abs[1])/25;
    ord=range(pca$x[,3]); ord=abs(ord[2]-ord[1])/25;
    plot(pca$x[,1], pca$x[,3],
         las = 1, cex = 2, pch = 16, col = col[as.integer(group)],
	     xlab = paste0("PC1 (",prp[1],"%)"), 
	     ylab = paste0("PC3 (",prp[3],"%)"), 
	     main = "Principal Component Analysis - Axes 1 and 3")
    abline(h=0,v=0,lty=2,col="lightgray")
    text(pca$x[,1] - ifelse(pca$x[,1]>0,abs,-abs), pca$x[,3] - ifelse(pca$x[,3]>0,ord,-ord),
         colnames(counts.trans), col=col[as.integer(group)])
  if (outfile) dev.off()

  return(invisible(pca$x))
}
