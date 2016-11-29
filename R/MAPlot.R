#' MA-plots
#'
#' MA-plot for each comparison: log2(FC) vs mean of normalized counts with one dot per feature (red dot for a differentially expressed feature, black dot otherwise)
#'
#' @param complete A \code{list} of \code{data.frame} containing features results (from \code{exportResults.DESeq2()} or \code{exportResults.edgeR()})
#' @param alpha cut-off to apply on each adjusted p-value
#' @param outfile TRUE to export the figure in a png file
#' @return A file named MAPlot.png in the figures directory containing one MA-plot per comparison
#' @author Marie-Agnes Dillies and Hugo Varet

MAPlot <- function(complete, alpha=0.05, outfile=TRUE){
  ncol <- ifelse(length(complete)<=4, ceiling(sqrt(length(complete))), 3)
  nrow <- ceiling(length(complete)/ncol)
  if (outfile) png(filename="figures/MAPlot.png", width=1800*ncol, height=1800*nrow, res=300)
    par(mfrow=c(nrow,ncol))
      for (name in names(complete)){
        complete.name <- complete[[name]]
        complete.name <- complete.name[complete.name$baseMean>0,]
        complete.name$padj <- ifelse(is.na(complete.name$padj),1,complete.name$padj)
        log2FC <- complete.name$log2FoldChange
        ylim <- 1.1 * c(-1,1) * quantile(abs(log2FC[is.finite(log2FC)]), probs=0.99)
        plot(complete.name$baseMean, pmax(ylim[1], pmin(ylim[2], log2FC)),
	         log = "x", cex=0.45, las = 1, ylim = ylim,
             col = ifelse(complete.name[,"padj"] < alpha, "red", "black"),
             pch = ifelse(log2FC<ylim[1], 6, ifelse(log2FC>ylim[2], 2, 20)),
		     xlab = "Mean of normalized counts",
		     ylab = expression(log[2]~fold~change),
  		     main = paste0("MA-plot - ",gsub("_"," ",name)))
        abline(h=0, lwd=1, col="lightgray")
      }
  if (outfile) dev.off()
}


