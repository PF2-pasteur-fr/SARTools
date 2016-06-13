#' Density plot of all samples
#'
#' Estimation the counts density for each sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the curves (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named densplot.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

densityPlot <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  if (outfile) png(filename="figures/densplot.png",width=1800,height=1800,res=300)
    counts <- removeNull(counts)
    plot(density(log2(counts[,1]+1)), las = 1, lwd = 2,
         main = "Density of counts distribution",
	     xlab = expression(log[2] ~ (raw ~ count + 1)),
	     ylim = c(0,max(apply(counts,2,function(x){max(density(log2(x+1))$y)}))*1.05),
         col = col[as.integer(group)[1]])
    for (i in 2:ncol(counts)){
      lines(density(log2(counts[,i]+1)),col=col[as.integer(group)[i]],lwd=2)
    }
  legend("topright", levels(group), lty=1, col=col[1:nlevels(group)], lwd=2, bty="n")
  if (outfile) dev.off()
}
