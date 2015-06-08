#' Total number of reads per sample
#'
#' Bar plot of the total number of reads per sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one color per biological condition)
#' @return A file named barplotTotal.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

barplotTotal <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen")){
  png(filename="figures/barplotTotal.png",width=min(800,400+200*ncol(counts)/10),height=400)
  barplot(colSums(counts),
          main = "Total read count per sample",
		  ylab = "Total read count",
		  ylim = c(0, max(colSums(counts))*1.2),
		  col = col[as.integer(group)],
		  las = 2)
  legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
  dev.off()
}
