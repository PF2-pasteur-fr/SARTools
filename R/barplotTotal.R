#' Total number of reads per sample
#'
#' Bar plot of the total number of reads per sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one color per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named barplotTotal.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

barplotTotal <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  if (outfile) png(filename="figures/barplotTotal.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
  libsize <- colSums(counts)/1e6
  barplot(libsize,
          main = "Total read count per sample (million)",
		  ylab = "Total read count (million)",
		  ylim = c(0, max(libsize)*1.2),
		  col = col[as.integer(group)],
		  las = 2)
  legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
  if (outfile) dev.off()
}
