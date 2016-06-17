#' Scatter plots for pairwise comparaisons of log counts
#'
#' Scatter plots for pairwise comparaisons of log counts
#'
#' @param counts \code{matrix} of raw counts
#' @param group factor vector of the condition from which each sample belongs
#' @param outfile TRUE to export the figure in a png file
#' @return A file named pairwiseScatter.png in the figures directory containing a pairwise scatter plot with the SERE statistics in the lower panel
#' @author Marie-Agnes Dillies and Hugo Varet

pairwiseScatterPlots <- function(counts, group, outfile=TRUE){
  ncol <- ncol(counts)
  if (outfile) png(filename="figures/pairwiseScatter.png",width=700*ncol,height=700*ncol,res=300)
    # defining panel and lower.panel functions
	panel <- function(x,y,...){points(x, y, pch=".");abline(a=0,b=1,lty=2);}
	lower.panel <- function(x,y,...){
	  horizontal <- (par("usr")[1] + par("usr")[2]) / 2;
      vertical <- (par("usr")[3] + par("usr")[4]) / 2; 
	  text(horizontal, vertical, round(SERE(2^cbind(x,y) - 1), digits=2), cex=ncol/2.5)
	}
	# use of the paris function
    pairs(log2(counts+1), panel=panel, lower.panel=lower.panel,
          las=1, labels=paste(colnames(counts),group,sep="\n"),
          main="Pairwise scatter plot",cex.labels=ncol/2,cex.main=ncol/4)
  if (outfile) dev.off()
}
