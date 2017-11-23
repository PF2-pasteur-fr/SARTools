#' MDS plot (for edgeR objects)
#'
#' Multi-Dimensional Scaling plot of samples based on the 500 most variant features (for edgeR analyses)
#'
#' @param dge a \code{DGEList} object
#' @param group vector of the condition from which each sample belongs
#' @param n number of features to keep among the most variant
#' @param gene.selection \code{"pairwise"} to choose the top features separately for each pairwise comparison between the samples or \code{"common"} to select the same features for all comparisons. Only used when \code{method="logFC"}
#' @param col colors to use (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named MDS.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

MDSPlot <- function(dge, group, n=min(500,nrow(dge$counts)), gene.selection=c("pairwise", "common"),
                    col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  if (outfile) png(filename="figures/MDS.png", width=1800, height=1800, res=300)
    coord <- plotMDS(dge, top=n, method="logFC", gene.selection=gene.selection[1], plot=FALSE)
    abs=range(coord$x); abs=abs(abs[2]-abs[1])/25;
    ord=range(coord$y); ord=abs(ord[2]-ord[1])/25;
    plot(coord$x,coord$y, col=col[as.integer(group)], las=1, main="Multi-Dimensional Scaling plot",
         xlab="Leading logFC dimension 1", ylab="Leading logFC dimension 2", cex=2, pch=16)
    abline(h=0,v=0,lty=2,col="lightgray")
    text(coord$x - ifelse(coord$x>0,abs,-abs), coord$y - ifelse(coord$y>0,ord,-ord),
         colnames(dge$counts), col=col[as.integer(group)])
  if (outfile) dev.off()      
}
