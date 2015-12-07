#' Clustering of the samples
#'
#' Clustering of the samples based on VST- or rlog-counts (if use of DESeq2) or cpm-counts (if use of edgeR)
#'
#' @param counts.trans a matrix a transformed counts (VST- or rlog-counts if use of DESeq2 or cpm-counts if use of edgeR)
#' @param group factor vector of the condition from which each sample belongs
#' @param outfile TRUE to export the figure in a png file
#' @return A file named cluster.png in the figures directory with the dendrogram of the clustering
#' @author Marie-Agnes Dillies and Hugo Varet

clusterPlot <- function(counts.trans, group, outfile=TRUE){
  hc <- hclust(dist(t(counts.trans)), method="ward.D")
  if (outfile) png(filename="figures/cluster.png",width=400,height=400) 
    plot(hc, hang=-1, ylab="Height", las=2, xlab="Method: Euclidean distance - Ward criterion", main="Cluster dendrogram")
  if (outfile) dev.off()
}
