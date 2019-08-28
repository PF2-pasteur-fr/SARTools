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
  if (outfile) png(filename="figures/cluster.png", width=1800, height=1800, res=300) 
  print(ggdendrogram(hc, theme_dendro=FALSE) +
          xlab("Samples") +
          ylab("Height") +
          ggtitle("Cluster dendrogram\nEuclidean distance, Ward criterion") +
          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
                axis.text.y=element_text(angle=0)) +
          scale_y_continuous(expand=expand_scale(mult=c(0.01, 0.05))))
  if (outfile) dev.off()
}
