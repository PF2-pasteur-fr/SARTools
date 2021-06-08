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
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return A file named MDS.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

MDSPlot <- function(dge, group, n=min(500, nrow(dge$counts)), gene.selection=c("pairwise", "common"),
                    col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE, ggplot_theme=theme_gray()){
  if (outfile) png(filename="figures/MDS.png", width=1800, height=1800, res=300)
    coord <- plotMDS(dge, top=n, method="logFC", gene.selection=gene.selection[1], plot=FALSE)
    d <- data.frame(x=coord$x, y=coord$y, group=group, 
                    sample=factor(colnames(coord$distance.matrix.squared), 
                           levels = colnames(coord$distance.matrix.squared)))
    print(ggplot(data=d, aes(x=.data$x, y=.data$y, color=group, label=sample)) + 
      geom_point(show.legend=TRUE, size=3) +
      labs(color="") +
      scale_colour_manual(values=col) +
      geom_text_repel(show.legend=FALSE, size=5, point.padding=0.2) +
      xlab("Leading logFC dimension 1") +
      ylab("Leading logFC dimension 2") +
      ggtitle("Multi-Dimensional Scaling plot") +
      ggplot_theme)
  if (outfile) dev.off()      
}
