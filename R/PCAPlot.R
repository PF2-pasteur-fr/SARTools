#' PCA of samples (if use of DESeq2)
#'
#' Principal Component Analysis of samples based on the 500 most variant features on VST- or rlog-counts (if use of DESeq2)
#'
#' @param counts.trans a matrix a transformed counts (VST- or rlog-counts)
#' @param group factor vector of the condition from which each sample belongs
#' @param n number of features to keep among the most variant
#' @param col colors to use (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return A file named PCA.png in the figures directory with a pairwise plot of the three first principal components
#' @author Marie-Agnes Dillies and Hugo Varet

PCAPlot <- function(counts.trans, group, n=min(500, nrow(counts.trans)), 
                    col=c("lightblue","orange","MediumVioletRed","SpringGreen"),
                    outfile=TRUE, ggplot_theme=theme_gray()){
  # PCA on the 500 most variables features
  rv = apply(counts.trans, 1, var, na.rm=TRUE)
  pca = prcomp(t(counts.trans[order(rv, decreasing = TRUE), ][1:n,]))
  prp <- pca$sdev^2 * 100 / sum(pca$sdev^2)
  prp <- round(prp[1:3],2)

  # create figure
  if (outfile) png(filename="figures/PCA.png", width=cairoSizeWrapper(1900*2), height=cairoSizeWrapper(1800), res=300)
  
  tmpFunction <- function(axes=c(1, 2)){
    index1 <- axes[1]
    index2 <- axes[2]
    d <- data.frame(x=pca$x[,index1], y=pca$x[,index2], 
                    group=group, sample=factor(rownames(pca$x), levels=rownames(pca$x)))
    ggplot(data=d, aes(x=.data$x, y=.data$y, color=group, label=sample)) + 
      geom_point(show.legend=TRUE, size=3) +
      labs(color="") +
      scale_colour_manual(values=col) +
      geom_text_repel(show.legend=FALSE, size=5, point.padding=0.2) +
      xlab(paste0("PC", index1, " (",prp[index1],"%)")) +
      ylab(paste0("PC", index2, " (",prp[index2],"%)")) +
      ggplot_theme
  }
  p1 <- tmpFunction(c(1, 2))
  p2 <- tmpFunction(c(1, 3))
  grid.arrange(p1, p2, nrow=1, ncol=2, 
               top=textGrob("Principal Component Analysis", x=0.01, hjust=0, gp=gpar(fontsize=15)))
  
  if (outfile) dev.off()

  return(invisible(pca$x))
}
