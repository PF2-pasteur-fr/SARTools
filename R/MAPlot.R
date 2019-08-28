#' MA-plots
#'
#' MA-plot for each comparison: log2(FC) vs mean of normalized counts with one dot per feature (red dot for a differentially expressed feature, black dot otherwise)
#'
#' @param complete A \code{list} of \code{data.frame} containing features results (from \code{exportResults.DESeq2()} or \code{exportResults.edgeR()})
#' @param alpha cut-off to apply on each adjusted p-value
#' @param outfile TRUE to export the figure in a png file
#' @param log2FClim numeric vector containing both upper and lower y-axis limits for all the MA-plots produced (NULL by default to set them automatically)
#' @return A file named MAPlot.png in the figures directory containing one MA-plot per comparison
#' @author Marie-Agnes Dillies and Hugo Varet

MAPlot <- function(complete, alpha=0.05, outfile=TRUE, log2FClim=NULL){
  ncol <- ifelse(length(complete)<=4, ceiling(sqrt(length(complete))), 3)
  nrow <- ceiling(length(complete)/ncol)
  if (outfile) png(filename="figures/MAPlot.png", width=cairoSizeWrapper(1800*ncol), height=cairoSizeWrapper(1800*nrow), res=300)
  p <- list()
  for (name in names(complete)){
    complete.name <- complete[[name]]
    complete.name <- complete.name[which(complete.name$baseMean>0),]
    complete.name$padj <- ifelse(is.na(complete.name$padj), 1, complete.name$padj)
    py <- complete.name$log2FoldChange
    if (is.null(log2FClim)) ymax <- quantile(abs(py[is.finite(py)]), probs=0.99) else ymax <- log2FClim
    complete.name$log2FoldChange[which(py > ymax)] <- ymax
    complete.name$log2FoldChange[which(py < -ymax)] <- -ymax
    complete.name$shape <- abs(py) > ymax
    p[[name]] <- ggplot(data=complete.name, aes(x=.data$baseMean, y=.data$log2FoldChange, color=.data$padj<=alpha, shape=.data$shape)) +
      scale_x_continuous(trans = log10_trans(),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(~10^.x))) +
      geom_point(show.legend=FALSE, alpha=0.5) +
      scale_colour_manual(values=c("black","red")) +
      xlab("") +
      ylab("") +
      ggtitle(gsub("_"," ",name))
  }
  tmpfun <- function(...) grid.arrange(..., nrow=nrow, ncol=ncol,
                                       top=textGrob("MA-plot", x=0.01, just="left", gp=gpar(fontsize=20)),
                                       bottom=textGrob("Mean of normalized counts", gp=gpar(fontsize=15)),
                                       left=textGrob(expression(log[2]~fold~change), rot=90, gp=gpar(fontsize=15)))
  do.call(tmpfun, p)
  if (outfile) dev.off()
}


