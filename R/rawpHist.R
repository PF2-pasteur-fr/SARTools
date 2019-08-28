#' Histograms of raw p-values
#'
#' Histogram of raw p-values for each comparison
#'
#' @param complete a list of \code{data.frames} created by \code{summaryResults.DESeq2()} or \code{summaryResults.edgeR()}
#' @param outfile TRUE to export the figure in a png file
#' @return A file named rawpHist.png in the figures directory with one histogram of raw p-values per comparison
#' @author Marie-Agnes Dillies and Hugo Varet

rawpHist <- function(complete, outfile=TRUE){
  ncol <- ifelse(length(complete)<=4, ceiling(sqrt(length(complete))), 3)
  nrow <- ceiling(length(complete)/ncol)
  if (outfile) png(filename="figures/rawpHist.png", width=cairoSizeWrapper(1800*ncol), height=cairoSizeWrapper(1800*nrow), res=300)
  p <- list()
  for (name in names(complete)){
    complete.name <- complete[[name]]
    complete.name <- complete.name[which(!is.na(complete.name$pvalue)),]
    p[[name]] <- ggplot(data=complete.name, aes(x=.data$pvalue)) +
      geom_histogram(binwidth=0.02) +
      scale_y_continuous(expand=expand_scale(mult=c(0.01, 0.05))) +
      xlab("") +
      ylab("") +
      ggtitle(gsub("_"," ",name))
  }
  tmpfun <- function(...) grid.arrange(..., nrow=nrow, ncol=ncol,
                                       top=textGrob("Distribution of raw p-values", x=0.01, just="left", gp=gpar(fontsize=20)),
                                       bottom=textGrob("Raw p-value", gp=gpar(fontsize=15)))
  do.call(tmpfun, p)
  if (outfile) dev.off()
}
