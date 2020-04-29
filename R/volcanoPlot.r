#' Volcano plots
#'
#' Volcano plot for each comparison: -log10(adjusted P value) vs log2(FC) with one dot per feature (red dot for a differentially expressed feature, black dot otherwise)
#'
#' @param complete A \code{list} of \code{data.frame} containing features results (from \code{exportResults.DESeq2()} or \code{exportResults.edgeR()})
#' @param alpha cut-off to apply on each adjusted p-value
#' @param outfile TRUE to export the figure in a png file
#' @param padjlim numeric value between 0 and 1 for the adjusted p-value upper limits for all the volcano plots produced (NULL by default to set them automatically)
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return A file named volcanoPlot.png in the figures directory containing one volcano plot per comparison
#' @author Hugo Varet

volcanoPlot <- function(complete, alpha=0.05, outfile=TRUE, padjlim=NULL, ggplot_theme=theme_gray()){
  ncol <- min(2, length(complete))
  nrow <- ceiling(length(complete)/ncol)
  if (outfile) png(filename="figures/volcanoPlot.png", width=cairoSizeWrapper(1800*ncol), height=cairoSizeWrapper(1800*nrow), res=300)
  p <- list()
  for (name in names(complete)){
    complete.name <- complete[[name]]
    complete.name$padj[which(complete.name$padj==0)] <- .Machine$double.xmin
    complete.name <- complete.name[which(!is.na(complete.name$padj)),]
    complete.name$DE <- factor(ifelse(complete.name$padj <= alpha, "yes", "no"), levels=c("no", "yes"))
    if (is.null(padjlim)) padjlim.name <- quantile(complete.name$padj, probs=0.01, na.rm=TRUE) else padjlim.name <- padjlim
    complete.name$outfield <- factor(ifelse(complete.name$padj < padjlim.name, "top", "in"), levels=c("in", "top"))
    complete.name$padj[which(complete.name$padj < padjlim.name)] <- padjlim.name
    reverselog_trans <- function(base = exp(1)) {
      trans <- function(x) -log(x, base)
      inv <- function(x) base^(-x)
      trans_new(paste0("reverselog-", format(base)), trans, inv,
                log_breaks(base = base),
                domain = c(.Machine$double.xmin, Inf))
    }
    p[[name]] <- ggplot(data=complete.name, 
                        aes(x=.data$log2FoldChange, y=.data$padj, color=.data$DE, shape=.data$outfield)) +
      geom_point(show.legend=FALSE, alpha=0.5) +
      scale_y_continuous(trans = reverselog_trans(10),
                         breaks = trans_breaks("log10", function(x) 10^x),
                         labels = trans_format("log10", math_format(~10^.x))) +
      scale_colour_manual(values=c("no"="black", "yes"="red"), drop=FALSE) +
      scale_shape_manual(values=c("in"=16, "top"=17), drop=FALSE) +
      xlab(expression(log[2]~fold~change)) +
      ylab("Adjusted P-value") +
      ggtitle(paste0("Volcano plot - ", gsub("_", " ", name))) +
      ggplot_theme
  }
  tmpfun <- function(...) grid.arrange(..., nrow=nrow, ncol=ncol)
  do.call(tmpfun, p)
  if (outfile) dev.off()
}
