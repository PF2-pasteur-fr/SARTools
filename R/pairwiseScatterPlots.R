#' Scatter plots for pairwise comparaisons of log counts
#'
#' Scatter plots for pairwise comparaisons of log counts
#'
#' @param counts \code{matrix} of raw counts
#' @param outfile TRUE to export the figure in a png file
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return A file named pairwiseScatter.png in the figures directory containing a pairwise scatter plot with the SERE statistics in the lower panel
#' @author Marie-Agnes Dillies and Hugo Varet

pairwiseScatterPlots <- function(counts, outfile=TRUE, ggplot_theme=theme_gray()){
  ncol <- ncol(counts)
  if (ncol <= 12){
    if (outfile) png(filename="figures/pairwiseScatter.png", width=cairoSizeWrapper(850*ncol), height=cairoSizeWrapper(700*ncol), res=300)
    d <- data.frame(counts+1)
    p <- list()
    layout_matrix = matrix(NA, ncol=ncol, nrow=ncol)
    k <- 1
    for (i in 1:ncol){
      for (j in 1:ncol){
        if (i==j) next
        if (i > j){
          p[[k]] <- ggplot(data=cbind(d, z=1), aes_string(x=names(d)[i], y=names(d)[j], z="z")) +
            stat_summary_2d(fun=function(z) log(sum(z)), bins=60, show.legend=FALSE) +
            scale_x_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(~10^.x))) +
            scale_y_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(~10^.x))) +
            xlab(names(d)[i]) + ylab(names(d)[j]) +
            geom_abline(slope=1, intercept=0, linetype="dashed", col="lightgrey") +
            ggplot_theme
        } else{
          text = as.character(round(SERE(counts[,c(i,j)]), digits=2))
          p[[k]] <- ggplot() + annotate("text", x=4, y=25, size=10, label=text) +
            theme_bw() + xlab(names(d)[i]) + ylab(names(d)[j]) +
            theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
                  axis.text.x=element_blank(), axis.ticks.x=element_blank(),
                  axis.text.y=element_blank(), axis.ticks.y=element_blank())
        }
        layout_matrix[j,i] <- k
        k <- k+1
      }
    }
    tmpfun <- function(...) grid.arrange(..., ncol=ncol, nrow=ncol, 
                                         layout_matrix=layout_matrix, 
                                         top=textGrob("Pairwise scatter plot", x=0.01, just="left", gp=gpar(fontsize=20)))
    do.call(tmpfun, p)
    if (outfile) dev.off()
  } else{
    warning("No pairwise scatter-plot produced because of a too high number of samples (>12).")
    if (outfile) png(filename="figures/pairwiseScatter.png", width=1900, height=1900, res=300)
    par(mar=c(1.5, 1.5, 1.5, 1.5))
    plot(0, 0, bty="o", pch=".", col="white", xaxt="n", yaxt="n", xlab="", ylab="")
    text(0, 0.3, "No pairwise scatter-plot produced because of\na too high number of samples (>12).", pos=1, cex=1.5)
    if (outfile) dev.off()
  }
}
