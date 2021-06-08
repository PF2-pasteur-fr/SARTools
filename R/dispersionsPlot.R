#' Plots about DESeq2 dispersions
#'
#' A plot of the mean-dispersion relationship and a diagnostic of log normality of the dispersions (if use of DESeq2)
#'
#' @param dds a \code{DESeqDataSet} object
#' @param outfile TRUE to export the figure in a png file
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return A file named dispersionsPlot.png in the figures directory containing the plot of the mean-dispersion relationship and a diagnostic of log normality of the dispersions
#' @author Marie-Agnes Dillies and Hugo Varet

dispersionsPlot <- function(dds, outfile=TRUE, ggplot_theme=theme_gray()){
  if (outfile) png(filename="figures/dispersionsPlot.png", width=3600, height=1800, res=300)
	
  # dispersions plot
  d <- as.data.frame(mcols(dds)[,c("baseMean", "dispGeneEst", "dispFit", "dispersion")])
  d <- d[which(d$baseMean > 0),]
  d <- data.frame(baseMean=rep(d$baseMean, 3),
                  value=c(d$dispGeneEst, d$dispersion, d$dispFit),
                  variable=factor(rep(c("dispGeneEst", "dispersion", "dispFit"), each=nrow(d)),
                                  levels=c("dispGeneEst", "dispersion", "dispFit")))
  p1 <- ggplot(d, aes(x=.data$baseMean, y=.data$value, colour=.data$variable)) + 
    geom_point(size=0.1) +
    scale_x_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format())) +
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format())) +
    ylab("Dispersion") + 
    xlab("Mean of normalized counts") +
    scale_colour_manual(
      values=c("Black", "#377eb8", "#e41a1c"),
      breaks=c("dispGeneEst", "dispersion", "dispFit"),
      labels=c("Estimate", "Final", "Fit"),
      name="") +
    guides(colour = guide_legend(override.aes = list(size=2))) +
    ggtitle("Dispersions") +
    ggplot_theme
  
  # diagnostic of log normality
  disp <- mcols(dds)$dispGeneEst
  disp <- disp[!is.na(disp)]
  disp <- disp[disp>1e-8]
  d <- data.frame(disp)
  p2 <- ggplot(data=d, aes(x=.data$disp)) +
    geom_histogram(bins=80, aes(y=.data$..density..)) +
    scale_x_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format())) +
    scale_y_continuous(expand=expansion(mult=c(0.01, 0.05))) +
    xlab("Feature dispersion estimate") +
    ylab("Density") +
    ggtitle("log-normality dispersion diagnostic") +
    stat_function(fun = function(x, ...){dnorm(log10(x), ...)}, 
                  args = list(mean = mean(log10(disp), na.rm=TRUE), 
                              sd = sd(log10(disp), na.rm=TRUE))) +
    ggplot_theme
  
  grid.arrange(p1, p2, layout_matrix=matrix(c(1, 1, 1, 1, 1, 2, 2, 2, 2), nrow=1))
  
  if (outfile) dev.off()
}
