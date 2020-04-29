#' Density plot of all samples
#'
#' Estimation the counts density for each sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the curves (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return A file named densplot.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

densityPlot <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE, ggplot_theme=theme_gray()){
  if (outfile) png(filename="figures/densplot.png", width=2000, height=1800, res=300)
    counts <- removeNull(counts)
    d <- stack(data.frame(counts))
    d$group <- rep(group, each=nrow(counts))
    print(ggplot(d, aes(x=.data$values+1)) +
            stat_density(aes(group=.data$ind, color=.data$group), position="identity", geom="line", show.legend=TRUE) +
            scale_x_continuous(trans = log10_trans(),
                               breaks = trans_breaks("log10", function(x) 10^x),
                               labels = trans_format("log10", math_format(~10^.x))) +
            labs(color="") +
            scale_colour_manual(values=col) +
            xlab("Raw counts") +
            ylab("Density") +
            ggtitle("Density of counts distribution") +
            ggplot_theme)
  if (outfile) dev.off()
}
