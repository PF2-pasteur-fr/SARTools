#' Percentage of null counts per sample
#'
#' Bar plot of the percentage of null counts per sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one color per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return A file named barplotNull.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

barplotNull <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE, ggplot_theme=theme_gray()){
  if (outfile) png(filename="figures/barplotNull.png", width=min(3600, 1800+800*ncol(counts)/10), height=1800, res=300)
    percentage <- apply(counts, 2, function(x){sum(x == 0)})*100/nrow(counts)
    percentage.allNull <- (nrow(counts) - nrow(removeNull(counts)))*100/nrow(counts)
    d <- data.frame(percentage=percentage, sample=factor(names(percentage), levels=names(percentage)), group)
    print(ggplot(d, aes(x=.data$sample, y=.data$percentage, fill=.data$group)) +
            geom_bar(stat="identity", show.legend=TRUE) +
            labs(fill="") +
            scale_fill_manual(values=col) +
            xlab("Samples") + 
            ylab("Percentage of null counts") +
            ggtitle("Percentage of null counts per sample") +
            ggplot_theme +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
            scale_y_continuous(expand=expand_scale(mult=c(0.01, 0.05))) +
            geom_hline(yintercept=percentage.allNull, linetype="dashed", color="black", size=1))
  if (outfile) dev.off()
}
