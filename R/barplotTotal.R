#' Total number of reads per sample
#'
#' Bar plot of the total number of reads per sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one color per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return A file named barplotTotal.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

barplotTotal <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE, ggplot_theme=theme_gray()){
  if (outfile) png(filename="figures/barplotTotal.png", width=min(3600, 1800+800*ncol(counts)/10), height=1800, res=300)
  d <- data.frame(tc=colSums(counts)/1e6, sample=factor(colnames(counts), colnames(counts)), group)
  print(ggplot(d, aes(x=.data$sample, y=.data$tc, fill=.data$group)) +
          geom_bar(stat="identity", show.legend=TRUE) +
          labs(fill="") +
          scale_fill_manual(values=col) +
          xlab("Samples") + 
          ylab("Total read count (million)") +
          ggtitle("Total read count per sample (million)") +
          ggplot_theme +
          theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
          scale_y_continuous(expand=expand_scale(mult=c(0.01, 0.05))))
  if (outfile) dev.off()
}
