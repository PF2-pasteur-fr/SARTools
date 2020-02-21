#' Most expressed sequences per sample
#'
#' Proportion of reads associated with the three most expressed sequences per sample
#'
#' @param counts \code{matrix} of counts
#' @param n number of most expressed sequences to return
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @param ggplot_theme ggplot2 theme function (\code{theme_gray()} by default)
#' @return A \code{matrix} with the percentage of reads of the three most expressed sequences and a file named majSeq.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

majSequences <- function(counts, n=3, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE, ggplot_theme=theme_gray()){

  seqnames <- apply(counts, 2, function(x){x <- sort(x, decreasing=TRUE); names(x)[1:n]})
  seqnames <- unique(unlist(as.character(seqnames)))

  sum <- apply(counts,2,sum)
  counts <- counts[seqnames,]
  sum <- matrix(sum,nrow(counts),ncol(counts),byrow=TRUE)
  p <- round(100*counts/sum,digits=3)

  if (outfile) png(filename="figures/majSeq.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
    maj <- apply(p, 2, max)
    seqname <- rownames(p)[apply(p, 2, which.max)]
    d <- data.frame(maj=maj, sample=factor(names(maj), levels=names(maj)), group, seqname=seqname)
    print(ggplot(d, aes(x=.data$sample, y=.data$maj, fill=.data$group)) +
            geom_bar(stat="identity", show.legend=TRUE) +
            labs(fill="") +
            scale_fill_manual(values=col) +
            xlab("Samples") + 
            ylab("Percentage of reads") +
            ggtitle("Percentage of reads from most expressed sequence") +
            ggplot_theme +
            theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)) +
            scale_y_continuous(expand=expand_scale(mult=c(0.01, 0.05))) +
            geom_text(aes(y=0.8*maj, label=seqname), color="black", size=2.5, angle=90, fontface="bold"))
  if (outfile) dev.off()
  
  return(invisible(p))
}
