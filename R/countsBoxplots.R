#' Box-plots of (normalized) counts distribution per sample
#'
#' Box-plots of raw and normalized counts distributions per sample to assess the effect of the normalization
#'
#' @param object a \code{DESeqDataSet} object from DESeq2 or a \code{DGEList} object from edgeR
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the boxplots (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A file named countsBoxplots.png in the figures directory containing boxplots of the raw and normalized counts
#' @author Marie-Agnes Dillies and Hugo Varet

countsBoxplots <- function(object, group, col = c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  if (class(object)=="DESeqDataSet"){
    counts <- counts(object)
    counts <- removeNull(counts)
    norm.counts <- counts(object, normalized=TRUE)
    norm.counts <- removeNull(norm.counts)  
  } else{
    counts <- object$counts
    counts <- removeNull(counts)
    tmm <- object$samples$norm.factors
    N <- colSums(object$counts)
    f <- tmm * N/mean(tmm * N)
    norm.counts <- scale(object$counts, center=FALSE, scale=f)
    norm.counts <- removeNull(norm.counts)    
  }

  if (outfile) png(filename="figures/countsBoxplots.png", width=2*min(2200, 1800+800*ncol(norm.counts)/10), height=1800, res=300)
  d <- stack(as.data.frame(counts))
  d$group <- rep(group, each=nrow(counts))
  p1 <- ggplot(d) + 
    geom_boxplot(aes(x=.data$ind, y=.data$values+1, fill=.data$group), show.legend=TRUE) +
    labs(fill="") +
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(~10^.x))) +
    scale_fill_manual(values=col) +
    xlab("Samples") +
    ylab("Raw counts") +
    ggtitle("Raw counts distribution") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  
  d <- stack(as.data.frame(norm.counts))
  d$group <- rep(group, each=nrow(norm.counts))
  p2 <- ggplot(d) + 
    geom_boxplot(aes(x=.data$ind, y=.data$values+1, fill=.data$group), show.legend=TRUE) +
    labs(fill="") +
    scale_y_continuous(trans = log10_trans(),
                       breaks = trans_breaks("log10", function(x) 10^x),
                       labels = trans_format("log10", math_format(~10^.x))) +
    scale_fill_manual(values=col) +
    xlab("Samples") +
    ylab("Normalized counts") +
    ggtitle("Normalized counts distribution") +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  
  grid.arrange(p1, p2, nrow=1, ncol=2)
  if (outfile) dev.off()
    
}
