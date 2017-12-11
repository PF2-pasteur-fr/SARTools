#' Most expressed sequences per sample
#'
#' Proportion of reads associated with the three most expressed sequences per sample
#'
#' @param counts \code{matrix} of counts
#' @param n number of most expressed sequences to return
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one per biological condition)
#' @param outfile TRUE to export the figure in a png file
#' @return A \code{matrix} with the percentage of reads of the three most expressed sequences and a file named majSeq.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

majSequences <- function(counts, n=3, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){

  seqnames <- apply(counts, 2, function(x){x <- sort(x, decreasing=TRUE); names(x)[1:n]})
  seqnames <- unique(unlist(as.character(seqnames)))

  sum <- apply(counts,2,sum)
  counts <- counts[seqnames,]
  sum <- matrix(sum,nrow(counts),ncol(counts),byrow=TRUE)
  p <- round(100*counts/sum,digits=3)

  if (outfile) png(filename="figures/majSeq.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
    maj <- apply(p, 2, max)
    seqname <- rownames(p)[apply(p, 2, which.max)]
    x <- barplot(maj, col=col[as.integer(group)], main="Percentage of reads from most expressed sequence",
	             ylim=c(0, max(maj)*1.2), las=2, ylab="Percentage of reads")
    legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
    for (i in 1:length(seqname)) text(x[i], maj[i]/2, seqname[i], cex=0.8, srt=90, adj=0)
  if (outfile) dev.off()
  
  return(invisible(p))
}
