#' Percentage of null counts per sample
#'
#' Bar plot of the percentage of null counts per sample
#'
#' @param counts \code{matrix} of counts
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors of the bars (one color per biological condition)
#' @return A file named barplotNull.png in the figures directory
#' @author Marie-Agnes Dillies and Hugo Varet

barplotNull <- function(counts, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen")){
  png(filename="figures/barplotNull.png",width=400,height=400)
    percentage <- apply(counts, 2, function(x){sum(x == 0)})*100/nrow(counts)
    percentage.allNull <- (nrow(counts) - nrow(removeNull(counts)))*100/nrow(counts)
    barplot(percentage, las = 2,
            col = col[as.integer(group)],
		    ylab = "Proportion of null counts",
		    main = "Proportion of null counts per sample", 
	  	    ylim = c(0,1.2*max(percentage)))
    abline(h = percentage.allNull, lty = 2, lwd = 2)
    legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
  dev.off()
}
