#' BCV plot (for edgeR dispersions)
#'
#' Biological Coefficient of Variation plot (for edgeR objects)
#'
#' @param dge a \code{DGEList} object
#' @param outfile TRUE to export the figure in a png file
#' @return A file named BCV.png in the figures directory with a BCV plot produced by the \code{plotBCV()} function of the edgeR package
#' @author Marie-Agnes Dillies and Hugo Varet

BCVPlot <- function(dge, outfile=TRUE){	
  if (outfile) png(filename="figures/BCV.png", width=2100, height=1800, res=300)
  A <- dge$AveLogCPM
  if (is.null(A)) A <- aveLogCPM(dge$counts, offset = getOffset(dge))
  disp <- getDispersion(dge)
  if (is.null(disp)) stop("No dispersions to plot")
  if (attr(disp, "type") == "common") disp <- rep(disp, length = length(A))
  d <- data.frame(A=A,
                  sqrtdisp=sqrt(disp),
                  sqrttagwise=sqrt(dge$tagwise.dispersion),
                  sqrttrended=sqrt(dge$trended.dispersion),
                  sqrtcommon=sqrt(dge$common.dispersion))
  print(ggplot() +
          scale_color_manual("", values=c("black", "blue", "red"), labels=c("Tagwise", "Trend", "Common")) +
          geom_point(data=d, mapping=aes(x=.data$A, y=.data$sqrttagwise, color="a"), size=0.5, alpha=0.5) +
          geom_line(data=d, mapping=aes(x=.data$A, y=.data$sqrttrended, color="b")) +
          geom_hline(data=d, aes(yintercept=.data$sqrtcommon, color="c")) +
          xlab("Average log CPM") +
          ylab("Biological coefficient of variation") +
          ggtitle("BCV plot"))
  if (outfile) dev.off()	
}
