#' Plots about DESeq2 dispersions
#'
#' A plot of the mean-dispersion relationship and a diagnostic of log normality of the dispersions (if use of DESeq2)
#'
#' @param dds a \code{DESeqDataSet} object
#' @return A file named dispersionsPlot.png in the figures directory containing the plot of the mean-dispersion relationship and a diagnostic of log normality of the dispersions
#' @author Marie-Agnes Dillies and Hugo Varet

dispersionsPlot <- function(dds){
  disp <- mcols(dds)$dispGeneEst
  disp <- disp[!is.na(disp)]
  disp <- disp[disp>1e-8]
  disp <- log(disp)
  mean.disp <- mean(disp,na.rm=TRUE)
  sd.disp <- sd(disp,na.rm=TRUE)
  png(filename="figures/dispersionsPlot.png",width=800,height=400)
	par(mfrow=c(1,2))
    # dispersions plot
	plotDispEsts(dds, main="Dispersions", las=1, xlab="Mean of normalized counts", ylab="Dispersion")
    # diagnostic of log normality
    hist(disp, freq=FALSE, nclass=50, xlab="Feature dispersion estimate", las=1,
         main = "log-normality dispersion diagnostic",col="skyblue")
    fun <- function(x){dnorm(x,mean=mean.disp,sd=sd.disp)}
    curve(fun,min(disp,na.rm=TRUE),max(disp,na.rm=TRUE),lwd=2,n=101,add=TRUE)
  dev.off()
}
