#' BCV plot (for edgeR dispersions)
#'
#' Biological Coefficient of Variation plot (for edgeR objects)
#'
#' @param dge a \code{DGEList} object
#' @return A file named BCV.png in the figures directory with a BCV plot produced by the \code{plotBCV()} function of the edgeR package
#' @author Marie-Agnes Dillies and Hugo Varet

BCVPlot <- function(dge){	
  png(filename="figures/BCV.png", width=400, height=400)
    plotBCV(dge, las = 1, main = "BCV plot")
  dev.off()	
}
