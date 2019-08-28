#' Assess the estimations of the size factors
#'
#' Plots to assess the estimations of the size factors
#'
#' @param dds a \code{DESeqDataSet} object
#' @param group factor vector of the condition from which each sample belongs
#' @param col colors for the plots
#' @param outfile TRUE to export the figure in a png file
#' @param plots vector of plots to generate
#' @return Two files in the figures directory: diagSizeFactorsHist.png containing one histogram per sample and diagSizeFactorsTC.png for a plot of the size factors vs the total number of reads
#' @author Marie-Agnes Dillies and Hugo Varet

diagSizeFactorsPlots <- function(dds, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), 
                                 outfile=TRUE, plots=c("diag","sf_libsize")){
  # histograms
  if ("diag" %in% plots){
    ncol <- ifelse(ncol(counts(dds))<=4, ceiling(sqrt(ncol(counts(dds)))), 3)
    nrow <- ceiling(ncol(counts(dds))/ncol)
    if (outfile) png(filename="figures/diagSizeFactorsHist.png", width=cairoSizeWrapper(1400*ncol), height=cairoSizeWrapper(1400*nrow), res=300)
    p <- list()
    geomeans <- exp(rowMeans(log(counts(dds))))
    samples <- colnames(counts(dds))
    counts.trans <- log2(counts(dds)/geomeans)
    counts.trans <- counts.trans[which(!is.na(apply(counts.trans, 1, sum))),]
    for (j in 1:ncol(dds)){
      d <- data.frame(x=counts.trans[,j])
      p[[j]] <- ggplot(data=d, aes(x=.data$x)) +
        geom_histogram(bins=100) +
        scale_y_continuous(expand=expand_scale(mult=c(0.01, 0.05))) +
        xlab("") +
        ylab("") +
        ggtitle(paste0("Sample ", samples[j])) +
        geom_vline(xintercept=log2(sizeFactors(dds)[j]), linetype="dashed", color="red", size=1)
    }
    tmpfun <- function(...) grid.arrange(..., nrow=nrow, ncol=ncol,
                                         bottom=textGrob(expression(log[2]~(counts/geometric~mean)), gp=gpar(fontsize=15)),
                                         top=textGrob("Size factors diagnostic", x=0.01, just="left", gp=gpar(fontsize=20)))
    do.call(tmpfun, p)
    if (outfile) dev.off()
  }
  
  # total read counts vs size factors
  if ("sf_libsize" %in% plots){
    if (outfile) png(filename="figures/diagSizeFactorsTC.png", width=1800, height=1800, res=300)
    d <- data.frame(sf=sizeFactors(dds), libsize=colSums(counts(dds))/1e6, group, sample=colnames(dds))
    print(ggplot(data=d, aes(x=.data$sf, y=.data$libsize, color=.data$group, label=.data$sample)) + 
            geom_point(show.legend=FALSE, size=3) +
            scale_colour_manual(values=col) +
            geom_text_repel(show.legend=FALSE, size=5, point.padding=0.2) +
            xlab("Size factors") +
            ylab("Total number of reads (millions)") +
            ggtitle("Diagnostic: size factors vs total number of reads") +
            geom_abline(slope=coefficients(lm(libsize ~ sf + 0, data=d)), intercept=0, show.legend = FALSE, linetype="dashed", color="grey"))
    if (outfile) dev.off()
  }
}


