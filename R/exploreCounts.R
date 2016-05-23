#' Explore counts structure
#'
#' Explore counts structure: PCA (DESeq2) or MDS (edgeR) and clustering
#'
#' @param object a \code{DESeqDataSet} from DESeq2 or \code{DGEList} object from edgeR
#' @param group factor vector of the condition from which each sample belongs
#' @param typeTrans transformation method for PCA/clustering with DESeq2: \code{"VST"} or \code{"rlog"}
#' @param gene.selection selection of the features in MDSPlot (\code{"pairwise"} by default)
#' @param col colors used for the PCA/MDS (one per biological condition)
#' @return A list containing the dds object and the results object
#' @author Hugo Varet

exploreCounts <- function(object, group, typeTrans="VST", gene.selection="pairwise",
                          col=c("lightblue","orange","MediumVioletRed","SpringGreen")){
  if (class(object)=="DESeqDataSet"){
    if (typeTrans == "VST") counts.trans <- assay(varianceStabilizingTransformation(object))
    else counts.trans <- assay(rlogTransformation(object))
    PCAPlot(counts.trans=counts.trans, group=group, col=col)
    clusterPlot(counts.trans=counts.trans, group=group)  
  } else if (class(object)=="DGEList"){
    MDSPlot(dge=object, group=group, col=col, gene.selection=gene.selection)
    clusterPlot(counts.trans=cpm(object, prior.count=2, log=TRUE), group=group)  
  } else{
    stop("The object is not a DESeqDataSet nor a DGEList")
  }
}
