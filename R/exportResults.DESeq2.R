#' Export results for DESeq2 analyses
#'
#' Export counts and DESeq2 results 
#'
#' @param out.DESeq2 the result of \code{run.DESeq2()}
#' @param group factor vector of the condition from which each sample belongs
#' @param alpha threshold to apply to adjusted p-values
#' @param export TRUE to export the results (up, down and complete table) in the tables directory
#' @return A list of \code{data.frame} containing counts, pvalues, FDR, log2FC...
#' @author Marie-Agnes Dillies and Hugo Varet

exportResults.DESeq2 <- function(out.DESeq2, group, alpha=0.05, export=TRUE){
  
  dds <- out.DESeq2$dds
  results <- out.DESeq2$results
  
  # comptages bruts et normalis?s
  counts <- data.frame(Id=rownames(counts(dds)), counts(dds), round(counts(dds, normalized=TRUE)))
  colnames(counts) <- c("Id", colnames(counts(dds)), paste0("norm.", colnames(counts(dds))))
  # baseMean avec identifiant
  bm <- data.frame(Id=rownames(results[[1]]),baseMean=round(results[[1]][,"baseMean"],2))
  # merge des info, comptages et baseMean selon l'Id
  base <- merge(counts, bm, by="Id", all=TRUE)
  tmp <- base[,paste("norm", colnames(counts(dds)), sep=".")]
  for (cond in levels(group)){
    base[,cond] <- round(apply(as.data.frame(tmp[,group==cond]),1,mean),0)
  }
  
  complete <- list()
  for (name in names(results)){
    complete.name <- base

    # ajout d'elements depuis results
    res.name <- data.frame(Id=rownames(results[[name]]),
                           FoldChange=round(2^(results[[name]][,"log2FoldChange"]), 3),
                           log2FoldChange=round(results[[name]][,"log2FoldChange"], 3),
                           stat=round(results[[name]][,"stat"], 3),
                           pvalue=results[[name]][,"pvalue"],
                           padj=results[[name]][,"padj"])
    complete.name <- merge(complete.name, res.name, by="Id", all=TRUE)
    # ajout d'elements depuis mcols(dds)
    mcols.add <- data.frame(Id=rownames(counts(dds)),dispGeneEst=round(mcols(dds)$dispGeneEst,4),
                            dispFit=round(mcols(dds)$dispFit,4),dispMAP=round(mcols(dds)$dispMAP,4),
                            dispersion=round(mcols(dds)$dispersion,4),betaConv=mcols(dds)$betaConv,
                            maxCooks=round(mcols(dds)$maxCooks,4))
    complete.name <- merge(complete.name, mcols.add, by="Id", all=TRUE)
    complete[[name]] <- complete.name
    
    if (export){
      # s?lection des up et down
      up.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange>=0),]
      up.name <- up.name[order(up.name$padj),]
      down.name <- complete.name[which(complete.name$padj <= alpha & complete.name$betaConv & complete.name$log2FoldChange<=0),]
      down.name <- down.name[order(down.name$padj),]
      
      # exports
      name <- gsub("_","",name)
      write.table(complete.name, file=paste0("tables/",name,".complete.txt"), sep="\t", row.names=FALSE, dec=".", quote=FALSE)
      write.table(up.name, file=paste0("tables/", name,".up.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
      write.table(down.name, file=paste0("tables/", name,".down.txt"), row.names=FALSE, sep="\t", dec=".", quote=FALSE)
    }
  }

  return(complete)
}
