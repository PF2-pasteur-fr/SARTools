#' Load count files
#'
#' Load one count file per sample thanks to the file names in the target file.
#'
#' @param target target \code{data.frame} of the project returned by \code{loadTargetFile()}
#' @param rawDir path to the directory containing the count files
#' @param header a logical value indicating whether the file contains the names of the variables as its first line
#' @param skip number of lines of the data file to skip before beginning to read data
#' @param featuresToRemove vector of feature Ids (or character string common to feature Ids) to remove from the counts
#' @return The \code{matrix} of raw counts with row names corresponding to the feature Ids and column names to the sample names as provided in the first column of the target.
#' @details If \code{featuresToRemove} is equal to \code{"rRNA"}, all the features containing the character string "rRNA" will be removed from the counts.
#' @author Marie-Agnes Dillies and Hugo Varet

loadCountData <- function(target, rawDir="raw", header=FALSE, skip=0,
                          featuresToRemove=c("alignment_not_unique", "ambiguous", "no_feature", "not_aligned", "too_low_aQual")){

  labels <- as.character(target[,1])
  files <- as.character(target[,2])

  rawCounts <- read.table(paste(rawDir,files[1],sep="/"), sep="\t", quote="\"", header=header, skip=skip)
  rawCounts <- rawCounts[,1:2]
  colnames(rawCounts) <- c("Id", labels[1])
  cat("Loading files:\n")
  cat(files[1],": ",length(rawCounts[,labels[1]])," rows and ",sum(rawCounts[,labels[1]]==0)," null count(s)\n",sep="")

  for (i in 2:length(files)){
  	tmp <- read.table(paste(rawDir,files[i],sep="/"), sep="\t", header=header, skip=skip)
    tmp <- tmp[,1:2]
  	colnames(tmp) <- c("Id", labels[i])
  	rawCounts <- merge(rawCounts, tmp, by="Id", all=TRUE)
    cat(files[i],": ",length(tmp[,labels[i]])," rows and ",sum(tmp[,labels[i]]==0)," null count(s)\n",sep="")
  }
  
  rawCounts[is.na(rawCounts)] <- 0
  counts <- as.matrix(rawCounts[,-1])
  rownames(counts) <- rawCounts[,1]
  counts <- counts[order(rownames(counts)),]
  
  cat("\nFeatures removed:\n")
  for (f in setdiff(featuresToRemove,"")){
    match <- grep(f, rownames(counts))
    if (length(match)>0){
	  cat(rownames(counts)[match],sep="\n")
	  counts <- counts[-match,]
	}
  }

  cat("\nTop of the counts matrix:\n")
  print(head(counts))
  cat("\nBottom of the counts matrix:\n")
  print(tail(counts))
  return(counts)
}
