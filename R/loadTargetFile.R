#' Load target file
#'
#' Load the target file containing sample information
#'
#' @param targetFile path to the target file
#' @param varInt variable on which sorting the target
#' @param condRef reference condition of \code{varInt}
#' @param batch batch effect to take into account
#' @return A \code{data.frame} containing the informations about the samples (name, file containing the counts and biological condition)
#' @details The \code{batch} parameter is used only to check if it is available in the target file before running the suite of the script.
#' @author Marie-Agnes Dillies and Hugo Varet

loadTargetFile <- function(targetFile, varInt, condRef, batch){
  target <- read.table(targetFile, header=TRUE, sep="\t", na.strings="")
  if (!I(varInt %in% names(target))) stop(paste("The factor of interest", varInt, "is not in the target file"))
  if (!is.null(batch) && !I(batch %in% names(target))) stop(paste("The batch effect", batch, "is not in the target file")) 
  target[,varInt] <- as.factor(target[,varInt])
  if (!I(condRef %in% as.character(target[,varInt]))) stop(paste("The reference level", condRef, "is not a level of the factor of interest"))
  target[,varInt] <- relevel(target[,varInt],ref=condRef)
  target <- target[order(target[,varInt]),]
  rownames(target) <- as.character(target[,1])
  # check if varInt contains replicates
  if (min(table(target[,varInt]))<2) stop(paste("The factor of interest", varInt, "has a level without replicates"))
  # check if NA in the target
  if (any(is.na(cbind(target[,c(varInt, batch)], target[,1:2])))) stop("NA are present in the target file")
  # warning message if batch is numeric
  if (!is.null(batch) && is.numeric(target[,batch])) warning(paste("The", batch, "variable is numeric. Use factor() or rename the levels with letters to convert it into a factor"))
  if (any(grepl("[[:punct:]]", as.character(target[,varInt])))) stop(paste("The", varInt, "variable contains punctuation characters, please remove them"))
  cat("Target file:\n")
  print(target)
  return(target)
}
