#' Check the parameters (when using edgeR)
#'
#' Check the format and the validity of the parameters which will be used for the analysis with edgeR.
#  For example, it is important that $alpha$ be a numeric of length 1 between 0 and 1. This function avoid
#  potential stupid bugs when running the suite of the script.
#'
#' @param projectName name of the project
#' @param author author of the statistical analysis/report
#' @param targetFile path to the design/target file
#' @param rawDir path to the directory containing raw counts files
#' @param featuresToRemove names of the features to be removed
#' @param varInt factor of interest
#' @param condRef reference biological condition
#' @param batch blocking factor in the design
#' @param alpha threshold of statistical significance
#' @param pAdjustMethod p-value adjustment method: \code{"BH"} (default) or \code{"BY"} for example
#' @param cpmCutoff counts-per-million cut-off to filter low counts
#' @param normalizationMethod normalization method: \code{"TMM"} (default), \code{"RLE"} (DESeq) or \code{"upperquartile"}
#' @param gene.selection selection of the features in MDSPlot
#' @param colors vector of colors of each biological condition on the plots
#' @return A boolean indicating if there is a problem in the parameters
#' @author Hugo Varet

checkParameters.edgeR <- function(projectName,author,targetFile,rawDir,
                                  featuresToRemove,varInt,condRef,batch,alpha,
                                  pAdjustMethod,cpmCutoff,gene.selection,
                                  normalizationMethod,colors){
  problem <- FALSE
  if (!is.character(projectName) | length(projectName)!=1){
    warning("projectName must be a character vector of length 1")
    problem <- TRUE
  }
  if (!is.character(author) | length(author)!=1){
    warning("author must be a character vector of length 1")
    problem <- TRUE
  }
  if (!is.character(targetFile) | length(targetFile)!=1 || !file.exists(targetFile)){
    warning("targetFile must be a character vector of length 1 specifying an accessible file")
    problem <- TRUE
  }
  if (!is.character(rawDir) | length(rawDir)!=1 || is.na(file.info(rawDir)[1,"isdir"]) | !file.info(rawDir)[1,"isdir"]){
    warning("rawDir must be a character vector of length 1 specifying an accessible directory")
    problem <- TRUE
  }  
  if (!is.null(featuresToRemove) && !is.character(featuresToRemove)){
    warning("featuresToRemove must be a character vector or equal to NULL")
    problem <- TRUE
  }
  if (!is.character(varInt) | length(varInt)!=1){
    warning("varInt must be a character vector of length 1")
    problem <- TRUE
  }
  if (!is.character(condRef) | length(condRef)!=1){
    warning("condRef must be a character vector of length 1")
    problem <- TRUE
  }
  if (!is.null(batch) && I(!is.character(batch) | length(batch)!=1)){
    warning("batch must be NULL or a character vector of length 1")
    problem <- TRUE
  }
  if (!is.numeric(alpha) | length(alpha)!=1 || I(alpha<=0 | alpha>=1)){
    warning("alpha must be a numeric vector of length 1 with a value between 0 and 1")
    problem <- TRUE
  }
  if (!is.character(pAdjustMethod) | length(pAdjustMethod)!=1 || !I(pAdjustMethod %in% p.adjust.methods)){
    warning(paste("pAdjustMethod must be a value in", paste(p.adjust.methods, collapse=", ")))
    problem <- TRUE
  }
  if (!is.numeric(cpmCutoff) | length(cpmCutoff)!=1 || cpmCutoff<0){
    warning("cpmCutoff must be a numeric vector of length 1 with a value equal to or greater than 0")
    problem <- TRUE
  }
  if (!is.character(normalizationMethod) | length(normalizationMethod)!=1 || !I(normalizationMethod %in% c("TMM","RLE","upperquartile"))){
    warning("gene.selection must be equal to 'TMM', 'RLE' or 'upperquartile'")
    problem <- TRUE
  }
  if (!is.character(gene.selection) | length(gene.selection)!=1 || !I(gene.selection %in% c("pairwise","common"))){
    warning("gene.selection must be equal to 'pairwise' or 'common'")
    problem <- TRUE
  }
  areColors <- function(col){
    sapply(col, function(X){tryCatch(is.matrix(col2rgb(X)), error=function(e){FALSE})})
  }
  if (!is.vector(colors) || !all(areColors(colors))){
    warning("colors must be a vector of colors")
    problem <- TRUE
  }
  
  if (!problem){
    print("All the parameters are correct")
  }
  return(invisible(problem))
}
