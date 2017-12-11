#' Check the parameters (when using DESeq2)
#'
#' Check the format and the validity of the parameters which will be used for the analysis with DESeq2.
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
#' @param fitType mean-variance relationship: "parametric" (default) or "local"
#' @param cooksCutoff outliers detection (TRUE or FALSE)
#' @param independentFiltering TRUE/FALSE to perform independent filtering
#' @param alpha threshold of statistical significance
#' @param pAdjustMethod p-value adjustment method: "BH" (default) or "BY" for example
#' @param typeTrans transformation for PCA/clustering: "VST" ou "rlog"
#' @param locfunc "median" (default) or "shorth" to estimate the size factors
#' @param colors vector of colors of each biological condition on the plots
#' @return A boolean indicating if there is a problem in the parameters
#' @author Hugo Varet

checkParameters.DESeq2 <- function(projectName,author,targetFile,rawDir,
                                   featuresToRemove,varInt,condRef,batch,fitType,
                                   cooksCutoff,independentFiltering,alpha,pAdjustMethod,
                                   typeTrans,locfunc,colors){
  problem <- FALSE
  if (!is.character(projectName) | length(projectName)!=1){
    message("projectName must be a character vector of length 1")
    problem <- TRUE
  }
  if (!is.character(author) | length(author)!=1){
    message("author must be a character vector of length 1")
    problem <- TRUE
  }
  if (!is.character(targetFile) | length(targetFile)!=1 || !file.exists(targetFile)){
    message("targetFile must be a character vector of length 1 specifying an accessible file")
    problem <- TRUE
  }
  if (!is.character(rawDir) | length(rawDir)!=1 || is.na(file.info(rawDir)[1,"isdir"]) | !file.info(rawDir)[1,"isdir"]){
    message("rawDir must be a character vector of length 1 specifying an accessible directory")
    problem <- TRUE
  }
  if (!is.null(featuresToRemove) && !is.character(featuresToRemove)){
    message("featuresToRemove must be a character vector or equal to NULL")
    problem <- TRUE
  }
  if (!is.character(varInt) | length(varInt)!=1){
    message("varInt must be a character vector of length 1")
    problem <- TRUE
  }
  if (!is.character(condRef) | length(condRef)!=1){
    message("condRef must be a character vector of length 1")
    problem <- TRUE
  }
  if (!is.null(batch) && I(!is.character(batch) | length(batch)!=1)){
    message("batch must be NULL or a character vector of length 1")
    problem <- TRUE
  }
  if (!is.character(fitType) | length(fitType)!=1 || !I(fitType %in% c("parametric","local"))){
    message("fitType must be equal to 'parametric' or 'local'")
    problem <- TRUE
  }
  if (!is.logical(cooksCutoff) | length(cooksCutoff)!=1){
    message("cooksCutoff must be a boolean vector of length 1")
    problem <- TRUE
  }
  if (!is.logical(independentFiltering) | length(independentFiltering)!=1){
    message("independentFiltering must be a boolean vector of length 1")
    problem <- TRUE
  }
  if (!is.numeric(alpha) | length(alpha)!=1 || I(alpha<=0 | alpha>=1)){
    message("alpha must be a numeric vector of length 1 with a value between 0 and 1")
    problem <- TRUE
  }
  if (!is.character(pAdjustMethod) | length(pAdjustMethod)!=1 || !I(pAdjustMethod %in% p.adjust.methods)){
    message(paste("pAdjustMethod must be a value in", paste(p.adjust.methods, collapse=", ")))
    problem <- TRUE
  }
  if (!is.character(typeTrans) | length(typeTrans)!=1 || !I(typeTrans %in% c("VST","rlog"))){
    message("typeTrans must be equal to 'VST' or 'rlog'")
    problem <- TRUE
  }
  if (!is.character(locfunc) | length(locfunc)!=1 || !I(locfunc %in% c("median","shorth"))){
    message("locfunc must be equal to 'median' or 'shorth'")
    problem <- TRUE
  } else{
    if (locfunc=="shorth" & !I("genefilter" %in% installed.packages()[,"Package"])){
      message("Package genefilter is needed if using locfunc='shorth'")
      problem <- TRUE
    }
  }
  areColors <- function(col){
    sapply(col, function(X){tryCatch(is.matrix(col2rgb(X)), error=function(e){FALSE})})
  }
  if (!is.vector(colors) || !all(areColors(colors))){
    message("colors must be a vector of colors")
    problem <- TRUE
  }
  
  if (!problem){
    print("All the parameters are correct")
  }
  return(invisible(problem))
}
