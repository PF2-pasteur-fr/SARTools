#' Write HTML report for DESeq2 analyses
#'
#' Write HTML report from graphs and tables created during the analysis with DESeq2
#'
#' @param target target \code{data.frame} of the project returned by \code{loadTargetFile()}
#' @param counts \code{matrix} of counts returned by \code{loadCountData()}
#' @param out.DESeq2 the result of \code{run.DESeq2()}
#' @param summaryResults the result of \code{summarizeResults.DESeq2()}
#' @param majSequences the result of \code{descriptionPlots()}
#' @param workDir working directory
#' @param projectName name of the project
#' @param author name of the author of the analysis
#' @param targetFile path to the target file
#' @param rawDir path to the directory containing the counts files
#' @param featuresToRemove vector of features to remove from the counts matrix
#' @param varInt factor of interest (biological condition)
#' @param condRef reference condition for the factor of interest
#' @param batch variable to take as a batch effect
#' @param fitType mean-variance relationship: \code{"parametric"} (default) or \code{"local"}
#' @param cooksCutoff outliers detection threshold
#' @param independentFiltering \code{TRUE} or \code{FALSE} to perform the independent filtering or not
#' @param alpha threshold of statistical significance
#' @param pAdjustMethod p-value adjustment method: \code{"BH"} or \code{"BY"} for instance
#' @param typeTrans transformation for PCA/clustering: \code{"VST"} or \code{"rlog"}
#' @param locfunc \code{"median"} (default) or \code{"shorth"} to estimate the size factors
#' @param colors vector of colors of each biological condition on the plots
#' @details This function generates the HTML report for a statistical analysis with DESeq2. It uses the tables and graphs created during the workflow as well as the parameters defined at the beginning of the script.
#' @author Hugo Varet

writeReport.DESeq2 <- function(target, counts, out.DESeq2, summaryResults, majSequences,
                               workDir, projectName, author, targetFile, rawDir,
                               featuresToRemove, varInt, condRef, batch, fitType,
                               cooksCutoff, independentFiltering, alpha, pAdjustMethod,
                               typeTrans, locfunc, colors){
  rmarkdown::render(input=system.file("report_DESeq2.rmd", package="SARTools"),
                    output_file=paste0(projectName, "_report.html"),
                    output_dir=workDir,
                    intermediates_dir=workDir,
                    knit_root_dir=workDir,
                    run_pandoc=TRUE,
                    quiet=TRUE,
                    clean=TRUE)
  # delete unwanted directory/file
  # unlink("cache",force=TRUE,recursive=TRUE)
  # unlink("report_DESeq2.md",force=TRUE)
  cat("HTML report created\n")
}
