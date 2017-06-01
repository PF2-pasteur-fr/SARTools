#' Write HTML report for edgeR analyses
#'
#' Write HTML report from graphs and tables created during the analysis with edgeR
#'
#' @param target target \code{data.frame} of the project returned by \code{loadTargetFile()}
#' @param counts \code{matrix} of counts returned by \code{loadCountData()}
#' @param out.edgeR the result of \code{run.edgeR()}
#' @param summaryResults the result of \code{summarizeResults.DESeq2()}
#' @param majSequences the result of \code{descriptionPlots()}
#' @param workDir path to the working directory
#' @param projectName name of the project
#' @param author name of the author of the analysis
#' @param targetFile path to the target file
#' @param rawDir path to the directory containing the counts files
#' @param featuresToRemove vector of features to remove from the counts matrix
#' @param varInt factor of interest (biological condition)
#' @param condRef reference condition for the factor of interest
#' @param batch variable to take as a batch effect
#' @param alpha threshold of statistical significance
#' @param pAdjustMethod p-value adjustment method: \code{"BH"} (default) or \code{"BY"}
#' @param colors vector of colors of each biological condition on the plots
#' @param gene.selection selection of the features in \code{MDSPlot()} (\code{"pairwise"} by default)
#' @param normalizationMethod normalization method: \code{"TMM"} (default), \code{"RLE"} (DESeq) or \code{"upperquartile"}
#' @details This function generates the HTML report for a statistical analysis with edgeR. It uses the tables and graphs created during the workflow as well as the parameters defined at the beginning of the script.
#' @author Hugo Varet

writeReport.edgeR <- function(target,counts,out.edgeR,summaryResults,majSequences,
                              workDir,projectName,author,targetFile,rawDir,
                              featuresToRemove,varInt,condRef,batch,
                              alpha,pAdjustMethod,colors,gene.selection,
                              normalizationMethod){
  rmarkdown::render(input=system.file("report_edgeR.rmd", package="SARTools"),
                    output_file=paste0(projectName, "_report.html"),
                    output_dir=workDir,
                    intermediates_dir=workDir,
                    knit_root_dir=workDir,
                    run_pandoc=TRUE,
                    quiet=TRUE,
                    clean=TRUE)
  # delete unwanted directory/file
  # unlink("cache",force=TRUE,recursive=TRUE)
  # unlink(paste0("report_edgeR.md"),force=TRUE)
  cat("HTML report created\n")
}
