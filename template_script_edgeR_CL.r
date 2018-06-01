################################################################################
### R script to compare several conditions with the SARTools and edgeR packages
### Hugo Varet
### May 16th, 2018
### designed to be executed with SARTools 1.6.3
### run "Rscript template_script_edgeR_CL.r --help" to get some help
################################################################################

rm(list=ls())                                        # remove all the objects from the R session
library(optparse)                                    # to run the script in command lines

# options list with associated default value.
option_list <- list( 
make_option(c("-P", "--projectName"),
			default=basename(getwd()),
			dest="projectName",
			help="name of the project used for the report [default: name of the current directory]."),

make_option(c("-A", "--author"),
			default=Sys.info()[7],
			dest="author",
			help="name of the report author [default: %default]."),

make_option(c("-t", "--targetFile"),
			default="target.txt",
			dest="targetFile",
			help="path to the design/target file [default: %default]."),

make_option(c("-r", "--rawDir"),
			default="raw",
			dest="rawDir",
			help="path to the directory containing the HTSeq files [default: %default]."),		

make_option(c("-F", "--featuresToRemove"),
			default="alignment_not_unique,ambiguous,no_feature,not_aligned,too_low_aQual",
			dest="FTR",
			help="names of the features to be removed, more than once can be specified [default: %default]"),
			
make_option(c("-v", "--varInt"),
			default="group",
			dest="varInt", 
			help="factor of interest [default: %default]"),

make_option(c("-c", "--condRef"),
			default="WT",
			dest="condRef",
			help="reference biological condition [default: %default]"),

make_option(c("-b", "--batch"),
			default=NULL,
			dest="batch",
			help="blocking factor [default: %default] or \"batch\" for example"),

make_option(c("-a", "--alpha"),
			default=0.05,
			dest="alpha", 
			help="threshold of statistical significance [default: %default]"),

make_option(c("-p", "--pAdjustMethod"),
			default="BH",
			dest="pAdjustMethod", 
			help="p-value adjustment method: \"BH\" or \"BY\" [default: %default]"),

make_option(c("-m", "--cpmCutoff"),
			default=1,
			dest="cpmCutoff", 
			help="counts-per-million cut-off to filter low counts"),
			
make_option(c("-g", "--gene.selection"),
			default="pairwise",
			dest="gene.selection", 
			help="selection of the features in MDSPlot [default: %default]"),
			
make_option(c("-n", "--normalizationMethod"),
			default="TMM",
			dest="normalizationMethod", 
			help="normalization method in calcNormFactors: \"TMM\", \"RLE\" or \"upperquartile\" [default: %default]"),

make_option(c("-C", "--colors"),
			default="dodgerblue,firebrick1,MediumVioletRed,SpringGreen,chartreuse,cyan,darkorchid,darkorange",
			dest="cols",
			help="colors of each biological condition on the plots\n\t\t\"col1,col2,col3,col4\"\n\t\t[default: %default]"),

make_option(c("-f", "--forceCairoGraph"),
            action="store_true",
            default=FALSE,
            dest="forceCairoGraph",
            help="activate cairo type")
)

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
					   option_list=option_list, 
					   description="Compare two or more biological conditions in a RNA-Seq framework with edgeR.",
					   epilogue="For comments, bug reports etc... please contact Hugo Varet <hugo.varet@pasteur.fr>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

# get options and arguments
workDir <- getwd()
projectName <- opt$projectName                       # name of the project
author <- opt$author	                             # author of the statistical analysis/report
targetFile <- opt$targetFile                         # path to the design/target file
rawDir <- opt$rawDir								 # path to the directory containing raw counts files
featuresToRemove <- unlist(strsplit(opt$FTR, ","))   # names of the features to be removed (specific HTSeq-count information and rRNA for example)
varInt <- opt$varInt                                 # factor of interest
condRef <- opt$condRef                               # reference biological condition
batch <- opt$batch                                   # blocking factor: NULL (default) or "batch" for example
alpha <- as.numeric(opt$alpha)                       # threshold of statistical significance
pAdjustMethod <- opt$pAdjustMethod                   # p-value adjustment method: "BH" (default) or "BY"
gene.selection <- opt$gene.selection                 # selection of the features in MDSPlot
normalizationMethod <- opt$normalizationMethod       # normalization method in calcNormFactors
cpmCutoff <- opt$cpmCutoff                           # counts-per-million cut-off to filter low counts
colors <- unlist(strsplit(opt$cols, ","))            # vector of colors of each biologicial condition on the plots
forceCairoGraph <- opt$forceCairoGraph				 # force cairo as plotting device if enabled
# print(paste("workDir", workDir))
# print(paste("projectName", projectName))
# print(paste("author", author))
# print(paste("targetFile", targetFile))
# print(paste("rawDir", rawDir))
# print(paste("varInt", varInt))
# print(paste("condRef", condRef))
# print(paste("batch", batch))
# print(paste("alpha", alpha))
# print(paste("pAdjustMethod", pAdjustMethod))
# print(paste("featuresToRemove", featuresToRemove))
# print(paste("colors", colors))
# print(paste("gene.selection", gene.selection))
# print(paste("normalizationMethod", normalizationMethod))
# print(paste("cpmCutoff", cpmCutoff))

################################################################################
###                             running script                               ###
################################################################################
# setwd(workDir)
library(SARTools)
if (forceCairoGraph) options(bitmapType="cairo")

# checking parameters
problem <- checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)
if (problem) quit(save="yes")
					  
# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# edgeR analysis
out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

# MDS + clustering
exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
                  majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                  targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                  condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, cpmCutoff=cpmCutoff,
                  colors=colors, gene.selection=gene.selection, normalizationMethod=normalizationMethod)
