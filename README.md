SARTools
========

SARTools is a R package dedicated to the differential analysis of RNA-seq data. It provides tools to generate descriptive and diagnostic graphs, to run the differential analysis with one of the well known DESeq2 or edgeR packages and to export the results into easily readable tab-delimited files. It also facilitates the generation of a HTML report which displays all the figures produced, explains the statistical methods and gives the results of the differential analysis. Note that SARTools does not intend to replace DESeq2 or edgeR: it simply provides an environment to go with them. For more details about the methodology behind DESeq2 or edgeR, the user should read their documentations and papers.

SARTools is distributed with two R script templates (`template_script_DESeq2.r` and `template_script_edgeR.r`) which use functions of the package. For a more fluid analysis and to avoid possible bugs when creating the final HTML report, the user is encouraged to use them rather than writing a new script. Two other scripts are available (`template_script_DESeq2_CL.r` and `template_script_edgeR_CL.r`) to run SARTools in a shell with the `Rscript` command. In that case, the `optparse` R package must be available to interpret the command line parameters.

How to install SARTools?
------------------------

### Within R

In addition to the SARTools package itself, the workflow requires the installation of several packages: DESeq2, edgeR, genefilter, xtable and knitr (all available online, see the dedicated webpages). SARTools needs R version 3.3.0 or higher, DESeq2 1.12.0 or higher and edgeR 3.12.0 or higher: old versions of DESeq2 or edgeR may be incompatible with SARTools.

To install the SARTools package from GitHub, open a R session and:
- Install devtools with `install.packages("devtools")` (if not installed yet)
- Notes:

	- Ubuntu users may have to install some libraries (libxml2-dev, libcurl4-openssl-dev and libssl-dev) to be able to install DESeq2 and devtools
	- Some users may have to install the pandoc and pandoc-citeproc libraries to be able to generate the final HTML reports

- For Windows users only, install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) or check that it is already installed (needed to build the package)
- Load the devtools R package with `library(devtools)`
- Run `install_github("PF2-pasteur-fr/SARTools", build_opts="--no-resave-data")`

### Using Conda

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/r-sartools/README.html)

1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Install the SARTools R library and its dependencies using conda `conda install r-sartools`

Note: if you want to set a dedicated conda environment for SARTools, use `conda create -n sartools r-sartools` and follow the instructions to activate it.

How to use SARTools?
--------------------

A HTML vignette is available within the vignettes folder on GitHub and provides extensive information on the use of SARTools. The user can also open it with `vignette("SARTools")` if it has been generated during the installation of the package. Note that it is not available when SARTools has been installed using conda.

Be careful to use the R script associated with the version of SARTools installed on your system.

Please read the NEWS file to see the latest improvements!

About SARTools
--------------
The SARTools package has been developped at PF2 - Institut Pasteur by M.-A. Dillies and H. Varet (hugo.varet@pasteur.fr). Thanks to cite H. Varet, L. Brillet-Guéguen, J.-Y. Coppee and M.-A. Dillies, _SARTools: A DESeq2- and EdgeR-Based R Pipeline for Comprehensive Differential Analysis of RNA-Seq Data_, PLoS One, 2016, doi: http://dx.doi.org/10.1371/journal.pone.0157022 when using this tool for any analysis published.
