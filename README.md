SARTools
========

SARTools is a R package dedicated to the differential analysis of RNA-seq data. It provides tools to generate descriptive and diagnostic graphs, to run the differential analysis with one of the well known DESeq2 or edgeR packages and to export the results into easily readable tab-delimited files. It also facilitates the generation of a HTML report which displays all the figures produced, explains the statistical methods and gives the results of the differential analysis. Note that SARTools does not intend to replace DESeq2 or edgeR: it simply provides an environment to go with them. For more details about the methodology behind DESeq2 or edgeR, the user should read their documentations and papers.

SARTools is distributed with two R script templates which use functions of the package. For a more fluid analysis and to avoid possible bugs when creating the final HTML report, the user is encouraged to use them rather than writing a new script.

How to install SARTools?
------------------------

In addition to the SARTools package itself, the workflow requires the installation of several packages: DESeq2, edgeR, genefilter, xtable and knitr (all available online, see the dedicated webpages). This current version of SARTools has been developed under R 3.1.1 and with DESeq2 1.6.1, edgeR 3.8.2, genefilter 1.48.1 and knitr 1.7. As a DESeq2 or edgeR update might make the workflow unusable due to modifications on the statistical models, care is recommended when updating these packages.

To install the SARTools package from GitHub, open a R session and:
- install DESeq2, edgeR, genefilter if not installed yet (see the dedicated webpages for Bioconductor packages)
- load the devtools R package: `library(devtools)` (after `install.packages("devtools")` if not installed yet)
- run `install_github("PF2-pasteur-fr/SARTools")`

How to use SARTools?
--------------------

A PDF vignette (tutorial.pdf) is available within the package and provides extensive information on the use of SARTools. To open it, run `vignette("tutorial",package="SARTools")`

About SARTools
--------------
The SARTools package has been developped at PF2 - Institut Pasteur by M.-A. Dillies and H. Varet (hugo.varet@pasteur.fr). Thanks to cite H. Varet, J.-Y. Coppee and M.-A. Dillies, _SARTools: a DESeq2- and edgeR-based R pipeline for comprehensive differential analysis of RNA-seq data_, 2014 (submitted) when using this tool for any analysis published.
