# ==========================================================================
# package initialization
# ==========================================================================
.onAttach = function(libname, pkgname) {
  msg <- c("----------------------------------------------",
          paste0("Welcome to SARTools version ", packageVersion("SARTools"),"."),
		  "R template scripts are available at the end of the vignette and on GitHub.")
  # checking DESeq2 version
  if (packageVersion("DESeq2") < "1.6.0" | packageVersion("DESeq2") >= "1.7.0"){
    msg <- c(msg,"warning: SARTools has been developped with DESeq2 1.6.X, your version of DESeq2 might be incompatible with the workflow.")
  }
  # checking edgeR version
  if (packageVersion("edgeR") < "3.8.0" | packageVersion("edgeR") >= "3.9.0"){
    msg <- c(msg,"warning: SARTools has been developped with edgeR 3.8.X, your version of edgeR might be incompatible with the workflow.")
  }
  msg <- c(msg,"----------------------------------------------")
  msg <- strwrap(msg, exdent=4, indent=4)
  packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}
