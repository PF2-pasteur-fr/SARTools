# ==========================================================================
# package initialization
# ==========================================================================
.onAttach = function(libname, pkgname) {
  msg <- c("----------------------------------------------",
          paste0("Welcome to SARTools version ", packageVersion("SARTools"),"."),
		  "R template scripts are available on GitHub.")
  # checking DESeq2 version
  if (packageVersion("DESeq2") < "1.6.0"){
    msg <- c(msg,"warning: SARTools needs DESeq2 1.6.X or higher, your version of DESeq2 might be incompatible with the workflow.")
  }
  # checking edgeR version
  if (packageVersion("edgeR") < "3.8.0"){
    msg <- c(msg,"warning: SARTools needs edgeR 3.8.X or higher, your version of edgeR might be incompatible with the workflow.")
  }
  msg <- c(msg,"----------------------------------------------")
  msg <- strwrap(msg, exdent=4, indent=4)
  packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}
