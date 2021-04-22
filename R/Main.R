#' General function to run the statistical test abstracting much of the individual steps
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' @param method the statistical test(s) to run on the sce object. Can be any of the following: Bayes, LRT.linear, LRT.multiple, ANOVA, KW, WRS, MAST, ChiSQ
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
DETest = function(sce, method = "All", verbose = FALSE){
  DETest.list = list()
  #Check/validate object
  if (method  == "Bayes" | method == "All"){
    if (verbose) {message("Running Bayes test...")}
    ###############
  }
  if (method  == "LRT.linear" | method == "All"){
    if (verbose) {message("Running LRT linear test...")}
    DETest.list[["LRTLin"]] = LRT_linearModel_reworked(sce)
  }
  if (method  == "LRT.multiple" | method == "All"){
    if (verbose) {message("Running LRT multiple test...")}
    #DETest needs to be fixed
  }
  if (method  == "CHISQ" | method == "All"){
    if (verbose) {message("Running Chi Squared test...")}
    DETest.list[["CHISQ"]] = runChi(sce)
  }
  if (method  == "ANOVA" | method == "All"){
    if (verbose) {message("Running ANOVA test...")}
    DETest.list[['ANOVA']] = batchANOVA(sce)
  }
  if (method  == "WRS" | method == "All"){
    if (verbose) {message("Running Wilcoxon rank sum test...")}
    DETest.list[["WRS"]] = batchWRS(sce)
  }
  if (method  == "KW" | method == "All"){
    if (verbose) {message("Running Kruskal Wallis test...")}
    DETest.list[['KW']] = batchKW(sce)
  }
  if (method  == "MAST" | method == "All"){
    if (verbose) {message("Running MAST test...")}
    #################
  }
  return(DETest.list)
}
