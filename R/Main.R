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


getDEGs = function(sce, DETestoutput, threshold = 0.05, bayes.threshold = 1/3, fc.threshold = 0, pct.expressed = 0, verbose = TRUE){
  fc = abs(calcFC(sce))
  fc.max = data.frame(apply(fc, 1, function(x) max(x)))
  pz = calcZeroP(sce)
  pz.min = data.frame(apply(pz, 1, function(x) min(x)))
  
  DEGenes.list = list()
  
  for (test in names(DETestoutput)){
    if (verbose) {message(paste0("Evaluating assay:", test))}
    filtered = NULL
    #How to deal with LRT_linear and Bayes?
    if (ncol(DETestoutput[[test]]) == 1){
      merged.df = Reduce(merge, lapply(list(DETestoutput[[test]], fc.max, pz.min), function(x) data.frame(x, rn = row.names(x))))
      colnames(merged.df) = c("Gene", "pvalue", "fc.max", "pz.min")
      filtered = merged.df %>% filter(pvalue <= threshold & fc.max >= fc.threshold & pz.min >= pct.expressed)
    } else if (test != 'LRTLin' & test != 'Bayes'){
      stopifnot(identical(dim(DETestoutput[[test]]), dim(fc)) == TRUE)
      DETestoutput[[test]][which(fc < fc.threshold)] = NA
      min.pval = data.frame(apply(DETestoutput[[test]], 1, function(x) min(x)))
      merged.df = Reduce(merge, lapply(list(min.pval, fc.max, pz.min), function(x) data.frame(x, rn = row.names(x))))
      colnames(merged.df) = c("Gene", "pvalue", "fc.max", "pz.min")
      filtered = merged.df %>% filter(pvalue <= threshold & fc.max >= fc.threshold & pz.min >= pct.expressed)
    }
    DEGenes.list[[test]] = filtered
  }
  return(DEGenes.list)
}