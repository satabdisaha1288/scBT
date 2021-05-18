#' Differential expression testing
#' 
#' The main function to run any of the statistical tests developed as part of
#' this package for the analysis of single-cell sequencing data.
#' 
#' @param sce SingleCellExperiment object with a logcounts assay and a "Dose" 
#' column in the rowData(sce)
#' @param method the statistical test(s) to run. Options are "BAYES", 
#' "LRT.linear", "LRT.multiple", "ANOVA", "KW", "WRS", "CHISQ", or "MAST". Leave
#' empty to run all tests on your dataset. 
#' 
#' @details 
#' Various statistical tests were adapted specifically for dose-response
#' single-cell/single-nuclei transcriptomic data. The \code{method} enable the 
#' following tests:
#' \enumerate{
#'   \item BAYES is our new Bayes implementation.. what happens when it goes to
#'   next line?
#'   \item LRT.linear
#' }
#' 
#' @return a list of data.frames containing the statistical output for each
#' individual tests performed. 
#' 
#' @references Authors... 2021
#' 
#' @examples
#' # Running all statistical tests
#' DEG.testing <- DETest(sim) 
#' 
#' @export
DETest <- function(sce, method = "All", verbose = FALSE){
  checkmate::assertClass(sce, "SingleCellExperiment")
  #Check if dose column is there
  #Check that method is valid
  stopifnot(validateDRsce(sce))
  
  if ("All" %in% method){
    method = c('BAYES', 'LRT.linear', 'LRT.multiple', 'CHISQ', 'ANOVA',
               'WRS', 'KW', 'MAST')
  }
  DETest.list = list()

  #if (method  == "BAYES" | method == "All"){
  if ("BAYES" %in% method){
    if (verbose) {message("Running Bayes test...")}
    priors <- sceCalcPriors(sce)
    DETest.list[["BAYES"]] <- bayesDETest(priors)
  }
  if ("LRT.linear" %in% method){
    if (verbose) {message("Running LRT linear test...")}
    DETest.list[["LRTLin"]] <- LRT_linearModel(sce)
  }
  if ("LRT.multiple" %in% method){  
    if (verbose) {message("Running LRT multiple test...")}
    priors <- sceCalcPriors(sce)
    DETest.list[["LRTMult"]] <- LRT_multipleModel(priors[[1]])
  }
  if ("CHISQ" %in% method){
    if (verbose) {message("Running Chi Squared test...")}
    DETest.list[["CHISQ"]] <- runChi(sce)
  }
  if ("ANOVA" %in% method){
    if (verbose) {message("Running ANOVA test...")}
    DETest.list[['ANOVA']] <- batchANOVA(sce)
  }
  if ("WRS" %in% method){
    if (verbose) {message("Running Wilcoxon rank sum test...")}
    DETest.list[["WRS"]] <- batchWRS(sce)
  }
  if ("KW" %in% method){
    if (verbose) {message("Running Kruskal Wallis test...")}
    DETest.list[['KW']] <- batchKW(sce)
  }
  if ("MAST" %in% method){
    if (verbose) {message("Running MAST test...")}
    #DETest.list[["MAST"]] <- runMASTDR(sce)
  }
  return(DETest.list)
}


#' General function to run the statistical test abstracting much of the individual steps
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' @param method the statistical test(s) to run on the sce object. Can be any of the following: Bayes, LRT.linear, LRT.multiple, ANOVA, KW, WRS, MAST, ChiSQ
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
getDEGs <- function(sce, DETestoutput, threshold = 0.05, bayes.threshold = 1/3, fc.threshold = 0, pct.expressed = 0, verbose = TRUE){
  fc <- abs(calcFC(sce))
  fc.max <- data.frame(apply(fc, 1, function(x) max(x)))
  pz <- calcZeroP(sce)
  pz.min <- data.frame(apply(pz, 1, function(x) min(x)))
  
  DEGenes.list <- list()
  
  for (test in names(DETestoutput)){
    if (verbose) {message(paste0("Evaluating assay:", test))}
    filtered <- NULL
    #How to deal with LRT_linear and Bayes?
    if (ncol(DETestoutput[[test]]) == 1){
      merged.df <- Reduce(merge, lapply(list(DETestoutput[[test]], fc.max, pz.min), function(x) data.frame(x, rn = row.names(x))))
      colnames(merged.df) <- c("Gene", "pvalue", "fc.max", "pz.min")
      filtered <- merged.df %>% filter(pvalue <= threshold & fc.max >= fc.threshold & pz.min <= (1-pct.expressed))
    } else if (test != 'LRTLin' & test != 'LRTMult' & test != 'BAYES'){
      stopifnot(identical(dim(DETestoutput[[test]]), dim(fc)) == TRUE)
      DETestoutput[[test]] <- DETestoutput[[test]][rownames(fc),]
      DETestoutput[[test]] <- as.matrix(DETestoutput[[test]])
      DETestoutput[[test]][which(fc < fc.threshold)] = NA
      min.pval <- data.frame(apply(DETestoutput[[test]], 1, function(x) min(x)))
      merged.df <- Reduce(merge, lapply(list(min.pval, fc.max, pz.min), function(x) data.frame(x, rn = row.names(x))))
      colnames(merged.df) <- c("Gene", "pvalue", "fc.max", "pz.min")
      filtered <- merged.df %>% filter(pvalue <= threshold & fc.max >= fc.threshold & pz.min <= (1-pct.expressed))
    } else if (test == 'BAYES') {
      bayes_exp_bf <- DETestoutput[[test]][,"exp_bf", drop = FALSE]
      merged.df <- Reduce(merge, lapply(list(bayes_exp_bf, fc.max, pz.min), function(x) data.frame(x, rn = row.names(x))))
      colnames(merged.df) <- c("Gene", "exp_bf", "fc.max", "pz.min")
      filtered <- merged.df %>% filter(exp_bf <= bayes.threshold & fc.max >= fc.threshold & pz.min <= (1-pct.expressed))
    } else if (test == 'LRTLin' | test == 'LRTMult') {
      FDR <- DETestoutput[[test]][,'FDR', drop = FALSE]
      merged.df <- Reduce(merge, lapply(list(FDR, fc.max, pz.min), function(x) data.frame(x, rn = row.names(x))))
      colnames(merged.df) <- c("Gene", "fdr", "fc.max", "pz.min")
      filtered <- merged.df %>% filter(fdr <= threshold & fc.max >= fc.threshold & pz.min <= (1-pct.expressed))
    }
    DEGenes.list[[test]] <- filtered
  }
  return(DEGenes.list)
}


#' General function to run the statistical test abstracting much of the individual steps
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' @param method the statistical test(s) to run on the sce object. Can be any of the following: Bayes, LRT.linear, LRT.multiple, ANOVA, KW, WRS, MAST, ChiSQ
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
TruthFromSim <- function(sim, fc.threshold = 0, pct.expressed = 0){
  fc <- abs(calcFC(sim))
  fc.max <- data.frame(fc.max = apply(fc, 1, function(x) max(x)))
  pz <- calcZeroP(sim)
  pz.min <- data.frame(pz.min = apply(pz, 1, function(x) min(x)))
  gene.meta <- data.frame(rowData(sim))
  merged <- Reduce(merge, lapply(list(gene.meta, fc.max, pz.min), function(x) data.frame(x, rn = row.names(x))))
  merged.deg <- merged %>% filter(DE_idx != 1 & fc.max >= fc.threshold & pz.min <= (1-pct.expressed))
  merged.deg <- merged.deg[,-1]
  return(merged.deg)
}

#' General function to run the statistical test abstracting much of the individual steps
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' @param method the statistical test(s) to run on the sce object. Can be any of the following: Bayes, LRT.linear, LRT.multiple, ANOVA, KW, WRS, MAST, ChiSQ
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
filterLow <- function(sce, pct.expressed = 0){
  pz <- calcZeroP(sce)
  pz.min <- data.frame(pz.min =  apply(pz, 1, function(x) min(x)))
  passFilt <- rownames(pz.min %>% filter(pz.min <= (1-pct.expressed)))
  sce <- sce[passFilt,]
  return(sce)
} 


#' General function to run the statistical test abstracting much of the individual steps
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' @param method the statistical test(s) to run on the sce object. Can be any of the following: Bayes, LRT.linear, LRT.multiple, ANOVA, KW, WRS, MAST, ChiSQ
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
validateDRsce <- function(sce){
  value <- FALSE
  checkmate::assertClass(sce, "SingleCellExperiment")
  if ("Dose" %in% colnames(colData(sce))){
    checkmate::assertClass(colData(sce)$Dose, "factor")
    value <- TRUE
  }
  return(value)
}


#' General function to run the statistical test abstracting much of the individual steps
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' @param method the statistical test(s) to run on the sce object. Can be any of the following: Bayes, LRT.linear, LRT.multiple, ANOVA, KW, WRS, MAST, ChiSQ
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
benchmarkDETests <- function(sce, dge.true, dge.list){
  all.genes <- rownames(sce)
  truth <- as.factor(all.genes %in% dge.true$Gene)
  names(truth) <- all.genes
  comparisons <- data.frame(truth)
  confusions.list <- list()
  confusions.results <- list()
  for (test in names(dge.list)){
    comparisons[,test] <- as.factor(all.genes %in% dge.list[[test]]$Gene)
    confusion <- caret::confusionMatrix(comparisons[,test], comparisons$truth)

    # Extract the classification table
    conf.mat <- confusion$table
    confusions.list[[test]] <- reshape2::melt(conf.mat)
    confusions.list[[test]]$tname <- test
    
    # Extract additional info
    confusions.results[[test]] <- rbind(as.matrix(confusion, what = "classes"), 
                                       as.matrix(confusion, what = "overall"))
    
  }
  df.confusionMat <- do.call(rbind, lapply(confusions.list, as.data.frame))
  df.confusionMat[which(df.confusionMat$Reference == TRUE & df.confusionMat$Prediction == TRUE), 'result'] <- 'True Positive'
  df.confusionMat[which(df.confusionMat$Reference == TRUE & df.confusionMat$Prediction == FALSE), 'result'] <- 'False Negative'
  df.confusionMat[which(df.confusionMat$Reference == FALSE & df.confusionMat$Prediction == TRUE), 'result'] <- 'False Positive'
  df.confusionMat[which(df.confusionMat$Reference == FALSE & df.confusionMat$Prediction == FALSE), 'result'] <- 'True Negative'
  
  df.results <- do.call(cbind, lapply(confusions.results, as.data.frame))
  colnames(df.results) <- names(confusions.results)
  df.results$value <- rownames(df.results)
  df.results <- reshape2::melt(df.results)
  colnames(df.results) <- c('metric', 'test', 'value')
  
  return(list(classification = df.confusionMat, results = df.results))
}
