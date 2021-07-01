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
DETest <- function(sce, method = "All", verbose = TRUE, fixed.priors = TRUE, return.time = TRUE){
  checkmate::assertClass(sce, "SingleCellExperiment")
  #Check if dose column is there
  #Check that method is valid
  stopifnot(validateDRsce(sce))
  
  if ("All" %in% method){
    method = c('BAYES', 'LRT.linear', 'LRT.multiple', 'CHISQ', 'ANOVA',
               'WRS', 'KW', 'MAST', 'LIMMA-TREND', 'SEURATBIMOD')
  }
  DETest.list = list()
  timing.summary = list()

  if ("BAYES" %in% method){
    if (verbose) {message("Running Bayes test...")}
    Y <- DoseMatrix2List(sce)
    
    timing <- system.time({
      optim_output_null <- optim(par = runif(6, -1.25, 1.25),   # Applying optim
                                 fn = log.lklh.marginal.null,
                                 Y = Y,
                                 control = list(maxit=1000), method="L-BFGS-B", upper=rep(5, 6),
                                 lower=rep(-5, 6))
      optim_output_alt <- optim(par = runif(6, -1.25, 1.25),   # Applying optim
                                fn = log.lklh.marginal.alt,
                                Y = Y,
                                control = list(maxit=1000), method="L-BFGS-B", upper=rep(5, 6),
                                lower=rep(-5, 6))
      opt_par <- function(x) {
        opt_par<-c(x[1], exp(x[-1]))
        return(opt_par)
      }
      
      par_null <- opt_par(optim_output_null$par)
      par_alt <- opt_par(optim_output_alt$par)
      prior.null <- 0.9
      prior.alt <- 0.1
      bayes.factor <- calculate_BF(Y, par_null , par_alt, prior.null, prior.alt)
      DETest.list[["ppH0"]] <- calculate_threshold_posterior_prob_null_bayes(bayes.factor$BF,seq(0.01,1,0.01),0.05)
      posterior_prob_null <- 1/(1+ (1/bayes.factor$BF))
      
      bayes.out <- data.frame(do.call('cbind', bayes.factor))
      bayes.out$adjusted.p <- posterior_prob_null
      DETest.list[["BAYES"]] <- bayes.out
    })
    timing.summary[["BAYES"]] <- timing
  }
  if ("LRT.linear" %in% method){
    if (verbose) {message("Running LRT linear test...")}
    timing <- system.time({
      DETest.list[["LRTLin"]] <- LRT_linearModel(sce)
    })
    timing.summary[["LRTLin"]] <- timing
  }
  if ("LRT.multiple" %in% method){  
    if (verbose) {message("Running LRT multiple test...")}
    timing <- system.time({
      priors <- sceCalcPriors(sce)
      DETest.list[["LRTMult"]] <- LRT_multipleModel(priors[[1]])
    })
    timing.summary[['LRTMult']] <- timing
  }
  if ("CHISQ" %in% method){
    if (verbose) {message("Running Chi Squared test...")}
    timing <- system.time({
      DETest.list[["CHISQ"]] <- batchChi(sce)
    })
    timing.summary[["CHISQ"]] <- timing
  }
  if ("ANOVA" %in% method){
    if (verbose) {message("Running ANOVA test...")}
    timing <- system.time({
      DETest.list[['ANOVA']] <- batchANOVA(sce)
    })
    timing.summary[['ANOVA']] <- timing
  }
  if ("WRS" %in% method){
    if (verbose) {message("Running Wilcoxon rank sum test...")}
    timing <- system.time({
      DETest.list[["WRS"]] <- batchWRS(sce)
    })
    timing.summary[["WRS"]] <- timing
  }
  if ("KW" %in% method){
    if (verbose) {message("Running Kruskal Wallis test...")}
    timing <- system.time({
      DETest.list[['KW']] <- batchKW(sce)
    })
    timing.summary[["KW"]] <- timing
  }
  if ("MAST" %in% method){
    if (verbose) {message("Running MAST test...")}
    #DETest.list[["MAST"]] <- runMASTDR(sce)
  }
  if ("LIMMA-TREND" %in% method){
    if (verbose) {message("Running limma-trend test...")}
    timing <- system.time({
      DETest.list[['LIMMA-TREND']] <- runLimmaTrend(sce)
    })
    timing.summary[['LIMMA-TREND']] <- timing
  }
  if ("SEURATBIMOD" %in% method){
    if (verbose) {message("Running Seurat Bimod test...")}
    timing <- system.time({
      DETest.list[['SEURATBIMOD']] <- runSeuratBimod(sce)
    })
    timing.summary[['SEURATBIMOD']] <- timing
  }
  if (return.time == TRUE){
    DETest.list[["timing"]] <- timing.summary
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


#' General function to run the statistical test abstracting much of the individual steps
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' @param method the statistical test(s) to run on the sce object. Can be any of the following: Bayes, LRT.linear, LRT.multiple, ANOVA, KW, WRS, MAST, ChiSQ
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
benchmarkSim <- function(sim, test, p.threshold, fc.threshold = 0, pct.expressed = 0){
  # Estimate base parameters from simulated data
  fc <- abs(calcFC(sim))
  pz <- calcZeroP(sim)
  
  test <- test[rownames(rowData(sim)), , drop = FALSE]
  indexed.df <- data.frame(
    adjusted.p = test$adjusted.p,
    truth = as.factor(ifelse(rowData(sim)$Model == 'Unchanged', 0, 1)),
    fc.max = data.frame(fc.max = apply(fc, 1, function(x) max(x))),
    pz.min = data.frame(pz.min = apply(pz, 1, function(x) min(x))) 
  )  
  
  pos.idx <- which(indexed.df$adjusted.p <= p.threshold &
                    indexed.df$fc.max >= fc.threshold &
                    indexed.df$pz.min <= (1-pct.expressed)
                  )
  
  pos <- indexed.df[pos.idx,]
  neg <- indexed.df[-pos.idx,]
  
  classification = data.frame(
    TP = sum(pos$truth == 1),
    FP = sum(pos$truth == 0),
    TN = sum(neg$truth == 0),
    FN = sum(neg$truth == 1)
  )
  
  return(classification)
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
benchmark <- function(indexed.df, test, p.threshold, fc.threshold, pct.expressed){
  # Estimate base parameters from simulated data
  
  pos.idx <- which(indexed.df$adjusted.p <= p.threshold &
                     indexed.df$fc.max >= fc.threshold &
                     indexed.df$pz.min <= (1-pct.expressed)
  )
  
  pos <- indexed.df[pos.idx,]
  neg <- indexed.df[-pos.idx,]
  
  classification = c(
    TP = sum(pos$truth == 1),
    FP = sum(pos$truth == 0),
    TN = sum(neg$truth == 0),
    FN = sum(neg$truth == 1),
    p.thresh = p.threshold
  )
  
  return(classification)
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
benchmarkSim_batch <- function(sim, test, fc.threshold = 0, pct.expressed = 0){
  # Estimate base parameters from simulated data
  fc <- abs(calcFC(sim))
  pz <- calcZeroP(sim)
  
  test <- test[rownames(rowData(sim)), , drop = FALSE]
  pval.list <- sort(unique(test$adjusted.p))
  
  indexed.df <- data.frame(
    adjusted.p = test$adjusted.p,
    truth = as.factor(ifelse(rowData(sim)$Model == 'Unchanged', 0, 1)),
    fc.max = data.frame(fc.max = apply(fc, 1, function(x) max(x))),
    pz.min = data.frame(pz.min = apply(pz, 1, function(x) min(x))) 
  )  

  classification.batch <- sapply(pval.list, 
                                 function(x) benchmark(indexed.df, 
                                                       test, x,
                                                       fc.threshold,
                                                       pct.expressed)
  )
  classification.batch <- data.frame(t(classification.batch))
  return(classification.batch)
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
summarizeBenchmark <- function(benchmark.out){
  TP = benchmark.out$TP
  FP = benchmark.out$FP
  TN = benchmark.out$TN
  FN = benchmark.out$FN
  
  df <- data.frame(
    FPR = FP/(FP + TN),
    TPR = TP/(TP + FN),
    FNR = FN/(FN + TP),
    TNR = TN/(TN + FP),
    precision = TP/(TP + FP),
    recall = TP/(TP + FN),
    ppv = TP/(TP + FP),
    npv = TN/(TN + FN),
    balancedAcc = ((TP/(TP + FN)) + (TN/(TN + FP)))/2 
  )
  df.out = cbind(benchmark.out, df)
  return(df.out)
}



