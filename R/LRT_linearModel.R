#' Title Reworked version of LRT_linearModel.R (no-loops)
#' @author Satabdi
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#'
#' @return A dataframe having the test statistic and p-values 
#' @export 
LRT_linearModel<-function(sce){
  data = as.matrix(logcounts(sce))
  dose = as.numeric(as.character(colData(sce)$Dose))
  
  output = apply(data, 1, function(x) runLRT_Linear(x, dose))
  output = data.frame(t(output))
  output$FDR = p.adjust(output$p.value.comb, 'fdr')
  return(output)
}


#' Title Reworked version of LRT_linearModel.R (no-loops)
#' @author Satabdi
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#'
#' @return A dataframe having the test statistic and p-values 
#' @export 
runLRT_Linear = function(data, dose){
  #fit.logistic <- glm(ifelse(data > 0, 1, 0) ~ 1, family = binomial)
  fit.logistic <- logistf(factor(ifelse(data > 0, 1, 0)) ~ 1)
  #fit.logistic.alt <- glm(ifelse(data > 0, 1, 0) ~ dose, family = binomial)
  fit.logistic.alt <- logistf(factor(ifelse(data > 0, 1, 0)) ~ dose)
  fit.mean <- lm(data[which(data > 0)] ~ 1)
  fit.mean.alt <- lm(data[which(data > 0)] ~ dose[which(data > 0)])
  loglik_normal <- (as.numeric(logLik(fit.mean)) - as.numeric(logLik(fit.mean.alt)))
  loglik_logistic <- (as.numeric(logLik(fit.logistic)) - as.numeric(logLik(fit.logistic.alt)))
  resultvec <- c(-2*loglik_normal, pchisq(-2*loglik_normal, df = 3-2, lower.tail = FALSE),
                 -2*loglik_logistic, pchisq(-2*loglik_logistic, df = 2-1, lower.tail = FALSE),
                 -2 * (loglik_normal + loglik_logistic), 
                 pchisq(-2 * (loglik_normal + loglik_logistic), df = 2, lower.tail = FALSE))
  names(resultvec) <- c("lrstat_norm",'p.value.norm',"lrstat_binom",'p.value.binom',"lrstat_comb",'p.value.comb')
  return(resultvec)
}