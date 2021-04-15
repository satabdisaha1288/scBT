#' Linear Model LRT for dose dependent differential gene expression
#' @author Satabdi Saha
#' @param sce A single cell object (Does not allow a logcount matrix with genes having
#' a dropout of 100%. Such genes should be filtered out before running the test)
#'
#' @return A data frame giving the test statistic and unadjusted p values
#' for the linear, logistic and combined regression models for p genes
#' @export
#'
#' @examples
batchLRT_linearModel <- function(sce){
  data = as.matrix(logcounts(sce))
  dose = colData(sce)$Dose
  output<- t(apply(data,1,runLRT_linearModel))
  colnames(output)<-c("teststat_norm", "pvalue_norm",
                      "teststat_logistic", "pvalue_logistic",
                      "teststat_comb", "pvalue_comb")
  return(output)
}
runLRT_linearModel<-function(x){
  fit.logistic<-glm(ifelse(x >0,1,0)~1,family = binomial)
  fit.logistic.alt<-glm(ifelse(x >0,1,0) ~ 
                          as.numeric(as.character(dose)),family = binomial)
  fit.mean<- lm(x[x>0]~1)
  fit.mean.alt<-lm(x[x>0] ~ as.numeric(as.character(dose)[x >0]))
  loglik_normal<-(as.numeric(logLik(fit.mean))-as.numeric(logLik(fit.mean.alt)))
  loglik_logistic<-(as.numeric(logLik(fit.logistic))-as.numeric(logLik(fit.logistic.alt)))
  resultvec <- c(-2*loglik_normal, pchisq(-2*loglik_normal, df = 3-2, lower.tail = FALSE),
                 -2*loglik_logistic, pchisq(-2*loglik_logistic, df = 2-1, lower.tail = FALSE),
                 -2 * (loglik_normal + loglik_logistic), 
                 pchisq(-2 * (loglik_normal + loglik_logistic), df = 2, lower.tail = FALSE))
  result <- matrix(resultvec, nrow=2, ncol=3, dimnames=list(metric=c('lrstat', 'p.value'), component=c('norm', 'binom', 'comb')))
  return(result)
}