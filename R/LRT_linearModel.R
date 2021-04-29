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
  mydf<-data.frame(t(data),dose)
  # Specify the columns that contain your predictor variables
  respIdx <- seq(1:(ncol(mydf)-1));
  # lm(y ~ x), for x being a single predictor
  output<-lapply(respIdx[-which(colSums(mydf)==0)], function(x) {
    fit.logistic<-glm(ifelse(mydf[,x] >0,1,0)~1,family = binomial)
    fit.logistic.alt<-glm(ifelse(mydf[,x] >0,1,0) ~ mydf[, ncol(mydf)],family = binomial)
    fit.mean<- lm(mydf[,x][mydf[,x]>0]~1)
    fit.mean.alt<-lm(mydf[,x][mydf[,x]>0] ~ mydf[, ncol(mydf)][mydf[,x]>0])
    loglik_normal<-(as.numeric(logLik(fit.mean))-as.numeric(logLik(fit.mean.alt)))
    loglik_logistic<-(as.numeric(logLik(fit.logistic))-as.numeric(logLik(fit.logistic.alt)))
    resultvec <- c(-2*loglik_normal, pchisq(-2*loglik_normal, df = 3-2, lower.tail = FALSE),
                   -2*loglik_logistic, pchisq(-2*loglik_logistic, df = 2-1, lower.tail = FALSE),
                   -2 * (loglik_normal + loglik_logistic), 
                   pchisq(-2 * (loglik_normal + loglik_logistic), df = 2, lower.tail = FALSE))
    names(resultvec) <- c("lrstat_norm",'p.value.norm',"lrstat_binom",'p.value.binom',"lrstat_comb",'p.value.comb')
    return(resultvec)
  })
  output<-data.frame(matrix(unlist(output), nrow=length(output), byrow=TRUE))
  colnames(output)<-c("lrstat_norm",'p.value.norm',"lrstat_binom",'p.value.binom',"lrstat_comb",'p.value.comb')
  rownames(output)<-rownames(logcounts(sce))[-which(rowSums(logcounts(sce))==0)]
  return(output) 
}
