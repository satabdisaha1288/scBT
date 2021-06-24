#Y<-readRDS("data_sim.rds")
#' Title Compute the marginal log-likehihood of the p-genes under the alternative hypotheesis
#'
#' @param par Parameters of the function
#' @author Satabdi Saha
#'
#' @return Marginal log-likelihood under the alternative hypothesis
#' @export
#'
#' @examples
log.lklh.marginal.alt <- function(par, Y)
{
  ## DO NOT CHANGE
  #K: number of dose groups
  K <- length(Y)
  # n: number of cells in each dose groups
  n <- sapply(Y, function(x) nrow(x))
  # Number of genes being tested in each dose group ( in case of independent testing p=1)
  p <- sapply(Y, function(x) ncol(x))
  # list of K dose level indicator matrices (1 for positive expression)
  R <- lapply(Y, function(x) ifelse(x>0,1,0))
  # vector of length dose_level giving the total number of positively expressed cells in each dose group
  R_colsum <- sapply(R, function(x) colSums(x))
  # constant denoting the total number of positively expressed cells across all dose groups
  R_grand_sum <- rowSums(R_colsum)
  #list of squared Y matrices where
  Y_squared <- Map("*",Y,Y)
  R_y_2 <- Map("*" ,R , Y_squared)
  #vector of length dose_level
  R_y_2_colsum <- sapply(R_y_2, function(x) colSums(x))
  R_y <- Map("*" ,R , Y)
  #vector of length dose_level
  R_y_colsum <- sapply(R_y, function(x) colSums(x))
  #constant
  R_y_grand_colsum <- rowSums(R_y_colsum)
  #constant
  R_y_2_grand_colsum <- rowSums(R_y_2_colsum)
  #vector of length dose_level
  A <- R_y_2_colsum + ((rep(par[1],K)^2)/rep(exp(par[2]),K)) -
    (((R_y_colsum + ((rep(par[1],K))/(rep(exp(par[2]),K))))^2)/ (R_colsum + (1/rep(exp(par[2]),K))))
  #constant
  A_tot <- R_y_2_grand_colsum + ((par[1]^2)/exp(par[2])) -
    (((R_y_grand_colsum + (par[1]/ exp(par[2])))^2)/ (R_grand_sum + (1/exp(par[2]))))
  ## END OF DO NOT CHANGE
  l0 <- (0.5 *R_grand_sum)* log(2*pi)
  l1 <- rowSums((0.5)* log(1 + (exp(par[2])*R_colsum)))
  l2 <-  lgamma(exp(par[3])) + (exp(par[3])*par[4])
  l3<- -lgamma(exp(par[3]) + R_grand_sum/2)
  l4<- (exp(par[3]) + R_grand_sum/2)*log (1/exp(par[4]) + (rowSums(A)/2))
  l5gr<-matrix(NA,nrow=nrow(R_colsum),ncol=K)
  for(k in 1:K)
  {
    l5gr[,k]<- -lbeta(exp(par[5]) + R_colsum[,k] , exp(par[6]) + n[k] - R_colsum[,k]) +
      lbeta(exp(par[5]),exp(par[6]))
  }
  dimnames(l5gr)<-dimnames(R_colsum)
  l5<- rowSums(l5gr)
  #l5<-rowSums(t(apply(l5gr, 1, function(x) x - (lbeta(rep(exp(par[5]),length(x)), rep(exp(par[6]),length(x)))))))
  log.lklh.marginal.alt<-l0+l1+l2+l3+l4+l5
  return(sum(log.lklh.marginal.alt))
}

