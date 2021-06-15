#' Title Compute Bayes factors for all the p-genes
#'
#' @param Y: Data list of K (number of dose groups) matrices having dimension n (cells) times p (genes)
#' @param par_null : Parameters obtained after optimizing the marginal likelihood under the null hypothesis
#' @param par_alt : Parameters obtained after optimizing the marginal likelihood under the alternative hypothesis
#' @param prior.null : Prior probability of the null hypothesis
#' @param prior.alt : Prior probability of the alternative hypothesis
#'
#' @return A list containing the log-likelihood under the null hypothesis,
#' log-likelihood under the alternative hypothesis, log Bayes factor and the Bayes Factor values
#' @export
#'
#' @examples
calculate_BF <- function(Y, par_null , par_alt, prior.null, prior.alt)
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
  A <- R_y_2_colsum + ((rep(par_alt[1],K)^2)/rep(par_alt[2],K)) -
    (((R_y_colsum + ((rep(par_alt[1],K))/(rep(par_alt[2],K))))^2)/ (R_colsum + (1/rep(par_alt[2],K))))
  #constant
  A_tot <- R_y_2_grand_colsum + ((par_null[1]^2)/par_null[2]) -
      (((R_y_grand_colsum + (par_null[1]/ par_null[2]))^2)/ (R_grand_sum + (1/par_null[2])))
  ## END OF DO NOT CHANGE
  l0 <- (0.5 *R_grand_sum)* log(2*pi)
  l1 <- (0.5)* log(1 + par_null[2]*R_grand_sum)
  l2 <-  lgamma(par_null[3]) + par_null[3]* log(par_null[4])
  l3<- -lgamma(par_null[3] + R_grand_sum/2)
  l4<- (par_null[3] + R_grand_sum/2)*log (1/par_null[4] + A_tot/2)
  l5<- -lbeta(par_null[5] + R_grand_sum , par_null[6] + sum(n) - R_grand_sum)
  l6<- lbeta(par_null[5], par_null[6])
  log.lklh.marginal.null<- -(l0+l1+l2+l3+l4+l5+l6)

  ## END OF DO NOT CHANGE
  l7 <- (0.5 *R_grand_sum)* log(2*pi)
  l8 <- rowSums((0.5)* log(1 + (par_alt[2]*R_colsum)))
  l9 <-  lgamma(par_alt[3]) + (par_alt[3]* log(par_alt[4]))
  l10<- -lgamma(par_alt[3] + R_grand_sum/2)
  l11<- (par_alt[3] + R_grand_sum/2)*log (1/par_alt[4] + (rowSums(A)/2))
  l12gr<-matrix(NA,nrow=nrow(R_colsum),ncol=K)
  for(k in 1:K)
  {
    l12gr[,k]<- -lbeta(par_alt[5] + R_colsum[,k] , par_alt[6] + n[k] - R_colsum[,k]) +
      lbeta(par_alt[5],par_alt[6])
  }
  dimnames(l12gr)<-dimnames(R_colsum)
  l12<- rowSums(l12gr)
  log.lklh.marginal.alt<- -(l7+l8+l9+l10+l11+l12)
  log.BF <- log.lklh.marginal.null - log.lklh.marginal.alt -
    rep(log(prior.null), p[1]) + rep(log(prior.alt), p[1])
  BF <- exp(log.BF)
  return(list(log.lklh.marginal.null = log.lklh.marginal.null,
              log.lklh.marginal.alt = log.lklh.marginal.alt,
              log.BF=log.BF,
              BF = BF))
  #sum(log.lklh.marginal.null)
}


