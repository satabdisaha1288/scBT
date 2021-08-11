#' Title Compute the marginal log-likehihood of the p-genes under the null hypotheesis
#'
#' @param par Parameters of the function
#' @author Satabdi Saha
#' @return Marginal log-likelihood under the null hypothesis
#' @export
#'
#' @examples
log.lklh.marginal.null <- function(par, Y)
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
  #constant
  A_tot <- (
    R_y_2_grand_colsum + (par[1]^2/exp(par[2])) -
      (R_y_grand_colsum + par[1]/ exp(par[2]))^2/ (R_grand_sum + 1/exp(par[2]))
  )

  #  A_tot <- R_y_2_grand_colsum + ((par[1]^2)/par[2]) -
  #    (((R_y_grand_colsum + (par[1]/ par[2]))^2)/ (R_grand_sum + (1/par[2])))
  ## END OF DO NOT CHANGE
  l0 <- (0.5 *R_grand_sum)* log(2*pi)
  l1 <- (0.5)* log(1 + exp(par[2])*R_grand_sum)
  l2 <-  lgamma(exp(par[3])) + exp(par[3])*par[4]
  l3<- -lgamma(exp(par[3]) + R_grand_sum/2)
  l4<- (exp(par[3]) + R_grand_sum/2)*log (1/exp(par[4]) + A_tot/2)
  l5<- -lbeta(exp(par[5]) + R_grand_sum , exp(par[6]) + sum(n) - R_grand_sum)
  l6<- lbeta(exp(par[5]), exp(par[6]))
  log.lklh.marginal.null<-l0+l1+l2+l3+l4+l5+l6
  return(sum(log.lklh.marginal.null))
  #sum(log.lklh.marginal.null
}


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


#' Title Calculate FDR adjusted threshold for posterior probabilities of the null distribution
#' Newton, Michael A., Amine Noueiry, Deepayan Sarkar, and Paul Ahlquist. "Detecting differential gene expression with a semiparametric hierarchical mixture method." Biostatistics 5, no. 2 (2004): 155-176.
#' Tadesse, Mahlet G., Joseph G. Ibrahim, Robert Gentleman, Sabina Chiaretti, Jerome Ritz, and Robin Foa. "Bayesian Error‐in‐Variable Survival Model for the Analysis of GeneChip Arrays." Biometrics 61, no. 2 (2005): 488-497.
#' @author Satabdi Saha
#' @param DETest_output_Bayes : named (gene-names) vector of Bayes Factor values obtained from the Bayes Test
#' @param kappa : sequence of thresholding values preferably; kappa<-seq(0.01,1,0.01)
#' @param alpha : Desired FDR (False Discovery Rate) level
#'
#' @return max_kappa : Threshold value: Reject all hypothesis with posterior null probabilies less than max_kappa
#' @export
#'
#' @examples
calculate_threshold_posterior_prob_null_bayes<-function(DETest_output_Bayes, kappa, alpha){
  #Posterior probabilities under the null
  posterior_prob_null<-1/(1+ (1/DETest_output_Bayes))
  names(posterior_prob_null)<-names(DETest_output_Bayes)
  jkappa<-list()
  ckappa<-vector()
  val<-vector()

  for(i in 1: length(kappa))
  {
    jkappa[[i]]<-names(which(posterior_prob_null< kappa[i]))
    ckappa[i]<-sum(posterior_prob_null[jkappa[[i]]])
    val[i]<-ckappa[i]/length(jkappa[[i]])
  }
  #Value of kappa for controlling FDR at alpha
  max_kappa<-max(kappa[which(val<= alpha)])
  return(max_kappa)
}

#' Title Calculate FDR adjusted threshold for posterior probabilities of the null distribution
#' Newton, Michael A., Amine Noueiry, Deepayan Sarkar, and Paul Ahlquist. "Detecting differential gene expression with a semiparametric hierarchical mixture method." Biostatistics 5, no. 2 (2004): 155-176.
#' Tadesse, Mahlet G., Joseph G. Ibrahim, Robert Gentleman, Sabina Chiaretti, Jerome Ritz, and Robin Foa. "Bayesian Error‐in‐Variable Survival Model for the Analysis of GeneChip Arrays." Biometrics 61, no. 2 (2005): 488-497.
#' @author Satabdi Saha
#' @param DETest_output_Bayes : named (gene-names) vector of Bayes Factor values obtained from the Bayes Test
#' @param kappa : sequence of thresholding values preferably; kappa<-seq(0.01,1,0.01)
#' @param alpha : Desired FDR (False Discovery Rate) level
#'
#' @return max_kappa : Threshold value: Reject all hypothesis with posterior null probabilies less than max_kappa
#' @export
#'
#' @examples
calcPosteriorProbNull <- function(bayes.out, kappa = seq(0.01, 1, 0.01), alpha = c(0.01, 0.05)){
  posterior_prob_null <- 1/(1+ (1/bayes.out$exp_bf))
  names(posterior_prob_null) <- rownames(bayes.out)
  jkappa <- list()
  ckappa <- vector()
  val <- vector()

  for(i in 1: length(kappa))
  {
    jkappa[[i]] <- names(which(posterior_prob_null < kappa[i]))
    ckappa[i] <- sum(posterior_prob_null[jkappa[[i]]])
    sum(posterior_prob_null[which(posterior_prob_null < kappa[i])])
    val[i] <- ckappa[i]/length(jkappa[[i]])
  }
  #Value of kappa for controlling FDR at alpha
  max_kappa <- sapply(alpha, function(x) max(kappa[which(val <= x)]))
  names(max_kappa) = alpha
  return(list(max_kappa = max_kappa, ppH0 = posterior_prob_null))
}

#' Title Calculate FDR adjusted threshold for posterior probabilities of the null distribution
#' Newton, Michael A., Amine Noueiry, Deepayan Sarkar, and Paul Ahlquist. "Detecting differential gene expression with a semiparametric hierarchical mixture method." Biostatistics 5, no. 2 (2004): 155-176.
#' Tadesse, Mahlet G., Joseph G. Ibrahim, Robert Gentleman, Sabina Chiaretti, Jerome Ritz, and Robin Foa. "Bayesian Error‐in‐Variable Survival Model for the Analysis of GeneChip Arrays." Biometrics 61, no. 2 (2005): 488-497.
#' @author Satabdi Saha
#' @param DETest_output_Bayes : named (gene-names) vector of Bayes Factor values obtained from the Bayes Test
#' @param kappa : sequence of thresholding values preferably; kappa<-seq(0.01,1,0.01)
#' @param alpha : Desired FDR (False Discovery Rate) level
#'
#' @return max_kappa : Threshold value: Reject all hypothesis with posterior null probabilies less than max_kappa
#' @export
#'
#' @examples
calcPosteriorProbNull <- function(bayes.out, kappa = seq(0.01, 1, 0.01), alpha = c(0.01, 0.05)){
  posterior_prob_null <- 1/(1+ (1/bayes.out$exp_bf))
  names(posterior_prob_null) <- rownames(bayes.out)
  jkappa <- list()
  ckappa <- vector()
  val <- vector()

  for(i in 1: length(kappa))
  {
    jkappa[[i]] <- names(which(posterior_prob_null < kappa[i]))
    ckappa[i] <- sum(posterior_prob_null[jkappa[[i]]])
    sum(posterior_prob_null[which(posterior_prob_null < kappa[i])])
    val[i] <- ckappa[i]/length(jkappa[[i]])
  }
  #Value of kappa for controlling FDR at alpha
  max_kappa <- sapply(alpha, function(x) max(kappa[which(val <= x)]))
  names(max_kappa) = alpha
  return(list(max_kappa = max_kappa, ppH0 = posterior_prob_null))
}


# Hyperparameter calculation functions

#' Calculate a_sigma, b_sigma
#' @author Satabdi Saha
#' @param sce Single Cell Object
#' @return a_sigma and b_sigma values
#
# calc_a_sigma_b_sigma <- function(sce){
#   data <- logcounts(sce)
#   sigma <- apply(data, 1, var)
#   mean_sigma <- mean(sigma[sigma>0])
#   var_sigma <- sd(sigma[sigma>0])^2
#   a_sigma <- (1/((mean_sigma^2)*var_sigma)) + 2
#   b_sigma <- 1/(mean_sigma* (a_sigma-1))
#   output <- c(a_sigma,b_sigma)
#   names(output) <- c("a_sigma","b_sigma")
#   return(output)
# }

#' Calculate a_w, b_w
#' @author Satabdi Saha
#' @param omega_mean Groupwise dropout means
#' @param omega_var Groupwise dropout variance
#' @return a_w and b_w values
#
# calc_a_w_b_w <- function(omega_mean,omega_var){
#   mean_omega_mean <- mean(omega_mean)
#   mean_omega_var <- mean(omega_var)
#   a_w <- round((((mean_omega_mean^2)*(1-mean_omega_mean)) / mean_omega_var) - mean_omega_mean ,2)
#   b_w <- round(a_w * ((1- mean_omega_mean)/mean_omega_mean) ,2 )
#   output <- c(a_w,b_w)
#   names(output) <- c("a_w","b_w")
#   return(output)
# }

#' Calculate a_w, b_w
#' @author Satabdi Saha
#' @param omega_mean Groupwise dropout means
#' @param omega_var Groupwise dropout variance
#' @return a_w and b_w values
#
# calc_alt_null <- function(sce, method = 'fixed', de.prob = 0.25){
#   if (method == 'fixed'){
#     de.prob = 0.25
#   } else {
#     ##Calcualte
#   }
#   prior_Alter <- c(1-((1-de.prob)^(1/nrow(sce))))
#   prior_Null <- 1-prior_Alter
#   p_alt <- c(prior_Alter = prior_Alter, prior_Null = prior_Null)
#   return(p_alt)
# }

# Bayes analysis functions

#' INSERT DESCRIPTION
#'
#' @param sce SingleCellExperiment object with a logcounts assay
#' and Dose column in the cell metadata
#'
#' @return INSERT DESCRIPTION
#'
#' @export
DoseMatrix2List <- function(sce) {
  data.list <- list()
  data <- t(as.matrix(logcounts(sce)))
  gene.metadata <- rowData(sce)
  cell.metadata <- colData(sce)
  dose_vec <- sort(unique(cell.metadata$Dose))

  for (dose in dose_vec){
    data.list[[dose]] <- data[which(cell.metadata$Dose == dose),]
  }

  return(data.list)
}

# Bayes analysis functions

#' INSERT DESCRIPTION
#'
#' @param sce SingleCellExperiment object with a logcounts assay
#' and Dose column in the cell metadata
#'
#' @return INSERT DESCRIPTION
#'
#' @export
sceCalcPriors <- function(sce, fixed.priors = TRUE){
  #Initialize list, tables, and vectors
  data.list <- list()
  m <- vector()
  priors <- data.frame(matrix(NA, nrow = 0, ncol = 5))
  colnames(priors) <- c('dose', 'sigma_mean', 'sigma_var', 'omega_mean', 'omega_var')

  data <- as.matrix(logcounts(sce))
  gene.metadata <- rowData(sce)
  cell.metadata <- colData(sce)
  dose_vec <- sort(unique(cell.metadata$Dose))

  for (dose in dose_vec){
    data.list[[dose]] <- data[,which(cell.metadata$Dose == dose)]

    # m[dose] <- mean(rowMeans(replace(data.list[[dose]], data.list[[dose]] == 0, NA), na.rm = TRUE), na.rm = TRUE)
    #
    # sigma_mean <- mean(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE), na.rm = TRUE)
    # sigma_var <- var(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE),na.rm = TRUE)
    # sigma <- calc_a_sigma_b_sigma(sce)
    #
    # #Dose group specific mean dropout proportions for all genes
    # omega_mean <- mean(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    # omega_var <- var(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    # omega <- calc_a_w_b_w(omega_mean, omega_var)
    #
    # alt <- calc_alt_null(sce)
    # priors <- c(sigma, omega, c(tau_k_mu = length(dose_vec)), alt)
    # if (fixed.priors){
    #   priors['a_w'] = 2
    #   priors['b_w'] = 2
    # }
    #TODO: Add the prior fixed or calculated to the returned values.
  }
  #return(list(split.simulated = data.list, priors = priors, m = m))
  return(list(split.simulated = data.list))
}

#' INSERT DESCRIPTION
#'
#' @param priors INSERT DESCRIPTION
#' @param detailed A boolean denoting how much info to show in result
#'
#' @return INSESRT DESCRIPTION - describe omega etc...
#'
# bayesDETest <- function(priors, detailed = FALSE){
#   library(dplyr)
#   data.list <- priors$split.simulated
#   #TODO: Look into estimating from real data to replace prior_Alter and prior_Null from KW/WRS/ANOVA. Add to priors step.
#
#   bf_multiple_01 <- list()
#   # For each gene
#   for(j in rownames(data.list[[1]])){
#     in.list <- data.list %>% purrr::map(~ as.matrix(.x[j,]))
#     names(in.list) <- paste0("Y_", 1:length(in.list))
#     bf_multiple_01[[j]] <- Bayes_factor_multiple(Y = in.list, priors, detailed)
#   }
#   #Converts to data frame
#   bf_multiple_01 <- do.call(rbind, lapply(bf_multiple_01, as.data.frame))
#   if (detailed){
#     bf_multiple_01$mu <- bf_multiple_01$l_D0
#     bf_multiple_01$omega <- bf_multiple_01$l_D1 + bf_multiple_01$l_D2 + bf_multiple_01$l_D3 + bf_multiple_01$l_D4 + bf_multiple_01$l_D4
#   }
#   #Extracts part of test that depends on dropout
#   bf_multiple_01$exp_bf <- exp(bf_multiple_01$l_Bayes_factor_01)
#
#   return(bf_multiple_01)
# }

#' Title Bayesian Test for dose-dependent differential gene expression analysis
#' @author Satabdi Saha
#' @param Y list of size K (dose levels) of vectors length n (number of cells).
#' @param prior priors object generated with sceCalcPriors.
#' @param detailed logical output detailed statistical estimates.
#'
#' @return output A list having the estimates of the log-likelihood function,
#' log prior odds, and the log Bayes factor test statistic.
#'
# Bayes_factor_multiple <- function(Y, prior, detailed = FALSE){
#   a_sigma <- prior$priors['a_sigma']
#   b_sigma <- prior$priors['b_sigma']
#   a_w <- prior$priors['a_w']
#   b_w <- prior$priors['b_w']
#   K <- prior$priors['tau_k_mu']
#   m <- prior$m
#   m_0 <- mean(m)
#   tau_k_mu <- rep(1,K)
#   tau_mu <- 1
#   prior_alt <- prior$priors['prior_Alter']
#   prior_null <- prior$priors['prior_Null']
#
#
#   #Creates empty vectors to put the output from the test. This will also change depending on the loop structure.
#   ind_D0 <- matrix(NA, ncol = ncol(Y[[1]]),nrow = K)
#   l_D0 <- vector(mode = "numeric",length = ncol(Y[[1]]))
#   l_D1 <- vector(mode = "numeric",length = ncol(Y[[1]]))
#   l_D2 <- vector(mode = "numeric",length = ncol(Y[[1]]))
#   ind_D3 <- matrix(NA, ncol = ncol(Y[[1]]),nrow = K)
#   l_D3 <- vector(mode = "numeric",length = ncol(Y[[1]]))
#   ind_D4 <- matrix(NA, ncol = ncol(Y[[1]]),nrow = K)
#   l_D4 <- vector(mode = "numeric",length = ncol(Y[[1]]))
#   l_D5 <- vector(mode = "numeric",length = ncol(Y[[1]]))
#   l_Bayes_factor_01 <- vector(mode = "numeric",length = ncol(Y[[1]]))
#
#   ## DO NOT CHANGE
#   # n: number of cells in each dose groups
#   n <- sapply(Y, function(x) nrow(x))
#   # Number of genes being tested in each dose group ( in case of independent testing p=1)
#   p <- sapply(Y, function(x) ncol(x))
#   # list of K dose level indicator matrices (1 for positive expression)
#   R <- lapply(Y, function(x) ifelse(x>0,1,0))
#   # vector of length dose_level giving the total number of positively expressed cells in each dose group
#   R_colsum <- sapply(R, function(x) colSums(x))
#   # constant denoting the total number of positively expressed cells across all dose groups
#   R_grand_sum <- sum(R_colsum)
#   #list of squared Y matrices where
#   Y_squared <- Map("*",Y,Y)
#   R_y_2 <- Map("*" ,R , Y_squared)
#   #vector of length dose_level
#   R_y_2_colsum <- sapply(R_y_2, function(x) colSums(x))
#   R_y <- Map("*" ,R , Y)
#   #vector of length dose_level
#   R_y_colsum <- sapply(R_y, function(x) colSums(x))
#   #constant
#   R_y_grand_colsum <- sum(R_y_colsum)
#   #constant
#   R_y_2_grand_colsum <- sum(R_y_2_colsum)
#   #vector of length dose_level
#   A <- R_y_2_colsum + ((m^2)/tau_k_mu) -
#        (((R_y_colsum + (m/tau_k_mu))^2)/ (R_colsum + (1/tau_k_mu)))
#   #constant
#   A_tot <- R_y_2_grand_colsum + ((m_0^2)/tau_mu) -
#      (((R_y_grand_colsum + (m_0/tau_mu))^2)/ (R_grand_sum + (1/tau_mu)))
#   ## END OF DO NOT CHANGE
#
#   #TODO: j = genes - remove the genes loop because the input is already 1 gene
#   for(j in 1:ncol(Y[[1]]))
#   {
#     #k - dose group. Try to change to an apply(matrix)/lapply(listfromlist)/tapply/sapply(list)
#     for(k in 1:K)
#     {
#       suppressWarnings(
#       if (is.nan(log(b_w + n[k] - R_colsum[k][j] - 1))){
#         a_w = 2
#         b_w = 5
#       }
#       )
#
#       ind_D0[k]<- (0.5)*log(1+ (tau_mu*R_colsum[k][j]))
#
#       ind_D3[k] <- (0.5) * (log(a_w + b_w + n[k] - 1) - log(a_w + R_colsum[k][j] - 1) -
#                           log(b_w + n[k] - R_colsum[k][j] - 1))
#
#       ind_D4[k] <- ((a_w + b_w + n[k] - 1)* log(a_w + b_w + n[k] - 1)) -
#         ((a_w + R_colsum[k][j] - 1) * log(a_w + R_colsum[k][j] - 1)) -
#         ((b_w + n[k] - R_colsum[k][j] - 1)*log(b_w + n[k] - R_colsum[k][j] - 1))
#     }
#     # We should then be able to remove this from being put into a list element [j]. We can also remove the [j] from all the inputs/variables
#     l_D0[j] <- (1-K) + ((0.5*(1-K))* log (2 * pi)) + colSums(ind_D0) - ((0.5)* log(1+ (tau_mu* R_grand_sum[j]))) +
#                (a_sigma + (0.5 *R_grand_sum[j])) * (- log ((1/b_sigma) + (0.5*A_tot[j])) +
#               (log ((1/b_sigma)+  (0.5* sum(A)[j]))))
#
#     l_D1[j] <- 0.5*(log(a_w + R_grand_sum[j] - 1) + log(b_w + sum(n) - R_grand_sum[j] - 1) -
#                      log(a_w + b_w + sum(n) - 1))
#     l_D2[j] <- ((a_w + R_grand_sum[j]-1)* log(a_w + R_grand_sum[j]-1)) -
#       ((a_w + b_w + sum(n) -1)* log(a_w + b_w + sum(n) -1)) +
#       ((b_w + sum(n) - R_grand_sum[j] -1)* log(b_w + sum(n) - R_grand_sum[j] -1))
#
#     #l is a sum of the ind
#     l_D3[j] <- colSums(ind_D3)
#     l_D4[j] <- colSums(ind_D4)
#     l_D5[j] <- (K-1)* lbeta(a_w,b_w)
#     l_likelihood <- l_D0[j] + l_D1[j] + l_D2[j] + l_D3[j] + l_D4[j] + l_D5[j]
#     l_prior_odds <- log(prior_alt) - log(prior_null)
#     l_Bayes_factor_01[j] <-  l_likelihood + l_prior_odds
#   }
#
#   output <- list(l_likelihood, l_prior_odds, l_Bayes_factor_01)
#   names(output) <- c("l_likelihood", "l_prior_odds", "l_Bayes_factor_01")
#   if (detailed){
#     output <- list(l_likelihood, l_prior_odds, l_Bayes_factor_01,
#                    l_D0, l_D1, l_D2, l_D3, l_D4, l_D5)
#     names(output) <- c("l_likelihood", "l_prior_odds", "l_Bayes_factor_01",
#                        "l_D0","l_D1","l_D2","l_D3","l_D4","l_D5")
#   }
#
#   return(output)
# }
