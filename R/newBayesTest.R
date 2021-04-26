#' INSERT DESCRIPTION
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return INSERT DESCRIPTION
#' 
#' @export
new_sceCalcPriors = function(sce){
  #Initialize list, tables, and vectors
  data.list = list()
  m = vector()
  priors = data.frame(matrix(NA, nrow = 0, ncol = 5))
  colnames(priors) = c('dose', 'sigma_mean', 'sigma_var', 'omega_mean', 'omega_var')
  
  data = as.matrix(logcounts(sce))
  gene.metadata = rowData(sce)
  cell.metadata = colData(sce)
  dose_vec = sort(unique(cell.metadata$Dose))
  
  for (dose in dose_vec){
    data.list[[dose]] = data[,which(cell.metadata$Dose == dose)]
    
    m[dose]<-mean(rowMeans(replace(data.list[[dose]], data.list[[dose]] == 0, NA), na.rm = TRUE), na.rm = TRUE)
    
    sigma_mean = mean(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE), na.rm = TRUE)
    sigma_var = var(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE),na.rm = TRUE)
    sigma = calc_a_sigma_b_sigma(sce)
    
    #Dose group specific mean dropout proportions for all genes
    omega_mean = mean(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    omega_var = var(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    omega = calc_a_w_b_w(omega_mean, omega_var)
    
    alt = calc_alt_null(sce)
    priors = c(sigma, omega, c(tau_k_mu = length(dose_vec)), alt)
    #TODO: Add the prior fixed or calculated to the returned values.
  }
  return(list(split.simulated = data.list, priors = priors, m = m))
}

#' INSERT DESCRIPTION
#' 
#' @param priors INSERT DESCRIPTION
#' @param detailed A boolean denoting how much info to show in result
#' 
#' @return INSESRT DESCRIPTION - describe omega etc...
#' 
#' @export
new_bayesDETest = function(priors, detailed = FALSE){
  library(dplyr)
  data.list = priors$split.simulated
  #TODO: Look into estimating from real data to replace prior_Alter and prior_Null from KW/WRS/ANOVA. Add to priors step.

  bf_multiple_01 = list()
  for(j in rownames(data.list[[1]])){ # For each gene
    in.list = data.list %>% purrr::map(~ as.matrix(.x[j,]))
    names(in.list) = paste0("Y_", 1:length(in.list))
    bf_multiple_01[[j]] <- new_Bayes_factor_multiple(Y = in.list, priors, detailed)
  }

  bf_multiple_01 = do.call(rbind, lapply(bf_multiple_01, as.data.frame)) #Converts to data frame
  if (detailed){
    bf_multiple_01$mu = bf_multiple_01$l_D0
    bf_multiple_01$omega = bf_multiple_01$l_D1 + bf_multiple_01$l_D2 + bf_multiple_01$l_D3 + bf_multiple_01$l_D4 + bf_multiple_01$l_D4 
  }
  bf_multiple_01$exp_bf = exp(bf_multiple_01$l_Bayes_factor_01) #Extracts part of test that depends on dropout
  
  return(bf_multiple_01)
}

#' Title Bayesian Test for dose-dependent differential gene expression analysis
#' @author Satabdi Saha 
#' @param Y list of size K (nlevels(dose_level)) of n (cells) times p(gene, in our case p=1) 
#' matrices of observed data having n cells and p genes in K treatment dose groups 
#' @param prior INSERT DESCRIPTION
#' @param detailed A boolean denoting how much info to show in result
#'
#' @return output: A list having the different parts of the log-likelihood function, log 
#' prior odds and the log Bayes factor test statistic
#' @export
new_Bayes_factor_multiple<-function(Y, prior, detailed = FALSE){
  a_sigma = prior$priors['a_sigma']
  b_sigma = prior$priors['b_sigma']
  a_w = prior$priors['a_w']
  b_w = prior$priors['b_w']
  K = prior$priors['tau_k_mu']
  m = prior$m
  m_0 = mean(m)
  tau_k_mu = rep(1:K)
  tau_mu = 1
  prior_alt = prior$priors['prior_Alter']
  prior_null = prior$priors['prior_Null']
  

  #Creates empty vectors to put the output from the test. This will also change depending on the loop structure.
  ind_D0<-matrix(NA, ncol = ncol(Y[[1]]),nrow = K)
  l_D0<-vector(mode = "numeric",length = ncol(Y[[1]]))
  l_D1<-vector(mode = "numeric",length = ncol(Y[[1]]))
  l_D2<-vector(mode = "numeric",length = ncol(Y[[1]]))
  ind_D3<-matrix(NA, ncol = ncol(Y[[1]]),nrow = K)
  l_D3<-vector(mode = "numeric",length = ncol(Y[[1]]))
  ind_D4<-matrix(NA, ncol = ncol(Y[[1]]),nrow = K)
  l_D4<-vector(mode = "numeric",length = ncol(Y[[1]]))
  l_D5<-vector(mode = "numeric",length = ncol(Y[[1]]))
  l_Bayes_factor_01<-vector(mode = "numeric",length = ncol(Y[[1]]))
  
  ## DO NOT CHANGE
  # n: number of cells in each dose groups
  n<- sapply(Y, function(x) nrow(x))
  # Number of genes being tested in each dose group ( in case of independent testing p=1)
  p<- sapply(Y, function(x) ncol(x))
  # list of K dose level indicator matrices (1 for positive expression)
  R <- lapply(Y, function(x) ifelse(x>0,1,0))
  # vector of length dose_level giving the total number of positively expressed cells in each dose group
  R_colsum<-sapply(R, function(x) colSums(x))
  # constant denoting the total number of positively expressed cells across all dose groups
  R_grand_sum<-sum(R_colsum)
  #list of squared Y matrices where
  Y_squared<- Map("*",Y,Y)
  R_y_2<-Map("*" ,R , Y_squared)
  #vector of length dose_level
  R_y_2_colsum<-sapply(R_y_2, function(x) colSums(x))
  R_y<-Map("*" ,R , Y)
  #vector of length dose_level
  R_y_colsum<-sapply(R_y, function(x) colSums(x))
  #constant
  R_y_grand_colsum<-sum(R_y_colsum)
  #constant
  R_y_2_grand_colsum<-sum(R_y_2_colsum)
  #vector of length dose_level
  A<- R_y_2_colsum + ((m^2)/tau_k_mu) - 
       (((R_y_colsum + (m/tau_k_mu))^2)/ (R_colsum + (1/tau_k_mu)))
  #constant
  A_tot<- R_y_2_grand_colsum + ((m_0^2)/tau_mu) -
     (((R_y_grand_colsum + (m_0/tau_mu))^2)/ (R_grand_sum + (1/tau_mu)))
  ## END OF DO NOT CHANGE
  
  #TODO: j = genes - remove the genes loop because the input is already 1 gene
  for(j in 1:ncol(Y[[1]]))
  {
    #k - dose group. Try to change to an apply(matrix)/lapply(listfromlist)/tapply/sapply(list)
    for(k in 1:K)
    {
      ind_D0[k]<- (0.5)*log(1+ (tau_mu*R_colsum[k][j]))
      
      ind_D3[k]<-0.5 * (log(a_w + b_w + n[k] - 1) - log(a_w + R_colsum[k][j] - 1) -
                          log(b_w + n[k] - R_colsum[k][j] - 1))
      
      ind_D4[k]<- ((a_w + b_w + n[k] - 1)* log(a_w + b_w + n[k] - 1)) -
        ((a_w + R_colsum[k][j] - 1) * log(a_w + R_colsum[k][j] - 1)) -
        ((b_w + n[k] - R_colsum[k][j] - 1)*log(b_w + n[k] - R_colsum[k][j] - 1))
    }
    
    # We should then be able to remove this from being put into a list element [j]. We can also remove the [j] from all the inputs/variables
    l_D0[j]<-(1-K) + ((0.5*(1-K))* log (2 * pi)) + colSums(ind_D0) - ((0.5)* log(1+ (tau_mu* R_grand_sum[j]))) +
               (a_sigma + (0.5 *R_grand_sum[j])) * (- log ((1/b_sigma) + (0.5*A_tot[j])) +
              (log ((1/b_sigma)+  (0.5* sum(A)[j]))))                
    l_D1[j]<- 0.5*(log(a_w + R_grand_sum[j] - 1) + log(b_w + sum(n) - R_grand_sum[j] - 1) -
                     log(a_w + b_w + sum(n) - 1))
    l_D2[j]<-((a_w + R_grand_sum[j]-1)* log(a_w + R_grand_sum[j]-1)) -
      ((a_w + b_w + sum(n) -1)* log(a_w + b_w + sum(n) -1)) +
      ((b_w + sum(n) - R_grand_sum[j] -1)* log(b_w + sum(n) - R_grand_sum[j] -1))
    
    #l is a sum of the ind
    l_D3[j]<- colSums(ind_D3)
    l_D4[j]<- colSums(ind_D4)
    l_D5[j]<- (K-1)* lbeta(a_w,b_w)
    l_likelihood<-  l_D0[j] + l_D1[j] + l_D2[j] + l_D3[j] + l_D4[j] +
    l_D5[j]
    l_prior_odds<- log(prior_alt) - log(prior_null)
    l_Bayes_factor_01[j]<-  l_likelihood + l_prior_odds
  }
  
  output <- list(l_likelihood, l_prior_odds, l_Bayes_factor_01)
  names(output) <- c("l_likelihood", "l_prior_odds", "l_Bayes_factor_01")
  if (detailed){
    output <- list(l_likelihood, l_prior_odds, l_Bayes_factor_01, 
                   l_D0, l_D1, l_D2, l_D3, l_D4, l_D5)
    names(output) <- c("l_likelihood", "l_prior_odds", "l_Bayes_factor_01",
                       "l_D0","l_D1","l_D2","l_D3","l_D4","l_D5")
  }

  return(output)
}
