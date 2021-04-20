#' Calculation of Bayes factor
#' 
#' @param Y_t list of matrices of observed data having n cells and p genes
#' @param Y_c list of matrices of observed data having n cells and p genes 
#' in the control group
#' @param m_t0 INSERT DESCR
#' @param m_c0 INSERT DESCR
#' @param m_0 INSERT DESCR
#' @param tau_t_mu INSERT DESCR
#' @param tau_c_mu INSERT DESCR
#' @param tau_mu INSERT DESCR
#' @param b_sigma INSERT DESCR
#' @param a_sigma INSERT DESCR
#' @param a_w INSERT DESCR
#' @param a_c_w INSERT DESCR
#' @param a_t_w INSERT DESCR
#' @param b_w INSERT DESCR
#' @param b_t_w INSERT DESCR
#' @param b_c_w INSERT DESCR
#' @param prior_alt INSERT DESCR
#' @param prior_null INSERT DESCR
#' 
#' @details
#'  
#' @return INSERT DESCR
#'  
Bayes_factor<-function(Y_t, Y_c, m_t0, m_c0, m_0, tau_t_mu, tau_c_mu,
                       tau_mu, b_sigma, a_sigma, a_w, a_c_w, a_t_w,
                       b_w, b_t_w, b_c_w, prior_alt, prior_null)
{
  #Bayes_factor_01<-vector(mode = "numeric",length = ncol(Y_t))
  l_D0<-vector(mode = "numeric",length = ncol(Y_t))
  l_D1<-vector(mode = "numeric",length = ncol(Y_t))
  l_D2<-vector(mode = "numeric",length = ncol(Y_t))
  l_D3<-vector(mode = "numeric",length = ncol(Y_t))
  l_D4<-vector(mode = "numeric",length = ncol(Y_t))
  l_D5<-vector(mode = "numeric",length = ncol(Y_t))
  l_D6<-vector(mode = "numeric",length = ncol(Y_t))
  l_D7<-vector(mode = "numeric",length = ncol(Y_t))
  l_Bayes_factor_01<-vector(mode = "numeric",length = ncol(Y_t))
  n_1<-nrow(Y_t)
  n_2<-nrow(Y_c)
  p<-ncol(Y_t)
  R_t<-ifelse( Y_t > 0, 1, 0)
  R_c<-ifelse( Y_c > 0, 1, 0)
  R_t_colsum<-colSums(R_t)
  R_c_colsum<-colSums(R_c)
  R_t_c_colsum<-colSums(R_t) + colSums(R_c)
  R_t_y_t_2_colsum<-colSums(R_t * (Y_t^2 ))
  R_c_y_c_2_colsum<-colSums(R_c * (Y_c^2 ))
  R_t_y_t_colsum<-colSums(R_t * Y_t)
  R_c_y_c_colsum<-colSums(R_c * Y_c)
  R_t_c_y_t_c_colsum<-colSums((R_t * Y_t)) + colSums((R_c * Y_c))
  R_t_c_y_t_c_2_colsum<-colSums((R_t * (Y_t^2))) + colSums((R_c * (Y_c^2)))
  A_t<- (R_t_y_t_2_colsum + ((m_t0^2)/tau_mu)) - 
                    (((R_t_y_t_colsum + (m_t0/tau_mu))^2)/
                           (R_t_colsum + (1/tau_mu)))
  A_c<- (R_c_y_c_2_colsum + ((m_c0^2)/tau_mu)) - 
                  (((R_c_y_c_colsum + (m_c0/tau_mu))^2)/
                      (R_c_colsum + (1/tau_mu)))
  A_t_c<- (R_t_c_y_t_c_2_colsum + ((m_0^2)/tau_mu))-
                  (((R_t_c_y_t_c_colsum +(m_0/tau_mu))^2)/
                  ( R_t_c_colsum + (1/tau_mu)))
  #eps <- 1e-14
  for(j in 1:p)
    {
      l_D0[j]<- (0.5 * log(R_t_colsum[j] + (1/tau_mu))) + (0.5* log(R_c_colsum[j] + (1/tau_mu)))-
              ( 0.5 * log(R_t_c_colsum[j] + (1/tau_mu))) + 
              ((a_sigma + (0.5*R_t_c_colsum[j]))*(log((1/b_sigma)+ (0.5* A_t[j])+(0.5* A_c[j]))))-
              ((a_sigma + (0.5*R_t_c_colsum[j]))* (log((1/b_sigma)+ (0.5* A_t_c[j]))))
      l_D1[j]<- 0.5*(log(a_w + R_t_c_colsum[j] - 1) + log(b_w + n_1 + n_2 - R_t_c_colsum[j] - 1) -
                       log(a_w + b_w + n_1 + n_2 - 1))
      l_D2[j]<-((a_w + R_t_c_colsum[j]-1)* log(a_w + R_t_c_colsum[j]-1)) -
               ((a_w + b_w + n_1 + n_2 -1)* log(a_w + b_w + n_1 + n_2 -1)) +
               ((b_w + n_1 + n_2 - R_t_c_colsum[j] -1)* log(b_w + n_1 + n_2 - R_t_c_colsum[j] -1))
      l_D3[j]<-0.5 * (log(a_w + b_w + n_1 - 1) - log(a_w + R_t_colsum[j] - 1) -
               log(b_w + n_1 - R_t_colsum[j] - 1))
      l_D4[j]<- ((a_w + b_w + n_1 - 1)* log(a_w + b_w + n_1 - 1)) -
                ((a_w + R_t_colsum[j] - 1) * log(a_w + R_t_colsum[j] - 1)) -
                ((b_w + n_1 - R_t_colsum[j] - 1)*log(b_w + n_1 - R_t_colsum[j] - 1))
      l_D5[j]<-0.5* (log(a_w - 1) + log(b_w - 1) - log(a_w + b_w - 1)) +
              ((a_w - 1)*log(a_w - 1)) + ((b_w - 1) * log(b_w - 1)) +
              ((a_w + b_w - 1)*log(a_w + b_w - 1))
      l_D6[j]<-0.5 * (log(a_w + b_w + n_2 - 1) - log(a_w + R_c_colsum[j] - 1) -
                        log(b_w + n_2 - R_c_colsum[j] - 1)) 
      l_D7[j]<-((a_w + b_w + n_2 - 1)*log(a_w + b_w + n_2 - 1)) -
               ((a_w + R_c_colsum[j] - 1)*log(a_w + R_c_colsum[j] - 1)) -
               ((b_w + n_2 - R_c_colsum[j] - 1)* log(b_w + n_2 - R_c_colsum[j] - 1))
      l_prior_odds<- log (prior_alt) - log(prior_null)
      l_likelihood<- (0.5 * log(tau_mu)) + l_D0[j] + l_D1[j] + l_D2[j] + l_D3[j] + l_D4[j] +
                          l_D5[j] + l_D6[j] + l_D7[j]
      l_Bayes_factor_01[j]<-  l_likelihood + l_prior_odds
  }
  output<-list(l_D0,l_D1,l_D2,l_D3,l_D4,l_D5,l_D6,l_D7, l_likelihood,
               l_prior_odds, l_Bayes_factor_01)
  names(output)<-c("l_D0","l_D1","l_D2","l_D3","l_D4","l_D5","l_D6","l_D7",
                  "l_likelihood", "l_prior_odds", "l_Bayes_factor_01")
  return(output)
}


#' Title Bayesian Test for dose-dependent differential gene expression analysis
#' @author Satabdi Saha 
#' @param Y : list of size K (nlevels(dose_level)) of n (cells) times p(gene, in our case p=1) 
#' matrices of observed data having n cells and p genes in K treatment dose groups 
#' @param m : vector of group means
#' @param m_0 : grand_mean of m
#' @param tau_k_mu 
#' @param K : number of dose groups
#' @param tau_mu 
#' @param b_sigma 
#' @param a_sigma 
#' @param a_w 
#' @param b_w 
#' @param prior_alt 
#' @param prior_null 
#'
#' @return output: A list having the different parts of the log-likelihood function, log 
#' prior odds and the log Bayes factor test statistic
#' @export
Bayes_factor_multiple<-function(Y, m, m_0, tau_k_mu = rep(1,9), K,
                                tau_mu = c(0.5, 1, 2, 3), b_sigma = 1, a_sigma = 6, a_w = 0.8,
                                b_w = 0.2, prior_alt, prior_null)
{
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
  
  for(j in 1:ncol(Y[[1]]))
  {
    for(k in 1:K)
    {
      ind_D0[k]<- (0.5)*log(1+ (tau_mu*R_colsum[k][j]))
    }
    l_D0[j]<-(1-K) + ((0.5*(1-K))* log (2 * pi)) + colSums(ind_D0) - ((0.5)* log(1+ (tau_mu* R_grand_sum[j]))) +
      (a_sigma + (0.5 *R_grand_sum[j])) * (- log ((1/b_sigma) + (0.5*A_tot[j])) +
                                             (log ((1/b_sigma)+  (0.5* sum(A)[j]))))                
    l_D1[j]<- 0.5*(log(a_w + R_grand_sum[j] - 1) + log(b_w + sum(n) - R_grand_sum[j] - 1) -
                     log(a_w + b_w + sum(n) - 1))
    l_D2[j]<-((a_w + R_grand_sum[j]-1)* log(a_w + R_grand_sum[j]-1)) -
      ((a_w + b_w + sum(n) -1)* log(a_w + b_w + sum(n) -1)) +
      ((b_w + sum(n) - R_grand_sum[j] -1)* log(b_w + sum(n) - R_grand_sum[j] -1))
    for( k in 1: K)
    {
      ind_D3[k]<-0.5 * (log(a_w + b_w + n[k] - 1) - log(a_w + R_colsum[k][j] - 1) -
                          log(b_w + n[k] - R_colsum[k][j] - 1))
      
      ind_D4[k]<- ((a_w + b_w + n[k] - 1)* log(a_w + b_w + n[k] - 1)) -
        ((a_w + R_colsum[k][j] - 1) * log(a_w + R_colsum[k][j] - 1)) -
        ((b_w + n[k] - R_colsum[k][j] - 1)*log(b_w + n[k] - R_colsum[k][j] - 1))
    }
    l_D3[j]<- colSums(ind_D3)
    l_D4[j]<- colSums(ind_D4)
    # l_D5[j]<- (K-1)*0.5* (log(a_w - 1) + log(b_w - 1) - log(a_w + b_w - 1)) +
    #   ((a_w - 1)*(K-1)*log(a_w - 1)) + ((b_w - 1)*(K-1) * log(b_w - 1)) -
    #   ((a_w + b_w - 1)*(K-1)*log(a_w + b_w - 1))
    l_D5[j]<- (K-1)* lbeta(a_w,b_w)
    l_likelihood<-  l_D0[j] + l_D1[j] + l_D2[j] + l_D3[j] + l_D4[j] +
      l_D5[j]
    l_prior_odds<- log (prior_alt) - log(prior_null)
    l_Bayes_factor_01[j]<-  l_likelihood + l_prior_odds
  }
  output<-list(l_D0,l_D1,l_D2,l_D3,l_D4,l_D5, l_likelihood,
               l_prior_odds, l_Bayes_factor_01)
  names(output)<-c("l_D0","l_D1","l_D2","l_D3","l_D4","l_D5",
                   "l_likelihood", "l_prior_odds", "l_Bayes_factor_01")
  
  return(output)
}

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
calc_priors = function(simulated.data, simulated.cell.metadata, simulated.gene.metadata){
  data.list = list()
  priors = data.frame(matrix(NA, nrow = 0, ncol = 5))
  m = vector()
  colnames(priors) = c('dose', 'sigma_mean', 'sigma_var', 'omega_mean', 'omega_var')
  
  for (dose in sort(unique(simulated.cell.metadata$Dose))){
    data.list[[dose]] = simulated.data[,which(simulated.cell.metadata$Dose == dose)]
    
    m[dose]<-mean(rowMeans(replace(data.list[[dose]], data.list[[dose]] == 0, NA), na.rm = TRUE),na.rm = TRUE)
    
    sigma_mean = mean(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE), na.rm = TRUE)
    sigma_var = var(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE),na.rm = TRUE)
    #Dose group specific mean dropout proportions for all genes
    omega_mean = mean(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    omega_var = var(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    #priors = rbind(priors, 
    #      c(dose, sigma_mean, sigma_var, omega_mean, omega_var))
    priors[dose,] = c(dose, sigma_mean, sigma_var, omega_mean, omega_var)
  }
  return(list(split.simulated = data.list, priors = priors, m = m))
}

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
bf_01 = function(data.list, m, tau_k_mu, tau_mu, prior_Alter, prior_Null){
  #bf_multiple_01<-rep(list(list()), nrow(data.list[["0"]]))
  bf_multiple_01 = list()
  for(j in rownames(data.list[["0"]])){
    bf_multiple_01[[j]]<-Bayes_factor_multiple(
      Y = list(Y_1=as.matrix(data.list[["0"]][j,]), 
               Y_2=as.matrix(data.list[["0.01"]][j,]),
               Y_3=as.matrix(data.list[["0.03"]][j,]),
               Y_4=as.matrix(data.list[["0.1"]][j,]),
               Y_5=as.matrix(data.list[["0.3"]][j,]),
               Y_6=as.matrix(data.list[["1"]][j,]),
               Y_7=as.matrix(data.list[["3"]][j,]),
               Y_8=as.matrix(data.list[["10"]][j,]),
               Y_9=as.matrix(data.list[["30"]][j,])),
      m_0=mean(m),m=m,K=9,
      tau_k_mu =tau_k_mu,tau_mu=tau_mu[2], 
      b_sigma=1,a_sigma=6,a_w=8,b_w=2,
      prior_alt =prior_Alter[4],
      prior_null = prior_Null[4])
  }
  bf_multiple_01 = do.call(rbind, lapply(bf_multiple_01, as.data.frame))
  bf_multiple_01$omega = bf_multiple_01$l_D1 + bf_multiple_01$l_D2 + bf_multiple_01$l_D3 + bf_multiple_01$l_D4 + bf_multiple_01$l_D4
  bf_multiple_01$exp_bf = exp(bf_multiple_01$l_Bayes_factor_01)
  
  return(bf_multiple_01)
}

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
bf_pair = function(data.list, bf_1, m, tau_k_mu, tau_mu, prior_Alter, prior_Null){
  #significant.genes = rownames(bf_1[which(bf_1$bf.t.0.1 == 'Positive'),])
  significant.genes = rownames(bf_1)
  bf_actual_01 = list()
  for (dose in names(data.list)[-1]){
    message(paste("performing pair-wise comparisons for dose ", dose, sep = ''))
    bf_actual_01.genes = list()
    for (g in significant.genes){
      bf_actual_01.genes[[g]] = Bayes_factor_multiple(Y = list(
        Y_1=as.matrix(data.list[["0"]][g,]), 
        Y_2=as.matrix(data.list[[dose]][g,])
      ),
      m_0=mean(m["0"] + m[dose]), 
      m=c(m["0"],m[dose]), 
      K=2,
      tau_k_mu =c(tau_k_mu["0"], 
                  tau_k_mu[dose]),
      tau_mu=tau_mu[2], 
      b_sigma=1,
      a_sigma=6,
      a_w=8,
      b_w=2,
      prior_alt =prior_Alter[4],
      prior_null = prior_Null[4]
      )
      bf_actual_01.genes[[g]]$Dose = dose
      bf_actual_01.genes[[g]]$Gene = g
    }
    bf_actual_01[[dose]] = do.call(rbind, lapply(bf_actual_01.genes, as.data.frame))
  }
  bf_pair = do.call(rbind, lapply(bf_actual_01, as.data.frame))
  bf_pair$exp_bf = exp(bf_pair$l_Bayes_factor_01)
  return(bf_pair)
}


#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
sceCalcPriors = function(sce){
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
    
    m[dose]<-mean(rowMeans(replace(data.list[[dose]], data.list[[dose]] == 0, NA), na.rm = TRUE),na.rm = TRUE)
    
    sigma_mean = mean(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE), na.rm = TRUE)
    sigma_var = var(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE),na.rm = TRUE)
    #Dose group specific mean dropout proportions for all genes
    omega_mean = mean(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    omega_var = var(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    priors[dose,] = c(dose, sigma_mean, sigma_var, omega_mean, omega_var)
  }
  return(list(split.simulated = data.list, priors = priors, m = m))
}

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
bayesDETest = function(data.list, m, tau_k_mu, tau_mu, prior_Alter, prior_Null){
  bf_multiple_01 = list()
  for(j in rownames(data.list[[1]])){
    in.list = data.list %>% map(as.matrix(~.x[j,]))
    names(in.list) = paste0("Y_", 1:length(in.list))
    bf_multiple_01[[j]]<-Bayes_factor_multiple(
      Y = in.list,
      m_0=mean(m),m=m,K=length(in.list),
      tau_k_mu =tau_k_mu,tau_mu=tau_mu[2], 
      b_sigma=1,a_sigma=6,a_w=8,b_w=2,
      prior_alt =prior_Alter[4],
      prior_null = prior_Null[4])
  }
  bf_multiple_01 = do.call(rbind, lapply(bf_multiple_01, as.data.frame))
  bf_multiple_01$omega = bf_multiple_01$l_D1 + bf_multiple_01$l_D2 + bf_multiple_01$l_D3 + bf_multiple_01$l_D4 + bf_multiple_01$l_D4
  bf_multiple_01$exp_bf = exp(bf_multiple_01$l_Bayes_factor_01)
  
  return(bf_multiple_01)
}
