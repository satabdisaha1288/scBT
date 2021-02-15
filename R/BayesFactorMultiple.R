#Calculation of Bayes factor
#Y : list of size K of n times p matrices of observed data having n cells and p genes in K treatment dose groups 
#in the treatment group
#K is the number of groups
#m is the vector of group means
Bayes_factor_multiple<-function(Y,m,m_0,tau_k_mu,K,
                       tau_mu, b_sigma,a_sigma,a_w,
                       b_w,prior_alt,prior_null)
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
  n<- sapply(Y, function(x) nrow(x))
  p<- sapply(Y, function(x) ncol(x))
  R <- lapply(Y, function(x) ifelse(x>0,1,0))
  R_colsum<-sapply(R, function(x) colSums(x))
  R_grand_sum<-sum(R_colsum)
  Y_squared<- Map("*",Y,Y)
  R_y_2<-Map("*" ,R , Y_squared)
  R_y_2_colsum<-sapply(R_y_2, function(x) colSums(x))
  R_y<-Map("*" ,R , Y)
  R_y_colsum<-sapply(R_y, function(x) colSums(x))
  R_y_grand_colsum<-sum(R_y_colsum)
  R_y_2_grand_colsum<-sum(R_y_2_colsum)
  A<- R_y_2_colsum + ((m^2)/tau_k_mu) - 
       (((R_y_colsum + (m/tau_k_mu))^2)/ (R_colsum + (1/tau_k_mu)))
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
