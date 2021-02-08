#Calculation of Bayes factor
#Y_t : n times p matrix of observed data having n cells and p genes 
#in the treatment group
#Y_c : n times p matrix of observed data having n cells and p genes 
#in the control group
Bayes_factor<-function(Y_t,Y_c,m_t0,m_c0,m_0,tau_t_mu,tau_c_mu,
                       tau_mu, b_sigma,a_sigma,a_w,a_c_w,a_t_w,
                       b_w,b_t_w,b_c_w,prior_alt,prior_null)
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
