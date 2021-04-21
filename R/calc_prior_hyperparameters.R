#Calculate a_sigma, b_sigma
calc_a_sigma_b_sigma<-function(sigma_mean,sigma_var){
  mean_sigma_mean<-mean(sigma_mean)
  mean_sigma_var<-mean(sigma_var)
  a_sigma<-(1/((mean_sigma_mean^2)*mean_sigma_var)) + 2
  b_sigma<-1/(mean_sigma_mean* (a_sigma-1))
  output<-c(a_sigma,b_sigma)
  names(output)<-c("a_sigma","b_sigma")
  return(output)
}
#Calculate a_w,b_w
calc_a_w_b_w<-function(omega_mean, omega_var){
  mean_omega_mean<-mean(omega_mean)
  mean_omega_var<-mean(omega_var)
  a_w = round((((mean_omega_mean^2)*(1-mean_omega_mean)) / mean_omega_var) - mean_omega_mean ,2)
  b_w = round(a_w * ((1- mean_omega_mean)/mean_omega_mean) ,2 )
  output<-c(a_w,b_w)
  names(output)<-c("a_w","b_w")
  return(output)
}




