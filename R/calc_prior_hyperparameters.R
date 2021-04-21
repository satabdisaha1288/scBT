#Calculate prior hyper-parameters

#Calculate a_sigma, b_sigma

calc_a_sigma_b_sigma<-function(sce){
  data<-logcounts(sce)
  sigma<- apply(data, 1, var)
  mean_sigma<-mean(sigma[sigma>0])
  var_sigma<-sd(sigma[sigma>0])^2
  a_sigma<-(1/((mean_sigma^2)*var_sigma)) + 2
  b_sigma<-1/(mean_sigma* (a_sigma-1))
  output<-c(a_sigma,b_sigma)
  names(output)<-c("a_sigma","b_sigma")
  return(output)
}
#Calculate a_w,b_w
calc_a_w_b_w<-function(omega_mean,omega_var)
  mean_omega_mean<-mean(omega_mean)
  mean_omega_var<-mean(omega_var)
  a_w = round((((mean_omega_mean^2)*(1-mean_omega_mean)) / mean_omega_var) - mean_omega_mean ,2)
  b_w = round(a_w * ((1- mean_omega_mean)/mean_omega_mean) ,2 )
  output<-c(a_w,b_w)
  names(output)<-c("a_w","b_w")
  return(output)
}




