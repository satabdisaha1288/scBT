#' Calculate a_sigma, b_sigma
#' @author Satabdi Saha
#' @param sce Single Cell Object
#' @return a_sigma and b_sigma values
#
#' @export
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

#' Calculate a_w, b_w
#' @author Satabdi Saha
#' @param omega_mean Groupwise dropout means
#' @param omega_var Groupwise dropout variance
#' @return a_w and b_w values
#
#' @export
calc_a_w_b_w<-function(omega_mean,omega_var){
  mean_omega_mean<-mean(omega_mean)
  mean_omega_var<-mean(omega_var)
  a_w = round((((mean_omega_mean^2)*(1-mean_omega_mean)) / mean_omega_var) - mean_omega_mean ,2)
  b_w = round(a_w * ((1- mean_omega_mean)/mean_omega_mean) ,2 )
  output<-c(a_w,b_w)
  names(output)<-c("a_w","b_w")
  return(output)
}

#' Calculate a_w, b_w
#' @author Satabdi Saha
#' @param omega_mean Groupwise dropout means
#' @param omega_var Groupwise dropout variance
#' @return a_w and b_w values
#
#' @export
calc_alt_null = function(sce, method = 'fixed', de.prob = 0.25){
  if (method == 'fixed'){
    de.prob = 0.25
  } else {
    ##Calcualte
  }
  prior_Alter<-c(1-((1-de.prob)^(1/nrow(sce))))
  prior_Null<- 1-prior_Alter
  p_alt = c(prior_Alter = prior_Alter, prior_Null = prior_Null)
  return(p_alt)
}
