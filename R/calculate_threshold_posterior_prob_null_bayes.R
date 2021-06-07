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
calculate_threshold_posterior_prob_null_bayes<-function(DETest_output_Bayes,kappa,alpha){
  #Posterior probabilities under the null
  posterior_prob_null<-1/(1+DETest_output_Bayes)
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


