#' Calculate concordance AUC Refer to: Soneson, Charlotte, and Mark D. Robinson.
#' "Bias, robustness and scalability in single-cell differential
#' expression analysis." Nature methods 15, no. 4 (2018): 255.
#'
#' @param DETestoutput
#' @author Satabdi Saha
#' @return Normed area under the concordance curve.
#' @export
#'
#' @examples
calculate_concordance_AUC<-function(DETestoutput)
{
  kw.data<-DETestoutput$KW
  kw.sort<-kw.data[order( kw.data[,1] ),,drop=FALSE]
  aov.data<-DETestoutput$ANOVA
  aov.sort<-aov.data[order( aov.data[,1] ),,drop=FALSE]
  h<-seq(1,100,1)
  common_dge<-vector(length = length(h))
  for(i in 1:length(h))
  {
    common_dge[i]<-length(intersect(rownames(kw.sort)[1:h[i]],
                                    rownames(aov.sort)[1:h[i]]))
  }
  library(zoo)
  id <- order(h)
  AUC <- sum(diff(h[id])*rollmean(common_dge[id],2))
  AUC_norm<-AUC/(h[length(h)]^2 /2)
  return(AUC_norm)
}
