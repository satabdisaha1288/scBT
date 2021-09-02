#' Title Computes ROC curve for DEG classification
#'
#' @param sim
#' @param DETestoutput
#' @author Satabdi Saha
#' @return ROC plot
#' @export
#'
#' @examples
compute_ROC_curve<-function(sim,DETestoutput){
  require(PRROC)
  true_model<-rowData(sim)[,"Model"]
  true_model<-ifelse(true_model== "Unchanged",0,1)
  wfg<- c(runif(300,min=0.5,max=1),runif(500,min=0,max=0.5))
  fg <- DETestoutput$KW$kw.pvalues[true_model == 1]
  bg <- DETestoutput$KW$kw.pvalues[true_model == 0]
  x_KW<-c(fg,bg)
  lab<-c(rep(1,length(fg)),rep(0,length(bg)))
  x_aov<-c(DETestoutput$ANOVA$aov.pvalues[true_model == 1],
       DETestoutput$ANOVA$aov.pvalues[true_model == 0])
  #ROC Curve
  wroc_KW<-roc.curve(scores.class0 = x_KW, weights.class0 = lab, curve = TRUE,
                  max.compute = T, min.compute = T, rand.compute = T)
  wroc_aov<-roc.curve(scores.class0 = x_aov, weights.class0 = lab, curve = TRUE,
                     max.compute = T, min.compute = T, rand.compute = T)
  wroc_plot<-plot(wroc_KW,max.plot = TRUE, min.plot = TRUE,
                  rand.plot = TRUE, fill.area = TRUE,color = 2,
                  scale.color = heat.colors(100))
  wroc_plot<-plot(wroc_aov, add = TRUE, color = 3);
  return(c(wroc_KW,wroc_aov,wroc_plot))
}



