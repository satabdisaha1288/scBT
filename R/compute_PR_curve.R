#' Title Computes PR curve for DEG classification
#'
#' @param sim
#' @param DETestoutput
#' @author Satabdi Saha
#' @return PR plot
#' @export
#'
#' @examples
compute_PR_curve<-function(sim,DETestoutput){
  require(PRROC)
  true_model<-rowData(sim)[,"Model"]
  true_model<-ifelse(true_model== "Unchanged",0,1)
  fg <- DETestoutput$KW$kw.pvalues[true_model == 1]
  bg <- DETestoutput$KW$kw.pvalues[true_model == 0]
  x_KW<-c(fg,bg)
  lab<-c(rep(1,length(fg)),rep(0,length(bg)))
  x_aov<-c(DETestoutput$ANOVA$aov.pvalues[true_model == 1],
           DETestoutput$ANOVA$aov.pvalues[true_model == 0])
  #ROC Curve
  wPR_KW<-pr.curve(scores.class0 = x_KW, weights.class0 = lab, curve = TRUE,
                     max.compute = T, min.compute = T, rand.compute = T)
  wPR_aov<-pr.curve(scores.class0 = x_aov, weights.class0 = lab, curve = TRUE,
                      max.compute = T, min.compute = T, rand.compute = T)
  wPR_plot<-plot(wPR_KW,max.plot = TRUE, min.plot = TRUE,
                  rand.plot = TRUE, fill.area = TRUE,color = 2,
                  scale.color = heat.colors(100))
  wPR_plot<-plot(wPR_aov, add = TRUE, color = 3);
  return(c(wPR_KW,wPR_aov,wPR_plot))
}
