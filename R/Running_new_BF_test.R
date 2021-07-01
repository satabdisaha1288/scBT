#Y is a Data list of K (number of dose groups) matrices having dimension n (cells) times p (genes)
# optim_output_null <- optim(par = runif(6, -1.25, 1.25),   # Applying optim
#                            fn = log.lklh.marginal.null,
#                            control = list(maxit=1000), method="L-BFGS-B", upper=rep(5, 6),
#                            lower=rep(-5, 6))
# optim_output_alt <- optim(par = runif(6, -1.25, 1.25),   # Applying optim
#                           fn = log.lklh.marginal.alt,
#                           control = list(maxit=1000), method="L-BFGS-B", upper=rep(5, 6),
#                           lower=rep(-5, 6))
# opt_par<- function(x)
# {
#   opt_par<-c(x[1],exp(x[-1]))
#   return(opt_par)
# }
# par_null<-opt_par(optim_output_null$par)
# par_alt<-opt_par(optim_output_alt$par)
# prior.null=0.75
# prior.alt=0.25
# bayes.factor<-calculate_BF(Y, par_null , par_alt, prior.null, prior.alt)
# thresh<-calculate_threshold_posterior_prob_null_bayes(bayes.factor$BF,seq(0.01,1,0.01),0.05)
# posterior_prob_null<-1/(1+ (1/bayes.factor$BF))
# #posterior_prob_null<-1/(1+ (1/DGEAnalysis$BAYES$exp_bf))
# bayes.factor.est<-factor(ifelse(posterior_prob_null <= thresh , "Pos","Neg"),
#                          levels = c("Pos", "Neg"))
# true_fac<-factor(ifelse(rowData(sim)$Model == "Unchanged","Neg","Pos"),
#                  levels = c("Pos", "Neg"))
# caret::confusionMatrix(bayes.factor.est,true_fac)
