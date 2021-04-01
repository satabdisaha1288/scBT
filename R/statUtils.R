calc_priors = function(simulated.data, simulated.cell.metadata, simulated.gene.metadata){
  data.list = list()
  priors = data.frame(matrix(NA, nrow = 0, ncol = 5))
  m = vector()
  colnames(priors) = c('dose', 'sigma_mean', 'sigma_var', 'omega_mean', 'omega_var')
  
  for (dose in sort(unique(simulated.cell.metadata$Dose))){
    data.list[[dose]] = simulated.data[,which(simulated.cell.metadata$Dose == dose)]
    
    m[dose]<-mean(rowMeans(replace(data.list[[dose]], data.list[[dose]] == 0, NA), na.rm = TRUE),na.rm = TRUE)
    
    sigma_mean = mean(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE), na.rm = TRUE)
    sigma_var = var(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE),na.rm = TRUE)
    #Dose group specific mean dropout proportions for all genes
    omega_mean = mean(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    omega_var = var(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    #priors = rbind(priors, 
    #      c(dose, sigma_mean, sigma_var, omega_mean, omega_var))
    priors[dose,] = c(dose, sigma_mean, sigma_var, omega_mean, omega_var)
  }
  return(list(split.simulated = data.list, priors = priors, m = m))
}

bf_01 = function(data.list, m, tau_k_mu, tau_mu, prior_Alter, prior_Null){
  #bf_multiple_01<-rep(list(list()), nrow(data.list[["0"]]))
  bf_multiple_01 = list()
  for(j in rownames(data.list[["0"]])){
    bf_multiple_01[[j]]<-Bayes_factor_multiple(
      Y = list(Y_1=as.matrix(data.list[["0"]][j,]), 
               Y_2=as.matrix(data.list[["0.01"]][j,]),
               Y_3=as.matrix(data.list[["0.03"]][j,]),
               Y_4=as.matrix(data.list[["0.1"]][j,]),
               Y_5=as.matrix(data.list[["0.3"]][j,]),
               Y_6=as.matrix(data.list[["1"]][j,]),
               Y_7=as.matrix(data.list[["3"]][j,]),
               Y_8=as.matrix(data.list[["10"]][j,]),
               Y_9=as.matrix(data.list[["30"]][j,])),
      m_0=mean(m),m=m,K=9,
      tau_k_mu =tau_k_mu,tau_mu=tau_mu[2], 
      b_sigma=1,a_sigma=6,a_w=8,b_w=2,
      prior_alt =prior_Alter[4],
      prior_null = prior_Null[4])
  }
  bf_multiple_01 = do.call(rbind, lapply(bf_multiple_01, as.data.frame))
  bf_multiple_01$omega = bf_multiple_01$l_D1 + bf_multiple_01$l_D2 + bf_multiple_01$l_D3 + bf_multiple_01$l_D4 + bf_multiple_01$l_D4
  bf_multiple_01$exp_bf = exp(bf_multiple_01$l_Bayes_factor_01)
  
  return(bf_multiple_01)
}

bf_pair = function(data.list, bf_1, m, tau_k_mu, tau_mu, prior_Alter, prior_Null){
  #significant.genes = rownames(bf_1[which(bf_1$bf.t.0.1 == 'Positive'),])
  significant.genes = rownames(bf_1)
  bf_actual_01 = list()
  for (dose in names(data.list)[-1]){
    message(paste("performing pair-wise comparisons for dose ", dose, sep = ''))
    bf_actual_01.genes = list()
    for (g in significant.genes){
      bf_actual_01.genes[[g]] = Bayes_factor_multiple(Y = list(
        Y_1=as.matrix(data.list[["0"]][g,]), 
        Y_2=as.matrix(data.list[[dose]][g,])
      ),
      m_0=mean(m["0"] + m[dose]), 
      m=c(m["0"],m[dose]), 
      K=2,
      tau_k_mu =c(tau_k_mu["0"], 
                  tau_k_mu[dose]),
      tau_mu=tau_mu[2], 
      b_sigma=1,
      a_sigma=6,
      a_w=8,
      b_w=2,
      prior_alt =prior_Alter[4],
      prior_null = prior_Null[4]
      )
      bf_actual_01.genes[[g]]$Dose = dose
      bf_actual_01.genes[[g]]$Gene = g
    }
    bf_actual_01[[dose]] = do.call(rbind, lapply(bf_actual_01.genes, as.data.frame))
  }
  bf_pair = do.call(rbind, lapply(bf_actual_01, as.data.frame))
  bf_pair$exp_bf = exp(bf_pair$l_Bayes_factor_01)
  return(bf_pair)
}

runANOVA_KW = function(logcounts, cell.meta){
  simulated.data.transposed = data.frame(t(as.matrix(logcounts)))
  simulated.data.transposed$Dose = as.numeric(cell.meta$Dose)
  kw.out = KW_test(simulated.data.transposed)
  kw.out = p.adjust(kw.out, 'fdr')
  anova.out = anova_test(simulated.data.transposed)
  anova.out = p.adjust(anova.out, 'fdr')
  
  ANOVA_KW.df = data.frame(KW.fdr = kw.out, anova.fdr = anova.out)
  rownames(ANOVA_KW.df) = colnames(simulated.data.transposed)[1:length(anova.out)]
  return(ANOVA_KW.df)
}

runWRS = function(logcounts, cell.meta){
  #simulated.data.transposed = data.frame(t(as.matrix(logcounts)))
  #simulated.data.transposed$Dose = as.numeric(cell.meta$Dose)
  WRS = data.frame(Gene = rownames(logcounts))
  rownames(WRS) = WRS$Gene
  for (dose in sort(unique(cell.meta$Dose))[-1]){
    message(paste('Working on dose ', dose, sep = ''))
    WRS[,paste('wrs.p.', dose, sep = '')] = NA
    
    for (g in rownames(logcounts)){
      WRS[g, paste('wrs.p.', dose, sep = '')] = wilcox.test(logcounts[g, which(cell.meta$Dose == 0)],
                                                            logcounts[g, which(cell.meta$Dose == dose)])$p.value
    }
    WRS[,paste('wrs.fdr.', dose, sep = '')] = p.adjust(WRS[ , paste('wrs.p.', dose, sep = '')], method = 'fdr')
  }
  message('Wilcox Rank Sum test complete')
  return(WRS)
}

runLRT = function(data.list){
  for (dose in names(data.list)){
    data.list[[dose]] = as.matrix(data.list[[dose]])
  }
  #LRT_mult<-rep(list(list()),nrow(data[["0"]]))
  LRT_mult = list()
  data_ind = lapply(data.list, function(x) ifelse(x > 0, 1, 0))
  for(g in 1: nrow(data.list[["0"]])){
    gene = rownames(data.list[["0"]])[g]
    LRT_mult[[gene]]<-LRT_multiple_groups(data = list(data.list[["0"]][g,][data_ind[["0"]][g,] > 0],
                                                      data.list[["0.01"]][g,][data_ind[["0.01"]][g,] > 0],
                                                      data.list[["0.03"]][g,][data_ind[["0.03"]][g,] > 0],
                                                      data.list[["0.1"]][g,][data_ind[["0.1"]][g,] > 0],
                                                      data.list[["0.3"]][g,][data_ind[["0.3"]][g,] > 0],
                                                      data.list[["1"]][g,][data_ind[["1"]][g,] > 0],
                                                      data.list[["3"]][g,][data_ind[["3"]][g,] > 0],
                                                      data.list[["10"]][g,][data_ind[["10"]][g,] > 0],
                                                      data.list[["30"]][g,][data_ind[["30"]][g,] > 0]),
                                          data_ind = list(data_ind[["0"]][g,],
                                                          data_ind[["0.01"]][g,],
                                                          data_ind[["0.03"]][g,],
                                                          data_ind[["0.1"]][g,],
                                                          data_ind[["0.3"]][g,],
                                                          data_ind[["1"]][g,],
                                                          data_ind[["3"]][g,],
                                                          data_ind[["10"]][g,],
                                                          data_ind[["30"]][g,]))
  }
  LRT_mult_p_value = sapply(LRT_mult,function(x) x[2,3])
  LRT_mult_p_value_adj = p.adjust(LRT_mult_p_value, method = "fdr")
  LRT.out = data.frame(pval = LRT_mult_p_value, pval.fdr = LRT_mult_p_value_adj)
  rownames(LRT.out) = names(LRT_mult_p_value_adj)
  return(LRT.out)
}


calcFC = function(sim){
  dose_vec = sort(unique(colData(sim)$Dose))
  m0 = rowMeans(as.matrix(logcounts(sim)[,which(colData(sim)$Dose == 0)]))
  fc = list()
  for (dose in dose_vec[-1]){
    temp.means = rowMeans(as.matrix(logcounts(sim)[,which(colData(sim)$Dose == dose)]))
    fc[[dose]] = temp.means-m0
  }
  fc.out = do.call(cbind, lapply(fc, as.data.frame))
  colnames(fc.out) = paste0('calculatedFC',names(fc))
  return(fc.out)
}

calcZeroP = function(sim){
  dose_vec = sort(unique(colData(sim)$Dose))
  pz = list()
  for (dose in dose_vec){
    temp.df = as.matrix(logcounts(sim)[,which(colData(sim)$Dose == dose)])
    pz[[dose]] = apply(temp.df, 1, function(x) sum(x == 0)/length(x))
  }
  pz.out = do.call(cbind, lapply(pz, as.data.frame))
  colnames(pz.out) = paste0('percent.zero',names(pz))
  return(pz.out)
}


runMAST = function(sce){
  scaRaw = FromMatrix(as.matrix(logcounts(sce)),
                      data.frame(colData(sce)),
                      data.frame(rowData(sce))
  )
  cdr <-colSums(assay(scaRaw) > 0)
  colData(scaRaw)$cngeneson = scale(cdr)
  filterCrit = with(colData(scaRaw),cngeneson > 1)
  rowData(scaRaw)$symbolid = rownames(rowData(sce))
  sca = subset(scaRaw,filterCrit)
  Dose = factor(colData(scaRaw)$Dose)
  Dose = relevel(Dose,"0")
  colData(scaRaw)$Dose = Dose
  zlmDose = zlm(~Dose + cngeneson, scaRaw)
  
  summaryDose<-rep(list(list()),times= nlevels(Dose)-1)
  names(summaryDose)<-levels(Dose)[2:9]
  for(i in levels(Dose)[2:9])
  {
    summaryDose[[i]]<-summary(zlmDose, doLRT=paste0("Dose",i))
  }
  
  summaryDT<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(summaryDT)<-levels(Dose)[2:9]
  fcHurdle<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(fcHurdle)<-levels(Dose)[2:9]
  fcHurdleSig<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(fcHurdleSig)<-levels(Dose)[2:9]
  fcHurdleNonSig<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(fcHurdleSig)<-levels(Dose)[2:9]
  TP_dose<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(TP_dose)<-levels(Dose)[2:9]
  FP_dose<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(FP_dose)<-levels(Dose)[2:9]
  TN_dose<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(TN_dose)<-levels(Dose)[2:9]
  FN_dose<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(FN_dose)<-levels(Dose)[2:9]
  
  for(i in levels(Dose)[2:9]){
    summaryDT[[i]] <-summaryDose[[i]]$datatable
    fcHurdle[[i]]<-summaryDT[[i]][contrast== paste0("Dose",i) & component=='H',.(primerid, `Pr(>Chisq)`)] #logFC coefficients
    fcHurdle[[i]][,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdleSig[[i]]<- fcHurdle[[i]][fdr<.05 ]
    setorder(fcHurdleSig[[i]], fdr)
    fcHurdleNonSig[[i]]<- fcHurdle[[i]][fdr >.05| fdr == 0.05]
    setorder(fcHurdleNonSig[[i]], fdr)
  }
  
  mast.out = do.call(cbind, lapply(fcHurdle, data.frame))
  
  if (sum(transform(mast.out[,grepl('primer', colnames(mast.out))], same = apply(mast.out[,grepl('primer', colnames(mast.out))], 1, function(x) length(unique(x)) == 1))$same == FALSE) == 0) {
    m = mast.out[,grepl('fdr', colnames(mast.out))]
    rownames(m) = mast.out[,1]
  } else {
    message('ERROR!!!!!!')
  }
  return(m)
}


sceCalcPriors = function(sce){
  if (class(sce) != "SingleCellExperiment"){
    message('ERROR: Input is not a SingleCellExperiment. Stopping...')
    break
  }
  
  #Initialize list, tables, and vectors
  data.list = list()
  m = vector()
  priors = data.frame(matrix(NA, nrow = 0, ncol = 5))
  colnames(priors) = c('dose', 'sigma_mean', 'sigma_var', 'omega_mean', 'omega_var')
  
  data = as.matrix(logcounts(sce))
  cell.metadata = rowData(sce)
  gene.metadata = colData(sce)
  dose_vec = sort(unique(gene.metadata$Dose))
  
  for (dose in dose_vec){
    data.list[[dose]] = data[,which(cell.metadata$Dose == dose)]
    
    m[dose]<-mean(rowMeans(replace(data.list[[dose]], data.list[[dose]] == 0, NA), na.rm = TRUE),na.rm = TRUE)
    
    sigma_mean = mean(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE), na.rm = TRUE)
    sigma_var = var(rowVars(replace(as.matrix(data.list[[dose]]), as.matrix(data.list[[dose]]) == 0, NA), na.rm = TRUE),na.rm = TRUE)
    #Dose group specific mean dropout proportions for all genes
    omega_mean = mean(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    omega_var = var(apply(as.matrix(data.list[[dose]]), 1, function(x) length(which(x==0))/length(x)))
    #priors = rbind(priors, 
    #      c(dose, sigma_mean, sigma_var, omega_mean, omega_var))
    priors[dose,] = c(dose, sigma_mean, sigma_var, omega_mean, omega_var)
  }
  return(list(split.simulated = data.list, priors = priors, m = m))
}

bf_01 = function(data.list, m, tau_k_mu, tau_mu, prior_Alter, prior_Null){
  #bf_multiple_01<-rep(list(list()), nrow(data.list[["0"]]))
  bf_multiple_01 = list()
  for(j in rownames(data.list[["0"]])){
    bf_multiple_01[[j]]<-Bayes_factor_multiple(
      Y = list(Y_1=as.matrix(data.list[["0"]][j,]), 
               Y_2=as.matrix(data.list[["0.01"]][j,]),
               Y_3=as.matrix(data.list[["0.03"]][j,]),
               Y_4=as.matrix(data.list[["0.1"]][j,]),
               Y_5=as.matrix(data.list[["0.3"]][j,]),
               Y_6=as.matrix(data.list[["1"]][j,]),
               Y_7=as.matrix(data.list[["3"]][j,]),
               Y_8=as.matrix(data.list[["10"]][j,]),
               Y_9=as.matrix(data.list[["30"]][j,])),
      m_0=mean(m),m=m,K=9,
      tau_k_mu =tau_k_mu,tau_mu=tau_mu[2], 
      b_sigma=1,a_sigma=6,a_w=8,b_w=2,
      prior_alt =prior_Alter[4],
      prior_null = prior_Null[4])
  }
  bf_multiple_01 = do.call(rbind, lapply(bf_multiple_01, as.data.frame))
  bf_multiple_01$omega = bf_multiple_01$l_D1 + bf_multiple_01$l_D2 + bf_multiple_01$l_D3 + bf_multiple_01$l_D4 + bf_multiple_01$l_D4
  bf_multiple_01$exp_bf = exp(bf_multiple_01$l_Bayes_factor_01)
  
  return(bf_multiple_01)
}
