#' Performed genewise Kruskal Wallis test on a cell x row matrix
#' 
#' @author Satabdi Saha
#' @param data matrix where cells are rows and genes are columns
#' 
#' @return a vector of p values from the Kruskal Wallis test
#' 
#' @export
#Function for carrying out the Kruskal Wallis Test
#data is a dataframe of n cells as rows and p+1 columns with first p columns as p genes and the last 
#column is the dose information
KW_test<-function(data){
  my_data_kruskal_p_value<-vector()
  for(i in 1: (ncol(data)-1))
  {
    my_data<-data.frame(data[,i],as.factor(data[,ncol(data)]))
    colnames(my_data)<-c("value","dose")
    res.kruskal<-kruskal.test(value ~ dose, data = my_data)
    my_data_kruskal_p_value[i]<-res.kruskal$p.value
  }
  return(my_data_kruskal_p_value)
}

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
runKW = function(logcounts, cell.meta){
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

#' Performed genewise Wilcoxon Rank Sum test on a cell x row matrix
#' 
#' @param data matrix where cells are rows and genes are columns
#' 
#' @return
#' 
#' @export
#Function for carrying out the Wicoxon Rank Sum test
#data is a dataframe of n cells as rows and p+1 columns with first p columns as p genes and the last 
#column is the dose information
WRS_test<-function(data,dose_level_ref,dose_level_test){
  my_data_wilcoxon_p_value<-vector()
  for(i in 1: (ncol(data)-1))
  {
    res.wrs<-wilcox.test(data[,i][data[,ncol(data)]== dose_level_ref]
                         ,data[,i][data[,ncol(data)]== dose_level_test])
    my_data_wilcoxon_p_value[i]<-res.wrs$p.value
  }
  return(my_data_wilcoxon_p_value)
}

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
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

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' @export
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

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
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

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
runANOVA = function(sce){
  logcounts = data.frame(t(as.matrix(logcounts(sce))))
  
  
  #logcounts, cell.meta
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

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
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

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
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

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
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