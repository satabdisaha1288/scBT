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