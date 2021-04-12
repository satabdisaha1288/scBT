#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @example 
#' 
#' @importFrom SingleCellExperiment logcounts colData
#' @export
batchANOVA = function(sce){
  data = as.matrix(logcounts(sce))
  dose = colData(sce)$Dose
  aov.pvalues = apply(data, 1, function(x) runAnova(x, dose))
  return(aov.pvalues)
}

#' Performs a ANOVA test on a logcounts vactor for a given dose
#' 
#' @param data The logcounts vector
#' @param dose The dose to analyze
#' 
#' @return A p value from the ANOVA test
#' 
#' @example 
#' 
#' @importFrom
#' @export
runAnova = function(data, dose){
  my_data = data.frame(value = data, dose = dose)
  res.aov = aov(value ~ dose, data = my_data)
  aov.pvalue = summary(res.aov)[[1]][["Pr(>F)"]][1]
  return(aov.pvalue)
}


#' Performed genewise ANOVA test on a cell x row matrix
#' 
#' @param data matrix where cells are rows and genes are columns
#' 
#' @return
#' 
#' @example 
#' 
#' @export
#Function for ANOVA test for each gene 
#data is a dataframe of n cells as rows and p+1 columns with first p columns as p genes and the last 
#column is the dose information
anova_test<-function(data){
  my_data_anova<-rep(list(list()),ncol(data)-1)
  my_data_anova_p_value<-vector()
  for(i in 1: (ncol(data)-1))
  {
    my_data<-data.frame(data[,i],as.factor(data[,ncol(data)]))
    colnames(my_data)<-c("value","dose")
    res.aov <- aov(value ~ dose, data = my_data)
    my_data_anova[[i]]<-summary.aov(res.aov)
    my_data_anova_p_value[i]<-my_data_anova[[i]][[1]][["Pr(>F)"]][1]
  }
  return(my_data_anova_p_value)
}

