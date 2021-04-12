# Author: Jack Dodson

#' Performs a genewise Kruskal Wallis test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the Kruskal Wallis test
#' 
#' @example 
#' 
#' @importFrom SingleCellExperiment logcounts colData
#' @export
batchKW = function(sce){
  data = as.matrix(logcounts(sce))
  dose = colData(sce)$Dose
  kw.pvalues = apply(data, 1, function(x) runKW(x, dose))
  return(kw.pvalues)
}

#' Performs a Kruskal Wallis test on a logcounts vactor for a given dose
#' 
#' @param data The logcounts vector
#' @param dose The dose to analyze
#' 
#' @return A p value from the Kruskal Wallis test
#' 
#' @example 
#' 
#' @importFrom
#' @export
runKW = function(data, dose){
  my_data <- data.frame(value = data, dose = dose)
  res.kw = kruskal.test(value ~ dose, data = my_data)
  kw.pvalue = summary(res.kw)[[1]][["Pr(>F)"]][1]
  return(kw.pvalue)
}

##########################################################################

# Author: Satabdi Saha

#' Performed genewise Kruskal Wallis test on a cell x row matrix
#' 
#' @param data matrix where cells are rows and genes are columns
#' 
#' @return
#' 
#' @example 
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

test_result <- KW_test(KW_input)
