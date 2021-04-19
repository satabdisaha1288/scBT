# Author: Jack Dodson

#' Performs a genewise Wilcoxon Rank Sum test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the Wilcoxon Rank Sum test
#' 
#' @example 
#' 
#' @importFrom SingleCellExperiment logcounts colData
#' @export
batchWRS = function(sce){
  data = as.matrix(logcounts(sce))
  dose = colData(sce)$Dose
  wrs.pvalues = apply(data, 1, function(x) runWRS(x, dose))
  return(wrs.pvalues)
}

#' Performs a Wilcoxon Rank Sum test on a logcounts vactor for a given dose
#'
#' @param data The logcounts vector
#' @param dose The dose to analyze
#'
#' @return A p value from the Wilcoxon Rank Sum test
#'
#' @example
#'
#' @importFrom
#' @export
runWRS = function(data, dose){
  my_data <- data.frame(value = data, dose = dose)
  res.wrs = wilcox.test(value ~ dose, data = my_data)
  wrs.pvalue = res.wrs[[3]]
  return(wrs.pvalue)
}

##########################################################################

# Author: Satabdi Saha

#' Performed genewise Wilcoxon Rank Sum test on a cell x row matrix
#' 
#' @param data matrix where cells are rows and genes are columns
#' 
#' @return
#' 
#' @example
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