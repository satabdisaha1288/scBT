#' Chi-square test of independence for binned gene expression data
#' @author Satabdi Saha
#' @param sce A Single Cell Object
#' @return FDR adjusted p-values
#
#' @export
runChi <- function(sce){
  # TODO: replace hardcoded values.
  data_logcounts <- as.matrix(logcounts(sce))
  dose <- colData(sce)$Dose
  data_dose <- data.frame(t(data_logcounts),dose)
  colnames(data_dose) <- c(rownames(data_logcounts),"Dose")
  #Break the dataframe to data lists for seperate dose groups
  data <- split_tibble(data_dose,"Dose")
  data <- lapply(data, function(x) t(x))
  max_value_data <- max(sapply(data, function(x) apply(x,1,max)))
  #Create a data list with genes binned into certain values
  data_bin <- lapply(data, function(x) t(apply(x,1,function(x){
    #Most data exists in these bins - fixed number of bins spanning range of
    #values but needs to be unequal to be more around where most of the data is.
    #Quantiles maybe, but fixed for all dose groups.
    new_breaks <- c(-0.01,0.5,1,2,3,4,ceiling(max_value_data))
    table(cut(x, new_breaks))
  })))
  #Create a list of data frames for applying chi-square test of independence
  df <- rep(list(data.frame()),length=nrow(data_bin[["0"]]))
  names(df) <- nrow(data_bin[["0"]])
  for(j in 1:nrow(data_bin[["0"]]))
  {
    a1 <- rapply(data_bin, classes = 'matrix', how = 'list', f = function(x) x[j, , drop = FALSE])
    # a second `lapply` is required to drop `NULL` entries
    a1Only <- lapply(a1, Filter, f = Negate(is.null))
    df[[j]] <- data.frame(matrix(unlist(a1Only), nrow=length(a1Only), byrow=TRUE))
    rownames(df[[j]]) <- levels(dose)
    colnames(df[[j]]) <- colnames(data_bin[["0.01"]])
  }
  #chisq_test<-lapply(df,function(x) chisq.test(x[,which(colSums(x)>0)]))
  chisq_test <- lapply(df,function(x) chisq.test(x[, which(colSums(x) > 0)],
                                                 simulate.p.value = TRUE))

  names(chisq_test)<-rownames(data[["0"]])
  adjusted.p <- p.adjust(sapply(chisq_test,function(x) x$p.value),method = "fdr")
  chisq_test_pvalue_adj = data.frame(adjusted.p)

  return(adjusted.p)
}

#' Chi-square test of independence for binned gene expression data
#' @author Satabdi Saha
#' @param sce A Single Cell Object
#' @return FDR adjusted p-values
#
#' @export
split_tibble <- function(tibble, column = 'col') {
  library(dplyr)
  tibble %>% split(., .[,column]) %>% lapply(., function(x) x[,setdiff(names(x),column)])
}
