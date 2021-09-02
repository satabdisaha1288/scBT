#' Performs a genewise ANOVA test on a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object with a logcounts assay
#' and "Dose" factor column in the cell metadata
#'
#' @return a vector of p values from the ANOVA test
#'
#' @export
batchANOVA <- function(sce){
  data <- as.matrix(logcounts(sce))
  dose <- colData(sce)$Dose
  aov.pvalues <- apply(data, 1, function(x) runAnova(x, dose))
  aov.out <- data.frame(aov.pvalues)
  aov.out$adjusted.p <- p.adjust(aov.out$aov.pvalues, 'fdr')
  aov.out[which(is.nan(aov.out$adjusted.p)),'adjusted.p'] <- 1
  return(aov.out)
}

#' Performs a ANOVA test on a log normalized vector for a given dose
#'
#' @param data The log normalized counts vector
#' @param dose The dose vector of the same length as data
#'
#' @return A p value from the ANOVA test
#'
#' @export
runAnova <- function(data, dose){
  my_data <- data.frame(value = data, dose = dose)
  res.aov <- aov(value ~ dose, data = my_data)
  aov.pvalue <- summary(res.aov)[[1]][["Pr(>F)"]][1]
  return(aov.pvalue)
}
