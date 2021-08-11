#' Performs a genewise Wilcoxon Rank Sum test on a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object with a logcounts assay
#' and Dose column in the cell metadata
#'
#' @return a vector of p values from the Wilcoxon Rank Sum test
#' @export
batchWRS <- function(sce){
  data <- as.matrix(logcounts(sce))
  dose <- colData(sce)$Dose
  wrs.pvalues <- apply(data, 1, function(x) runWRS(x, dose))
  rownames(wrs.pvalues) <- levels(dose)[-1]
  wrs.pvalues <- data.frame(t(wrs.pvalues))
  wrs.pvalues[is.na(wrs.pvalues)] <- 1
  wrs.pvalues$adjusted.p <- apply(wrs.pvalues, 1, function(x) p.adjust(min(x), 'fdr'))
  #wrs.adj <- data.frame(apply(wrs.pvalues, 2, function(x) p.adjust(x, 'fdr')))
  #colnames(wrs.adj) <- paste0('adjusted.p.', colnames(wrs.adj))
  #wrs.adj$adjusted.p <- apply(wrs.adj, 1, function(x) min(x))
  #wrs = cbind(wrs.pvalues, wrs.adj)
  #wrs <- wrs.pvalues
  #wrs[which(is.nan(wrs$adjusted.p)), 'adjusted.p'] <- 1

  return(wrs.pvalues)
}

#' Performs a Wilcoxon Rank Sum test on a logcounts vector for a given dose
#'
#' @param data The logcounts vector
#' @param dose The dose to analyze
#'
#' @return A p value from the Wilcoxon Rank Sum test
#' @export
runWRS <- function(data, dose){
  my_data <- data.frame(value = data, dose = dose)
  my_control <- my_data %>% filter(dose == levels(dose)[1])
  res.wrs <- my_data %>%
    filter(dose != levels(dose)[1]) %>%
    group_by(dose) %>%
    summarise(p_value = wilcox.test(my_control$value, value)$p.value, .groups = 'drop')
  wrs.pvalue <- t(data.frame(p.value = res.wrs$p_value))
  return(wrs.pvalue)
}
