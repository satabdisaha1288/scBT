#' Performs a genewise Wilcoxon Rank Sum test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the Wilcoxon Rank Sum test
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
#' @export
runWRS = function(data, dose){
  my_data <- data.frame(value = data, dose = dose)
  res.wrs = wilcox.test(value ~ dose, data = my_data)
  wrs.pvalue = res.wrs[[3]]
  return(wrs.pvalue)
}

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
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
