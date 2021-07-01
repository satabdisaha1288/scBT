#' Performs a genewise Wilcoxon Rank Sum test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the Wilcoxon Rank Sum test
#' 
#' @import limma
#' @export
runLimmaTrend <- function(sce){
  data <- as.matrix(logcounts(sce))
  dose <- colData(sce)$Dose
  zeroSums <- which(rowSums(data) == 0)
  expr <- data[-zeroSums,]
  
  design <- model.matrix(~dose)
  fit <- limma::lmFit(expr, design = design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, n = Inf, adjust.method = 'BH')
  
  missing <- data.frame(matrix(ncol = ncol(tt), nrow = length(zeroSums)))
  colnames(missing) <- colnames(tt)
  rownames(missing) <- names(zeroSums)
  missing$adj.P.Val <- 1
  limma.out <- rbind(tt, missing)
  limma.out <- limma.out[rownames(data),]
  colnames(limma.out) <- gsub('adj.P.Val', 'adjusted.p', colnames(limma.out))
  
  return(limma.out)
}
