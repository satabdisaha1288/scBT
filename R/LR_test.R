#' INSERT DESCRIPTION HERE
#' 
#' @author
#' @param w.x vector of zeros/ones for expressed or not in each group
#' @param w.y vector of zeros/ones for expressed or not in each group
#' @param x vector of the positive observations (must be of length sum(w.x) and sum(w.y))
#' @param y vector of the positive observations (must be of length sum(w.x) and sum(w.y))
#' 
#' @return INSERT RETURN DESCR HERE
#' 
#' @export
lrtest <- function(w.x, w.y, x, y){
  e.x <- sum(w.x)
  e.y <- sum(w.y)
  n.x <- length(w.x)
  n.y <- length(w.y)
  stopifnot(e.x == length(x) && e.y == length(y))
  
  
  p.0 <- (e.x+e.y)/(n.x + n.y)
  p.x <- e.x/n.x
  p.y <- e.y/n.y
  
  m0 <- (sum(x)+sum(y))/(e.x+e.y)
  mu.x <- mean(x)
  mu.y <- mean(y)
  
  Tstar <- 1+e.x*e.y/(e.x+e.y) * (mu.x - mu.y)^2/(sum((mu.x - x)^2) + sum((mu.y-y)^2))
  
  if(!is.finite(Tstar)){
    Tstar <- 1
  }
  
  binom <- logProd(e.x, p.0/p.x) +
    logProd(e.y, p.0/p.y) +
    logProd(n.x-e.x, (1-p.0)/(1-p.x)) +
    logProd(n.y-e.y, (1-p.0)/(1-p.y))
  binomsign <- (p.y>p.x)*2 -1
  
  norm <- -(e.x+e.y)/2 * log(Tstar)
  normsign <- (mu.y>mu.x)*2-1
  
  logLR <- binom+norm
  
  maxsign <- c(binomsign, normsign)[which.min(c(binom, norm))]
  resultvec <- c(-2*binom, binomsign, pchisq(-2*binom, 1, lower.tail = FALSE),
                 -2*norm, normsign, pchisq(-2*norm, 1, lower.tail = FALSE),
                 -2*logLR, maxsign, pchisq(-2*logLR, 2, lower.tail = FALSE))
  result <- matrix(resultvec, nrow=3, ncol=3, dimnames=list(metric = c('lrstat', 'direction', 'p.value'), component=c('binom', 'norm', 'comb')))
  result_comb <- cbind(rownames(data_0_filt_binary[i,]),-2*logLR, maxsign, pchisq(-2*logLR, 2, lower.tail=FALSE))
  colnames(result_comb) <- c("Gene_name","lrstat","Direction","p-value")
  return(result_comb)
}

#' INSERT DESCRIPTION HERE
#' 
#' @author
#' @param prod INSERT DESCR
#' @param logand INSERT DESCR
#' 
#' @export
logProd <- function(prod, logand){
  ifelse(prod == 0, 0, prod * log(logand))
}

#' Performs a genewise ANOVA test on a SingleCellExperiment object
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
runLRT <- function(data.list){
  for (dose in names(data.list)){
    data.list[[dose]] <- as.matrix(data.list[[dose]])
  }

  LRT_mult <- list()
  data_ind <- lapply(data.list, function(x) ifelse(x > 0, 1, 0))
  for(g in 1: nrow(data.list[["0"]])){
    gene <- rownames(data.list[["0"]])[g]
    LRT_mult[[gene]] <- LRT_multiple_groups(data = list(data.list[["0"]][g,][data_ind[["0"]][g,] > 0],
                                                      data.list[["0.01"]][g,][data_ind[["0.01"]][g,] > 0],
                                                      data.list[["0.03"]][g,][data_ind[["0.03"]][g,] > 0],
                                                      data.list[["0.1"]][g,][data_ind[["0.1"]][g,] > 0],
                                                      data.list[["0.3"]][g,][data_ind[["0.3"]][g,] > 0],
                                                      data.list[["1"]][g,][data_ind[["1"]][g,] > 0],
                                                      data.list[["3"]][g,][data_ind[["3"]][g,] > 0],
                                                      data.list[["10"]][g,][data_ind[["10"]][g,] > 0],
                                                      data.list[["30"]][g,][data_ind[["30"]][g,] > 0]),
                                          data_ind = list(data_ind[["0"]][g,],
                                                          data_ind[["0.01"]][g,],
                                                          data_ind[["0.03"]][g,],
                                                          data_ind[["0.1"]][g,],
                                                          data_ind[["0.3"]][g,],
                                                          data_ind[["1"]][g,],
                                                          data_ind[["3"]][g,],
                                                          data_ind[["10"]][g,],
                                                          data_ind[["30"]][g,]))
  }
  LRT_mult_p_value <- sapply(LRT_mult,function(x) x[2,3])
  LRT_mult_p_value_adj <- p.adjust(LRT_mult_p_value, method = "fdr")
  LRT.out <- data.frame(pval = LRT_mult_p_value, pval.fdr = LRT_mult_p_value_adj)
  rownames(LRT.out) <- names(LRT_mult_p_value_adj)
  return(LRT.out)
}
