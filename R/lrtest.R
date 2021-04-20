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
#' @example 
#' 
#' @export
lrtest <- function(w.x, w.y, x, y){
  e.x <- sum(w.x)
  e.y <-  sum(w.y)
  n.x <-  length(w.x)
  n.y <-  length(w.y)
  stopifnot(e.x == length(x) && e.y == length(y))
  
  
  p.0 <- (e.x+e.y)/(n.x + n.y)
  p.x <- e.x/n.x
  p.y <- e.y/n.y
  
  m0 <-  (sum(x)+sum(y))/(e.x+e.y)
  mu.x <-  mean(x)
  mu.y <-  mean(y)
  
  Tstar <-  1+e.x*e.y/(e.x+e.y)* (mu.x - mu.y)^2/(sum((mu.x - x)^2) + sum((mu.y-y)^2))
  
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
  resultvec <- c(-2*binom, binomsign, pchisq(-2*binom, 1, lower.tail=FALSE),
                 -2*norm, normsign, pchisq(-2*norm, 1, lower.tail=FALSE),
                 -2*logLR, maxsign, pchisq(-2*logLR, 2, lower.tail=FALSE))
  result <- matrix(resultvec, nrow=3, ncol=3, dimnames=list(metric=c('lrstat', 'direction', 'p.value'), component=c('binom', 'norm', 'comb')))
  result_comb<-cbind(rownames(data_0_filt_binary[i,]),-2*logLR, maxsign, pchisq(-2*logLR, 2, lower.tail=FALSE))
  colnames(result_comb)<-c("Gene_name","lrstat","Direction","p-value")
  return(result_comb)
}

#' INSERT DESCRIPTION HERE
#' 
#' @author
#' @param prod INSERT DESCR
#' @param logand INSERT DESCR
#' 
#' @example 
#' 
#' @export
logProd <- function(prod, logand){
  ifelse(prod==0, 0, prod*log(logand))
}


