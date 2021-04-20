#' INSERT DESCRIPTION HERE
#' 
#' @author
#' @param data the list of vectors of the positive observations( must of length sum(each vector in data_ind))
#' @param data_ind a list having K vectors of the zeroes/ones for expressed or not in each of the K -groups
#' 
#' @return INSERT RETURN DESCR HERE
#' 
#' @example 
#' 
#' @export
LRT_multiple_groups <- function(data, data_ind){
  data_pos<-lapply(data_ind,function(x) sum(x))
  data_size<-lapply(data_ind,function(x) length(x))
  dim.check<-Map("==",data_pos, lapply(data,function(x) length(x)))
  allSame <- function(x) length(unique(x)) == 1
  stopifnot(allSame(dim.check == TRUE))
 
  w_null <- unlist(Map("/",Reduce("+",data_pos), Reduce("+",data_size)))
  w_group<- Map("/",data_pos,data_size) 
  
  library(dplyr)
  
  mu_null<- as.numeric(bind_rows(lapply(data, as.data.frame)) %>% colSums(na.rm=TRUE)/
                bind_rows(lapply(data_pos, as.data.frame)) %>% colSums(na.rm=TRUE))
  mu_group<-lapply(data, function(x) mean(x))
 
ssg<-vector()
sse<-data
 for(k in 1: length(data))
 {
   ssg[k]<-data_pos[[k]]*((mu_group[[k]] - mu_null)^2)
    for(i in 1: length(data[[k]]))
    {
     sse[[k]][i]<- (data[[k]][i]-mu_group[[k]])^2
    }
 }
ss_between<-sum(ssg)
ss_residuals<-sum(sapply(sse, function(x) sum(x)))
Tstar <-  1+(ss_between/ss_residuals)

if(!is.finite(Tstar)){
  Tstar <- 1
}
norm <- -((Reduce("+",data_pos))/2) * log(Tstar)
a<-vector()
b<-vector()
for(k in 1: length(data))
{
  a[k]<-logProd(data_pos[[k]], w_null/w_group[[k]])
  b[k]<-logProd(data_size[[k]]-data_pos[[k]], (1 - w_null)/(1- w_group[[k]]))
}
binom <- sum(a)+sum(b)

logLR <- binom + norm

resultvec <- c(-2*binom,  pchisq(-2*binom, length(data) - 1 , lower.tail=FALSE),
               -2*norm,  pchisq(-2*norm, length(data) - 1, lower.tail=FALSE),
               -2*logLR, pchisq(-2*logLR,2*length(data)-2 , lower.tail=FALSE))
result <- matrix(resultvec, nrow=2, ncol=3, dimnames=list(metric=c('lrstat', 'p.value'), component=c('binom', 'norm', 'comb')))
return(result)

logProd <- function(prod, logand){
  ifelse(prod==0, 0, prod*log(logand))
}
}