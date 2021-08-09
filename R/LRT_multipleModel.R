#' INSERT DESCRIPTION
#'
#' @param data.list INSERT DESCRIPTION
#'
#' @return INSERT DESCRIPTION
#'
#' @export
LRT_multipleModelNew <- function(data.list){
  lrt_multiple_01 <- list()
  for(j in rownames(data.list[[1]])){
    in.list <- data.list %>% purrr::map(~ as.matrix(.x[j,]))
    lrt_multiple_01[[j]] <- lrt.mult(in.list)
  }

  lrt_multiple.out <- data.frame(t(do.call(cbind, lapply(lrt_multiple_01, data.frame))))
  rownames(lrt_multiple.out) <- names(lrt_multiple_01)
  lrt_multiple.out$adjusted.p <- p.adjust(lrt_multiple.out$p.value.comb, 'fdr')
  return(lrt_multiple.out)
}


#' INSERT DESCRIPTION HERE
#'
#' @author Satabdi Saha
#' @param data the list of vectors of the observations
#' @param data_ind a list having K vectors of the zeroes/ones for expressed or not in each of the K -groups
#'
#' @return INSERT RETURN DESCR HERE
#'
#' @export
lrt.mult <- function(data){
  # data_ind <- lapply(data, function(x) ifelse(x > 0, 1, 0))
  # mu_null <- sum(unlist(lapply(Map("*",data_ind,data),sum)))/ sum(unlist(lapply(data_ind, sum)))
  # mu_group <- Map("/",lapply(Map("*",data_ind,data),sum), (lapply(data_ind, sum)))
  # a <- sum(unlist(lapply(data_ind, sum))) #constant
  # b <- lapply(data_ind, sum)
  # c <- lapply(Map("*",lapply(Map("-",data,mu_group),function(x) x^2),data_ind), sum)
  # d <- sum(unlist(lapply(Map("*",lapply(data, function(x) (x-mu_null)^2),data_ind),sum))) #constant
  # f <- Map("-",lapply(c, log),lapply(b, log))
  # g <- lapply(f, function(x) x + log(a) - log(d))
  # l.norm <- sum(unlist(Map("*",lapply(b, function(x) x/2),g)))

  l.norm <- lrt_norm(data)
  l.binom <- lrt_binom(data)
  # a0 <- lapply(data_ind,sum)
  # b0 <- sum(unlist(a0))#constant
  # c0 <- lapply(data_ind,length)
  # d0 <- sum(unlist(c0))#constant
  # e0 <- b0/d0 #constant
  # #f0 <- lapply(Map("-",lapply(a0, log) ,lapply(c0, log)),function(x) x-log(e0))
  # f0 <- lapply(Map("-",lapply(a0, log) ,lapply(c0, log)),function(x) log(e0)-x)
  # g0 <- sum(unlist(Map("*",f0,a0)))
  # h0 <- log(1-e0) #constant
  # i0 <- lapply(lapply(Map("/",a0,c0), function(x) 1-x),log)
  # j0 <- Map("-",c0,a0)
  # k0 <- Map("*",j0,lapply(i0, function(x) h0-x))
  # l.binom <- sum(unlist(k0)) + g0

  logLR = l.norm + l.binom

  resultvec <- c(-2*l.binom, pchisq(-2*l.binom, length(data) - 1, lower.tail=FALSE),
                 -2* l.norm,  pchisq(-2*l.norm, length(data) - 1, lower.tail=FALSE),
                 -2*logLR, pchisq(-2*logLR,2*length(data)-2 , lower.tail=FALSE))
  names(resultvec) <- c("lrstat.binom", "p.value.binom", "lrstat.norm", "p.value.norm", "lrstat.comb", "p.value.comb")

  if (is.nan(resultvec['lrstat.binom']) & !is.nan(resultvec['lrtstat.norm'])){
    resultvec['lrstat.comb'] <- resultvec['lrstat.norm']
    resultvec['p.value.comb'] <- resultvec['p.value.norm']
  }

  if (!is.nan(resultvec['lrstat.binom']) & is.nan(resultvec['lrtstat.norm'])){
    resultvec['lrstat.comb'] <- resultvec['lrstat.binom']
    resultvec['p.value.comb'] <- resultvec['p.value.binom']
  }

  if (is.infinite(resultvec['lrstat.norm'])){
    resultvec['lrstat.comb'] <- resultvec['lrstat.binom']
    resultvec['p.value.comb'] <- resultvec['p.value.binom']
  }
  if (is.nan(resultvec['lrstat.binom']) & is.nan(resultvec['lrstat.norm'])){
    resultvec['lrstat.comb'] <- NaN
    resultvec['p.value.comb'] <- 1
  }
  return(resultvec)
}


#' INSERT DESCRIPTION HERE
#'
#' @author Satabdi Saha
#' @param data the list of vectors of the observations
#' @param data_ind a list having K vectors of the zeroes/ones for expressed or not in each of the K -groups
#'
#' @return INSERT RETURN DESCR HERE
#'
#' @export
lrt_binom <- function(data){
  data_ind<-lapply(data, function(x) ifelse(x>0,1,0))
  data_pos<-lapply(data_ind,function(x) sum(x))
  data_size<-lapply(data_ind,function(x) length(x))

  w_null <- unlist(Map("/",Reduce("+",data_pos), Reduce("+",data_size)))
  w_group<- Map("/",data_pos,data_size)
  a<-unlist(Map("logProd",data_pos,lapply(w_group, function(x) w_null/x)))
  b<-unlist(Map("logProd", Map("-",data_size,data_pos),lapply(lapply(w_group,
                                                                     function(x) 1-x),function(x) (1-w_null)/x )))

  binom <- sum(a)+sum(b)
  resultvec <- c(-2*binom,  pchisq(-2*binom, length(data) - 1 , lower.tail=FALSE))
  names(resultvec)<-c("lrstat.binom","p.value.binom")
  return(binom)
}

#' INSERT DESCRIPTION HERE
#'
#' @author Satabdi Saha
#' @param data the list of vectors of the observations
#' @param data_ind a list having K vectors of the zeroes/ones for expressed or not in each of the K -groups
#'
#' @return INSERT RETURN DESCR HERE
#'
#' @export
lrt_norm <- function(data){
  data_ind<-lapply(data, function(x) ifelse(x>0,1,0))
  mu_null<- sum(unlist(lapply(Map("*",data_ind,data),sum)))/   sum(unlist(lapply(data_ind, sum)))
  mu_group<- Map("/",lapply(Map("*",data_ind,data),sum), (lapply(data_ind, sum)))
  a = sum(unlist(lapply(data_ind, sum))) #constant
  b=lapply(data_ind, sum)
  c=lapply(Map("*",lapply(Map("-",data,mu_group),function(x) x^2),data_ind), sum)
  d = sum(unlist(lapply(Map("*",lapply(data, function(x) (x-mu_null)^2),data_ind),sum))) #constant
  e = lapply(b, function(x) a/x)
  f = lapply(c, function(x) x/d)
  Tstar<-Map("*",e,f)
  Tstar<-lapply(Tstar, function(x) ifelse(!is.finite(x),1,x))
  norm = sum(unlist(Map("logProd", lapply(b, function(x) x/2), Tstar)))
  resultvec <- c(-2*norm,  pchisq(-2*norm, length(data) - 1, lower.tail=FALSE))
  names(resultvec)<-c("lrstat.norm","p.value.norm")
  return(norm)
}

