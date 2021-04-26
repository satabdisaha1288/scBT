#' INSERT DESCRIPTION
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return a vector of p values from the ANOVA test
#' 
#' @export
runLRT = function(data.list){
  lrt_multiple_01 = list()
  for(j in rownames(data.list[[1]])){
    in.list = data.list %>% purrr::map(as.matrix(~.x[j,]))
    names(in.list) = paste0("Y_", 1:length(in.list))
    lrt_multiple_01 = LRT_multiple_groups(in.list)
  }
  
  for (dose in names(data.list)){
    data.list[[dose]] = as.matrix(data.list[[dose]])
  }
  #LRT_mult<-rep(list(list()),nrow(data[["0"]]))
  # LRT_mult = list()
  # data_ind = lapply(data.list, function(x) ifelse(x > 0, 1, 0))
  # for(g in 1: nrow(data.list[["0"]])){
  #   gene = rownames(data.list[["0"]])[g]
  #   LRT_mult[[gene]]<-LRT_multiple_groups(data = list(data.list[["0"]][g,][data_ind[["0"]][g,] > 0],
  #                                                     data.list[["0.01"]][g,][data_ind[["0.01"]][g,] > 0],
  #                                                     data.list[["0.03"]][g,][data_ind[["0.03"]][g,] > 0],
  #                                                     data.list[["0.1"]][g,][data_ind[["0.1"]][g,] > 0],
  #                                                     data.list[["0.3"]][g,][data_ind[["0.3"]][g,] > 0],
  #                                                     data.list[["1"]][g,][data_ind[["1"]][g,] > 0],
  #                                                     data.list[["3"]][g,][data_ind[["3"]][g,] > 0],
  #                                                     data.list[["10"]][g,][data_ind[["10"]][g,] > 0],
  #                                                     data.list[["30"]][g,][data_ind[["30"]][g,] > 0]),
  #                                         data_ind = list(data_ind[["0"]][g,],
  #                                                         data_ind[["0.01"]][g,],
  #                                                         data_ind[["0.03"]][g,],
  #                                                         data_ind[["0.1"]][g,],
  #                                                         data_ind[["0.3"]][g,],
  #                                                         data_ind[["1"]][g,],
  #                                                         data_ind[["3"]][g,],
  #                                                         data_ind[["10"]][g,],
  #                                                         data_ind[["30"]][g,]))
  # }
  # LRT_mult_p_value = sapply(LRT_mult,function(x) x[2,3])
  # LRT_mult_p_value_adj = p.adjust(LRT_mult_p_value, method = "fdr")
  # LRT.out = data.frame(pval = LRT_mult_p_value, pval.fdr = LRT_mult_p_value_adj)
  # rownames(LRT.out) = names(LRT_mult_p_value_adj)
  # return(LRT.out)
}


#' INSERT DESCRIPTION HERE
#' 
#' @author
#' @param data the list of vectors of the positive observations( must of length sum(each vector in data_ind))
#' @param data_ind a list having K vectors of the zeroes/ones for expressed or not in each of the K -groups
#' 
#' @return INSERT RETURN DESCR HERE
#' 
#' @export
LRT_multiple_groups <- function(data.list){
  
  data = lapply(data.list, function(x) x[which(x > 0)])
  data_pos <- lapply(data.list, function(x) length(x[which(x != 0)])) #Number of non-zero genes
  data_size<-lapply(data.list, function(x) length(x)) #Both zero's and ones
  dim.check<-Map("==", data_pos, lapply(data, function(x) length(x)))
  allSame <- function(x) length(unique(x)) == 1
  stopifnot(allSame(dim.check == TRUE))
  
  w_null <- unlist(Map("/",Reduce("+",data_pos), Reduce("+",data_size)))
  w_group<- Map("/",data_pos,data_size) 

  mu_null<- as.numeric(bind_rows(lapply(data, as.data.frame)) %>% colSums(na.rm=TRUE)/
                         bind_rows(lapply(data_pos, as.data.frame)) %>% colSums(na.rm=TRUE))
  mu_group<-lapply(data, function(x) mean(x))
  
  ssg<-vector()
  sse<-data
  # TODO: Look into possible using an apply function
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
  # TODO: Look into changing to apply also
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
  names(resultvec) <- c("lrstat.binom", "p.value.binom", "lrstat.norm", "p.value.norm", "lrstat.comb", "p.value.comb")
  return(resultvec)
  #TODO: result is a 4 x 4 matrix but we want it a as a data frame. Gene = rows, columns = values 
}