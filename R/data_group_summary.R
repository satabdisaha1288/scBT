#' Summarises data by groups
#' 
#' @param data a data object
#' 
#' @return INSETRT DESCRIPTION
#' 
#' @export
data_group_summary<-function(data){
  my_data_summary<-rep(list(data.frame()),ncol(data)-1)
  #Look at summary statistics for each gene by group
  for(i in 1: (ncol(data)-1))
  {
    my_data<-data.frame(data[,i],as.factor(data[,ncol(data)]))
    colnames(my_data)<-c("value","dose")
    library(dplyr)
    my_data_summary[[i]]<-group_by(my_data, dose) %>%
      summarise(
        count = n(),
        mean = mean(value, na.rm = TRUE),
        mean_pos = mean(value[value>0], na.rm = TRUE),
        quantile_25 = quantile(value,probs=0.25,na.rm=TRUE),
        quantile_50 = quantile(value,probs=0.50,na.rm=TRUE),
        quantile_75 = quantile(value,probs=0.75,na.rm=TRUE),
        sd = sd(value, na.rm = TRUE),
        drop_prop = length(which(value == 0))/ length(value),
      )
  }
  return(my_data_summary)
}