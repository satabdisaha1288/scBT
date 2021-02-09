#Function for getting the summary statistics of each gene by group
#data is a dataframe of n cells as rows and p+1 columns with first p columns as p genes and the last 
#column is the dose information
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
        sd = sd(value, na.rm = TRUE)
      )
  }
  return(my_data_summary)
}
