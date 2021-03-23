#Function for carrying out the Wicoxon Rank Sum test
#data is a dataframe of n cells as rows and p+1 columns with first p columns as p genes and the last 
#column is the dose information
WRS_test<-function(data,dose_level_ref,dose_level_test){
  my_data_wilcoxon_p_value<-vector()
  for(i in 1: (ncol(data)-1))
  {
    res.wrs<-wilcox.test(data[,i][data[,ncol(data)]== dose_level_ref]
                         ,data[,i][data[,ncol(data)]== dose_level_test])
    my_data_wilcoxon_p_value[i]<-res.wrs$p.value
  }
  return(my_data_wilcoxon_p_value)
}