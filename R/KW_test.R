#Function for carrying out the Kruskal Wallis Test
#data is a dataframe of n cells as rows and p+1 columns with first p columns as p genes and the last 
#column is the dose information
KW_test<-function(data){
  my_data_kruskal_p_value<-vector()
  for(i in 1: (ncol(data)-1))
  {
    my_data<-data.frame(data[,i],as.factor(data[,ncol(data)]))
    colnames(my_data)<-c("value","dose")
    res.kruskal<-kruskal.test(value ~ dose, data = my_data)
    my_data_kruskal_p_value[i]<-res.kruskal$p.value
  }
  return(my_data_kruskal_p_value)
}