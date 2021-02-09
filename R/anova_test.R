#Function for ANOVA test for each gene 
#data is a dataframe of n cells as rows and p+1 columns with first p columns as p genes and the last 
#column is the dose information
anova_test<-function(data){
  my_data_anova<-rep(list(list()),ncol(data)-1)
  my_data_anova_p_value<-vector()
  for(i in 1: (ncol(data)-1))
  {
    my_data<-data.frame(data[,i],as.factor(data[,ncol(data)]))
    colnames(my_data)<-c("value","dose")
    res.aov <- aov(value ~ dose, data = my_data)
    my_data_anova[[i]]<-summary.aov(res.aov)
    my_data_anova_p_value[i]<-my_data_anova[[i]][[1]][["Pr(>F)"]][1]
  }
  return(my_data_anova_p_value)
}
