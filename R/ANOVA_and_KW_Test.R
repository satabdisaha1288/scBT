my_data_summary<-rep(list(data.frame()),ncol(simulated.data_dose)-1)
#Look at summary statistics for each gene by group
for(i in 1: (ncol(simulated.data_dose)-1))
{
  my_data<-data.frame(simulated.data_dose[,i],as.factor(simulated.data_dose[,ncol(simulated.data_dose)]))
  colnames(my_data)<-c("value","dose")
  library(dplyr)
  my_data_summary[[i]]<-group_by(my_data, dose) %>%
    summarise(
      count = n(),
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE)
    )
}
#Conduct the ANOVA and Kruskal-Wallis Test for each gene
my_data_anova<-rep(list(list()),ncol(simulated.data_dose)-1)
my_data_anova_p_value<-vector()
my_data_kruskal_p_value<-vector()
for(i in 1: (ncol(simulated.data_dose)-1))
{
  my_data<-data.frame(simulated.data_dose[,i],as.factor(simulated.data_dose[,ncol(simulated.data_dose)]))
  colnames(my_data)<-c("value","dose")
  res.aov <- aov(value ~ dose, data = my_data)
  my_data_anova[[i]]<-summary.aov(res.aov)
  my_data_anova_p_value[i]<-my_data_anova[[i]][[1]][["Pr(>F)"]][1]
  res.kruskal<-kruskal.test(value ~ dose, data = my_data)
  my_data_kruskal_p_value[i]<-res.kruskal$p.value
}
my_data_anova_pred_dec<-ifelse(my_data_anova_p_value < 0.01, "Positive","Negative")
my_data_kruskal_pred_dec<-ifelse(my_data_kruskal_p_value < 0.01, "Positive","Negative")
length(which(my_data_anova_pred_dec == "Negative"))
length(which(my_data_kruskal_pred_dec == "Negative"))
my_data_true_dec<-ifelse(simulated.gene.metadata$Direction == "Upregulated" | 
                           simulated.gene.metadata$Direction == "Downregulated" , "Positive","Negative")


