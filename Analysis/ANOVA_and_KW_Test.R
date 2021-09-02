
#Conduct the ANOVA and Kruskal-Wallis Test for each gene
source('~/Documents/Documents/Summer_2019/Mouse_Liver/Rance/snSeq-DGE/R/anova_test.R')
source('~/Documents/Documents/Summer_2019/Mouse_Liver/Rance/snSeq-DGE/R/KW_test.R')
#Here data is a dataframe of n cells as rows and p+1 columns with first p columns as p genes and the last 
#column is the dose information
mydata<-simulated.data_dose
my_data_kruskal_p_value<-KW_test(simulated.data_dose)
my_data_anova_p_value<-anova_test(simulated.data_dose)
my_data_anova_pred_dec<-ifelse(my_data_anova_p_value < 0.05, "Positive","Negative")
my_data_kruskal_pred_dec<-ifelse(my_data_kruskal_p_value < 0.05, "Positive","Negative")
my_data_true_dec<-ifelse(simulated.gene.metadata$Direction == "Upregulated" | 
                           simulated.gene.metadata$Direction == "Downregulated" , "Positive","Negative")

#Generate the confusion matrix for the ANOVA and KW test wrt the reference data
caret::confusionMatrix(factor(my_data_anova_pred_dec,levels=c("Positive","Negative")),
                       factor(my_data_true_dec,levels=c("Positive","Negative")))
caret::confusionMatrix(factor(my_data_kruskal_pred_dec,levels=c("Positive","Negative")),
                       factor(my_data_true_dec,levels=c("Positive","Negative")))


