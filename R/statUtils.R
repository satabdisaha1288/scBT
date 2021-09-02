#' INSERT DESCRIPTION
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return INSERT DESCRIPTION
#' 
#' @export
runMAST <- function(sce){
  scaRaw <- FromMatrix(as.matrix(logcounts(sce)),
                      data.frame(colData(sce)),
                      data.frame(rowData(sce))
  )
  cdr <- colSums(assay(scaRaw) > 0)
  colData(scaRaw)$cngeneson <- scale(cdr)
  #filterCrit = with(colData(scaRaw), cngeneson > 1)
  rowData(scaRaw)$symbolid <- rownames(rowData(sce))
  Dose <- factor(colData(scaRaw)$Dose)
  Dose <- relevel(Dose,"0")
  colData(scaRaw)$Dose <- Dose
  zlmDose <- zlm(~Dose + cngeneson, sca = scaRaw)
  
  summaryDose <- rep(list(list()), times = nlevels(Dose)-1)
  names(summaryDose) <- levels(Dose)[-1]
  for(i in levels(Dose)[-1])
  {
    summaryDose[[i]] <- summary(zlmDose, doLRT=paste0("Dose",i))
  }
  
  summaryDT <- rep(list(data.table::data.table()),times=nlevels(Dose)-1)
  names(summaryDT) <- levels(Dose)[-1]
  fcHurdle <- rep(list(data.table::data.table()),times=nlevels(Dose)-1)
  names(fcHurdle) <- levels(Dose)[-1]
  fcHurdleSig <- rep(list(data.table::data.table()),times=nlevels(Dose)-1)
  names(fcHurdleSig) <- levels(Dose)[-1]
  fcHurdleNonSig <- rep(list(data.table::data.table()),times=nlevels(Dose)-1)
  names(fcHurdleSig) <- levels(Dose)[-1]
  TP_dose <- rep(list(data.table::data.table()),times=nlevels(Dose)-1)
  names(TP_dose) <- levels(Dose)[-1]
  FP_dose <- rep(list(data.table::data.table()),times=nlevels(Dose)-1)
  names(FP_dose) <- levels(Dose)[-1]
  TN_dose <- rep(list(data.table::data.table()),times=nlevels(Dose)-1)
  names(TN_dose) <- levels(Dose)[-1]
  FN_dose <- rep(list(data.table::data.table()),times=nlevels(Dose)-1)
  names(FN_dose) <- levels(Dose)[-1]
  
  for(i in levels(Dose)[-1]){
    summaryDT[[i]] <- summaryDose[[i]]['datatable']
    fcHurdle[[i]] <- summaryDT[[i]][summaryDT[[i]]$contrast == paste0("Dose",i) & summaryDT[[i]]$component=='H',.(primerid, `Pr(>Chisq)`)] #logFC coefficients
    fcHurdle[[i]][,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdleSig[[i]] <- fcHurdle[[i]][fdr<.05 ]
    setorder(fcHurdleSig[[i]], fdr)
    fcHurdleNonSig[[i]] <- fcHurdle[[i]][fdr >.05| fdr == 0.05]
    setorder(fcHurdleNonSig[[i]], fdr)
  }
  
  mast.out <- do.call(cbind, lapply(fcHurdle, data.frame))
  
  if (sum(transform(mast.out[,grepl('primer', colnames(mast.out))], same = apply(mast.out[,grepl('primer', colnames(mast.out))], 1, function(x) length(unique(x)) == 1))$same == FALSE) == 0) {
    m <- mast.out[,grepl('fdr', colnames(mast.out))]
    rownames(m) <- mast.out[,1]
  } else {
    message('ERROR!!!!!!')
  }
  return(m)
}


#' INSERT DESCRIPTION
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return INSERT DESCRIPTION
#' 
#' @export
runMASTDR <- function(sce){
  scaRaw <- FromMatrix(as.matrix(logcounts(sce)),
                       data.frame(colData(sce)),
                       data.frame(rowData(sce))
  )
  cdr <- colSums(assay(scaRaw) > 0)
  colData(scaRaw)$cngeneson <- scale(cdr)
  rowData(scaRaw)$symbolid <- rownames(rowData(sce))
  Dose <- factor(colData(scaRaw)$Dose)
  Dose <- relevel(Dose,"0")
  colData(scaRaw)$Dose = Dose
  cngeneson <- colData(scaRaw)$cngeneson
  zlmDose = MAST::zlm(~Dose + cngeneson, scaRaw)
  print(zlmDose)
  print(summary(zlmDose))
  print(names(zlmDose))
  m <- rep(list(list()), times = nlevels(Dose)-1)
  names(m) <- levels(Dose)[-1]
  for(i in levels(Dose)[1])
  {
    group <- paste0("Dose", i)
    summaryDose <- summary(zlmDose, doLRT="Dose0.01")
    print(names(summaryDose[[1]]))
    print(head(summaryDose[[1]]))
    print(class(summaryDose[[1]]))
    fcHurdle <- data.frame(
      summaryDose[[1]] %>% 
      dplyr::filter(contrast == group & component == "H")
      )
    fcHurdle$fdr <- p.adjust(fcHurdle$`Pr..Chisq`, 'fdr')
    m[[i]] <- data.frame(FDR = fcHurdle$fdr)
    rownames(m[[1]]) <- fcHurdle$primerid
  }
  
  mast.out <- do.call(cbind, lapply(m, data.frame))
  colnames(mast.out) <- paste0(colnames(mast.out), levels(Dose)[-1])
  return(mast.out)
}

#' Calculate concordance AUC Refer to: Soneson, Charlotte, and Mark D. Robinson.
#' "Bias, robustness and scalability in single-cell differential
#' expression analysis." Nature methods 15, no. 4 (2018): 255.
#'
#' @param DETestoutput
#' @author Satabdi Saha
#' @return Normed area under the concordance curve.
#' @export
#'
#' @examples
calculate_concordance_AUC<-function(DETestoutput)
{
  kw.data<-DETestoutput$KW
  kw.sort<-kw.data[order( kw.data[,1] ),,drop=FALSE]
  aov.data<-DETestoutput$ANOVA
  aov.sort<-aov.data[order( aov.data[,1] ),,drop=FALSE]
  h<-seq(1,100,1)
  common_dge<-vector(length = length(h))
  for(i in 1:length(h))
  {
    common_dge[i]<-length(intersect(rownames(kw.sort)[1:h[i]],
                                    rownames(aov.sort)[1:h[i]]))
  }
  library(zoo)
  id <- order(h)
  AUC <- sum(diff(h[id])*rollmean(common_dge[id],2))
  AUC_norm<-AUC/(h[length(h)]^2 /2)
  return(AUC_norm)
}

#' Title Computes ROC curve for DEG classification
#'
#' @param sim
#' @param DETestoutput
#' @author Satabdi Saha
#' @return ROC plot
#' @export
#'
#' @examples
compute_ROC_curve<-function(sim,DETestoutput){
  require(PRROC)
  true_model<-rowData(sim)[,"Model"]
  true_model<-ifelse(true_model== "Unchanged",0,1)
  wfg<- c(runif(300,min=0.5,max=1),runif(500,min=0,max=0.5))
  fg <- DETestoutput$KW$kw.pvalues[true_model == 1]
  bg <- DETestoutput$KW$kw.pvalues[true_model == 0]
  x_KW<-c(fg,bg)
  lab<-c(rep(1,length(fg)),rep(0,length(bg)))
  x_aov<-c(DETestoutput$ANOVA$aov.pvalues[true_model == 1],
           DETestoutput$ANOVA$aov.pvalues[true_model == 0])
  #ROC Curve
  wroc_KW<-roc.curve(scores.class0 = x_KW, weights.class0 = lab, curve = TRUE,
                     max.compute = T, min.compute = T, rand.compute = T)
  wroc_aov<-roc.curve(scores.class0 = x_aov, weights.class0 = lab, curve = TRUE,
                      max.compute = T, min.compute = T, rand.compute = T)
  wroc_plot<-plot(wroc_KW,max.plot = TRUE, min.plot = TRUE,
                  rand.plot = TRUE, fill.area = TRUE,color = 2,
                  scale.color = heat.colors(100))
  wroc_plot<-plot(wroc_aov, add = TRUE, color = 3);
  return(c(wroc_KW,wroc_aov,wroc_plot))
}

#' Summarises data by groups
#' 
#' @param data a data object
#' 
#' @return INSERT DESCRIPTION
#' 
#' @export
data_group_summary <- function(data){
  my_data_summary <- rep(list(data.frame()),ncol(data)-1)
  #Look at summary statistics for each gene by group
  for(i in 1: (ncol(data)-1))
  {
    my_data <- data.frame(data[,i],as.factor(data[,ncol(data)]))
    colnames(my_data) <- c("value","dose")
    library(dplyr)
    my_data_summary[[i]] <- group_by(my_data, dose) %>%
      summarise(
        count <- n(),
        mean <- mean(value, na.rm = TRUE),
        mean_pos <- mean(value[value>0], na.rm = TRUE),
        quantile_25 <- quantile(value,probs=0.25,na.rm=TRUE),
        quantile_50 <- quantile(value,probs=0.50,na.rm=TRUE),
        quantile_75 <- quantile(value,probs=0.75,na.rm=TRUE),
        sd <- sd(value, na.rm = TRUE),
        drop_prop <- length(which(value == 0))/ length(value),
      )
  }
  return(my_data_summary)
}

#' Title Computes PR curve for DEG classification
#'
#' @param sim
#' @param DETestoutput
#' @author Satabdi Saha
#' @return PR plot
#' @export
#'
#' @examples
compute_PR_curve<-function(sim,DETestoutput){
  require(PRROC)
  true_model<-rowData(sim)[,"Model"]
  true_model<-ifelse(true_model== "Unchanged",0,1)
  fg <- DETestoutput$KW$kw.pvalues[true_model == 1]
  bg <- DETestoutput$KW$kw.pvalues[true_model == 0]
  x_KW<-c(fg,bg)
  lab<-c(rep(1,length(fg)),rep(0,length(bg)))
  x_aov<-c(DETestoutput$ANOVA$aov.pvalues[true_model == 1],
           DETestoutput$ANOVA$aov.pvalues[true_model == 0])
  #ROC Curve
  wPR_KW<-pr.curve(scores.class0 = x_KW, weights.class0 = lab, curve = TRUE,
                   max.compute = T, min.compute = T, rand.compute = T)
  wPR_aov<-pr.curve(scores.class0 = x_aov, weights.class0 = lab, curve = TRUE,
                    max.compute = T, min.compute = T, rand.compute = T)
  wPR_plot<-plot(wPR_KW,max.plot = TRUE, min.plot = TRUE,
                 rand.plot = TRUE, fill.area = TRUE,color = 2,
                 scale.color = heat.colors(100))
  wPR_plot<-plot(wPR_aov, add = TRUE, color = 3);
  return(c(wPR_KW,wPR_aov,wPR_plot))
}

