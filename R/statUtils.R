#' INSERT DESCRIPTION
#' 
#' @param sce SingleCellExperiment object with a logcounts assay 
#' and Dose column in the cell metadata
#' 
#' @return INSERT DESCRIPTION
#' 
#' @export
runMAST = function(sce){
  library(MAST)
  scaRaw = FromMatrix(as.matrix(logcounts(sce)),
                      data.frame(colData(sce)),
                      data.frame(rowData(sce))
  )
  cdr <-colSums(assay(scaRaw) > 0)
  colData(scaRaw)$cngeneson = scale(cdr)
  filterCrit = with(colData(scaRaw),cngeneson > 1)
  rowData(scaRaw)$symbolid = rownames(rowData(sce))
  sca = subset(scaRaw,filterCrit)
  Dose = factor(colData(scaRaw)$Dose)
  Dose = relevel(Dose,"0")
  colData(scaRaw)$Dose = Dose
  zlmDose = zlm(~Dose + cngeneson, sca = scaRaw)
  
  summaryDose<-rep(list(list()),times= nlevels(Dose)-1)
  names(summaryDose)<-levels(Dose)[2:9]
  for(i in levels(Dose)[2:9])
  {
    summaryDose[[i]]<-summary(zlmDose, doLRT=paste0("Dose",i))
  }
  
  summaryDT<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(summaryDT)<-levels(Dose)[2:9]
  fcHurdle<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(fcHurdle)<-levels(Dose)[2:9]
  fcHurdleSig<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(fcHurdleSig)<-levels(Dose)[2:9]
  fcHurdleNonSig<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(fcHurdleSig)<-levels(Dose)[2:9]
  TP_dose<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(TP_dose)<-levels(Dose)[2:9]
  FP_dose<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(FP_dose)<-levels(Dose)[2:9]
  TN_dose<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(TN_dose)<-levels(Dose)[2:9]
  FN_dose<-rep(list(data.table()),times=nlevels(Dose)-1)
  names(FN_dose)<-levels(Dose)[2:9]
  
  for(i in levels(Dose)[2:9]){
    summaryDT[[i]] <-summaryDose[[i]]$datatable
    fcHurdle[[i]]<-summaryDT[[i]][contrast== paste0("Dose",i) & component=='H',.(primerid, `Pr(>Chisq)`)] #logFC coefficients
    fcHurdle[[i]][,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdleSig[[i]]<- fcHurdle[[i]][fdr<.05 ]
    setorder(fcHurdleSig[[i]], fdr)
    fcHurdleNonSig[[i]]<- fcHurdle[[i]][fdr >.05| fdr == 0.05]
    setorder(fcHurdleNonSig[[i]], fdr)
  }
  
  mast.out = do.call(cbind, lapply(fcHurdle, data.frame))
  
  if (sum(transform(mast.out[,grepl('primer', colnames(mast.out))], same = apply(mast.out[,grepl('primer', colnames(mast.out))], 1, function(x) length(unique(x)) == 1))$same == FALSE) == 0) {
    m = mast.out[,grepl('fdr', colnames(mast.out))]
    rownames(m) = mast.out[,1]
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
clean_sce = function(sce){
  library(Matrix)
  assays(sce)$BatchCellMeans = NULL
  assays(sce)$BaseCellMeans = NULL
  assays(sce)$BCV = NULL
  assays(sce)$CellMeans = NULL
  assays(sce)$TrueCounts = NULL	
  assays(sce)$counts = Matrix(assays(simDR[[vec_elem]])$counts, sparse = TRUE)
  assays(sce)$logcounts = Matrix(assays(simDR[[vec_elem]])$logcounts, sparse = TRUE)
  return(sce)
}
