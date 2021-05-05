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
  cdr <-colSums(assay(scaRaw) > 0)
  colData(scaRaw)$cngeneson <- scale(cdr)
  #filterCrit = with(colData(scaRaw), cngeneson > 1)
  rowData(scaRaw)$symbolid <- rownames(rowData(sce))
  Dose <- factor(colData(scaRaw)$Dose)
  Dose <- relevel(Dose,"0")
  colData(scaRaw)$Dose = Dose
  zlmDose = zlm(~Dose + cngeneson, sca = scaRaw)
  
  summaryDose <- rep(list(list()), times= nlevels(Dose)-1)
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
runMASTDR = function(sce){
  scaRaw <- FromMatrix(as.matrix(logcounts(sce)),
                       data.frame(colData(sce)),
                       data.frame(rowData(sce))
  )
  cdr <-colSums(assay(scaRaw) > 0)
  colData(scaRaw)$cngeneson <- scale(cdr)
  rowData(scaRaw)$symbolid <- rownames(rowData(sce))
  Dose <- factor(colData(scaRaw)$Dose)
  Dose <- relevel(Dose,"0")
  colData(scaRaw)$Dose = Dose
  cngeneson = colData(scaRaw)$cngeneson
  zlmDose = MAST::zlm(~Dose + cngeneson, scaRaw)
  print(zlmDose)
  print(summary(zlmDose))
  print(names(zlmDose))
  m <- rep(list(list()), times= nlevels(Dose)-1)
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
  
  mast.out = do.call(cbind, lapply(m, data.frame))
  colnames(mast.out) = paste0(colnames(mast.out), levels(Dose)[-1])
  return(mast.out)
}




