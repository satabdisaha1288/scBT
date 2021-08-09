#' Performs a genewise Wilcoxon Rank Sum test on a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object with a logcounts assay
#' and Dose column in the cell metadata
#'
#' @return a vector of p values from the Wilcoxon Rank Sum test
#'
#' @import Seurat
#' @export
runSeuratBimod <- function(sce){
  seurat.data <- Seurat::as.Seurat(sce, counts = 'counts', data = 'logcounts')
  Seurat::Idents(seurat.data) <- seurat.data@meta.data$Dose
  factorDose <- as.factor(seurat.data@meta.data$Dose)
  res.list <- list()
  for (d in 2:length(levels(factorDose))){
    doseName <- levels(factorDose)[d]
    res.list[[doseName]] <- Seurat::FindMarkers(seurat.data,
                                                ident.1 = levels(factorDose)[1],
                                                ident.2 = levels(factorDose)[d],
                                                logfc.threshold = -Inf,
                                                test.use = 'bimod',
                                                min.pct = -Inf,
                                                min.cell.features = -Inf,
                                                min.cell.groups = -Inf
                                                )
    colnames(res.list[[doseName]]) <- paste0(doseName,
                                             '_',
                                             colnames(res.list[[doseName]])
                                             )



  }
  seur.bimod <- Reduce(merge, lapply(res.list, function(x) data.frame(x, rn = row.names(x))))
  rownames(seur.bimod) <- seur.bimod$rn
  seur.bimod <- seur.bimod[rownames(seurat.data), ]
  p.val <- seur.bimod[,grepl('p_val$', colnames(seur.bimod))]
  seur.bimod$adjusted.p <- apply(p.val,
                                 1, function(x) p.adjust(min(x), 'fdr')
  )

  #p.val.adj <- seur.bimod[,grepl('p_val_adj', colnames(seur.bimod))]
  #seur.bimod$adjusted.p <- apply(p.val.adj,
  #                               1, function(x) min(x))

  return(seur.bimod)
}
