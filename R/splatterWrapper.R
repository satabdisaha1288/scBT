splatSimulateDR = function(params = newSplatParams(), 
                           method = c('doseresp'), 
                           sparsify = TRUE, 
                           verbose = TRUE,
                           dose.names = c(0), 
                           dose.prob = c(1),
                           lib.scale = 1.3, ...){
  checkmate::assertClass(params, "SplatParams")
  seed = getParam(params, "seed")
  set.seed(seed)

  # Get the parameters we are going to use
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  nBatches <- getParam(params, "nBatches")
  batch.cells <- getParam(params, "batchCells")
  
  # Set up name vectors
  if (verbose) {message("Creating simulation object...")}
  cell.names <- paste0("Cell", seq_len(nCells))
  gene.names <- paste0("Gene", seq_len(nGenes))
  batch.names <- paste0("Batch", seq_len(nBatches))
  
  # Create SingleCellExperiment to store simulation
  cells =  data.frame(Cell = cell.names)
  rownames(cells) = cell.names
  features = data.frame(Gene = gene.names)
  rownames(features) = gene.names
  sim = SingleCellExperiment(rowData = features, colData = cells,
                              metadata = list(Params = params))
  
  doses = sample(dose.names, nCells, prob = dose.prob, replace = TRUE)
  colData(sim)$Dose = factor(doses, levels = dose.names)
  
  if (verbose) {message("Simulating library sizes...")}
  sim = splatSimLibSizes(sim, params)
  sim$ExpLibSize = sim$ExpLibSize * lib.scale
  
  if (verbose) {message("Simulating Gene means...")}
  sim = splatSimGeneMeans(sim, params)
  
  if (verbose) {message("Simulating batch means...")}
  sim = splatSimBatchCellMeans(sim, params)
  
  ##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (verbose) {message("Simulating dose-response models...")}
  sim = splatSimDoseResponse(sim, params)
  
  if (verbose) {message("Simulating dose-response models, part 2...")}
  sim = splatSimDoseResponseModel(sim, params)

  if (verbose) {message("Simulating BCV...")}
  sim <- splatSimBCVMeans(sim, params)
  
  if (verbose) {message("Simulating counts...")}
  sim <- splatSimTrueCounts(sim, params)
  
  if (verbose) {message("Simulating dropout (if needed)...")}
  sim <- splatSimDropout(sim, params)
  assays(sim)$counts[is.na(assays(sim)$counts)] = 0
  assays(sim)$logcounts = log1p(t(t(assays(sim)$counts)/colSums(assays(sim)$counts))*10000)
  

  return(sim)
}




#EDIT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#' Simulate group differential expression
#'
#' Simulate differential expression. Differential expression factors for each
#' group are produced using \code{\link{getLNormFactors}} and these are added
#' along with updated means for each group. For paths care is taken to make sure
#' they are simulated in the correct order.
#'
#' @param sim SingleCellExperiment to add differential expression to.
#' @param params splatParams object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated differential expression.
#'
#' @name splatSimDE
NULL

#' @rdname splatSimDE
#' @importFrom SummarizedExperiment rowData
splatSimDoseResponse <- function(sim, params) {
  
  nGenes <- getParam(params, "nGenes")
  nGroups <- getParam(params, "nGroups")
  de.prob <- getParam(params, "de.prob")
  de.downProb <- getParam(params, "de.downProb")
  de.facLoc <- getParam(params, "de.facLoc")
  de.facScale <- getParam(params, "de.facScale")
  means.gene <- rowData(sim)$GeneMean
  
  
  de.facs <- getLNormFactors(nGenes, de.prob, de.downProb, de.facLoc, de.facScale)
  rowData(sim)$DE_idx = de.facs
  
  return(sim)
}


splatSimDoseResponseModel = function(sim, params, models.prob = rep(1/6, 6)){
  
  nCells <- getParam(params, "nCells")
  nGenes <- getParam(params, "nGenes")
  nDoses = length(unique(colData(sim)$Dose))
  doses <- colData(sim)$Dose
  dose.names <- c(0,0.01,0.03,0.1,0.3,1,3,10,30)
  exp.lib.sizes <- colData(sim)$ExpLibSize
  batch.means.cell <- assays(sim)$BatchCellMeans
  cell.names <- paste0("Cell", seq_len(nCells))
  gene.names <- paste0("Gene", seq_len(nGenes))
  
  de.ans = rowData(sim)$DE_idx
  de.idx = which(de.ans != 1)
  models.num = floor(models.prob*rep(length(de.idx)))
  models = rep(c('Hill','Power','Linear', 'Exp', 'Exp2', 'ExpB'), models.num)
  
  rowData(sim)$Model = 'Unchanged'
  rowData(sim)[sample(de.idx, sum(models.num), replace = FALSE),'Model'] = models
  
  #Create mat
  dose.means = data.frame(matrix(ncol = length(dose.names), nrow = nGenes))
  colnames(dose.means) = paste0("DEFac",dose.names)
  
  m.list = list()
  #Change to apply
  for (ix in 1:nrow(rowData(sim))){
    if (ix%%100 == 0){
      it.percent = (ix/nrow(rowData(sim)))*100
      message(paste(it.percent, '%----'), appendLF = TRUE)
    }
    if (rowData(sim)[ix, 'Model'] == 'Hill'){
      o = splatsimHill(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'])
      dose.means[ix,] = o$resp
      m.list[[rowData(sim)[ix,'Gene']]] = o
      
    } else if (rowData(sim)[ix, 'Model'] == 'Exp'){
      o = splatsimExp(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'])
      dose.means[ix,] = o$resp
      m.list[[rowData(sim)[ix,'Gene']]] = o
      
    } else if (rowData(sim)[ix, 'Model'] == 'Exp2'){
      o = splatsimExp(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'], power = TRUE)
      dose.means[ix,] = o$resp
      m.list[[rowData(sim)[ix,'Gene']]] = o
      
    } else if (rowData(sim)[ix, 'Model'] == 'ExpB'){
      o = splatsimExpB(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'])
      dose.means[ix,] = o$resp
      m.list[[rowData(sim)[ix,'Gene']]] = o
      
    } else if (rowData(sim)[ix, 'Model'] == 'Power'){
      o = splatsimPower(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'])
      if (is.null(o$resp)){
        print(paste('power_',ix))
        dose.means[ix,] = rep(rowData(sim)[ix,'GeneMean'], 9)
      } else {
        dose.means[ix,] = o$resp
        m.list[[rowData(sim)[ix,'Gene']]] = o
      }
    } else if (rowData(sim)[ix, 'Model'] == 'Linear'){
      o = splatsimLinear(dose.names, mean = rowData(sim)[ix, 'GeneMean'], fc = rowData(sim)[ix, 'DE_idx'])
      dose.means[ix,] = o$resp
      m.list[[rowData(sim)[ix,'Gene']]] = o
    } else {
      dose.means[ix,] = rep(rowData(sim)[ix,'GeneMean'], 9)
    }
  }
  metadata(sim)$modelFits = m.list
  
  rowData(sim) = cbind(rowData(sim), dose.means)
  dose.facs.gene = dose.means/rowData(sim)$GeneMean
  dose.facs.gene[which(is.na(dose.facs.gene)),'DEFac0'] = 1.0000
  cell.facs.gene = as.matrix(dose.facs.gene[, paste0('DEFac', doses)])
  cell.means.gene <- batch.means.cell * cell.facs.gene
  cell.props.gene <- t(t(cell.means.gene) / colSums(cell.means.gene))
  base.means.cell <- t(t(cell.props.gene) * exp.lib.sizes)
  
  colnames(base.means.cell) <- cell.names
  rownames(base.means.cell) <- gene.names
  assays(sim)$BaseCellMeans <- base.means.cell
  
  return(sim)
}

