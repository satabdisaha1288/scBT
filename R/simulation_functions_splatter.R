splatsimHill = function(doses, mean = 1, fc = 1.5){
  #' Simulate expression that follows a Hill model
  #' response = gamma + (V * doses^n)/(k^n + doses^n)
  #' 
  #' @param doses A vector of doses to model
  #' @param mean.range A vector of means obtain from realData
  #' @param fc.range a vector with the minimum and maximum |fold-change| (e.g., c(1.5, 5))
  #' @param downregulated set a TRUE to model repression instead of induction.
  #' @return a list of the model fit parameters including
  #' @examples 
  #' simHill(c(0,1,3,10,30), realData$mean)
  #' @export
  
  k.vals = c()
  for (d in 2:length(doses)){
    k.vals = c(k.vals, seq(doses[d-1], doses[d], length.out = 100))
  }
  gamma = mean
  V = (fc*mean)-mean
  n = runif(1, 1, 100)
  k = sample(k.vals, 1)
    
  resp = modelHill(doses, gamma, V, n, k)
  
  return(list(resp = resp, gamma = gamma, fc = fc, V = V, n = n, k = k))
}


splatsimExp = function(doses, mean, fc, power = FALSE){
  #' Simulate expression that follows an exponential model
  #' response = response = a * exp(sign * (b * dose)^d)
  #' 
  #' @param doses A vector of doses to model
  #' @param mean.range A vector of means obtain from realData
  #' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
  #' @param downregulated set a TRUE to model repression instead of induction
  #' @param power set as TRUE to randomize variable d (default = 1)
  #' @return a list of the model fit parameters including
  #' @examples 
  #' simExp(c(0,1,3,10,30), realData$mean)
  #' @export
  #' 
  #' 
  fc.max = 1E100
  mult.fac = 0.01
  iter = 0
  while (fc.max > fc + fc*mult.fac | fc.max < fc - fc*mult.fac){
    iter = iter + 1
    a = mean
    if (fc < 1){
      b = runif(1, -0.1, 0)
    } else {
      b = runif(1, 0, 0.1)
    }
    if (power){
      d = runif(1, 0,10)
    } else {
      d = 1
    }
    
    resp = modelExp(doses, a, b, d)
    fc.max = max(resp)/min(resp)
    if (fc < 1){
      fc.max = 1/fc.max
    }
    if (iter%%10000 == 0){
      {message("Relaxing fold-change range criteria...")}
      mult.fac = mult.fac * 10
    }
  }
  return(list(resp = resp, a = a, fc = fc, b = b, d = d))
}

splatsimExpB = function(doses, mean, fc = 1.5){
  #' Simulate expression that follows an exponential model (2 or 3)
  #' response = a(c-(c-1)exp???(-1 (bdose)^d ))
  #' 
  #' @param doses A vector of doses to model
  #' @param mean.range A vector of means obtain from realData
  #' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
  #' @param downregulated set a TRUE to model repression instead of induction
  #' @return a list of the model fit parameters including
  #' @examples 
  #' simExpB(c(0,1,3,10,30), realData$mean)
  #' @export

  a = mean
  b = runif(1, 0, 1) #CHANGE
  c = fc
  d = runif(1, 0, 4)
    
    
  resp = modelExpB(doses, a, b, c, d)
  return(list(resp = resp, a = a, fc = fc, b = b, c = c, d = d))
}


splatsimPower = function(doses, mean, fc, downregulated = FALSE){
  #' Simulate expression that follows an power model
  #' response = gamma + beta * doses^delta
  #' 
  #' @param doses A vector of doses to model
  #' @param mean.range A vector of means obtain from realData
  #' @param fc.range a vector with the minimum and maximim |fold-change| (e.g., c(1.5, 5))
  #' @param downregulated set a TRUE to model repression instead of induction
  #' @return a list of the model fit parameters including
  #' @examples 
  #' simPower(c(0,1,3,10,30), realData$mean)
  #' @export
  
  fc.max = 1E100
  mult.fac = 0.01
  iter = 0
  while (fc.max > fc + fc*mult.fac | fc.max < fc - fc*mult.fac){
    iter = iter + 1
    gamma = mean
    if (fc < 1){
      beta = runif(1, -10^-floor(abs(log(mean))), 0)
    } else {
      beta = runif(1, 0, 10^-floor(abs(log(mean))))
    }
    
    delta = runif(1, 0, 5)
    
    resp = modelPower(doses, gamma, beta, delta)
    fc.max = max(resp)/min(resp)
    if (fc < 1){
      fc.max = 1/fc.max
    }
    if (iter%%10000 == 0){
      {message("Relaxing fold-change range criteria...")}
      mult.fac = mult.fac * 10
    }
  }
  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta, delta = delta))
}

splatsimLinear = function(doses, mean, fc){
  gamma = mean
  beta = ((mean * fc) - mean)/max(doses)
  resp = modelPolynomial(doses, gamma, beta)

  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta))
}















simPolynomial = function(doses, mean.range, fc.range = c(1.3, 5), downregulated = FALSE){
  resp = max(mean.range) * 2 # To initiate while loop
  polyN = sample(c(2,3,4), 1)
  
  while (max(resp) > max(mean.range) | min(resp) < min(mean.range)){
    gamma = sample(mean.range, 1)
    if (downregulated){
      fc = 1/runif(1, fc.range[1], fc.range[2])  
    } else {
      fc = runif(1, fc.range[1], fc.range[2])
    }
    
    beta = runif(polyN,0,0.001)
    resp = modelPolynomial(doses, gamma, beta)
  }
  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta))
}


# Unchanged genes
simUnchanged = function(doses, realVec){
  n.doses = length(doses)
  resp = sample(realVec, n.doses, replace = FALSE)
  return(list(resp = resp))
}

simUnchanged2 = function(doses, realVec, sd){
  resp = rnorm(length(doses), mean = sample(realVec, 1), sd = sd)
  return(list(resp = resp))
}


runSimulation = function(n = 50, start_id = 1, class = 'Hill', realData = NULL, downregulated = FALSE, power = FALSE, fc.range = c(1.5, 5), sd = 1){
  num.sim = n
  
  df = data.frame(item = c(), doses = c(), resp = c())
  par = list()
  
  items = seq(start_id + 1, start_id + num.sim, by = 1)
  for (i in items){
    if (class == 'Hill'){ par[[paste('gene',i,sep='')]] = simHill(sim.doses, realData$mean, downregulated = downregulated, fc.range = fc.range) }
    if (class == 'Exp'){ par[[paste('gene',i,sep='')]] = simExp(sim.doses, realData$mean, downregulated = downregulated, power = power, fc.range = fc.range) }
    if (class == 'ExpB'){ par[[paste('gene',i,sep='')]] = simExpB(sim.doses, realData$mean, downregulated = downregulated, fc.range = fc.range) }
    if (class == 'Power'){ par[[paste('gene',i,sep='')]] = simPower(sim.doses, realData$mean, downregulated = downregulated, fc.range = fc.range) }
    if (class == 'Linear'){ par[[paste('gene',i,sep='')]] = simLinear(sim.doses, realData$mean, downregulated = downregulated) }
    if (class == 'Unchanged'){ par[[paste('gene',i,sep='')]] = simUnchanged2(sim.doses, realData$mean, sd) }
    
    temp.df = data.frame(item = rep(i, num.doses), doses = sim.doses, resp = par[[paste('gene',i,sep='')]]$resp)
    par[[paste('gene',i,sep='')]]$resp = NULL
    df = rbind(df, temp.df)
    start_id = i
  }
  return(list(data = df, parameters = data.table::rbindlist(par, fill = TRUE), stop_id = i))
}

convert2snseq = function(reference, modeledData, ncells, sd.scale = 1, p.zero = 2, by.doses = FALSE){
  #' 
  #' 
  #' @param reference = realData
  #' @param simulatedMeans
  #' @param..
  
  # Find the data point with the closest mean in the real data and return its index.
  # The out data frame joins the simulated mean with estimates from the real dataset. 
  diff.table = vapply(modeledData$resp, function(x) x - reference$mean, numeric(length(reference$mean)))
  indx = apply(abs(diff.table), 2, which.min)
  out = cbind(modeledData, reference[indx, ])
  out = out[, c('item','doses','resp','mean', 'sd','percent_zero')]
  out$resp = as.numeric(out$resp)
  out$sd = as.numeric(out$sd) * sd.scale
  if (p.zero > 1){
    out$percent_zero = as.numeric(out$percent_zero)
  } else {
    out$percent_zero = p.zero
  }
  
  modNorm = t(apply(out, 1, function(x) adjustSampling(as.numeric(x[3]), as.numeric(x[5]), as.numeric(x[6]), ncells)))  
  
  newout = cbind(out[,c('item', 'doses')], data.frame(modNorm))
  newout2 = reshape2::melt(newout, id.vars = c('item','doses'))
  return(list(diftab = diff.table, norm = modNorm, out1 = newout, out2 = newout2))
}


convert2snseq2 = function(reference, modeledData, ncells, sd.scale = 1, p.zero = 2, fitPar = c(5, 0.5)){
  #' 
  #' 
  #' @param reference = realData
  #' @param simulatedMeans
  #' @param..
  
  # Find the data point with the closest mean in the real data and return its index.
  # The out data frame joins the simulated mean with estimates from the real dataset. 
  diff.table = vapply(modeledData$resp, function(x) x - reference$mean, numeric(length(reference$mean)))
  indx = apply(abs(diff.table), 2, which.min)
  out = cbind(modeledData, reference[indx, ])
  #out = out[, c('item','doses','resp','mean', 'sd','percent_zero')]
  out$resp = as.numeric(out$resp)
  out$sd = as.numeric(out$median.sd) * sd.scale
  if (p.zero > 1){
    out$median.zero = as.numeric(out$median.zero)
  } else {
    out$median.zero = p.zero
  }
  #print(head(out, n=50))
  modNorm = t(apply(out, 1, function(x) adjustSampling2(as.numeric(x[3]),
                                                       c(as.numeric(x[11]), as.numeric(x[12]), as.numeric(x[14])),
                                                       c(as.numeric(x[7]), as.numeric(x[8])),
                                                       ncells,
                                                       fitPar)))
  #print(head(modNorm, n = 50))
    
  newout = cbind(out[,c('item', 'doses')], data.frame(modNorm))
  newout2 = reshape2::melt(newout, id.vars = c('item','doses'))
  return(list(diftab = diff.table, norm = modNorm, out1 = newout, out2 = newout2))
}



adjustSampling = function(mean, sd, percent.zero, ncells){
  #Sample count values from a normal distribution
  vec = rnorm(ncells, mean, sd)
  vec[vec < 0] = 0
  num.zero = sum(vec == 0)
  add.zero = round((ncells * percent.zero) - num.zero)
  replace.ind = sample(which(vec > 0), add.zero)
  vec[replace.ind] = 0
  return(vec)
}

adjustSampling2 = function(mean, sd.range, percent.zero.range, ncells, fitPar){
  k = -1 * fitPar[1]
  x = fitPar[2]
  i = 0
  #Sample count values from a normal distribution
  #vec = rnorm(ncells, mean, runif(1, min = sd.range[1]*3, max = sd.range[2]*3))
  #m.sd = 0.476169*exp(-1/2*(mean-1.117201)^2/0.473428^2)
  #vec = rnorm(ncells, mean, runif(1, min = sd.range[1], max = sd.range[2] + (2*1/exp(mean)^3.2)))
  vec = rnorm(ncells, mean, runif(1, min = sd.range[1], max = sd.range[2]))
  vec[vec < 0] = 0
  
  D = 1/(1 + exp(-k*(log(vec)-x)))
  D[which(D>1)] = 1
  D[which(D<0)] = 0
  dropout.indicator = rbinom(ncells, 1, D)
  dropout.indicator[is.na(dropout.indicator)] = 0  
  new.vec = vec * dropout.indicator

  #num.zero = sum(vec == 0)
  #add.zero = round((ncells * percent.zero) - num.zero)
  #replace.ind = sample(which(vec > 0), add.zero)
  #vec[replace.ind] = 0
  return(new.vec)
}


importSeu = function(path) {
  seu = readRDS(path)
  metadata = seu@meta.data
  norm = seu@assays$RNA@data
  sct = seu@assays$SCT@data
  return(list(meta = metadata, norm = norm, sct = sct))
}


fitStats = function(expression, name) {
  ks = 'NA'
  mean = 'NA'
  sd = 'NA'
  category = 'fail'
  df.all = data.frame(count = expression)
  df.nonzero = df.all %>% filter(count > 0)
  mean_nonzero  = mean(df.nonzero$count)
  var_nonzero = var(df.nonzero$count)
  mean_real = mean(df.all$count)
  var_real = var(df.all$count)
  num_zeros = NROW(df.all) - NROW(df.nonzero)
  percent.zero = num_zeros/NROW(df.all)
  
  if (length(unique(df.nonzero$count)) > 2){
    tryCatch({
      fit_n = fitdist(df.nonzero %>% pull(count), 'norm')
      fit_n.stat = gofstat(fit_n)
      ks = fit_n.stat$ks
      mean = fit_n$estimate[1]
      sd = fit_n$estimate[2]
      category = 'pass'},
      error = function(c){
        print('fail')
        print('error: continuing with remaining genes')
      }
      )
  } 
  return(c(name, mean_nonzero, var_nonzero, mean_real, var_real, percent.zero, ks, mean, sd, category))
}


realParams = function(dataset, ncores, assay = 'Norm', metacol = NULL, metacriteria = '', funpath = './simulation_functions.R') {
  print(paste('Working on dataset ', metacriteria, sep = ''))
  if (!is.null(metacol)){
    dataset$norm = dataset$norm[, which(dataset$meta[metacol] == metacriteria)]
    dataset$sct = dataset$sct[, which(dataset$meta[metacol] == metacriteria)]
  }
  
  c1 = makeCluster(ncores)
  registerDoParallel(c1)
  
  if (assay == 'Norm'){
    tab = foreach(i = 1:NROW(dataset$norm), .packages = c("Seurat","dplyr", "fitdistrplus")) %dopar% {
      source(funpath)
      fitStats(dataset$norm[i,], rownames(dataset$norm)[i])
    }
  }
  
  if (assay == 'SCT'){
    tab = foreach(i = 1:NROW(dataset$sct), .packages = c("Seurat","dplyr", "fitdistrplus")) %dopar% {
      source(funpath)
      fitStats(dataset$sct[i,], rownames(dataset$sct)[i])
    }
  }
  stopCluster(c1)
  ans = data.frame(Reduce("rbind", tab))
  colnames(ans) = c('Symbol', 'mean_nonzero', 'var_nonzero', 'mean_real', 'var_real', 'percent_zero', 'ks', 'mean', 'sd', 'category')
  ans$mean = suppressWarnings(as.numeric(ans$mean))
  ans$percent_zero = suppressWarnings(as.numeric(ans$percent_zero))
  ans$ks = suppressWarnings(as.numeric(ans$ks))
  ans$mean_nonzero = as.numeric(ans$mean_nonzero)
  ans$mean_real = as.numeric(ans$mean_real)
  ans$var_nonzero = as.numeric(ans$var_nonzero)
  ans$var_real = as.numeric(ans$var_real)
  ans$sd = as.numeric(ans$sd)
  ans$meta_var = metacriteria
  return(ans)
}

reformatForExport = function(simData, GeneMeta){
  tempData = simData$out2[order(simData$out2$doses),]
  tempData$barcode = paste(tempData$variable, tempData$doses, sep = '_')
  simCellMeta = unique(tempData[,c('barcode', 'doses')])
  rownames(simCellMeta) = simCellMeta$barcode
  colnames(simCellMeta) = c('barcode', 'Dose')
  tempData = tempData[,c('item', 'barcode', 'value')]
  tempData = reshape2::dcast(tempData, item ~ barcode)
  rownames(tempData) = paste('gene', tempData$item, sep = '')
  tempData = tempData[,-1]
  simDataFinal = tempData[,simCellMeta$barcode]
  
  SimGeneMeta = data.frame(Gene = rownames(simDataFinal))
  SimGeneMeta = cbind(SimGeneMeta, GeneMeta)
  
  return(list(norm = simDataFinal, meta = simCellMeta, geneMeta = SimGeneMeta))
}