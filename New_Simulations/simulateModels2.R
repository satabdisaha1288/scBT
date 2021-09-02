# Simulate the data
simHill = function(doses, mean.range, fc.range = c(1.3, 5), downregulated = FALSE){
  resp = max(mean.range) * 2 # To initiate while loop
  
  while (max(resp) > max(mean.range) | min(resp) < min(mean.range)){
    gamma = sample(mean.range, 1)
    if (downregulated){
      fc = 1/runif(1, fc.range[1], fc.range[2])  
    } else {
      fc = runif(1, fc.range[1], fc.range[2])
    }
    V = (log(fc) + gamma) - gamma
    n = runif(1, 1, 100)
    k = runif(1, min(doses), max(doses))
    
    resp = modelHill(doses, gamma, V, n, k)
  }
  return(list(resp = resp, gamma = gamma, fc = fc, V = V, n = n, k = k))
}

simExp = function(doses, mean.range, fc.range = c(1.3, 5), downregulated = FALSE, power = FALSE){
  resp = max(mean.range) * 2 # To initiate while loop
  fc = fc.range[2] * 2
  
  while (max(resp) > max(mean.range) | min(resp) < min(mean.range) | fc < fc.range[1] | fc > fc.range[2]){
    a = sample(mean.range, 1)
    if (downregulated){
      b = runif(1, -1, 0)
    } else {
      b = runif(1, 0, 1)
    }
    if (power){
      d = runif(1, 0,10)
    } else {
      d = 1
    }
    
    resp = modelExp(doses, a, b, d)
    fc = exp(max(resp) - min(resp))
  }
  return(list(resp = resp, a = a, fc = fc, b = b, d = d))
}


simExpB = function(doses, mean.range, fc.range = c(1.3, 5), downregulated = FALSE){
  resp = max(mean.range) * 2 # To initiate while loop
  
  while (max(resp) > max(mean.range) | min(resp) < min(mean.range)){
    a = sample(mean.range, 1)
    
    if (downregulated){
      fc = 1/runif(1, fc.range[1], fc.range[2])  
    } else {
      fc = runif(1, fc.range[1], fc.range[2])
    }
    
    b = runif(1, 0, 1) #CHANGE
    c = log((exp(a) -1) * fc)/a
    d = runif(1, 0, 4)
    
    
    resp = modelExpB(doses, a, b, c, d)
    resp[is.nan(resp)] = max(mean.range) * 2
  }
  return(list(resp = resp, a = a, fc = fc, b = b, c = c, d = d))
}


simPower = function(doses, mean.range, fc.range = c(1.3, 5), downregulated = FALSE){
  resp = max(mean.range) * 2 # To initiate while loop
  fc = fc.range[2] * 2
  
  while (max(resp) > max(mean.range) | min(resp) < min(mean.range) | fc < fc.range[1] | fc > fc.range[2]){
    gamma = runif(1, mean.range[1], mean.range[2])
    if (downregulated){
      beta = runif(1, -1, 0)
    } else {
      beta = runif(1, 0, 1)
    }
    
    delta = runif(1, 0, 5)

    resp = modelPower(doses, gamma, beta, delta)
    fc = exp(max(resp) - min(resp))
  }
  return(list(resp = resp, gamma = gamma, fc = fc, beta = beta, delta = delta))
}


simLinear = function(doses, mean.range, fc.range = c(1.3, 5), downregulated = FALSE){
  resp = max(mean.range) * 2 # To initiate while loop
  
  while (max(resp) > max(mean.range) | min(resp) < min(mean.range)){
    gamma = sample(mean.range, 1)
    if (downregulated){
      fc = 1/runif(1, fc.range[1], fc.range[2])  
    } else {
      fc = runif(1, fc.range[1], fc.range[2])
    }
    
    beta = (log1p(exp((log(fc) + log(exp(gamma)-1))))-gamma)/(max(doses)-min(doses))
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


simLinear = function(doses, mean.range, fc.range = c(1.3, 5), downregulated = FALSE){
  resp = max(mean.range) * 2 # To initiate while loop
  
  while (max(resp) > max(mean.range) | min(resp) < min(mean.range)){
    gamma = sample(mean.range, 1)
    if (downregulated){
      fc = 1/runif(1, fc.range[1], fc.range[2])  
    } else {
      fc = runif(1, fc.range[1], fc.range[2])
    }
    
    beta = (log1p(exp((log(fc) + log(exp(gamma)-1))))-gamma)/(max(doses)-min(doses))
    resp = modelPolynomial(doses, gamma, beta)
  }
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


runSimulation = function(n = 50, start_id = 1, class = 'Hill', realData = NULL, downregulated = FALSE, power = FALSE, fc.range = c(1.5, 5)){
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
    if (class == 'Unchanged'){ par[[paste('gene',i,sep='')]] = simUnchanged(sim.doses, realData$mean) }
    
    temp.df = data.frame(item = rep(i, num.doses), doses = sim.doses, resp = par[[paste('gene',i,sep='')]]$resp)
    par[[paste('gene',i,sep='')]]$resp = NULL
    df = rbind(df, temp.df)
    start_id = i
  }
  return(list(data = df, parameters = data.table::rbindlist(par, fill = TRUE), stop_id = i))
}

convert2snseq = function(reference, modeledData, ncells, sd.scale = 1, p.zero = 2){
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
  out = out[, c('item','doses','resp','mean', 'sd','percent.zero')]
  out$resp = as.numeric(out$resp)
  out$sd = as.numeric(out$sd) * sd.scale
  if (p.zero > 1){
    out$percent.zero = as.numeric(out$percent.zero)
  } else {
    out$percent.zero = p.zero
  }
  
  modNorm = t(apply(out, 1, function(x) adjustSampling(as.numeric(x[3]), as.numeric(x[5]), as.numeric(x[6]), ncells)))  
  
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