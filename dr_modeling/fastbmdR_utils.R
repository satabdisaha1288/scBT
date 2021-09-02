##### Other functions for fastbmdR ######

### Hill model and starting values
formHill <- as.formula(signal ~ c + (d - c) / (1 + (dose/e)^b ) )
startvalHillnls2 <- function(x, y, xm, ym, increase) 
  # requires the definition of increase from min and max values
  # which is the first one
  # inputs
  # - x values of the dose
  # - y values the corresponding signal
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
  # 
{
  maxi <- max(y, na.rm = TRUE)
  mini <- min(y, na.rm = TRUE)
  ampl <- maxi - mini
  
  # inflate maxi and mini so as all values are strictly inside the interval [mini; maxi]
  maxi <- maxi + 0.001 * ampl
  mini <- mini - 0.001 * ampl
  
  # initial value of c
  c <- ifelse(increase, maxi, mini) 
  # initial value of d
  d <-ifelse(increase, mini, maxi) 
  # initial value of e and b from regression
  yreg <- log((d - c) / (y[x!=0] - c) - 1)
  xreg <- log(x[x!=0])
  reg <- lm(yreg ~ xreg)
  b <- reg$coefficients[2]
  e <- reg$coefficients[1] / (-b)
  startval <- list(b = b, c = c, d = d, e = e)
}


### Exp2 model
# define the model formula
formExp2 <- as.formula(signal ~ e*exp(b*dose))
# get starting values
startvalExp2 <- function(xm, ym)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
{
  # initial value of a
  e <- ym[1]
  
  # transform y for regression
  yreg <- log(ym[xm != 0])
  reg <- lm(yreg ~ xm[xm != 0])
  
  # estimate slope from regression
  b <- coef(reg)[2]
  
  startval <- list(e = e, b = b)
}


### Exp3 model
# define the model formula
formExp3 <- as.formula(signal ~ e*(exp(sign(b)*(abs(b)*dose)^d)))

# get starting values
startvalExp3 <- function(xm, ym)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
{
  # initial value of e
  e <- ym[1]
  
  # transform y for regression
  yreg <- log(ym[xm != 0])
  reg <- lm(yreg ~ xm[xm != 0])
  
  # estimate b and d from regression
  b <- coef(reg)[2]
  
  d <- (exp(coef(reg)[1]))/e
  
  startval <- list(e = e, b = b, d = d)
}


### Exp4 model
formExp4 <- as.formula(signal ~ e*(c - (c-1)*exp((-1)*b*dose)))

# get starting values
startvalExp4 <- function(xm, ym, ad.dir)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
{
  # initial value of e
  e <- ym[1]
  
  # initial value of c (since the asymptote is always to the right, we can calculate c based on a and the max/min val)
  if(ad.dir == TRUE){
    c <- max(ym)/e + 0.001*(max(ym)/e)
  }else{
    c <- min(ym)/e - 0.001*(min(ym)/e)
  }
  
  # initial value of b
  yreg <- log((ym - e*c)/(e-e*c))
  reg <- lm(yreg ~ xm)
  
  b <- abs(coef(reg)[2])
  
  startval <- list(e = e, b = b, c = c)
}


#### Exp5 ####
formExp5 <- as.formula(signal ~ e*(c - (c-1)*exp((-1)*(b*dose)^d)))

# get starting values
startvalExp5 <- function(xm, ym, ad.dir)
  # inputs
  # - xm unique values of the dose (sorted by dose)
  # - ym means of the signal at each value of xm (sorted by dose)
{
  # initial value of e
  e <- ym[1]
  
  # initial value of c
  if(ad.dir){
    c <- max(ym)/e + 0.001*(max(ym)/e)
  }else{
    c <- min(ym)/e - 0.001*(min(ym)/e)
  }
  
  # initial value of b and d
  yreg <- log((ym - e*c)/(e-e*c))
  reg <- lm(yreg ~ xm)
  
  b <- abs(coef(reg)[2])
  d <- exp(coef(reg)[1])
  
  startval <- list(e = e, b = b, c = c, d = d)
}


#### power ####
formPow <- as.formula(signal ~ e + b*(dose^c))

# get starting values
# Power model has trouble converging, so we run nls() twice (once here, once in main function)
startvalPow <- function(xm, ym, ad.dir, dset){
  
  require(dplyr)
  
  if(ad.dir){
    
    e <- min(ym) - 0.001*min(ym)
    
    yreg <- log(ym[xm!=0] - e)[-1]
    xreg <- log(xm[xm!=0])[-1]
    
    reg <- lm(yreg ~ xreg)
    
    c <- max(c(coef(reg)[2],1))
    b <- exp(coef(reg)[1])
    
  } else {
    
    e <- max(ym) + 0.001*max(ym)
    
    yreg <- log((ym[xm!=0] - e)*(-1))
    xreg <- log(xm[xm!=0])
    
    reg <- lm(yreg ~ xreg)
    
    c <- max(c(coef(reg)[2],1))
    b <- exp(coef(reg)[1])*(-1)
  }
  
  start.1 <- list(e = e, b = b, c = c)
  
  Pow <- suppressWarnings(try(nls(formPow, start = start.1, data = dset,
                                  lower = c(1, -Inf, 0.999), control = nls.control(maxiter = 500),
                                  upper = c(Inf, Inf, 18), algorithm = "port",
                                  weights = 1/(signal^2)), silent = TRUE))
  
  if (!inherits(Pow, "try-error")){
    startval <- coef(Pow) %>% as.list()
  } else {
    startval <- list(e = e, b = b, c = c)
  }
  
}

#### function to get bmd results
bmdres <- function(fit){
  
  bmd.ci <- suppressMessages(confint(fit, "bmd"))
  bmd.mean <- coef(fit)[1] # get bmd estimate from fit
  c(bmd.mean, bmd.ci)
  
}


#### I use the pureErrorAnova function from alr3. alr3 is now deprecated, so extracted these 
#### lines from the alr3 source code in the CRAN archive
pureErrorAnova <- function(mod){UseMethod("pureErrorAnova")}
pureErrorAnova.lm <- function(mod) {
  if (is.null(mod$model)) mod <- update(mod, model=TRUE)
  p <- dim(mod$model)[2] -1
  mod$model$Lack.of.Fit <-
    factor(randomLinComb(model.matrix(mod), 101319853))
  aov1 <- anova(mod)
  #set.seed(save.seed) # restore random number seed
  if (length(levels(mod$model$Lack.of.Fit)) == length(mod$model$Lack.of.Fit))
    aov1 else {
      aov2 <- anova(lm(mod$model[ , 1]~mod$model$Lack.of.Fit, weights=weights(mod)))
      rrow <- dim(aov1)[1]
      aov2[1, 1] <- aov1[rrow, 1]-aov2[2, 1]
      aov2[1, 2] <- aov1[rrow, 2]-aov2[2, 2]
      aov2[1, 3] <- aov2[1, 2]/aov2[1, 1]
      aov1[1:(rrow-1), 4] <- aov1[1:(rrow-1), 3]/aov2[2, 3]
      aov2[1, 4] <- aov2[1, 3]/aov2[2, 3]
      row.names(aov2) <- c(" Lack of fit", " Pure Error")
      aov <- rbind(aov1, aov2)
      aov[ , 5] <- pf(aov[ , 4], aov[ , 1], aov2[2, 1], lower.tail=FALSE)
      aov
    }}


randomLinComb <- function(X, seed=NULL) {UseMethod("randomLinComb")}

randomLinComb.default <- function(X, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  std <- function(x){
    s <- sd(x)
    if( s > 0) (x-mean(x))/s else x}
  as.vector(apply(X, 2, std)%*% as.vector(2*rnorm(dim(X)[2])-1) )
}

randomLinComb.lm <- function(X, ...) {
  randomLinComb(model.matrix(X), ...)}

randomLinComb.lm <- function(X, seed=NULL) {
  if(is.null(X$model)) X <- update(X, model=TRUE)
  randomLinComb(X$model[ , -1], seed=seed)}

