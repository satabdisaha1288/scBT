# Example BMD analysis
# Jessica Ewald
# March 10, 2021

# The curve fitting functions assume that the data matrix has already been pre-filtered
# Thus, curves will be fit to every row in the data

source("C:\\Users\\15177\\OneDrive - Michigan State University\\Documents\\CodeSnippets\\fastbmdR\\fastbmdR\\fastbmdR_main.R")
source("C:\\Users\\15177\\OneDrive - Michigan State University\\Documents\\CodeSnippets\\fastbmdR\\fastbmdR\\fastbmdR_utils.R")

# read in example data object
data <- readRDS("C:\\Users\\15177\\OneDrive - Michigan State University\\Documents\\CodeSnippets\\fastbmdR\\fastbmdR\\data.rds")

# set up parameters (model choice and associated doses)
models = c("Exp2","Exp3","Exp4","Exp5","Poly2","Lin","Power","Hill")
dose = c(rep(0,5), rep(25,5), rep(100,5), rep(200,5), rep(300,5), rep(400,5))

# set the number of threads for parallel computing (defaults to 1)
ncpus = 1

# perform curve fitting
res <- PerformCurveFitting(data = data, dose = dose, ncpus = ncpus, models = models)

# filter to select best fit model for each gene
res <- FilterDRFit(res, lof.pval = 0.1)

# perform benchmark dose calculation
bmd.res <- PerformBMDCalc(res, ncpus = ncpus, num.sds = 1, sample.mean = TRUE)

# Results are returned for each gene in res. See the last few columns of bmd.res 
# for whether a gene BMD passed various quality criteria.

