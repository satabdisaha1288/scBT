% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Main.R
\name{DETest}
\alias{DETest}
\title{Differential expression testing}
\usage{
DETest(
  sce,
  method = "All",
  verbose = TRUE,
  fixed.priors = TRUE,
  return.time = TRUE
)
}
\arguments{
\item{sce}{SingleCellExperiment object with a logcounts assay and a "Dose"
column in the rowData(sce)}

\item{method}{the statistical test(s) to run. Options are "BAYES",
"LRT.linear", "LRT.multiple", "ANOVA", "KW", "WRS", "CHISQ", or "MAST". Leave
empty to run all tests on your dataset.}
}
\value{
a list of data.frames containing the statistical output for each
individual tests performed.
}
\description{
The main function to run any of the statistical tests developed as part of
this package for the analysis of single-cell sequencing data.
}
\details{
Various statistical tests were adapted specifically for dose-response
single-cell/single-nuclei transcriptomic data. The \code{method} enable the
following tests:
\enumerate{
\item BAYES is our new Bayes implementation.. what happens when it goes to
next line?
\item LRT.linear
}
}
\examples{
# Running all statistical tests
DEG.testing <- DETest(sim)

}
\references{
Authors... 2021
}
