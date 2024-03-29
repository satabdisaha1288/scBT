% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bayes_test.R
\name{calcPosteriorProbNull}
\alias{calcPosteriorProbNull}
\title{Title Calculate FDR adjusted threshold for posterior probabilities of the null distribution
Newton, Michael A., Amine Noueiry, Deepayan Sarkar, and Paul Ahlquist. "Detecting differential gene expression with a semiparametric hierarchical mixture method." Biostatistics 5, no. 2 (2004): 155-176.
Tadesse, Mahlet G., Joseph G. Ibrahim, Robert Gentleman, Sabina Chiaretti, Jerome Ritz, and Robin Foa. "Bayesian Error‐in‐Variable Survival Model for the Analysis of GeneChip Arrays." Biometrics 61, no. 2 (2005): 488-497.}
\usage{
calcPosteriorProbNull(
  bayes.out,
  kappa = seq(0.01, 1, 0.01),
  alpha = c(0.01, 0.05)
)

calcPosteriorProbNull(
  bayes.out,
  kappa = seq(0.01, 1, 0.01),
  alpha = c(0.01, 0.05)
)
}
\arguments{
\item{kappa}{: sequence of thresholding values preferably; kappa<-seq(0.01,1,0.01)}

\item{alpha}{: Desired FDR (False Discovery Rate) level}

\item{DETest_output_Bayes}{: named (gene-names) vector of Bayes Factor values obtained from the Bayes Test}
}
\value{
max_kappa : Threshold value: Reject all hypothesis with posterior null probabilies less than max_kappa

max_kappa : Threshold value: Reject all hypothesis with posterior null probabilies less than max_kappa
}
\description{
Title Calculate FDR adjusted threshold for posterior probabilities of the null distribution
Newton, Michael A., Amine Noueiry, Deepayan Sarkar, and Paul Ahlquist. "Detecting differential gene expression with a semiparametric hierarchical mixture method." Biostatistics 5, no. 2 (2004): 155-176.
Tadesse, Mahlet G., Joseph G. Ibrahim, Robert Gentleman, Sabina Chiaretti, Jerome Ritz, and Robin Foa. "Bayesian Error‐in‐Variable Survival Model for the Analysis of GeneChip Arrays." Biometrics 61, no. 2 (2005): 488-497.

Title Calculate FDR adjusted threshold for posterior probabilities of the null distribution
Newton, Michael A., Amine Noueiry, Deepayan Sarkar, and Paul Ahlquist. "Detecting differential gene expression with a semiparametric hierarchical mixture method." Biostatistics 5, no. 2 (2004): 155-176.
Tadesse, Mahlet G., Joseph G. Ibrahim, Robert Gentleman, Sabina Chiaretti, Jerome Ritz, and Robin Foa. "Bayesian Error‐in‐Variable Survival Model for the Analysis of GeneChip Arrays." Biometrics 61, no. 2 (2005): 488-497.
}
\author{
Satabdi Saha

Satabdi Saha
}
