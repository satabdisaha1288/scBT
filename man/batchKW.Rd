% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KW_test.R
\name{batchKW}
\alias{batchKW}
\title{Performs a genewise Kruskal Wallis test on a SingleCellExperiment object}
\usage{
batchKW(sce)
}
\arguments{
\item{sce}{SingleCellExperiment object with a logcounts assay
and Dose column in the cell metadata}
}
\value{
a vector of p values from the Kruskal Wallis test
}
\description{
Performs a genewise Kruskal Wallis test on a SingleCellExperiment object
}
\author{
Jack Dodson
}
