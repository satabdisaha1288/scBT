# snSeq-DGE
This package implements a new Bayesian test for detecting differential gene expression over multiple dose groups in single cell gene 
expression studies. 

## DGE testing
scBT is an R package for differential gene expression (DGE) analysis in multiple group study designs for single-cell RNA sequencing data. scBT contains a new Bayesian test of the same name designed along with 9 other benchmarking algorithms frequently used for the DGE analysis in multiple group experimental designs. The tests present in scBT are:

* Seurat Bimod ( Two sample test of mean for zero inflated continuous data)
* Wilcoxon Rank Sum Test (Non-parametric two sample test for mean)
* ANOVA ( Parametric K-sample test of mean for samples from a normal distribution)
* KW (Non-parametric K-sample test of mean)
* limma-trend
* LRT-multiple(k-sample test of mean for zero inflated continuous data)
* LRT-Linear( Regression model based test of DGE for zero inflated continuous data)
* MAST (Regression model based test of DGE for zero inflated continuous data with Bayesian estimation)

## Installation
The developmental version of scBT can be installed from Github:
```{r}
library("devtools")
devtools::install_github("zacharewskilab/scBT")
```

## Getting Started
Once installed the best place to get started is the vignette. The Quickstart vignette can be accessed as:

```
library(scBT)
DETest(sce, method = 'BAYES')
```

## Citing scBT
Please cite ["Nault, R., Saha, S., Bhattacharya, S., Dodson, J., Sinha, S., Maiti, T. and Zacharewski, T. (2021). Benchmarking of a Bayesian single cell RNAseq differential gene expression test for dose-response study designs. bioRxiv; doi.org/10.1101/2021.09.08.459475"][paper]

```
@article{,
   author = {Nault, Rance and Saha, Satabdi and Bhattacharya, Sudin and Dodson, Jack and Sinha, Samiran and Maiti, Tapabrata and Zacharewski, Tim},
   title = {Benchmarking of a Bayesian single cell RNAseq differential gene expression test for dose-response study designs},
   journal = {bioRxiv},
   pages = {2021.09.08.459475},
   DOI = {10.1101/2021.09.08.459475},
   url = {http://biorxiv.org/content/early/2021/09/10/2021.09.08.459475.abstract},
   year = {2021},
   type = {Journal Article}
}
```

[paper]: https://www.biorxiv.org/content/10.1101/2021.09.08.459475v1.full
