# scBT
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
Install dependencies
```{r}
# brglm
install.packages('brglm')

# Seurat
install.packages('remotes')
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')

# limma
BiocManager::install("limma")
```

The developmental version of scBT can be installed from Github:
```{r}
library("devtools")
devtools::install_github("satabdisaha1288/scBT")
```

## Getting Started
Once installed the best place to get started is the vignette. The Quickstart vignette can be accessed as:

```
library(scBT)
DETest(sce, method = 'BAYES')
```

## Citing scBT
Please cite ["Nault, R., Saha, S., Bhattacharya, S., Dodson, J., Sinha, S., Maiti, T. and Zacharewski, T., 2022. Benchmarking of a Bayesian single cell RNAseq differential gene expression test for dose–response study designs. Nucleic acids research, 50(8), pp.e48-e48."][paper]

```
@article{nault2022benchmarking,
  title={Benchmarking of a Bayesian single cell RNAseq differential gene expression test for dose--response study designs},
  author={Nault, Rance and Saha, Satabdi and Bhattacharya, Sudin and Dodson, Jack and Sinha, Samiran and Maiti, Tapabrata and Zacharewski, Tim},
  journal={Nucleic acids research},
  volume={50},
  number={8},
  pages={e48--e48},
  year={2022},
  publisher={Oxford University Press}
}
```

[paper]: https://academic.oup.com/nar/article/50/8/e48/6513570
