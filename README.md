# snSeq-DGE
This package implements a new Bayesian test for detecting differential gene expression over multiple dose groups in single cell gene 
expression studies. 

## Simulating dose-response data
In order to develop new testing methods for differential expression analysis we need data for which the ground truth is known.
Here we simulate dose-response data based on known dose-response models and using parameters derived from real data.

![simulation](/Analysis/SimulationWorkflow.png)

_dose-response models derived from [BMD Express](https://bmdexpress-2.readthedocs.io/en/feature-readthedocs/) and the [US EPA BMDS Documentation](https://www.epa.gov/bmds/benchmark-dose-software-bmds-version-27-user-manual)_

* Hill
* Exponential 2 - 5
* Power
* Linear
* Polynomial 2-4

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
The development version of scBT can be installed from Github:

```

```

## Getting Started
Once installed the best place to get started is the vignette. The Quickstart vignette can be accessed as:

```
library(scBT)
browseVignettes("sCBT")
```

## Citing scBT
If you use scBT please cite our paper
