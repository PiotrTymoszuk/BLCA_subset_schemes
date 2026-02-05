# BLCA_subset_schemes

_Transcriptomic similarity and dissimilarity between molecular subset schemes in bladder cancer_

## Summary

This repository contains R code for a comparison of transcriptomes differentially regulated in our novel bladder cancer clusters (see: [repository with with clustering and phenotyping pipelines](https://github.com/PiotrTymoszuk/BLCA-cluster-paper) and [repository with interactive supplementary material](https://github.com/PiotrTymoszuk/BLCA_cluster_supplements/edit/main/README.md)), UROMOL classes of non-muscle invasive bladder cancer (NMIBC) by Hedegard et al. (DOI: 10.1016/j.ccell.2016.05.004), and consensus classes of muscle invasive bladder cancer (MIBC) by Kamoun et al. (DOI: 10.1016/J.EURURO.2019.09.006). 

## Terms of use 

The repository is open-source under a [GPL-3 license](https://github.com/PiotrTymoszuk/BLCA_subset_schemes/blob/main/LICENSE). 
Source data dn code requests should be addressed to [Prof. Renate Pichler](mailto:Renate.Pichler@i-med.ac.at) or [Prof. Zoran Culig](mailto:Zoran.Culig@i-med.ac.at). The repository maintainer is [Piotr Tymoszuk](mailto:piotr.s.tymoszuk@gmail.com).

## Usage

The repository's code used several packages available from GitHub, which may be installed with `devtools` R package: 

```r

devtools::install_github("PiotrTymoszuk/soucer") ## execution of R scripts in the pipeline
devtools::install_github("PiotrTymoszuk/trafo") ## extra transformation tools for tables and lists

devtools::install_github("PiotrTymoszuk/microViz") ## visualization of multi-omics data, differential  gene expression
devtools::install_github("PiotrTymoszuk/fastTest") ## statistical tests
devtools::install_github("PiotrTymoszuk/exDA") ## exploratory and plotting tools

devtools::install_github("PiotrTymoszuk/graphExtra") ## utilities of graph/igraph objects


devtools::install_github("PiotrTymoszuk/figur") ## utilities for figure and table objects

```

You can launch the entire pipeline by sourcing the `exec.R` file: 

```r

source("exec.R")

```
