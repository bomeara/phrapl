<!-- badges: start -->
  [![Travis build status](https://travis-ci.com/bomeara/phrapl.svg?branch=master)](https://travis-ci.com/bomeara/phrapl) [![Codecov test coverage](https://codecov.io/gh/bomeara/phrapl/branch/master/graph/badge.svg)](https://codecov.io/gh/bomeara/phrapl?branch=master)
  <!-- badges: end -->

<img src="https://cloud.githubusercontent.com/assets/7992349/15196013/661bb536-1798-11e6-964f-100ae6aa1b7a.jpg" alt "phraplLogo" width= "600", height="200">

PHRAPL (Phylogeography using approximate likelihoods) is a method for phylogeographic model selection. This method can estimate the probability of observing a set of gene trees across a set of models that include population genetic parameters, e.g., population coalescence times and gene flow, in order to infer the model(s) in the set that most likely gave rise to the observed data.

Note that PHRAPL has recently moved from R-forge to Github for development. Moving to Github allows users to install pre-CRAN versions using devtools. To do this, within R, type

```
install.packages("devtools")
devtools::install_github("bomeara/phrapl")
```
