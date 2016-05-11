PHRAPL (Phylogeography using approximate likelihoods) is a method for phylogeographic model selection. This method can estimate the probability of observing a set of gene trees across a set of models that include population genetic parameters, e.g., population coalescence times and gene flow, in order to infer the model(s) in the set that most likely underlie the observed data.

Note that PHRAPL has recently moved from R-forge to Github for development. Moving to Github allows users to install pre-CRAN versions using devtools. To do this, within R, type

````
install.packages("devtools")
devtools::install_github("bomeara/phrapl")
````

