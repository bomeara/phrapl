phrapl (phylogeography using approximated likelihoods) has moved from R-forge to github for development. Building the package on R-forge was failing; moving to github will let users install pre-CRAN versions using devtools. To do this,

```r
install.packages("devtools")
devtools::install_github("phrapl", "bomeara")
```
