---
layout: page
title: Installation
description: ~
---

`mFABIO` is implemented as a R package, which can be installed from GitHub.

### Dependencies 
* R libraries: ggplot2, irlba, matrixStats, progress

#### 1. Install `mFABIO`
```r
install.packages("devtools")
devtools::install_github("superggbond/mFABIO")
library(mFABIO)
```
#### 2. Check the input options for the main function to perform TWAS fine-mapping on binary outcomes
```r
?run_mfabio
```
