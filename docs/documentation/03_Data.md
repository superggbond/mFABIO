---
layout: page
title: Data Input
description: ~
---
There are two mandatory and one optional inputs for mFABIO to perform multi-tissue TWAS fine-mapping on a binary trait of interest. Example data can be loaded within the mFABIO R package:
```r
    library(mFABIO)
    data(example_data)
 ```

#### 1. Predicted GReX of the TWAS cohort at gene-tissue pair level
  * An example R matrix can be loaded like this:
 ```r
    G <- example_data$G
 ```
  
#### 2. Binary phenotype of the TWAS cohort
  * An example R vector coding case as 1 and control as 0 can be loaded like this:
 ```r
    y <- example_data$y
 ```

#### 3. Corresponding genotypes of the cis-SNPs in the region of interest (optional)
  * This input is optional, and you can include it when considering the pleiotropic effects. An example R matrix can be loaded like this:
 ```r
    X <- example_data$X
 ```
