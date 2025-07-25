---
layout: page
title: Tutorial
description: ~
---
This page provides a tutorial for multi-tissue TWAS fine-mapping using mFABIO. Before runing the example code, make sure that the mFABIO software is installed successfully. For instructions on installation, please see the [Installation section](https://superggbond.github.io/mFABIO/documentation/02_Installation.html).

### mFABIO
The example data for mFABIO tutorial can be accessed following this [page](https://superggbond.github.io/mFABIO/documentation/03_Data.html). Here are the details about the input data formats and how to run mFABIO. 

### 1.1 Formats of input files for mFABIO
* Predicted GReX: We require the predicted GReX of the TWAS cohort built up using standard softwares like [SuSiE](https://github.com/stephenslab/susieR). The input GReX is expected to be a R matrix, which has column names in the format of [Tissue Name]-[Gene Name/ID]. No requirement on the row names of this R matrix. An example data matrix can be found here:
```r
library(FABIO)
data(example_data)
> example_data$G[1:5,1:5]
     Tissue1-Gene001 Tissue2-Gene001 Tissue3-Gene001 Tissue4-Gene001 Tissue5-Gene001
[1,]      -0.7550567     -1.42914073     -0.88960517     -0.11741247      -0.4213615
[2,]       0.8011494      0.32952949     -0.14682934     -0.07053631       0.1506835
[3,]       0.5685453      0.33475624      0.01785558      0.17683326       0.1383405
[4,]       0.6705940      0.99009064     -0.47588513      0.33221315       0.5847007
[5,]      -0.4316600     -0.07432518      0.62807354     -1.01307066      -0.5182504
```

* Binary phenotypes: We also require the observed binary phenotypes of the TWAS cohort. The input phenotypes are expected as a R vector, with 1 coding for case and 0 for control. The order of the individuals here should be consistent with the order of rows in the predicted GReX file. An example file can be downloaded [here](https://www.dropbox.com/scl/fo/fxynm8uvedgvy7ni6hcbt/AAfTQVo89s78DsRNwpBH3lU?dl=0&e=2&preview=example_pheno.txt&rlkey=nbqwrdi2r5y1bbojzf7z8ev7h&st=yz28n4nj).

### 1.2 Preperation of the predicted GReX file
FABIO provides a function to help generate the predicted GReX file in the FABIO-requested format, with the inputs:
* PLINK 1 binary files (.bed+.bim+.fam) of genotypes sorted by chromosome. Please refer to this [webpage](https://www.cog-genomics.org/plink/2.0/input#bed) for the details of the format.
* Pre-trained eQTL weights first sorted by chromosome, then grouped by gene. The eQTL weights we used in the manuscript in the preferred file structures can be downloaded [here](https://www.dropbox.com/scl/fo/fxynm8uvedgvy7ni6hcbt/AAfTQVo89s78DsRNwpBH3lU?dl=0&e=1&preview=GEUVADIS_BSLMM_weights.zip&rlkey=nbqwrdi2r5y1bbojzf7z8ev7h&st=yz28n4nj). These weights were trained by BSLMM using GEUVADIS data.

This helper function can be applied like this:
```r
library(FABIO)

chr <- 22
geno_dir <- "/path/to/plink/files"
weight_dir <- './GEUVADIS_BSLMM_weights/chr22'
# The eQTL weight files can be customized, but should be organized in the same format and file structures as our example files.
prepGReX(chr, geno_dir, weight_dir)
# The results will be saved as a file named "grex_for_fabio.txt.gz", in the format illustrated in the section 1.1 above.
```

### 2. Running FABIO
The TWAS fine-mapping can be performed using the following scripts with our example data:
```r
library(FABIO)

grex <- data.table::fread('./example_grex.txt.gz')
pheno <- scan('./example_pheno.txt',numeric())
fabio(grex, pheno, beta_a=0, beta_b=0, w_step=100, s_step=1000)
# The output results will be saved as a file named "FABIO_out.csv"
# w_step and s_step are set to be small for a shorter running time as an example here
```
The inputs are:
- grex: the predicted GReX
- pheno: the TWAS phenotype vector
- beta_a: alpha of the prior beta distribution on pi, both beta_a and beta_b are 0 by default, leading to the default uniform prior on log pi
- beta_b: beta of the prior beta distribution on pi, both beta_a and beta_b are 0 by default, leading to the default uniform prior on log pi
- w-step: the number of warm-up steps in MCMC, default = 6000
- s-step: the number of sampling steps in MCMC, default = 20000

### 3. FABIO output
FABIO will output a summary table with three columns, and save it as a .csv file:

|Gene|PIP|estFDR|
|---|---|---|
|DLEU2L|1|0|
|DNTTIP2|1|0|
|FAM73A|1|0|
|FHL3|1|0|
|PSRC1|1|0|
|TCTEX1D4|1|0|
|UBE4B|1|0|
|GCLM|0.84|0.16|
|PSMA5|0.41|0.75|
|FBXO42|0.4|1|
|...|...|...|

- Gene: name of each input gene; the genes are ordered decreasingly by PIP
- PIP: corresponding PIP of each input gene
- estFDR: estimated FDR using that gene as cutoff
