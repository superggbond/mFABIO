---
layout: full
homepage: true
disable_anchors: true
description: multi-tissue Fine-mApping of causal genes for BInary Outcomes
---
## mFABIO
![mfabio\_pipeline](mFABIO_scheme.png)
FABIO is a TWAS fine-mapping method that relies on a probit model to directly relate multiple genetically regulated gene expression (GReX) to binary outcome in TWAS fine-mapping. Additionally, it jointly models all genes located on the same chromosome to account for the correlation among GReX arising from cis-SNP LD and expression correlation across genomic regions. Through a Markov chain Monte Carlo (MCMC) algorithm, it obtains the posterior probability of having a non-zero effect for each gene, which is also known as the posterior inclusion probability (PIP). PIP serves as an important measure of evidence for the gene’s association with the outcome trait, and FABIO nominates signal genes based PIP. FABIO is implemented as a R package freely available at [https://github.com/superggbond/mFABIO/](https://github.com/superggbond/mFABIO/).

### Example Analysis with mFABIO: [here](https://superggbond.github.io/mFABIO/documentation/04_mFABIO_Example.html).
