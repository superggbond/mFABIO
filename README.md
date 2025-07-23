# mFABIO

## multi-tissue Fine-mApping of causal genes for BInary Outcomes (mFABIO)
![mfabio](mFABIO_scheme.png)
mFABIO is a TWAS fine-mapping method that utilizes a probit model to directly analyze binary outcomes and jointly models gene-tissue pairs across multiple tissues within a given locus. This approach accounts for correlations in genetically regulated expression (GReX) among genes and tissues. Through a Bayesian variational inference algorithm, it obtains the posterior probability of having a non-zero effect for each gene-tissue pair, which is also known as the posterior inclusion probability (PIP). PIP serves as an important measure of evidence for the gene’s association with the outcome trait, and mFABIO generates PIPs at both the gene-tissue pair level and the gene level to nominate signal gene-tissue pairs and genes.

## How to use mFABIO

See [Tutorial](https://superggbond.github.io/mFABIO/) for detailed documentation and examples.

## Issues

For any questions or feedbacks on mFABIO software, please email to Haihan Zhang (hhzhang@umich.edu).

## Citation

Haihan Zhang, Kevin He, Lam C. Tsoi, and Xiang Zhou. "mFABIO: Multi-tissue TWAS fine-mapping to prioritize causal genes and tissues for binary traits."


