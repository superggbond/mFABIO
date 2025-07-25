% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_mfabio.R
\name{run_mfabio}
\alias{run_mfabio}
\title{Multi-tissue TWAS fine-mapping for binary outcomes (mFABIO)}
\usage{
run_mfabio(
  G,
  y,
  X = NULL,
  L = 10,
  K = 10,
  prior_var_mu = 100,
  prior_a_b = 1,
  prior_d_b = 1,
  prior_a_alpha = 1,
  prior_d_alpha = 1,
  base_prior_delta_gene = 1,
  base_prior_delta_tissue = 1,
  prior_update_dampening = 0.5,
  standardize = TRUE,
  max_iter = 100,
  tol = 1e-05,
  verbose = TRUE,
  use_adaptive_dampening = TRUE,
  min_dampening = 0.1,
  use_momentum = TRUE,
  momentum_rate = 0.9,
  check_elbo = TRUE,
  elbo_frequency = 5
)
}
\arguments{
\item{G}{A numeric matrix of predicted gene expression values (n individuals x m gene-tissue pairs).
Column names must be in the format 'Tissue-GeneID'.}

\item{y}{A numeric vector of binary outcomes (0 or 1) for n individuals.}

\item{X}{An optional numeric matrix of genotypes (n individuals x p SNPs). Defaults to \code{NULL}.}

\item{L}{An integer specifying the maximum number of causal effects to fit.}

\item{K}{An integer specifying the number of principal components to extract from \code{X}. Defaults to 10.}

\item{base_prior_delta_gene}{A numeric base hyperparameter for the Dirichlet prior on gene selection.}

\item{base_prior_delta_tissue}{A numeric base hyperparameter for the Dirichlet prior on tissue selection.}

\item{standardize}{Logical. If TRUE, G and PCs will be standardized. Defaults to TRUE.}

\item{max_iter}{An integer for the maximum number of VB iterations.}

\item{tol}{A numeric tolerance for the convergence criterion.}

\item{verbose}{Logical. If TRUE, prints progress messages.}

\item{use_adaptive_dampening}{Logical. If TRUE, uses adaptive dampening for better convergence.}

\item{min_dampening}{Minimum dampening factor for adaptive scheme.}
}
\value{
A list containing the model results.
}
\description{
This main function applies the mFABIO algorithm
}
