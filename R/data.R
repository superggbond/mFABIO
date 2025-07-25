#' Example gene expression, phenotype, and genotype Data for mFABIO
#'
#' A simulated dataset containing gene expression values across multiple tissues,
#' corresponding genotypes, and binary phenotype data for demonstration of the mFABIO package.
#'
#' @format A list containing:
#' \describe{
#'   \item{G}{A numeric matrix of predicted gene expression values (3,000 individuals x 300 gene-tissue pairs).
#'             Column names are in format 'Tissue-GeneID'.}
#'   \item{y}{A numeric vector of binary outcomes (0 or 1) for 3,000 individuals.}
#'   \item{X}{A numeric matrix of genotypes (3,000 individuals x 1,000 SNPs).}
#' }
#'
#' @examples
#' data(example_data)
#'
#' # Check dimensions
#' dim(example_data$G)
#' table(example_data$y)
#'
#' # Run mfabio on example data
#' \dontrun{
#' results <- run_mfabio(
#'   G = example_data$G,
#'   y = example_data$y,
#'   X = example_data$X
#' )
#' }
"example_data"
