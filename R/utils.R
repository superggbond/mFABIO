#' Calculate Expected Value of a Truncated Normal Distribution
#'
#' An efficient, vectorized function to compute the mean of a Normal
#' distribution that is truncated to a given interval [lower, upper]. It
#' includes robust handling of edge cases like zero variance, infinite bounds,
#' and truncation intervals with near-zero probability mass.
#'
#' @param mean A numeric vector of means for the underlying Normal distributions.
#' @param sd A numeric vector of standard deviations. Must be positive.
#' @param lower A numeric vector of lower truncation bounds.
#' @param upper A numeric vector of upper truncation bounds.
#'
#' @return A numeric vector of the same length as `mean`, containing the
#'   expected values of the truncated normal distributions.
#'
#' @keywords internal
#' @noRd
.get_truncated_normal_mean <- function(mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  
  # Vectorize inputs and handle edge cases where sd <= 0
  if (any(sd <= 0)) {
    # Initialize output with the mean
    mean_trunc <- mean
    # For invalid sd, the result is the mean, unless it's outside the bounds,
    # in which case it's clamped to the nearest bound.
    idx_invalid_sd <- which(sd <= 0)
    mean_trunc[idx_invalid_sd] <- pmax(lower[idx_invalid_sd], pmin(upper[idx_invalid_sd], mean[idx_invalid_sd]))
    
    # If all sd values are invalid, return early
    if (length(idx_invalid_sd) == length(sd)) return(mean_trunc)
    
    # If only some are invalid, proceed with the valid ones
    idx_valid_sd <- which(sd > 0)
    mean <- mean[idx_valid_sd]
    sd <- sd[idx_valid_sd]
    lower <- lower[idx_valid_sd]
    upper <- upper[idx_valid_sd]
  } else {
    # If all are valid, create the output vector
    mean_trunc <- numeric(length(mean))
    idx_valid_sd <- seq_along(mean)
  }
  
  # Standardize bounds for valid sd cases
  alpha <- (lower - mean) / sd
  beta <- (upper - mean) / sd
  
  # Calculate the normalization constant Z = Phi(beta) - Phi(alpha)
  Z <- pnorm(beta) - pnorm(alpha)
  
  # Calculate numerator of the adjustment term: pdf(alpha) - pdf(beta)
  pdf_alpha <- dnorm(alpha)
  pdf_beta <- dnorm(beta)
  
  # The main formula for the mean
  mean_adj <- mean + sd * (pdf_alpha - pdf_beta) / Z
  mean_trunc[idx_valid_sd] <- mean_adj
  
  # --- Handle Numerical Stability & Special Cases ---
  
  # Case 1: Z is effectively zero (truncation far in the tails)
  # Fallback to the nearest bound as an approximation.
  idx_zero_z <- idx_valid_sd[Z < 1e-16]
  if (length(idx_zero_z) > 0) {
    mean_at_zero_z <- mean[Z < 1e-16]
    lower_at_zero_z <- lower[Z < 1e-16]
    upper_at_zero_z <- upper[Z < 1e-16]
    closer_bound <- ifelse(abs(lower_at_zero_z - mean_at_zero_z) < abs(upper_at_zero_z - mean_at_zero_z),
                           lower_at_zero_z, upper_at_zero_z)
    mean_trunc[idx_zero_z] <- closer_bound
  }
  
  # Case 2: lower and upper bounds are identical (point mass)
  idx_point_mass <- idx_valid_sd[lower == upper]
  if (length(idx_point_mass) > 0) {
    mean_trunc[idx_point_mass] <- lower[lower == upper]
  }
  
  # Final safety check for any remaining NaNs
  if (any(is.nan(mean_trunc))) {
    nan_idx <- which(is.nan(mean_trunc))
    # Fallback to the original mean for any case that slipped through
    mean_trunc[nan_idx] <- mean[nan_idx]
  }
  
  return(mean_trunc)
}

#' Efficient Principal Component Analysis Wrapper
#'
#' This function computes the first K principal components of a matrix. It acts as a
#' wrapper that preferentially uses the fast, partial PCA implementation from the
#' 'irlba' package. If 'irlba' is not installed, it falls back to the standard
#' base R `prcomp` function.
#'
#' @param X A numeric matrix (n x p) on which to perform PCA.
#' @param K The number of principal components to compute.
#' @param center Logical or a numeric vector of length p for centering.
#' @param scale. Logical or a numeric vector of length p for scaling.
#'
#' @return A list with a structure similar to `prcomp`, containing:
#'   \item{x}{An n x K matrix of principal components (scores).}
#'   \item{rotation}{A p x K matrix of variable loadings (eigenvectors).}
#'   \item{sdev}{A vector of standard deviations of the principal components.}
#'   \item{center}{The centering vector used.}
#'   \item{scale}{The scaling vector used.}
#'
#' @keywords internal
#' @noRd
.run_pca_wrapper <- function(X, K, center = TRUE, scale. = TRUE) {
  
  # Preferentially use the fast irlba package if available
  if (requireNamespace("irlba", quietly = TRUE)) {
    # --- Input Validation for K ---
    if (K > ncol(X)) {
      warning(paste("Requested K=", K, "PCs, but X only has", ncol(X), "columns. Using K=", ncol(X), "."), call. = FALSE)
      K <- ncol(X)
    }
    # Handle edge case where K=0
    if (K == 0) {
      return(list(x = matrix(0.0, nrow = nrow(X), ncol = 0),
                  rotation = matrix(0.0, nrow = ncol(X), ncol = 0),
                  sdev = numeric(0), center = NA, scale = NA))
    }
    
    # --- Pre-computation for Centering and Scaling ---
    # Ensure X is in matrix format
    if (!is.matrix(X)) X <- as.matrix(X)
    
    center_vals <- if (isTRUE(center)) colMeans(X, na.rm = TRUE) else if (is.numeric(center)) center else rep(0, ncol(X))
    scale_vals <- if (isTRUE(scale.)) {
      # Use matrixStats::colSds for efficiency if available
      sds <- if (requireNamespace("matrixStats", quietly = TRUE)) {
        matrixStats::colSds(X)
      } else {
        apply(X, 2, sd)
      }
      sds[sds == 0] <- 1 # Avoid division by zero
      sds
    } else if (is.numeric(scale.)) {
      scale.
    } else {
      rep(1, ncol(X))
    }
    
    # --- Run IRLBA ---
    # suppressWarnings to hide potential non-fatal convergence warnings
    pca_res <- suppressWarnings(
      irlba::irlba(X, nv = K, center = center_vals, scale. = scale_vals)
    )
    
    # Reconstruct PC scores
    pcs <- pca_res$u %*% diag(pca_res$d[1:K], nrow = K, ncol = K)
    
    return(list(
      x = pcs,
      rotation = pca_res$v,
      sdev = pca_res$d[1:K] / sqrt(max(1, nrow(X) - 1)),
      center = center_vals,
      scale = scale_vals
    ))
    
  } else {
    # --- Fallback to base R prcomp ---
    warning("Package 'irlba' not found. Falling back to slower 'prcomp'. For large datasets, install 'irlba'.", call. = FALSE)
    # Adjust K if necessary for prcomp
    K <- min(K, nrow(X), ncol(X))
    if (K == 0) {
      return(list(x = matrix(0, nrow = nrow(X), ncol = 0), rotation = matrix(0, nrow = ncol(X), ncol = 0), sdev = numeric(0), center=NA, scale=NA))
    }
    
    pca_res <- prcomp(X, rank. = K, center = center, scale. = scale.)
    
    return(list(
      x = pca_res$x[, 1:K, drop = FALSE],
      rotation = pca_res$rotation[, 1:K, drop = FALSE],
      sdev = pca_res$sdev[1:K],
      center = pca_res$center,
      scale = pca_res$scale
    ))
  }
}

#' Default Value for NULL
#'
#' This infix operator returns the right-hand side value if the left-hand side
#' is NULL. It's a convenient way to provide default values.
#'
#' @param a The object to check for NULL.
#' @param b The default value to return if `a` is NULL.
#'
#' @return `a` if `a` is not NULL, otherwise `b`.
#'
#' @keywords internal
#' @noRd
#' @name op-null-default
#'
#' @examples
#' \dontrun{
#'   x <- NULL
#'   y <- 5
#'   x %||% y  # Returns 5
#'
#'   z <- 10
#'   z %||% y  # Returns 10
#' }
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

#' Parse 'Condition-GeneID' Column Names
#'
#' An efficient, vectorized function to split column names into condition and
#' gene ID parts. It assumes the last hyphen in the string is the separator,
#' which allows condition names themselves to contain hyphens.
#'
#' @param cnames A character vector of column names (e.g., "Heart-Left-Ventricle-GENE123").
#'
#' @return A data frame with columns: `col_index`, `col_name`, `condition_name`, `gene_id`.
#'   If a name cannot be parsed, the full name is used as the condition and a
#'   placeholder "unknown_gene_..." is assigned.
#'
#' @keywords internal
#' @noRd
.parse_colnames <- function(cnames) {
  if (is.null(cnames)) {
    warning("Column names are NULL. Cannot parse.", call. = FALSE)
    return(NULL)
  }
  # Handle empty input
  if (length(cnames) == 0) {
    return(data.frame(col_index = integer(), col_name = character(),
                      condition_name = character(), gene_id = character(),
                      stringsAsFactors = FALSE))
  }
  
  # --- Vectorized Parsing using Base R `sub` ---
  
  # 1. Extract Gene ID: The greedy '.*-' matches everything up to the LAST hyphen
  # and replaces it with nothing, leaving only the gene ID.
  gene_id <- sub(".*-", "", cnames)
  
  # 2. Extract Condition Name: The '-[^-]*$' matches the LAST hyphen and any
  # subsequent characters that are NOT hyphens, and replaces it with nothing.
  condition_name <- sub("-[^-]*$", "", cnames)
  
  
  # --- Handle Cases with No Hyphen ---
  # If a name has no hyphen, the substitutions above will fail, and
  # gene_id and condition_name will both equal the original column name.
  no_hyphen_idx <- (gene_id == cnames)
  if (any(no_hyphen_idx)) {
    # The condition_name is already correct (it's the full name).
    # We just need to assign a placeholder gene_id for these cases.
    # Using seq_along ensures unique placeholders even if names are duplicated.
    gene_id[no_hyphen_idx] <- paste0("unknown_gene_", seq_along(cnames)[no_hyphen_idx])
  }
  
  # --- Create the final data frame ---
  parsed_df <- data.frame(
    col_index = seq_along(cnames),
    col_name = cnames,
    condition_name = trimws(condition_name),
    gene_id = trimws(gene_id),
    stringsAsFactors = FALSE
  )
  
  # --- Final checks and warnings (same as your original, good practice) ---
  num_unknown <- sum(grepl("^unknown_gene_", parsed_df$gene_id))
  if (num_unknown > 0 && num_unknown == length(cnames)) {
    warning("Could not parse GeneIDs from any column names based on hyphen. Using full column names as GeneIDs for summaries.", call. = FALSE)
    # Use full name as ID if all failed, so gene summary is still possible
    parsed_df$gene_id <- parsed_df$col_name
  } else if (num_unknown > 0) {
    warning(paste("Could not parse GeneID for", num_unknown, "columns using hyphen. Placeholders assigned."), call. = FALSE)
  }
  
  return(parsed_df)
}
