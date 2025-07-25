#' Multi-tissue TWAS fine-mapping for binary outcomes (mFABIO)
#'
#' This main function applies the mFABIO algorithm
#'
#' @param G A numeric matrix of predicted gene expression values (n individuals x m gene-tissue pairs).
#'   Column names must be in the format 'Tissue-GeneID'.
#' @param y A numeric vector of binary outcomes (0 or 1) for n individuals.
#' @param X An optional numeric matrix of genotypes (n individuals x p SNPs). Defaults to `NULL`.
#' @param L An integer specifying the maximum number of causal effects to fit.
#' @param K An integer specifying the number of principal components to extract from `X`. Defaults to 10.
#' @param base_prior_delta_gene A numeric base hyperparameter for the Dirichlet prior on gene selection.
#' @param base_prior_delta_tissue A numeric base hyperparameter for the Dirichlet prior on tissue selection.
#' @param standardize Logical. If TRUE, G and PCs will be standardized. Defaults to TRUE.
#' @param max_iter An integer for the maximum number of VB iterations.
#' @param tol A numeric tolerance for the convergence criterion.
#' @param verbose Logical. If TRUE, prints progress messages.
#' @param use_adaptive_dampening Logical. If TRUE, uses adaptive dampening for better convergence.
#' @param min_dampening Minimum dampening factor for adaptive scheme.
#'
#' @return A list containing the model results.
#' @importFrom graphics barplot hist par
#' @importFrom stats dnorm kmeans pnorm prcomp qnorm rnorm sd setNames var
#' @export
run_mfabio <- function(
    # --- Core Inputs ---
  G, y, X = NULL,
  # --- Model Parameters ---
  L = 10, K = 10,
  # --- Priors ---
  prior_var_mu = 100,
  prior_a_b = 1.0, prior_d_b = 1.0,
  prior_a_alpha = 1.0, prior_d_alpha = 1.0,
  base_prior_delta_gene = 1.0,
  base_prior_delta_tissue = 1.0,
  prior_update_dampening = 0.5,
  # --- Control Parameters ---
  standardize = TRUE,
  max_iter = 100,
  tol = 1e-5,
  verbose = TRUE,
  # --- New optimization parameters ---
  use_adaptive_dampening = TRUE,
  min_dampening = 0.1,
  use_momentum = TRUE,
  momentum_rate = 0.9,
  check_elbo = TRUE,
  elbo_frequency = 5
) {

  # --- Check for progress package ---
  use_progress_bar <- FALSE
  if (verbose && requireNamespace("progress", quietly = TRUE)) {
    use_progress_bar <- TRUE
  } else if (verbose) {
    warning("Package 'progress' not found. Progress bar not displayed.", call. = FALSE)
  }

  # --- 0. Setup & Input Checks ---
  if (verbose) cat("--- Initializing mFABIO ---\n")
  if (!is.matrix(G)) G <- as.matrix(G)
  n <- nrow(G); m <- ncol(G)

  # Validate inputs
  if (length(y) != n || !all(y %in% c(0, 1))) stop("Input 'y' is invalid.")
  if (m == 0 && L > 0) { L <- 0; warning("G has zero columns, setting L=0.", call. = FALSE) }
  if(prior_update_dampening <= 0 || prior_update_dampening > 1){
    warning("'prior_update_dampening' must be in (0, 1]. Setting to 1.0.", call.=FALSE)
    prior_update_dampening <- 1.0
  }

  # Parse column names
  column_info <- .parse_colnames(colnames(G))
  if(is.null(column_info) || nrow(column_info) != m) {
    stop("Colnames of G are missing or could not be parsed.", call.=FALSE)
  }
  names(column_info)[names(column_info) == "condition_name"] <- "tissue_name"

  # Create gene and tissue mappings
  gene_ids <- column_info$gene_id
  unique_gene_ids <- unique(gene_ids)
  Ng <- length(unique_gene_ids)

  # Optimized mapping creation
  map_col_to_gene_idx <- match(column_info$gene_id, unique_gene_ids)
  tissue_gene_maps <- vector("list", Ng)
  names(tissue_gene_maps) <- unique_gene_ids
  map_col_to_tissue_idx_in_gene <- integer(m)

  for (g_idx in 1:Ng) {
    gene_id_g <- unique_gene_ids[g_idx]
    cols_for_this_gene <- which(column_info$gene_id == gene_id_g)
    tissues_for_this_gene <- column_info$tissue_name[cols_for_this_gene]
    tissue_gene_maps[[g_idx]] <- list(
      original_indices = cols_for_this_gene,
      tissue_names = tissues_for_this_gene,
      n_tissue = length(cols_for_this_gene)
    )
    map_col_to_tissue_idx_in_gene[cols_for_this_gene] <- seq_along(cols_for_this_gene)
  }

  # Handle PC term
  include_pleio_term <- !is.null(X)
  if (include_pleio_term) {
    if (!is.matrix(X)) X <- as.matrix(X)
    if (nrow(X) != n) stop("Dim X/y mismatch.")
    p_input <- ncol(X)
    if (p_input < K && p_input > 0) { K <- p_input }
    else if (p_input == 0 || K <= 0) { include_pleio_term <- FALSE; X <- NULL; K <- 0 }
  } else { K <- 0 }

  if (verbose) {
    cat(sprintf("Model setup: n=%d, m=%d gene-tissue pairs (%d unique genes), L=%d, Include X effects=%s \n",
                n, m, Ng, L, include_pleio_term))
  }

  # --- 1. Preprocessing with improved numerical stability ---
  if (verbose) cat("Preprocessing data...\n")

  # Standardize G with numerical stability checks
  G_attr <- list(center = rep(0, m), scale = rep(1, m))
  if (standardize && m > 0) {
    G_means <- colMeans(G, na.rm = TRUE)

    # Use more stable SD calculation
    G_sds <- if(requireNamespace("matrixStats", quietly = TRUE)) {
      matrixStats::colSds(G, na.rm = TRUE)
    } else {
      apply(G, 2, function(x) sqrt(var(x, na.rm = TRUE)))
    }

    # Handle zero variance columns more carefully
    zero_var_cols <- which(G_sds < 1e-10)
    if(length(zero_var_cols) > 0) {
      warning(sprintf("Found %d columns with near-zero variance. These will be excluded from the model.",
                      length(zero_var_cols)), call. = FALSE)
      G_sds[zero_var_cols] <- 1
    }

    G <- sweep(G, 2, G_means, "-")
    G <- sweep(G, 2, G_sds, "/")
    G[is.nan(G) | is.infinite(G)] <- 0

    G_attr$center <- G_means
    G_attr$scale <- G_sds
  }

  # Pre-compute column squared sums for efficiency
  G_col_sq_sum <- colSums(G^2, na.rm = TRUE)

  # PCA handling with improved stability
  X_pc_used_in_model <- matrix(0.0, nrow = n, ncol = 0)
  Xpc_attr <- list()
  pca_ran_successfully <- FALSE
  input_was_pcs <- FALSE
  Xpc_col_sq_sum <- numeric(0)

  if (include_pleio_term && K > 0) {
    p_input <- ncol(X)
    if (p_input == K) {
      if (verbose) cat("Input X has K columns, assuming pre-computed PCs. Skipping PCA step.\n")
      X_pc <- X
      Xpc_attr <- list(note = "Input assumed PCs")
      input_was_pcs <- TRUE
      pca_ran_successfully <- TRUE
    } else if (p_input > K) {
      if (verbose) cat("Loading genotype matrix X...\n")
      tryCatch({
        pca_out <- .run_pca_wrapper(X, K)
        X_pc <- pca_out$x
        Xpc_attr <- list(
          center = pca_out$center,
          scale = pca_out$scale,
          rotation = pca_out$rotation
        )
        pca_ran_successfully <- TRUE
      }, error = function(e){
        warning("PCA calculation failed: ", e$message, ". Proceeding without X covariates.", call. = FALSE)
        include_pleio_term <<- FALSE
        K <<- 0
        X_pc <<- matrix(0.0, nrow = n, ncol = 0)
        Xpc_attr <<- list(error = "PCA failed")
        pca_ran_successfully <<- FALSE
      })
    }

    if(pca_ran_successfully && standardize) {
      pc_sds <- apply(X_pc, 2, sd, na.rm = TRUE)
      zero_var_idx <- which(pc_sds < 1e-10)
      if(length(zero_var_idx) > 0){
        warning(paste("PC(s)", paste(zero_var_idx, collapse=","), "have near-zero variance."), call. = FALSE)
        pc_sds[zero_var_idx] <- 1
      }
      if (K > 0) {
        X_pc <- sweep(X_pc, 2, pc_sds, "/")
        X_pc[is.nan(X_pc) | is.infinite(X_pc)] <- 0
        Xpc_attr$pc_scale <- pc_sds
      }
    }

    X_pc_used_in_model <- X_pc
    Xpc_col_sq_sum <- if(K > 0 && pca_ran_successfully) {
      colSums(X_pc_used_in_model^2, na.rm = TRUE)
    } else {
      numeric(0)
    }
  }

  # --- 2. Enhanced Initialization ---
  if (verbose) cat("Initializing VB parameters...\n")

  # Better initialization for intercept
  y_prop <- mean(y)
  mu_mu <- qnorm(max(0.01, min(0.99, y_prop)))
  s2_mu <- 1.0

  # Initialize PC effects
  if (include_pleio_term && K > 0) {
    # Initialize with small random values for better exploration
    mu_alpha <- rnorm(K, 0, 0.1)
    s2_alpha <- rep(1.0, K)
    E_alpha_pc <- mu_alpha
    E_alpha_pc_sq <- mu_alpha^2 + s2_alpha
    a_alpha_tilde <- prior_a_alpha + K / 2.0
    d_alpha_tilde <- prior_d_alpha + 0.5 * sum(E_alpha_pc_sq)
    E_inv_s2_alpha <- a_alpha_tilde / max(d_alpha_tilde, 1e-10)
  } else {
    E_alpha_pc <- numeric(0)
  }

  # Initialize latent variables with better starting values
  E_z <- qnorm(pmax(0.01, pmin(0.99, (y + 0.5) / 2)))

  # Initialize variational parameters
  if (m > 0 && L > 0) {
    # Use K-means initialization for alpha if L is small
    if (L <= 10 && m > L && L > 1) {
      tryCatch({
        set.seed(123)  # For reproducibility
        init_clusters <- kmeans(t(G), centers = L, nstart = 5, iter.max = 20)
        # Ensure we get a proper matrix
        if (!is.null(init_clusters$cluster)) {
          alpha_l <- matrix(0, nrow = L, ncol = m)
          for (j in 1:m) {
            alpha_l[init_clusters$cluster[j], j] <- 1
          }
          # Add small noise and normalize
          alpha_l <- alpha_l + 0.01
          alpha_l <- alpha_l / rowSums(alpha_l)
        } else {
          alpha_l <- matrix(1/m, nrow = L, ncol = m)
        }
      }, error = function(e) {
        # Fallback to uniform initialization
        alpha_l <<- matrix(1/m, nrow = L, ncol = m)
      })
    } else {
      alpha_l <- matrix(1/m, nrow = L, ncol = m)
    }
  } else if (L == 0 || m == 0) {
    alpha_l <- matrix(0, nrow = max(L, 0), ncol = max(m, 0))
  } else {
    alpha_l <- matrix(1/m, nrow = L, ncol = m)
  }

  # Initialize effect sizes
  mu_b <- rnorm(L, 0, 0.1)
  s2_b <- rep(1.0, L)
  E_beta <- numeric(m)
  E_b <- mu_b
  E_b2 <- mu_b^2 + s2_b

  # Ensure dimensions are correct
  if (L > 0 && m > 0) {
    if (length(mu_b) != L) stop("Dimension mismatch in mu_b initialization")
    if (nrow(alpha_l) != L || ncol(alpha_l) != m) {
      stop(sprintf("Dimension mismatch in alpha_l: expected %d x %d, got %d x %d",
                   L, m, nrow(alpha_l), ncol(alpha_l)))
    }
  }

  # Initialize hyperparameters
  a_b_tilde <- prior_a_b + L / 2.0
  d_b_tilde <- prior_d_b + 0.5 * sum(E_b2)
  E_inv_s2_b <- a_b_tilde / max(d_b_tilde, 1e-10)

  # Initialize Dirichlet parameters
  tilde_delta_G <- if(Ng > 0) rep(base_prior_delta_gene, Ng) else numeric(0)
  tilde_delta_Tg <- vector("list", Ng)
  if(Ng > 0) {
    names(tilde_delta_Tg) <- unique_gene_ids
    for (g_idx in 1:Ng) {
      n_tissue_g <- tissue_gene_maps[[g_idx]]$n_tissue %||% 0
      tilde_delta_Tg[[g_idx]] <- if(n_tissue_g > 0) {
        rep(base_prior_delta_tissue, n_tissue_g)
      } else {
        numeric(0)
      }
    }
  }

  # Initialize tracking variables
  param_proxy_history <- rep(NA_real_, max_iter)
  elbo_history <- rep(NA_real_, max_iter)

  # Momentum variables
  if (use_momentum) {
    momentum_E_beta <- numeric(m)
    momentum_alpha <- matrix(0, nrow = L, ncol = m)
  }

  # Adaptive dampening
  if (use_adaptive_dampening) {
    current_dampening <- prior_update_dampening
    dampening_history <- rep(NA_real_, max_iter)
  }

  # Setup Progress Bar
  if (use_progress_bar) {
    pb <- progress::progress_bar$new(
      format = " Running mFABIO [:bar] :percent | Iter: :current/:total | Change: :change",
      total = max_iter,
      width = 80,
      clear = FALSE
    )
    pb$tick(0, tokens = list(change = "---"))
  } else if (verbose) {
    cat("Starting VB iterations...\n")
  }

  # --- 3. Optimized VB Iterations ---
  iter <- 0
  converged <- FALSE
  divergence_count <- 0
  max_divergence <- 5

  for (iter in 1:max_iter) {
    E_beta_old <- E_beta
    alpha_l_old <- alpha_l

    # --- Update intercept ---
    Xpc_E_alpha_term <- if(include_pleio_term && K > 0) {
      X_pc_used_in_model %*% E_alpha_pc
    } else {
      0
    }

    G_E_beta_term <- if(m > 0) G %*% E_beta else 0

    residual_mu <- E_z - G_E_beta_term - Xpc_E_alpha_term
    s2_mu_inv <- n + (1 / prior_var_mu)
    s2_mu <- 1 / s2_mu_inv
    mu_mu <- s2_mu * sum(residual_mu)
    E_mu <- mu_mu

    # --- Update PC effects ---
    if (include_pleio_term && K > 0) {
      residual_alpha <- E_z - E_mu - G_E_beta_term
      s2_alpha_inv <- Xpc_col_sq_sum + E_inv_s2_alpha
      s2_alpha <- 1 / (s2_alpha_inv + 1e-10)
      mu_alpha <- s2_alpha * as.numeric(crossprod(X_pc_used_in_model, residual_alpha))
      E_alpha_pc <- mu_alpha
      E_alpha_pc_sq <- mu_alpha^2 + s2_alpha
      d_alpha_tilde <- prior_d_alpha + 0.5 * sum(E_alpha_pc_sq)
      E_inv_s2_alpha <- a_alpha_tilde / max(d_alpha_tilde, 1e-10)
      Xpc_E_alpha_term <- X_pc_used_in_model %*% E_alpha_pc
    } else {
      Xpc_E_alpha_term <- 0
    }

    # --- Compute expected log priors with numerical stability ---
    if(Ng > 0 && m > 0 && length(tilde_delta_G) == Ng){
      # Stabilized digamma calculations
      sum_tilde_delta_G <- sum(tilde_delta_G)
      psi_sum_tilde_delta_G <- digamma(max(sum_tilde_delta_G, 1e-10))
      E_log_pi_G <- digamma(pmax(tilde_delta_G, 1e-10)) - psi_sum_tilde_delta_G
      E_log_pi_G_full <- E_log_pi_G[map_col_to_gene_idx]

      E_log_pi_Cg_full <- numeric(m)
      for(g_idx in 1:Ng){
        delta_g <- tilde_delta_Tg[[g_idx]]
        if(length(delta_g) > 0){
          sum_delta_g <- sum(delta_g)
          psi_sum_delta_g <- digamma(max(sum_delta_g, 1e-10))
          E_log_pi_Cg_for_g <- digamma(pmax(delta_g, 1e-10)) - psi_sum_delta_g
          cols_for_this_gene <- tissue_gene_maps[[g_idx]]$original_indices
          E_log_pi_Cg_full[cols_for_this_gene] <- E_log_pi_Cg_for_g
        }
      }

      # Ensure numerical stability
      E_log_pi_G_full[!is.finite(E_log_pi_G_full)] <- -50
      E_log_pi_Cg_full[!is.finite(E_log_pi_Cg_full)] <- -50
    } else {
      E_log_pi_G_full <- numeric(m)
      E_log_pi_Cg_full <- numeric(m)
    }

    # --- Update variational parameters for effects (optimized) ---
    E_fit_minus_beta <- E_mu + Xpc_E_alpha_term
    E_b2_sum <- 0

    if (m > 0 && L > 0) {
      # Pre-compute E_beta from current parameters if first iteration
      if (iter == 1) {
        E_beta <- colSums(alpha_l * mu_b, na.rm = TRUE)
      }

      # Pre-compute for efficiency
      E_beta_l_terms <- matrix(0, nrow = L, ncol = m)
      for (l in 1:L) {
        E_beta_l_terms[l, ] <- alpha_l[l, ] * mu_b[l]
      }

      # Store old value for comparison
      E_beta_old_inner <- E_beta

      for (l in 1:L) {
        # Efficient residual calculation
        E_beta_minus_l <- E_beta - E_beta_l_terms[l, ]
        G_beta_minus_l <- G %*% E_beta_minus_l
        R_l <- E_z - E_fit_minus_beta - G_beta_minus_l

        # Single-effect regression
        SER_l <- as.numeric(crossprod(G, R_l))

        # Compute log posterior weights with improved stability
        Likelihood_term_lj <- SER_l * E_b[l] - 0.5 * E_b2[l] * G_col_sq_sum
        log_prior_term_lj <- E_log_pi_G_full + E_log_pi_Cg_full

        log_tilde_alpha_lj <- Likelihood_term_lj + log_prior_term_lj

        # Improved log-sum-exp trick
        max_log_alpha <- max(log_tilde_alpha_lj[is.finite(log_tilde_alpha_lj)], na.rm = TRUE)
        if(!is.finite(max_log_alpha)) max_log_alpha <- 0

        # Compute weights with underflow protection
        log_w_l <- log_tilde_alpha_lj - max_log_alpha
        log_w_l[!is.finite(log_w_l)] <- -50

        # Stable normalization
        log_sum_w_l <- log(sum(exp(pmin(log_w_l, 50))))
        log_alpha_l <- log_w_l - log_sum_w_l

        # Apply momentum if enabled
        if (use_momentum && iter > 1) {
          new_alpha_l <- exp(log_alpha_l)
          alpha_l[l, ] <- momentum_rate * momentum_alpha[l, ] + (1 - momentum_rate) * new_alpha_l
          momentum_alpha[l, ] <- alpha_l[l, ]
        } else {
          alpha_l[l, ] <- exp(log_alpha_l)
        }

        # Ensure normalization
        alpha_l[l, ] <- alpha_l[l, ] / sum(alpha_l[l, ])

        # Update effect size parameters
        s2_blj_inv <- G_col_sq_sum + E_inv_s2_b
        s2_blj <- 1 / (s2_blj_inv + 1e-10)
        mu_blj <- s2_blj * SER_l

        mu_b[l] <- sum(alpha_l[l, ] * mu_blj, na.rm = TRUE)
        s2_b[l] <- sum(alpha_l[l, ] * (mu_blj^2 + s2_blj), na.rm = TRUE) - mu_b[l]^2
        s2_b[l] <- max(s2_b[l], 1e-12)

        E_b[l] <- mu_b[l]
        E_b2[l] <- mu_b[l]^2 + s2_b[l]

        # Update E_beta_l_terms for next iteration
        E_beta_l_terms[l, ] <- alpha_l[l, ] * mu_b[l]
      }

      # Update combined effects
      E_beta_new <- colSums(E_beta_l_terms, na.rm = TRUE)

      # Apply momentum to E_beta if enabled
      if (use_momentum && iter > 1) {
        E_beta <- momentum_rate * momentum_E_beta + (1 - momentum_rate) * E_beta_new
        momentum_E_beta <- E_beta
      } else {
        E_beta <- E_beta_new
      }

      # Update variance hyperparameter
      d_b_tilde <- prior_d_b + 0.5 * sum(E_b2, na.rm = TRUE)
      E_inv_s2_b <- a_b_tilde / max(d_b_tilde, 1e-10)
    } else {
      E_beta <- numeric(0)
    }

    # --- Update Dirichlet parameters with adaptive dampening ---
    if (m > 0 && L > 0 && Ng > 0 && !anyNA(alpha_l)){
      tilde_delta_G_old <- tilde_delta_G
      tilde_delta_Tg_old <- tilde_delta_Tg

      # Efficient aggregation
      alpha_sums_j <- colSums(alpha_l, na.rm = TRUE)
      alpha_sums_j[!is.finite(alpha_sums_j)] <- 0

      # Update gene-level Dirichlet
      tilde_delta_G_new <- rep(base_prior_delta_gene, Ng)
      for (g_idx in 1:Ng) {
        cols_g <- tissue_gene_maps[[g_idx]]$original_indices
        tilde_delta_G_new[g_idx] <- base_prior_delta_gene + sum(alpha_sums_j[cols_g])
      }

      # Apply dampening (adaptive if enabled)
      if (use_adaptive_dampening) {
        # Compute change in parameters
        param_change_delta <- mean(abs(tilde_delta_G_new - tilde_delta_G_old))

        # Adjust dampening based on convergence behavior
        if (iter > 5) {
          if (param_change_delta < tol * 10) {
            current_dampening <- max(min_dampening, current_dampening * 0.95)
          } else if (param_change_delta > 0.1) {
            current_dampening <- min(prior_update_dampening, current_dampening * 1.05)
          }
        }
        dampening_history[iter] <- current_dampening
      } else {
        current_dampening <- prior_update_dampening
      }

      tilde_delta_G <- (1 - current_dampening) * tilde_delta_G_old + current_dampening * tilde_delta_G_new
      tilde_delta_G[tilde_delta_G <= 0 | !is.finite(tilde_delta_G)] <- 1e-10

      # Update tissue-level Dirichlet
      for(g_idx in 1:Ng){
        cols_for_this_gene <- tissue_gene_maps[[g_idx]]$original_indices
        if(length(cols_for_this_gene) > 0){
          alpha_sums_subset <- alpha_sums_j[cols_for_this_gene]
          alpha_sums_subset[!is.finite(alpha_sums_subset)] <- 0
          tilde_delta_Tg_new_g <- base_prior_delta_tissue + alpha_sums_subset

          tilde_delta_Tg[[g_idx]] <- (1 - current_dampening) * tilde_delta_Tg_old[[g_idx]] +
            current_dampening * tilde_delta_Tg_new_g
          tilde_delta_Tg[[g_idx]][tilde_delta_Tg[[g_idx]] <= 0 | !is.finite(tilde_delta_Tg[[g_idx]])] <- 1e-10
        } else {
          tilde_delta_Tg[[g_idx]] <- numeric(0)
        }
      }
    }

    # --- Update latent variables with truncated normal ---
    E_beta_final <- E_beta
    E_alpha_pc_final <- if(include_pleio_term && K > 0) E_alpha_pc else numeric(0)

    E_linpred <- E_mu +
      (if(m > 0) G %*% E_beta_final else 0) +
      (if(include_pleio_term && K > 0) X_pc_used_in_model %*% E_alpha_pc_final else 0)

    E_z <- .get_truncated_normal_mean(
      mean = E_linpred,
      sd = 1,
      lower = ifelse(y == 1, 0, -Inf),
      upper = ifelse(y == 1, Inf, 0)
    )

    # --- Compute ELBO if requested ---
    if (check_elbo && iter %% elbo_frequency == 0) {
      # Simplified ELBO computation for monitoring
      log_lik <- sum(dnorm(E_z, mean = E_linpred, sd = 1, log = TRUE))
      kl_terms <- 0

      # KL for effects
      if (L > 0) {
        kl_terms <- kl_terms + 0.5 * sum(log(s2_b) - s2_b - mu_b^2 + 1)
      }

      # KL for PC effects
      if (include_pleio_term && K > 0) {
        kl_terms <- kl_terms + 0.5 * sum(log(s2_alpha) - s2_alpha - mu_alpha^2 + 1)
      }

      elbo_history[iter] <- log_lik - kl_terms

      # Check for divergence
      if (iter > elbo_frequency && !is.na(elbo_history[iter - elbo_frequency])) {
        if (elbo_history[iter] < elbo_history[iter - elbo_frequency] - 1e-3) {
          divergence_count <- divergence_count + 1
          if (divergence_count > max_divergence) {
            warning("ELBO decreasing consistently. Stopping early.", call. = FALSE)
            break
          }
        } else {
          divergence_count <- 0
        }
      }
    }

    # --- Check convergence ---
    param_change <- if (m > 0) {
      mean(abs(E_beta - E_beta_old), na.rm = TRUE)
    } else {
      0
    }

    # Also check alpha convergence
    alpha_change <- mean(abs(alpha_l - alpha_l_old), na.rm = TRUE)
    total_change <- max(param_change, alpha_change)

    if(!is.finite(total_change)) total_change <- Inf
    param_proxy_history[iter] <- total_change

    # Update progress bar
    if (use_progress_bar) {
      pb$tick(tokens = list(change = sprintf("%.2e", total_change)))
    }

    # Check for numerical issues
    if (!is.finite(total_change)){
      converged <- FALSE
      if(use_progress_bar && !pb$finished) pb$update(1.0)
      warning("Numerical instability detected. Stopping.", call. = FALSE)
      break
    }

    # Check convergence with stricter criteria
    if (iter > 5) {
      # Check both parameter change and relative change
      rel_change <- total_change / (mean(abs(E_beta)) + 1e-10)
      if (total_change < tol && rel_change < tol * 10) {
        converged <- TRUE
        if(use_progress_bar && !pb$finished) pb$update(1.0)
        break
      }
    }
  }

  if (!converged && iter == max_iter) {
    if (verbose && !use_progress_bar) {
      cat(sprintf("Warning: Max iterations (%d) reached. Final change: %.2e\n", max_iter, total_change))
    }
    if(use_progress_bar && !pb$finished) pb$update(1.0)
  } else if (converged && verbose) {
    cat(sprintf("Converged after %d iterations. Final change: %.2e\n", iter, total_change))
  }

  # --- 4. Compute final quantities (focused on PIPs) ---
  if (verbose) cat("\n");cat("Computing final posterior quantities...\n")

  # === PRIMARY OUTPUT 1: Gene-tissue pair PIPs ===
  pip_pair <- rep(0, m)
  if (m > 0 && L > 0 && !is.null(alpha_l) && all(is.finite(alpha_l))) {
    pip_pair <- 1 - apply(1 - alpha_l, 2, prod, na.rm = TRUE)
    names(pip_pair) <- column_info$col_name
  }

  # === PRIMARY OUTPUT 2: Gene-level PIPs ===
  pip_gene <- rep(0, Ng)
  if(Ng > 0 && !is.null(tilde_delta_G) && length(tilde_delta_G) == Ng && all(is.finite(tilde_delta_G))){
    sum_delta_G <- sum(tilde_delta_G, na.rm = TRUE)
    if(sum_delta_G > 1e-10){
      pip_gene <- tilde_delta_G / sum_delta_G
      names(pip_gene) <- unique_gene_ids
    }
  }

  # === Tissue-given-gene PIPs (conditional probabilities) ===
  pip_tissue_given_gene <- vector("list", Ng)
  if(Ng > 0) names(pip_tissue_given_gene) <- unique_gene_ids

  if(Ng > 0 && !is.null(tilde_delta_Tg) && length(tilde_delta_Tg) == Ng){
    for(g_idx in 1:Ng){
      delta_g <- tilde_delta_Tg[[g_idx]]
      tissue_names_g <- tissue_gene_maps[[g_idx]]$tissue_names
      n_tissue_g <- length(tissue_names_g)

      pip_tg <- rep(0, n_tissue_g)
      if(n_tissue_g > 0 && !is.null(delta_g) && length(delta_g) == n_tissue_g && all(is.finite(delta_g))){
        sum_delta_g <- sum(delta_g, na.rm = TRUE)
        if(sum_delta_g > 1e-10){
          pip_tg <- delta_g / sum_delta_g
          names(pip_tg) <- tissue_names_g
        }
      }
      pip_tissue_given_gene[[g_idx]] <- pip_tg
    }
  }

  # === Gene summary with focus on PIPs ===
  gene_summary <- NULL
  if (Ng > 0 && !is.null(column_info)) {
    tryCatch({
      # Aggregate PIPs by gene
      pip_by_gene <- tapply(pip_pair, column_info$gene_id, sum, na.rm = TRUE)
      max_pip_by_gene <- tapply(pip_pair, column_info$gene_id, max, na.rm = TRUE)
      n_tissues_by_gene <- tapply(column_info$gene_id, column_info$gene_id, length)

      gene_summary <- data.frame(
        gene_id = unique_gene_ids,
        gene_pip = pip_gene[unique_gene_ids],
        sum_pair_pip = as.numeric(pip_by_gene[unique_gene_ids]),
        max_pair_pip = as.numeric(max_pip_by_gene[unique_gene_ids]),
        n_tissues = as.numeric(n_tissues_by_gene[unique_gene_ids]),
        avg_pair_pip = as.numeric(pip_by_gene[unique_gene_ids]) / as.numeric(n_tissues_by_gene[unique_gene_ids]),
        stringsAsFactors = FALSE
      )

      # Sort by gene PIP (descending)
      gene_summary <- gene_summary[order(-gene_summary$gene_pip), ]
      rownames(gene_summary) <- NULL

    }, error = function(e){
      warning("Gene summary computation failed: ", e$message, call. = FALSE)
      gene_summary <<- NULL
    })
  }

  # Compute posterior means of variance components
  post_mean_sigma2_b <- if(a_b_tilde > 1 && is.finite(d_b_tilde)) {
    d_b_tilde / (a_b_tilde - 1)
  } else {
    NA
  }

  post_mean_sigma2_alpha <- if(include_pleio_term && K > 0 && a_alpha_tilde > 1 && is.finite(d_alpha_tilde)) {
    d_alpha_tilde / (a_alpha_tilde - 1)
  } else {
    NA
  }

  # Prepare final results
  results <- list(
    # === PRIMARY OUTPUTS ===
    pip_gene = pip_gene,                      # Gene-level PIPs
    pip_pair = pip_pair,                      # Gene-tissue pair PIPs
    pip_tissue_given_gene = pip_tissue_given_gene,  # Conditional tissue PIPs
    gene_summary = gene_summary,              # Summary table for genes

    # === MODEL PARAMETERS ===
    model_params = list(
      alpha = alpha_l,                       # Allocation probabilities
      mu_b = mu_b,                           # Effect means
      s2_b = s2_b,                           # Effect variances
      E_sigma2_b = post_mean_sigma2_b,       # Posterior mean of variance
      intercept = list(mu = mu_mu, s2 = s2_mu),
      posterior_delta_G = setNames(tilde_delta_G, unique_gene_ids),
      posterior_delta_Tg = setNames(tilde_delta_Tg, unique_gene_ids)
    ),

    # === CONVERGENCE DIAGNOSTICS ===
    convergence = list(
      converged = converged,
      n_iter = iter,
      param_change_history = param_proxy_history[1:iter],
      elbo_history = if(check_elbo) elbo_history[!is.na(elbo_history)] else NULL,
      dampening_history = if(use_adaptive_dampening) dampening_history[1:iter] else NULL,
      final_change = total_change
    ),

    # === MODEL SETTINGS ===
    settings = list(
      L = L,
      K = K,
      n = n,
      m = m,
      n_genes = Ng,
      standardize = standardize,
      use_adaptive_dampening = use_adaptive_dampening,
      use_momentum = use_momentum,
      check_elbo = check_elbo
    )
  )

  # Add pleiotropic effects info if applicable
  if (include_pleio_term && K > 0) {
    results$pleio_effects <- list(
      mu = mu_alpha,
      s2 = s2_alpha,
      E_sigma2_alpha = post_mean_sigma2_alpha
    )

    # Add PCA information
    pca_info_list <- list()
    pca_info_list$X_pc <- X_pc_used_in_model

    if(pca_ran_successfully && !input_was_pcs){
      pca_info_list$rotation = Xpc_attr$rotation
      pca_info_list$center = Xpc_attr$center
      pca_info_list$scale = Xpc_attr$scale
      pca_info_list$pc_scale = Xpc_attr$pc_scale
      pca_info_list$type = "Calculated from X"
    } else if (input_was_pcs){
      pca_info_list$pc_scale = Xpc_attr$pc_scale
      pca_info_list$type = "Input assumed PCs"
    } else {
      pca_info_list$error = "PCA not run or failed"
      pca_info_list$type = "None"
    }

    results$pca_info <- pca_info_list
  }

  # Add standardization info
  if (standardize) {
    results$standardization <- list()
    if (m > 0 && !is.null(G_attr$center)) {
      results$standardization$G_center = G_attr$center
      results$standardization$G_scale = G_attr$scale
    }
    if(include_pleio_term && K > 0 && !is.null(Xpc_attr$pc_scale)){
      results$standardization$Xpc_scale = Xpc_attr$pc_scale
    }
  }

  # Add diagnostic plots function
  results$plot_diagnostics <- function() {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      warning("Package 'ggplot2' needed for diagnostics plots. Please install it.", call. = FALSE)
      return(NULL)
    }

    # Create diagnostic plots
    par(mfrow = c(2, 2))

    # 1. Convergence plot
    plot(results$convergence$param_change_history, type = "l",
         xlab = "Iteration", ylab = "Parameter Change",
         main = "Convergence History", log = "y")

    # 2. PIP distribution
    hist(results$pip_pair, breaks = 30,
         xlab = "PIP", main = "Distribution of PIPs",
         col = "lightblue")

    # 3. ELBO if available
    if (!is.null(results$convergence$elbo_history)) {
      plot(results$convergence$elbo_history, type = "l",
           xlab = "Iteration", ylab = "ELBO",
           main = "ELBO History")
    }

    # 4. Top gene-tissue pairs by PIP
    if (!is.null(results$pip_pair) && length(results$pip_pair) > 0) {
      # Get top 10 pairs
      n_pairs_to_show <- min(10, length(results$pip_pair))
      top_pairs_idx <- order(results$pip_pair, decreasing = TRUE)[1:n_pairs_to_show]
      top_pairs_pip <- results$pip_pair[top_pairs_idx]
      top_pairs_names <- names(results$pip_pair)[top_pairs_idx]

      # Shorten names if too long
      short_names <- sapply(top_pairs_names, function(x) {
        if (nchar(x) > 20) {
          paste0(substr(x, 1, 17), "...")
        } else {
          x
        }
      })

      # Create barplot
      par(mar = c(8, 4, 2, 2))  # Increase bottom margin for labels
      barplot(top_pairs_pip,
              names.arg = short_names,
              las = 2,  # Vertical labels
              main = "Top 10 Gene-Tissue Pairs by PIP",
              ylab = "PIP",
              col = "skyblue",
              cex.names = 0.8)
    }

    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))  # Reset to defaults
  }

  if (verbose) cat("--- mFABIO analysis complete ---\n")

  return(results)
}
