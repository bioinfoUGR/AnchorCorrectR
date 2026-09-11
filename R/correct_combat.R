# Anchor-aware ComBat batch correction
# Uses anchors to estimate batch parameters, then applies empirical Bayes shrinkage
# Standardizes data (ComBat-style) before parameter estimation
# Applies both location and scale correction
# Works on log scale; returns in input format.

#' Anchor-aware ComBat batch correction
#' 
#' Estimates batch location and scale parameters using anchor samples (replicates
#' across batches), then applies empirical Bayes shrinkage to stabilize estimates
#' across genes. This combines the reliability of anchor-based estimation with
#' ComBat's shrinkage approach.
#' 
#' The implementation follows ComBat's standardization approach: data are standardized
#' (mean=0, variance=1 per gene) before parameter estimation, and both location and
#' scale corrections are applied. This matches the standard ComBat methodology while
#' using anchor samples for more reliable batch effect estimation.
#' 
#' @param x matrix genes x samples (counts or log)
#' @param batch factor, length ncol(x)
#' @param sample_id factor of replicate IDs (IDs in >=2 batches act as anchors)
#' @param input_type "auto","counts","log"
#' @param ref_batch optional single batch label to use as uncorrected reference.
#'   When specified, this batch will have zero correction and other batches
#'   will be adjusted relative to it. If NULL, parameters are centered to preserve
#'   global means.
#' @param clip_counts logical; clip negatives to zero when returning counts
#' @param round_counts logical; round counts on return
#' @param verbose logical; print diagnostic messages
#' @return corrected matrix in the same format (counts or log) as the input
#' @export
#' @examples
#' \dontrun{
#' # Simulate counts with batch effects
#' set.seed(42)
#' counts <- matrix(rpois(5000 * 20, 100), nrow = 500)
#' colnames(counts) <- paste0("S", 1:20)
#' rownames(counts) <- paste0("Gene", 1:500)
#' 
#' # Add batch effect to second batch
#' batch <- factor(rep(c("A", "B"), each = 10))
#' counts[, batch == "B"] <- counts[, batch == "B"] * 1.5
#' 
#' # Create sample IDs (replicates across batches)
#' sample_id <- factor(rep(1:10, 2))
#' 
#' # Correct using anchor-aware ComBat
#' corrected <- correct_combat_anchor(counts, batch, sample_id)
#' }
correct_combat_anchor <- function(
  x, batch, sample_id,
  input_type = c("auto","counts","log"),
  ref_batch = NULL,
  clip_counts = TRUE,
  round_counts = FALSE,
  verbose = TRUE
) {
  input_type <- match.arg(input_type)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  batch <- droplevels(factor(batch))
  sample_id <- droplevels(factor(sample_id))
  B <- levels(batch)
  G <- nrow(x)
  N <- ncol(x)
  n_batches <- length(B)

  # Auto-detect and transform
  if (input_type == "auto") {
    input_type <- detect_input_type(x)
    if (verbose) message(sprintf("ComBat-anchor: auto-detected input type: %s", input_type))
  }
  
  tf <- make_transform_cpm_log1p(x, input_type = input_type)
  X <- tf$forward(x) # log-scale matrix

  # Validate reference batch
  ref_batch <- validate_ref_batch(ref_batch, batch, sample_id)

  # Identify anchor sample_ids observed in >= 2 batches
  anchor_info <- identify_anchors(batch, sample_id, min_batches = 2)
  anchor_ids <- anchor_info$anchor_ids
  
  if (verbose) {
    message(sprintf(
      "ComBat-anchor: using %d anchor sample IDs across %d observations",
      length(anchor_ids),
      length(anchor_info$anchor_indices)
    ))
  }

  # Step 0: Standardize data (ComBat-style standardization)
  # Compute grand mean per gene (overall mean across all samples)
  grand_mean <- rowMeans(X, na.rm = TRUE)
  
  # Compute pooled variance per gene (variance across all samples)
  # Handle missing values
  X_centered <- X - grand_mean
  var_pooled <- apply(X_centered, 1, function(r) {
    r_finite <- r[is.finite(r)]
    if (length(r_finite) < 2) return(1)
    var_r <- stats::var(r_finite, na.rm = TRUE)
    if (!is.finite(var_r) || var_r <= 0) return(1)
    var_r
  })
  
  # Standardize: (X - grand_mean) / sqrt(var_pooled)
  # Add small epsilon to avoid division by zero
  sqrt_var_pooled <- sqrt(pmax(var_pooled, 1e-8))
  stand_mean <- matrix(grand_mean, nrow = G, ncol = N)
  s_data <- (X - stand_mean) / sqrt_var_pooled
  
  # Replace any non-finite values with 0
  s_data[!is.finite(s_data)] <- 0
  
  if (verbose) {
    message("ComBat-anchor: standardized data (mean=0, var=1 per gene)")
  }

  # Step 1: Estimate batch location parameters (gamma) using anchors on STANDARDIZED data
  # For each gene and batch, compute location parameter from anchors
  gamma_gb <- matrix(0, nrow = G, ncol = n_batches, dimnames = list(rownames(X), B))
  n_anchors_gb <- matrix(0, nrow = G, ncol = n_batches, dimnames = list(rownames(X), B))
  
  for (b in B) {
    cols_b <- which(batch == b)
    if (!length(cols_b)) next
    
    # Which sample_ids are anchors and present in this batch?
    s_in_b <- intersect(anchor_ids, as.character(sample_id[cols_b]))
    if (!length(s_in_b)) {
      if (verbose) {
        warning(sprintf("Batch '%s' has no anchor samples - will use non-anchor samples for estimation", b))
        # Fall back to all samples in this batch
        s_in_b <- unique(as.character(sample_id[cols_b]))
      } else {
        next
      }
    }

    # For each anchor sample_id, compute its mean across all batches it appears in
    for (sid in s_in_b) {
      idxs_all <- which(sample_id == sid)
      if (length(idxs_all) == 0) next
      
      # Get this sample's value in batch b (using STANDARDIZED data)
      idxs_b <- idxs_all[batch[idxs_all] == b]
      if (length(idxs_b) == 0) next
      
      # Compute cross-batch mean on STANDARDIZED data
      if (length(idxs_all) > 1) {
        mu_sid_std <- rowMeans(s_data[, idxs_all, drop = FALSE], na.rm = TRUE)
      } else {
        mu_sid_std <- s_data[, idxs_all, drop = FALSE]
      }
      
      # Contribution to batch effect: difference from cross-batch mean (on standardized scale)
      contrib <- s_data[, idxs_b, drop = FALSE] - mu_sid_std
      gamma_gb[, b] <- gamma_gb[, b] + rowMeans(contrib, na.rm = TRUE)
      n_anchors_gb[, b] <- n_anchors_gb[, b] + 1
    }
  }
  
  # Average across anchors for each batch
  for (b in B) {
    mask <- n_anchors_gb[, b] > 0
    gamma_gb[mask, b] <- gamma_gb[mask, b] / n_anchors_gb[mask, b]
  }
  
  # Step 2: Estimate batch scale parameters (delta) using anchors
  # For each gene and batch, compute scale parameter
  delta_gb <- matrix(1, nrow = G, ncol = n_batches, dimnames = list(rownames(X), B))
  
  for (b in B) {
    cols_b <- which(batch == b)
    if (!length(cols_b)) next
    
    # Use anchors if available, otherwise all samples in batch
    s_in_b <- intersect(anchor_ids, as.character(sample_id[cols_b]))
    if (length(s_in_b) == 0) {
      s_in_b <- unique(as.character(sample_id[cols_b]))
    }
    
    # For each anchor, compute residual variance after removing location effect
    residuals_list <- list()
    for (sid in s_in_b) {
      idxs_all <- which(sample_id == sid)
      idxs_b <- idxs_all[batch[idxs_all] == b]
      if (length(idxs_b) == 0) next
      
      # Get cross-batch mean for this sample (on STANDARDIZED data)
      if (length(idxs_all) > 1) {
        mu_sid_std <- rowMeans(s_data[, idxs_all, drop = FALSE], na.rm = TRUE)
      } else {
        mu_sid_std <- s_data[, idxs_all, drop = FALSE]
      }
      
      # Residuals after removing location effect (on standardized scale)
      res <- s_data[, idxs_b, drop = FALSE] - mu_sid_std
      residuals_list[[length(residuals_list) + 1]] <- res
    }
    
    if (length(residuals_list) > 0) {
      # Combine residuals and compute variance
      all_res <- do.call(cbind, residuals_list)
      delta_gb[, b] <- apply(all_res, 1, function(r) {
        r_finite <- r[is.finite(r)]
        if (length(r_finite) < 2) return(1)
        var_r <- stats::var(r_finite, na.rm = TRUE)
        if (!is.finite(var_r) || var_r <= 0) return(1)
        sqrt(var_r)
      })
    } else {
      # Fallback: use variance of all samples in batch (on STANDARDIZED data)
      s_data_b <- s_data[, cols_b, drop = FALSE]
      delta_gb[, b] <- apply(s_data_b, 1, function(r) {
        r_finite <- r[is.finite(r)]
        if (length(r_finite) < 2) return(1)
        var_r <- stats::var(r_finite, na.rm = TRUE)
        if (!is.finite(var_r) || var_r <= 0) return(1)
        sqrt(var_r)
      })
    }
  }
  
  # Step 3: Apply empirical Bayes shrinkage to location parameters (gamma)
  if (verbose) message("ComBat-anchor: applying empirical Bayes shrinkage to location parameters...")
  gamma_star <- .ebayes_shrink_location(gamma_gb, n_anchors_gb, verbose)
  
  # Step 4: Apply empirical Bayes shrinkage to scale parameters (delta)
  if (verbose) message("ComBat-anchor: applying empirical Bayes shrinkage to scale parameters...")
  delta_star <- .ebayes_shrink_scale(delta_gb, n_anchors_gb, verbose)
  
  # Step 5: Recenter parameters based on reference batch or sum-to-zero
  if (!is.null(ref_batch)) {
    jref <- match(ref_batch, B)
    # Subtract reference batch location from all batches
    gamma_star <- sweep(gamma_star, 1, gamma_star[, jref], "-")
    gamma_star[, jref] <- 0
    # Scale parameters relative to reference (divide by reference)
    delta_star <- sweep(delta_star, 1, pmax(delta_star[, jref], 1e-6), "/")
    delta_star[, jref] <- 1
    if (verbose) {
      message(sprintf("ComBat-anchor: reference batch '%s' will remain uncorrected", ref_batch))
    }
  } else {
    # Sum-to-zero for location (preserve global mean per gene)
    gamma_mean <- rowMeans(gamma_star, na.rm = TRUE)
    gamma_star <- sweep(gamma_star, 1, gamma_mean, "-")
    # Scale parameters: normalize to geometric mean = 1
    delta_geom_mean <- exp(rowMeans(log(pmax(delta_star, 1e-6)), na.rm = TRUE))
    delta_star <- sweep(delta_star, 1, pmax(delta_geom_mean, 1e-6), "/")
    if (verbose) {
      message("ComBat-anchor: centering parameters to preserve global gene means and scales")
    }
  }
  
  # Step 6: Apply corrections gene by gene on STANDARDIZED data
  # ComBat correction: (s_data - gamma*) / sqrt(delta*)
  bayesdata <- s_data
  
  for (g in 1:G) {
    for (j in seq_along(B)) {
      cols <- which(batch == B[j])
      if (!length(cols)) next
      
      gamma_gj <- gamma_star[g, j]
      delta_gj <- delta_star[g, j]
      
      if (!is.finite(gamma_gj)) gamma_gj <- 0
      if (!is.finite(delta_gj) || delta_gj <= 0) delta_gj <- 1
      
      # ComBat correction on standardized data: (s_data - gamma*) / sqrt(delta*)
      # Note: delta_star is already on variance scale, so we use sqrt
      bayesdata[g, cols] <- (bayesdata[g, cols] - gamma_gj) / sqrt(pmax(delta_gj, 1e-8))
    }
  }
  
  # Step 7: Unstandardize: X_corrected = bayesdata * sqrt(var_pooled) + grand_mean
  # Broadcast sqrt_var_pooled and grand_mean across columns (samples)
  Xc <- sweep(bayesdata, 1, sqrt_var_pooled, "*")
  Xc <- sweep(Xc, 1, grand_mean, "+")
  
  # Handle reference batch: don't change reference batch data
  if (!is.null(ref_batch)) {
    jref <- match(ref_batch, B)
    cols_ref <- which(batch == B[jref])
    if (length(cols_ref) > 0) {
      Xc[, cols_ref] <- X[, cols_ref]
      if (verbose) {
        message(sprintf("ComBat-anchor: reference batch '%s' data unchanged", ref_batch))
      }
    }
  }
  
  if (tf$input_type == "counts") {
    out <- tf$inverse_counts_from_log(Xc)
    out <- .post_counts_guard(out, clip_counts, round_counts)
    return(out)
  } else {
    return(Xc)
  }
}

# ---- Empirical Bayes shrinkage functions ----

#' Empirical Bayes shrinkage for location parameters (gamma)
#' @keywords internal
.ebayes_shrink_location <- function(gamma_gb, n_anchors_gb, verbose = TRUE) {
  G <- nrow(gamma_gb)
  n_batches <- ncol(gamma_gb)
  
  # Estimate hyperparameters from data
  # Prior mean: assume zero (centered)
  # Prior variance: estimate from data
  
  # For each batch, estimate prior variance from genes
  gamma_star <- gamma_gb
  
  for (j in 1:n_batches) {
    gamma_j <- gamma_gb[, j]
    n_j <- n_anchors_gb[, j]
    
    # Only shrink genes with sufficient anchor support
    has_support <- n_j >= 2
    if (sum(has_support) < 10) {
      # Not enough data for shrinkage, use raw estimates
      if (verbose && j == 1) {
        message("  Insufficient anchor support for shrinkage, using raw estimates")
      }
      next
    }
    
    # Estimate prior variance from genes with support
    gamma_supported <- gamma_j[has_support]
    prior_var <- stats::var(gamma_supported, na.rm = TRUE)
    
    if (!is.finite(prior_var) || prior_var <= 0) {
      prior_var <- 1
    }
    
    # Estimate sampling variance for each gene (inverse of number of anchors)
    # More anchors = lower sampling variance = less shrinkage
    sampling_var <- 1 / pmax(n_j, 1)
    sampling_var[!has_support] <- Inf  # Don't shrink genes without support
    
    # Empirical Bayes estimate: weighted average of prior (0) and data
    # Shrinkage factor: prior_var / (prior_var + sampling_var)
    shrinkage <- prior_var / (prior_var + sampling_var)
    shrinkage[!is.finite(shrinkage)] <- 0
    
    # Shrink toward zero (centered)
    gamma_star[, j] <- gamma_j * (1 - shrinkage)
  }
  
  gamma_star
}

#' Empirical Bayes shrinkage for scale parameters (delta)
#' @keywords internal
.ebayes_shrink_scale <- function(delta_gb, n_anchors_gb, verbose = TRUE) {
  G <- nrow(delta_gb)
  n_batches <- ncol(delta_gb)
  
  # Work on log scale for scale parameters
  log_delta_gb <- log(pmax(delta_gb, 1e-6))
  
  delta_star <- delta_gb
  
  for (j in 1:n_batches) {
    log_delta_j <- log_delta_gb[, j]
    n_j <- n_anchors_gb[, j]
    
    # Only shrink genes with sufficient anchor support
    has_support <- n_j >= 2
    if (sum(has_support) < 10) {
      if (verbose && j == 1) {
        message("  Insufficient anchor support for scale shrinkage, using raw estimates")
      }
      next
    }
    
    # Estimate prior mean and variance on log scale
    log_delta_supported <- log_delta_j[has_support]
    prior_mean <- mean(log_delta_supported, na.rm = TRUE)
    prior_var <- stats::var(log_delta_supported, na.rm = TRUE)
    
    if (!is.finite(prior_mean)) prior_mean <- 0
    if (!is.finite(prior_var) || prior_var <= 0) prior_var <- 1
    
    # Sampling variance (inverse of number of anchors)
    sampling_var <- 1 / pmax(n_j, 1)
    sampling_var[!has_support] <- Inf
    
    # Empirical Bayes on log scale
    shrinkage <- prior_var / (prior_var + sampling_var)
    shrinkage[!is.finite(shrinkage)] <- 0
    
    # Shrink toward prior mean on log scale
    log_delta_star <- log_delta_j * (1 - shrinkage) + prior_mean * shrinkage
    
    # Transform back
    delta_star[, j] <- exp(log_delta_star)
    delta_star[!has_support, j] <- delta_gb[!has_support, j]  # Keep raw for unsupported
  }
  
  delta_star
}

