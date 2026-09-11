# Ridge implementation using glmnet on log data.
# Estimates per-gene batch effects via shared lambda (median from CV subset).
# Standardizes data (mean=0, variance=1 per gene) before fitting for consistent
# lambda selection and regularization across genes.

# Catch glmnet failure when y is constant (e.g. zero-variance genes after standardization)
.ridge_stop_if_constant_response <- function(e) {
  msg <- conditionMessage(e)
  if (grepl("y is constant", msg)) {
    stop(
      "Ridge: glmnet failed because the response is constant for at least one gene. ",
      "Remove genes that are all zero or have very low variance across samples, then retry.",
      call. = FALSE
    )
  }
  stop(e)
}

#' Fit ridge regression to estimate batch effects from anchor samples
#' 
#' Uses glmnet with alpha=0 (ridge penalty) to estimate per-gene batch effects.
#' A common lambda is chosen via cross-validation on a random subset of genes,
#' then applied to all genes for consistency.
#' 
#' Data are standardized (mean=0, variance=1 per gene) before fitting to ensure
#' consistent lambda selection and balanced regularization across genes with
#' different scales. Standardization parameters are stored and used during
#' correction application.
#' 
#' @param Xlog numeric matrix (genes x samples) in log scale
#' @param batch factor of batch labels
#' @param sample_id factor of sample identifiers (anchors span multiple batches)
#' @param lambda optional pre-specified lambda; if NULL, estimated via CV
#' @param seed random seed for gene sampling in CV
#' @param n_genes_cv number of genes to use for lambda CV
#' @param verbose logical; print diagnostic messages
#' @return object of class "anchor_ridge_fit" with components:
#'   - beta: matrix (genes x batches) of batch effect estimates (on standardized scale)
#'   - batches: batch levels
#'   - lambda: regularization parameter used
#'   - gene_mean: vector of gene means (for unstandardization)
#'   - gene_sd: vector of gene standard deviations (for unstandardization)
#' @export
#' @examples
#' \dontrun{
#' # Simulate log-scale data
#' set.seed(123)
#' Xlog <- matrix(rnorm(1000 * 30, mean = 5, sd = 2), nrow = 1000)
#' batch <- factor(rep(c("A", "B", "C"), each = 10))
#' sample_id <- factor(rep(1:10, 3))
#' 
#' # Fit ridge model
#' fit <- fit_anchor_ridge(Xlog, batch, sample_id)
#' print(dim(fit$beta))
#' }
fit_anchor_ridge <- function(Xlog, batch, sample_id, lambda = NULL, 
                             seed = 12345, n_genes_cv = 1000, verbose = TRUE) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet package is required for ridge method. Install it with: install.packages('glmnet')")
  }
  
  Xlog <- as.matrix(Xlog)
  storage.mode(Xlog) <- "double"
  batch <- droplevels(factor(batch))
  
  # One-hot encode batches
  Xb <- stats::model.matrix(~ 0 + batch)
  colnames(Xb) <- sub("^batch", "", colnames(Xb))
  # Center columns so effects sum to zero across batches
  Xb <- scale(Xb, center = TRUE, scale = FALSE)

  set.seed(seed)
  G <- nrow(Xlog)
  
  # Standardize data: mean=0, variance=1 per gene
  # Compute gene means and standard deviations
  gene_mean <- rowMeans(Xlog, na.rm = TRUE)
  gene_sd <- apply(Xlog, 1, function(r) {
    r_finite <- r[is.finite(r)]
    if (length(r_finite) < 2) return(1)
    sd_r <- sd(r_finite, na.rm = TRUE)
    if (!is.finite(sd_r) || sd_r <= 0) return(1)
    sd_r
  })
  
  # Standardize: (X - mean) / sd
  Xlog_std <- sweep(Xlog, 1, gene_mean, "-")
  Xlog_std <- sweep(Xlog_std, 1, pmax(gene_sd, 1e-8), "/")
  
  # Replace any non-finite values with 0
  Xlog_std[!is.finite(Xlog_std)] <- 0
  
  if (verbose) {
    message("Ridge: standardized data (mean=0, variance=1 per gene)")
  }
  
  # Choose lambda if not given via CV on gene subset (using STANDARDIZED data)
  if (is.null(lambda)) {
    if (verbose) message("Ridge: estimating lambda via cross-validation...")
    n_cv <- min(G, n_genes_cv)
    idx <- sample.int(G, n_cv)
    lam_vals <- numeric(length(idx))
    for (k in seq_along(idx)) {
      y <- Xlog_std[idx[k], ]
      fitcv <- tryCatch(
        glmnet::cv.glmnet(x = Xb, y = y, alpha = 0, intercept = TRUE),
        error = .ridge_stop_if_constant_response
      )
      lam_vals[k] <- fitcv$lambda.min
    }
    lambda <- median(lam_vals, na.rm = TRUE)
    if (verbose) message(sprintf("Ridge: selected lambda = %.4f (median from %d genes)", lambda, n_cv))
  } else {
    if (verbose) message(sprintf("Ridge: using pre-specified lambda = %.4f", lambda))
  }

  # Fit per-gene ridge at chosen lambda on STANDARDIZED data
  if (verbose) message("Ridge: fitting per-gene models on standardized data...")
  beta <- matrix(0, nrow = G, ncol = ncol(Xb))
  for (g in 1:G) {
    y <- Xlog_std[g, ]
    fitg <- tryCatch(
      glmnet::glmnet(x = Xb, y = y, alpha = 0, lambda = lambda, intercept = TRUE),
      error = .ridge_stop_if_constant_response
    )
    coefs <- as.matrix(stats::coef(fitg, s = lambda))
    # Drop intercept; keep batch columns (match colnames in Xb)
    nm <- rownames(coefs)
    keep <- match(colnames(Xb), nm)
    beta[g, ] <- as.numeric(coefs[keep, , drop = TRUE])
  }
  colnames(beta) <- colnames(Xb)
  rownames(beta) <- rownames(Xlog)

  structure(list(
    beta = beta,
    batches = colnames(Xb),
    lambda = lambda,
    gene_mean = gene_mean,
    gene_sd = gene_sd
  ), class = "anchor_ridge_fit")
}

#' Apply ridge correction to data
#' 
#' Subtracts estimated batch effects from samples. Data are standardized before
#' applying corrections (using stored standardization parameters), then unstandardized
#' after correction. If ref_batch is specified, effects are recentered so the
#' reference batch has zero correction. Otherwise, global means are preserved.
#' 
#' @param fit object of class "anchor_ridge_fit" from fit_anchor_ridge
#' @param Xlog numeric matrix (genes x samples) in log scale
#' @param batch factor of batch labels
#' @param ref_batch optional reference batch label to keep uncorrected
#' @return corrected matrix in log scale
#' @export
#' @examples
#' \dontrun{
#' # After fitting
#' Xc <- apply_correction_ridge(fit, Xlog, batch)
#' 
#' # With reference batch
#' Xc_ref <- apply_correction_ridge(fit, Xlog, batch, ref_batch = "A")
#' }
apply_correction_ridge <- function(fit, Xlog, batch, ref_batch = NULL) {
  if (!inherits(fit, "anchor_ridge_fit")) {
    stop("fit must be an object of class 'anchor_ridge_fit' from fit_anchor_ridge()")
  }
  
  Xlog <- as.matrix(Xlog)
  storage.mode(Xlog) <- "double"
  batch <- factor(batch, levels = fit$batches)
  if (any(is.na(batch))) {
    unknown_batches <- unique(as.character(batch[is.na(batch)]))
    stop(sprintf(
      "Unknown batches in data: %s. Batches must match those used in fit_anchor_ridge: %s",
      paste(unknown_batches, collapse = ", "),
      paste(fit$batches, collapse = ", ")
    ))
  }

  if (is.null(fit$gene_mean) || is.null(fit$gene_sd)) {
    stop(
      "Fit object is missing gene_mean/gene_sd. Refit with fit_anchor_ridge().",
      call. = FALSE
    )
  }

  # Standardize data using stored parameters
  gene_mean <- fit$gene_mean
  gene_sd <- fit$gene_sd
  
  # Ensure gene_mean and gene_sd match Xlog dimensions
  if (length(gene_mean) != nrow(Xlog)) {
    stop(sprintf(
      "Mismatch: fit was trained on %d genes, but Xlog has %d genes",
      length(gene_mean), nrow(Xlog)
    ))
  }
  
  # Standardize: (X - mean) / sd
  Xlog_std <- sweep(Xlog, 1, gene_mean, "-")
  Xlog_std <- sweep(Xlog_std, 1, pmax(gene_sd, 1e-8), "/")
  Xlog_std[!is.finite(Xlog_std)] <- 0

  beta <- fit$beta
  B <- fit$batches

  # Handle reference batch: recenter effects relative to ref_batch
  if (!is.null(ref_batch)) {
    ref_batch <- as.character(ref_batch)
    if (!ref_batch %in% B) {
      stop(sprintf(
        "ref_batch '%s' not found in fitted batches: %s",
        ref_batch,
        paste(B, collapse = ", ")
      ))
    }
    jref <- match(ref_batch, B)
    # Subtract reference batch effect from all batches
    beta <- sweep(beta, 1, beta[, jref], "-")
    beta[, jref] <- 0
  }

  # Apply corrections on STANDARDIZED data: subtract batch effects
  Xc_std <- Xlog_std
  for (j in seq_along(B)) {
    cols <- which(batch == B[j])
    if (!length(cols)) next
    off <- beta[, j]
    off[!is.finite(off)] <- 0
    Xc_std[, cols] <- sweep(Xc_std[, cols, drop = FALSE], 1, off, "-")
  }

  # Unstandardize: X_corrected = Xc_std * sd + mean
  Xc <- sweep(Xc_std, 1, pmax(gene_sd, 1e-8), "*")
  Xc <- sweep(Xc, 1, gene_mean, "+")
  
  # Handle reference batch: restore original data for reference batch
  if (!is.null(ref_batch)) {
    jref <- match(ref_batch, B)
    cols_ref <- which(batch == B[jref])
    if (length(cols_ref) > 0) {
      Xc[, cols_ref] <- Xlog[, cols_ref]
    }
  } else {
    # Preserve global means when no reference is specified
    mu_before <- rowMeans(Xlog, na.rm = TRUE)
    mu_after  <- rowMeans(Xc,   na.rm = TRUE)
    Xc <- sweep(Xc, 1, (mu_after - mu_before), "-")
  }
  
  Xc
}

#' Print method for anchor_ridge_fit
#'
#' @param x Object to print.
#' @param ... Additional arguments (unused).
#' @export
print.anchor_ridge_fit <- function(x, ...) {
  cat("== anchorCorrectR: Ridge regression fit ==\n")
  cat(sprintf("Genes: %d\n", nrow(x$beta)))
  cat(sprintf("Batches: %s\n", paste(x$batches, collapse = ", ")))
  cat(sprintf("Lambda: %.4f\n", x$lambda))
  cat("\nBatch effect summary (mean absolute coefficient):\n")
  mean_abs_coef <- colMeans(abs(x$beta), na.rm = TRUE)
  print(mean_abs_coef)
  invisible(x)
}
