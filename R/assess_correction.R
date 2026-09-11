#' Assess batch-correction performance
#'
#' A collection of metrics to quantify (1) how much batch variability was removed
#' and (2) how much biological variability and structure were retained, comparing
#' a matrix before and after correction.
#'
#' Matrices can be counts or log-scale. Internally, counts are transformed to
#' log1p-CPM using \code{make_transform_cpm_log1p()} and results are computed
#' on principal components (PCs) of the \emph{samples} unless stated otherwise.
#'
#' @param x_before Matrix: genes x samples, before correction (counts or log).
#' @param x_after  Matrix: genes x samples, after correction (same scale as x_before).
#' @param batch    Factor/character, length ncol(x): batch labels.
#' @param biology  Optional factor/character/numeric, length ncol(x): biological labels or values.
#' @param sample_id Optional factor/character, length ncol(x): sample IDs for replicate RMSE calculation.
#' @param exclude_replicates Logical; controls whether replicated samples are excluded from calculations.
#'   If TRUE (default) and sample_id is provided, excludes ALL samples that have replicates 
#'   from all metrics except RMSE. If FALSE, includes ALL samples (including replicates) in 
#'   all metric calculations. When sample_id is NULL, this parameter has no effect.
#' @param input_type "auto","counts","log". Applied to both matrices.
#' @param n_pcs    Number of PCs to use for PCA-based metrics (default 30).
#' @param k        Neighborhood size for kNN-based metrics (default 20).
#' @param top_n_hvg Number of top variable genes to consider for HVG overlap (default 2000).
#' @param verbose  Logical; message progress.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{knn_jaccard}: average kNN Jaccard before-vs-after.
#'   \item \code{hvg_overlap}: Jaccard of top variable genes before/after.
#'   \item \code{biology_concordance}: consistency metric for biological signal (varies by biology type).
#'   \item \code{biology_type}: type of biological variable detected.
#'   \item \code{n_samples_used}: number of unique samples used for metrics.
#'   \item \code{n_samples_excluded}: number of replicated samples excluded.
#'   \item \code{pvca}: list with PVCA results (batch_variance, biology_variance, residual_variance) before/after.
#'   \item \code{adjusted_rand_index}: ARI comparing clustering structure before vs after correction.
#'   \item \code{ari_n_clusters}: number of clusters used for ARI computation.
#'   \item \code{cluster_compactness_separation}: list with within/between cluster distances before/after.
#'   \item \code{rmse_replicates}: list with \code{before} and \code{after} scalar summaries (see Details).
#' }
#'
#' @details
#' Replicate RMSE (\code{rmse_replicates}): For each \code{sample_id} that appears in at least two
#' batches, the function builds per-batch mean expression (log1p-CPM) across samples for that ID,
#' then for each unordered pair of batches computes \code{sqrt(mean_g (a_g - b_g)^2)} over genes.
#' Pairwise RMSEs are averaged within each anchor, then those anchor-level values are averaged
#' across anchors so each anchor counts equally. Values are on the log-expression scale and are
#' not bounded to a small interval-substantial batch effects can yield RMSE of 1 or larger.
#'
#' @export
assess_correction <- function(
    x_before, x_after, batch, biology = NULL, sample_id = NULL,
    exclude_replicates = TRUE,
    input_type = c("auto","counts","log"),
    n_pcs = 30, k = 20, top_n_hvg = 2000, verbose = TRUE
) {
  input_type <- match.arg(input_type)
  stopifnot(ncol(x_before) == ncol(x_after))
  stopifnot(length(batch) == ncol(x_before))
  if (!is.null(biology)) stopifnot(length(biology) == ncol(x_before))
  if (!is.null(sample_id)) stopifnot(length(sample_id) == ncol(x_before))
  
  # Transform to log scale consistently
  tr <- make_transform_cpm_log1p(x_before, input_type = input_type)
  Xb_log <- tr$forward(x_before)
  Xa_log <- tr$forward(x_after)
  
  batch <- as.factor(batch)
  
  # Detect biology type and handle appropriately
  biology_type <- NULL
  if (!is.null(biology)) {
    if (is.numeric(biology)) {
      biology_type <- "continuous"
      if (verbose) message("Detected continuous biological variable")
    } else {
      biology <- as.factor(biology)
      n_levels <- nlevels(biology)
      if (n_levels == 2) {
        biology_type <- "binary"
        if (verbose) message("Detected binary biological variable")
      } else if (n_levels > 2) {
        biology_type <- "multi_category"
        if (verbose) message(sprintf("Detected multi-category biological variable (%d levels)", n_levels))
      }
    }
  }
  
  sample_id <- if (is.null(sample_id)) NULL else as.factor(sample_id)
  
  # Determine which samples to use for most metrics (excluding ALL replicates if requested)
  if (exclude_replicates && !is.null(sample_id)) {
    if (verbose) message("Excluding all replicated samples: keeping only unique sample_ids...")
    
    # Count occurrences of each sample_id
    sample_counts <- table(sample_id)
    
    # Keep only samples where sample_id appears exactly once (non-replicates)
    unique_sample_ids <- names(sample_counts)[sample_counts == 1]
    keep_idx <- sample_id %in% unique_sample_ids
    
    n_excluded <- sum(!keep_idx)
    n_unique_samples <- length(unique_sample_ids)
    n_replicated_samples <- length(sample_counts) - n_unique_samples
    
    if (verbose) {
      message(sprintf("  Found %d unique sample_ids and %d replicated sample_ids", 
                      n_unique_samples, n_replicated_samples))
      message(sprintf("  Excluded %d samples (all samples with replicates)", n_excluded))
      message(sprintf("  Using %d unique samples for metrics", sum(keep_idx)))
    }
  } else {
    # Include all samples (including replicates) in calculations
    keep_idx <- rep(TRUE, ncol(Xb_log))
    n_excluded <- 0
    if (verbose && !is.null(sample_id)) {
      message("Including all samples (including replicates) in metric calculations")
    }
  }
  
  # Subset data for most metrics
  Xb_log_filtered <- Xb_log[, keep_idx, drop = FALSE]
  Xa_log_filtered <- Xa_log[, keep_idx, drop = FALSE]
  batch_filtered <- batch[keep_idx]
  biology_filtered <- if (is.null(biology)) NULL else {
    if (is.numeric(biology)) biology[keep_idx] else biology[keep_idx]
  }
  
  # PCA on samples (using filtered data)
  if (verbose) message("Running PCA...")
  pcs_before <- .pca_scores(Xb_log_filtered, n_pcs = n_pcs)
  pcs_after  <- .pca_scores(Xa_log_filtered, n_pcs = n_pcs)
  
  # === 1) kNN graph conservation (before vs after) on PCs
  if (verbose) message("Computing kNN Jaccard...")
  knn_jacc <- .knn_jaccard(pcs_before, pcs_after, k = k)
  
  # === 2) HVG overlap
  if (verbose) message("Computing HVG overlap...")
  hvg_jacc <- .hvg_overlap(Xb_log_filtered, Xa_log_filtered, top_n = top_n_hvg)
  
  # === 3) Biology concordance (varies by type)
  if (verbose) message("Computing biology concordance (if applicable)...")
  biology_concordance <- list()
  
  if (!is.null(biology_filtered)) {
    if (biology_type == "binary") {
      # Original two-group logFC correlation
      de_conc <- .two_group_logfc_cor(Xb_log_filtered, Xa_log_filtered, biology_filtered)
      biology_concordance <- list(
        metric = "logFC_correlation",
        value = de_conc,
        interpretation = "Higher correlation indicates better preservation of differential expression"
      )
    } else if (biology_type == "multi_category") {
      # ANOVA F-statistic correlation
      f_conc <- .multi_group_fstat_cor(Xb_log_filtered, Xa_log_filtered, biology_filtered)
      eta_conc <- .multi_group_eta_squared_cor(Xb_log_filtered, Xa_log_filtered, biology_filtered)
      biology_concordance <- list(
        metric = "ANOVA_concordance",
        f_stat_correlation = f_conc$f_stat_cor,
        eta_squared_correlation = eta_conc$eta_cor,
        mean_eta_squared_before = eta_conc$mean_eta_before,
        mean_eta_squared_after = eta_conc$mean_eta_after,
        interpretation = "Higher F-stat/eta2 correlation indicates better preservation of group differences"
      )
    } else if (biology_type == "continuous") {
      # Correlation and R2 preservation
      cont_conc <- .continuous_concordance(Xb_log_filtered, Xa_log_filtered, biology_filtered)
      biology_concordance <- list(
        metric = "continuous_concordance",
        correlation_preservation = cont_conc$cor_preservation,
        r2_preservation = cont_conc$r2_preservation,
        mean_abs_cor_before = cont_conc$mean_abs_cor_before,
        mean_abs_cor_after = cont_conc$mean_abs_cor_after,
        interpretation = "Higher preservation values indicate better retention of biological associations"
      )
    }
  } else {
    biology_concordance <- list(metric = "none", value = NA_real_)
  }
  
  # === 4) RMSE between technical replicates
  rmse_before <- NA_real_
  rmse_after <- NA_real_
  if (!is.null(sample_id)) {
    if (verbose) message("Computing RMSE for technical replicates (using all samples including replicates)...")
    rmse_before <- .replicate_rmse(Xb_log, sample_id, batch)  # Use original unfiltered data
    rmse_after  <- .replicate_rmse(Xa_log, sample_id, batch)  # Use original unfiltered data
  } else {
    if (verbose) message("Skipping replicate RMSE (sample_id not provided)")
  }
  
  # === 5) PVCA - Principal Variance Component Analysis
  if (verbose) message("Computing PVCA...")
  pvca_before <- .pvca(pcs_before, batch_filtered, biology_filtered)
  pvca_after <- .pvca(pcs_after, batch_filtered, biology_filtered)
  
  # === 6) Adjusted Rand Index (ARI) - clustering agreement
  if (verbose) message("Computing Adjusted Rand Index...")
  ari_result <- .adjusted_rand_index(pcs_before, pcs_after, n_clusters = NULL)
  
  # === 7) Cluster compactness/separation (using unsupervised clustering)
  if (verbose) message("Computing cluster compactness/separation...")
  # Use the clustering labels from ARI computation
  compact_sep_before <- .cluster_compactness_separation(pcs_before, ari_result$labels_before)
  compact_sep_after <- .cluster_compactness_separation(pcs_after, ari_result$labels_after)
  
  res <- list(
    knn_jaccard = knn_jacc,
    hvg_overlap = hvg_jacc,
    biology_concordance = biology_concordance,
    biology_type = biology_type,
    n_samples_used = sum(keep_idx),
    n_samples_excluded = n_excluded,
    pvca = list(
      before = pvca_before,
      after = pvca_after
    ),
    adjusted_rand_index = ari_result$ari,
    ari_n_clusters = ari_result$n_clusters,
    cluster_compactness_separation = list(
      before = compact_sep_before,
      after = compact_sep_after
    ),
    rmse_replicates = list(
      before = rmse_before,
      after = rmse_after
    )
  )
  class(res) <- "anchor_assessment"
  res
}

# ---- Enhanced helper functions ----

.pca_scores <- function(Xlog, n_pcs = 30) {
  Xlog <- as.matrix(Xlog)
  storage.mode(Xlog) <- "double"
  Xlog[!is.finite(Xlog)] <- 0
  pr <- stats::prcomp(t(Xlog), center = TRUE, scale. = TRUE, rank. = n_pcs)
  pr$x  # samples x PCs
}

.compute_r2 <- function(y, design) {
  fit <- stats::lm.fit(design, y)
  pred <- fit$fitted.values
  rss <- colSums((y - pred)^2, na.rm = TRUE)
  tss <- colSums(scale(y, scale = FALSE)^2, na.rm = TRUE)
  1 - rss / pmax(tss, 1e-8)
}

.per_gene_r2 <- function(Xlog, fac) {
  fac <- as.factor(fac)
  design <- stats::model.matrix(~ 0 + fac)
  y <- t(Xlog)  # samples x genes
  as.numeric(.compute_r2(y, design))
}

# R2 for continuous variables
.per_gene_r2_continuous <- function(Xlog, cont_var) {
  # For each gene: regress expression ~ continuous variable
  r2_values <- apply(Xlog, 1, function(gene_expr) {
    fit <- stats::lm(gene_expr ~ cont_var)
    summary(fit)$r.squared
  })
  r2_values
}

.pc_r2 <- function(pcs, fac) {
  fac <- as.factor(fac)
  design <- stats::model.matrix(~ 0 + fac)
  y <- as.matrix(pcs)  # samples x PCs
  r2_pc <- .compute_r2(y, design)
  w <- apply(y, 2, stats::var, na.rm = TRUE)
  sum(r2_pc * w, na.rm = TRUE) / pmax(sum(w, na.rm = TRUE), 1e-8)
}

# PC R2 for continuous variables
.pc_r2_continuous <- function(pcs, cont_var) {
  # For each PC: regress PC ~ continuous variable
  y <- as.matrix(pcs)
  r2_pc <- apply(y, 2, function(pc) {
    fit <- stats::lm(pc ~ cont_var)
    summary(fit)$r.squared
  })
  w <- apply(y, 2, stats::var, na.rm = TRUE)
  sum(r2_pc * w, na.rm = TRUE) / pmax(sum(w, na.rm = TRUE), 1e-8)
}

.silhouette_mean <- function(pcs, labels) {
  labels <- as.factor(labels)
  if (nlevels(labels) < 2) return(NA_real_)
  d <- stats::dist(pcs)
  cl <- as.integer(labels)
  if (!requireNamespace("cluster", quietly = TRUE)) return(NA_real_)
  sil <- cluster::silhouette(cl, dmatrix = as.matrix(d))
  mean(sil[, "sil_width"], na.rm = TRUE)
}

.knn_jaccard <- function(pcs_before, pcs_after, k = 20) {
  if (!requireNamespace("FNN", quietly = TRUE)) {
    d1 <- as.matrix(stats::dist(pcs_before))
    d2 <- as.matrix(stats::dist(pcs_after))
    knn_idx1 <- apply(d1, 1, function(row) order(row)[2:(k+1)])
    knn_idx2 <- apply(d2, 1, function(row) order(row)[2:(k+1)])
  } else {
    knn_idx1 <- FNN::get.knn(pcs_before, k = k)$nn.index
    knn_idx2 <- FNN::get.knn(pcs_after,  k = k)$nn.index
  }
  if (is.vector(knn_idx1)) knn_idx1 <- matrix(knn_idx1, ncol = k)
  if (is.vector(knn_idx2)) knn_idx2 <- matrix(knn_idx2, ncol = k)
  n <- nrow(knn_idx1)
  jacc <- numeric(n)
  for (i in seq_len(n)) {
    a <- unique(knn_idx1[i, ])
    b <- unique(knn_idx2[i, ])
    inter <- length(intersect(a, b))
    uni <- length(union(a, b))
    jacc[i] <- if (uni == 0) NA_real_ else inter / uni
  }
  mean(jacc, na.rm = TRUE)
}

.hvg_overlap <- function(Xlog_before, Xlog_after, top_n = 2000) {
  v1 <- apply(Xlog_before, 1, stats::var, na.rm = TRUE)
  v2 <- apply(Xlog_after,  1, stats::var, na.rm = TRUE)
  g1 <- order(v1, decreasing = TRUE)[seq_len(min(top_n, length(v1)))]
  g2 <- order(v2, decreasing = TRUE)[seq_len(min(top_n, length(v2)))]
  inter <- length(intersect(g1, g2))
  uni <- length(union(g1, g2))
  inter / pmax(uni, 1L)
}

.replicate_rmse <- function(Xlog, sample_id, batch) {
  sample_id <- as.factor(sample_id)
  batch <- as.factor(batch)
  ids <- levels(sample_id)
  # One summary per anchor (sample_id spanning batches); equal weight per anchor
  rmse_by_anchor <- numeric(0)
  for (id in ids) {
    cols_id <- which(sample_id == id)
    if (length(cols_id) < 2) next
    batches_id <- split(cols_id, batch[cols_id])
    batches_id <- batches_id[lengths(batches_id) > 0]
    if (length(batches_id) < 2) next
    means_by_batch <- lapply(batches_id, function(idx) {
      if (length(idx) == 1) Xlog[, idx, drop = FALSE] else rowMeans(Xlog[, idx, drop = FALSE], na.rm = TRUE)
    })
    M <- do.call(cbind, means_by_batch)
    if (is.null(dim(M)) || ncol(M) < 2) next
    pair_idx <- utils::combn(ncol(M), 2)
    rmse_pairs <- numeric(ncol(pair_idx))
    for (p in seq_len(ncol(pair_idx))) {
      a <- M[, pair_idx[1, p], drop = FALSE]
      b <- M[, pair_idx[2, p], drop = FALSE]
      diff <- as.numeric(a - b)
      rmse_pairs[p] <- sqrt(mean(diff^2, na.rm = TRUE))
    }
    rmse_by_anchor <- c(rmse_by_anchor, mean(rmse_pairs, na.rm = TRUE))
  }
  if (length(rmse_by_anchor) == 0) return(NA_real_)
  mean(rmse_by_anchor, na.rm = TRUE)
}

.two_group_logfc_cor <- function(Xlog_before, Xlog_after, biology) {
  biology <- as.factor(biology)
  if (nlevels(biology) != 2) return(NA_real_)
  idx1 <- which(biology == levels(biology)[1])
  idx2 <- which(biology == levels(biology)[2])
  mean1_b <- rowMeans(Xlog_before[, idx1, drop = FALSE], na.rm = TRUE)
  mean2_b <- rowMeans(Xlog_before[, idx2, drop = FALSE], na.rm = TRUE)
  mean1_a <- rowMeans(Xlog_after[,  idx1, drop = FALSE], na.rm = TRUE)
  mean2_a <- rowMeans(Xlog_after[,  idx2, drop = FALSE], na.rm = TRUE)
  lfc_b <- (mean2_b - mean1_b)
  lfc_a <- (mean2_a - mean1_a)
  stats::cor(lfc_b, lfc_a, use = "pairwise.complete.obs", method = "spearman")
}

#' Multi-group F-statistic correlation
#' @keywords internal
.multi_group_fstat_cor <- function(Xlog_before, Xlog_after, biology) {
  biology <- as.factor(biology)
  
  # Calculate F-statistics for each gene
  f_stats_before <- apply(Xlog_before, 1, function(gene_expr) {
    fit <- stats::aov(gene_expr ~ biology)
    f_stat <- summary(fit)[[1]]$F[1]
    if (is.na(f_stat)) 0 else f_stat
  })
  
  f_stats_after <- apply(Xlog_after, 1, function(gene_expr) {
    fit <- stats::aov(gene_expr ~ biology)
    f_stat <- summary(fit)[[1]]$F[1]
    if (is.na(f_stat)) 0 else f_stat
  })
  
  # Correlation between F-statistics
  f_cor <- stats::cor(f_stats_before, f_stats_after, 
                      use = "pairwise.complete.obs", method = "spearman")
  
  list(f_stat_cor = f_cor, 
       f_stats_before = f_stats_before, 
       f_stats_after = f_stats_after)
}

#' Multi-group eta-squared correlation
#' @keywords internal
.multi_group_eta_squared_cor <- function(Xlog_before, Xlog_after, biology) {
  biology <- as.factor(biology)
  
  # Calculate eta-squared for each gene
  eta_squared_before <- apply(Xlog_before, 1, function(gene_expr) {
    ss_total <- sum((gene_expr - mean(gene_expr))^2, na.rm = TRUE)
    if (ss_total == 0) return(0)
    group_means <- tapply(gene_expr, biology, mean, na.rm = TRUE)
    group_sizes <- tapply(gene_expr, biology, length)
    grand_mean <- mean(gene_expr, na.rm = TRUE)
    ss_between <- sum(group_sizes * (group_means - grand_mean)^2, na.rm = TRUE)
    ss_between / ss_total
  })
  
  eta_squared_after <- apply(Xlog_after, 1, function(gene_expr) {
    ss_total <- sum((gene_expr - mean(gene_expr))^2, na.rm = TRUE)
    if (ss_total == 0) return(0)
    group_means <- tapply(gene_expr, biology, mean, na.rm = TRUE)
    group_sizes <- tapply(gene_expr, biology, length)
    grand_mean <- mean(gene_expr, na.rm = TRUE)
    ss_between <- sum(group_sizes * (group_means - grand_mean)^2, na.rm = TRUE)
    ss_between / ss_total
  })
  
  # Correlation between eta-squared values
  eta_cor <- stats::cor(eta_squared_before, eta_squared_after,
                        use = "pairwise.complete.obs", method = "spearman")
  
  list(eta_cor = eta_cor,
       mean_eta_before = mean(eta_squared_before, na.rm = TRUE),
       mean_eta_after = mean(eta_squared_after, na.rm = TRUE))
}

#' Continuous variable concordance
#' @keywords internal
.continuous_concordance <- function(Xlog_before, Xlog_after, cont_var) {
  # Calculate correlation between each gene and continuous variable
  cors_before <- apply(Xlog_before, 1, function(gene_expr) {
    stats::cor(gene_expr, cont_var, use = "pairwise.complete.obs", method = "spearman")
  })
  
  cors_after <- apply(Xlog_after, 1, function(gene_expr) {
    stats::cor(gene_expr, cont_var, use = "pairwise.complete.obs", method = "spearman")
  })
  
  # Calculate R2 for linear models
  r2_before <- apply(Xlog_before, 1, function(gene_expr) {
    fit <- stats::lm(gene_expr ~ cont_var)
    summary(fit)$r.squared
  })
  
  r2_after <- apply(Xlog_after, 1, function(gene_expr) {
    fit <- stats::lm(gene_expr ~ cont_var)
    summary(fit)$r.squared
  })
  
  # Preservation metrics
  cor_preservation <- stats::cor(cors_before, cors_after, 
                                 use = "pairwise.complete.obs", method = "spearman")
  r2_preservation <- stats::cor(r2_before, r2_after,
                                use = "pairwise.complete.obs", method = "spearman")
  
  list(
    cor_preservation = cor_preservation,
    r2_preservation = r2_preservation,
    mean_abs_cor_before = mean(abs(cors_before), na.rm = TRUE),
    mean_abs_cor_after = mean(abs(cors_after), na.rm = TRUE),
    mean_r2_before = mean(r2_before, na.rm = TRUE),
    mean_r2_after = mean(r2_after, na.rm = TRUE)
  )
}

#' Principal Variance Component Analysis (PVCA)
#' Quantifies variance explained by batch vs biological factors
#' @keywords internal
.pvca <- function(pcs, batch, biology = NULL) {
  # Use top PCs (weighted by variance explained)
  pcs <- as.matrix(pcs)
  
  # Calculate variance of each PC
  pc_var <- apply(pcs, 2, stats::var, na.rm = TRUE)
  total_var <- sum(pc_var, na.rm = TRUE)
  
  if (total_var == 0) {
    return(list(
      batch_variance = NA_real_,
      biology_variance = NA_real_,
      residual_variance = NA_real_
    ))
  }
  
  # Normalize PCs by their variance contribution
  weights <- pc_var / total_var
  
  # Compute R2 for batch on each PC
  batch_r2 <- numeric(ncol(pcs))
  for (i in seq_len(ncol(pcs))) {
    fit <- stats::lm(pcs[, i] ~ batch)
    batch_r2[i] <- summary(fit)$r.squared
  }
  batch_variance <- sum(batch_r2 * weights, na.rm = TRUE)
  
  # Compute R2 for biology on each PC (if provided)
  biology_variance <- NA_real_
  if (!is.null(biology)) {
    biology_r2 <- numeric(ncol(pcs))
    for (i in seq_len(ncol(pcs))) {
      if (is.numeric(biology)) {
        fit <- stats::lm(pcs[, i] ~ biology)
      } else {
        fit <- stats::lm(pcs[, i] ~ as.factor(biology))
      }
      biology_r2[i] <- summary(fit)$r.squared
    }
    biology_variance <- sum(biology_r2 * weights, na.rm = TRUE)
  }
  
  # Residual variance (1 - batch - biology, but account for overlap)
  # For simplicity, we compute residual as 1 - max(batch, biology) if biology is NULL
  # Otherwise, we use a model with both factors
  if (!is.null(biology)) {
    combined_r2 <- numeric(ncol(pcs))
    for (i in seq_len(ncol(pcs))) {
      if (is.numeric(biology)) {
        fit <- stats::lm(pcs[, i] ~ batch + biology)
      } else {
        fit <- stats::lm(pcs[, i] ~ batch + as.factor(biology))
      }
      combined_r2[i] <- summary(fit)$r.squared
    }
    combined_variance <- sum(combined_r2 * weights, na.rm = TRUE)
    residual_variance <- 1 - combined_variance
  } else {
    residual_variance <- 1 - batch_variance
  }
  
  list(
    batch_variance = batch_variance,
    biology_variance = biology_variance,
    residual_variance = residual_variance
  )
}

#' Adjusted Rand Index (ARI)
#' Compares clustering agreement between before and after correction
#' @keywords internal
.adjusted_rand_index <- function(pcs_before, pcs_after, n_clusters = NULL) {
  # Determine number of clusters if not provided
  # Use elbow method or default to sqrt(n/2)
  n_samples <- nrow(pcs_before)
  if (is.null(n_clusters)) {
    n_clusters <- max(2, floor(sqrt(n_samples / 2)))
    n_clusters <- min(n_clusters, n_samples - 1)
  }
  
  # Perform k-means clustering
  set.seed(123)  # For reproducibility
  km_before <- stats::kmeans(pcs_before, centers = n_clusters, iter.max = 100, nstart = 10)
  km_after <- stats::kmeans(pcs_after, centers = n_clusters, iter.max = 100, nstart = 10)
  
  labels_before <- km_before$cluster
  labels_after <- km_after$cluster
  
  # Compute ARI manually
  n <- length(labels_before)
  contingency <- table(labels_before, labels_after)
  
  # Sum over rows and columns
  sum_comb_rows <- sum(choose(table(labels_before), 2))
  sum_comb_cols <- sum(choose(table(labels_after), 2))
  sum_comb_table <- sum(choose(contingency, 2))
  
  # Expected index
  expected_index <- sum_comb_rows * sum_comb_cols / choose(n, 2)
  
  # Maximum index
  max_index <- (sum_comb_rows + sum_comb_cols) / 2
  
  # ARI
  if (max_index == expected_index) {
    ari <- 0
  } else {
    ari <- (sum_comb_table - expected_index) / (max_index - expected_index)
  }
  
  list(
    ari = ari,
    n_clusters = n_clusters,
    labels_before = labels_before,
    labels_after = labels_after
  )
}

#' Factor variance preservation
#' Assesses how well biological variance is preserved after correction
#' @keywords internal
.factor_variance_preservation <- function(pcs_before, pcs_after, biology) {
  if (is.null(biology)) {
    return(list(
      variance_preservation = NA_real_,
      variance_before = NA_real_,
      variance_after = NA_real_
    ))
  }
  
  # Calculate variance explained by biology factor on each PC
  biology_r2_before <- numeric(ncol(pcs_before))
  biology_r2_after <- numeric(ncol(pcs_after))
  
  pc_var_before <- apply(pcs_before, 2, stats::var, na.rm = TRUE)
  pc_var_after <- apply(pcs_after, 2, stats::var, na.rm = TRUE)
  
  total_var_before <- sum(pc_var_before, na.rm = TRUE)
  total_var_after <- sum(pc_var_after, na.rm = TRUE)
  
  weights_before <- pc_var_before / pmax(total_var_before, 1e-8)
  weights_after <- pc_var_after / pmax(total_var_after, 1e-8)
  
  for (i in seq_len(ncol(pcs_before))) {
    if (is.numeric(biology)) {
      fit_before <- stats::lm(pcs_before[, i] ~ biology)
      fit_after <- stats::lm(pcs_after[, i] ~ biology)
    } else {
      fit_before <- stats::lm(pcs_before[, i] ~ as.factor(biology))
      fit_after <- stats::lm(pcs_after[, i] ~ as.factor(biology))
    }
    biology_r2_before[i] <- summary(fit_before)$r.squared
    biology_r2_after[i] <- summary(fit_after)$r.squared
  }
  
  variance_before <- sum(biology_r2_before * weights_before, na.rm = TRUE)
  variance_after <- sum(biology_r2_after * weights_after, na.rm = TRUE)
  
  # Preservation ratio (how much variance is retained)
  variance_preservation <- if (variance_before > 0) {
    variance_after / variance_before
  } else {
    NA_real_
  }
  
  list(
    variance_preservation = variance_preservation,
    variance_before = variance_before,
    variance_after = variance_after
  )
}

#' Cluster compactness and separation
#' Computes within-cluster and between-cluster distances
#' @keywords internal
.cluster_compactness_separation <- function(pcs, labels) {
  if (length(unique(labels)) < 2) {
    return(list(
      within_cluster_dist = NA_real_,
      between_cluster_dist = NA_real_,
      compactness_ratio = NA_real_
    ))
  }
  
  # Compute distance matrix
  dist_matrix <- as.matrix(stats::dist(pcs))
  
  # Within-cluster distances
  within_dists <- numeric()
  for (cluster_id in unique(labels)) {
    cluster_samples <- which(labels == cluster_id)
    if (length(cluster_samples) > 1) {
      cluster_dist <- dist_matrix[cluster_samples, cluster_samples]
      # Get upper triangle (excluding diagonal)
      within_dists <- c(within_dists, cluster_dist[upper.tri(cluster_dist)])
    }
  }
  
  # Between-cluster distances
  between_dists <- numeric()
  unique_clusters <- unique(labels)
  n_clusters <- length(unique_clusters)
  
  if (n_clusters > 1) {
    for (i in seq_len(n_clusters - 1)) {
      for (j in (i + 1):n_clusters) {
        cluster_i_samples <- which(labels == unique_clusters[i])
        cluster_j_samples <- which(labels == unique_clusters[j])
        between_dist <- dist_matrix[cluster_i_samples, cluster_j_samples]
        between_dists <- c(between_dists, as.numeric(between_dist))
      }
    }
  }
  
  within_mean <- mean(within_dists, na.rm = TRUE)
  between_mean <- mean(between_dists, na.rm = TRUE)
  
  # Compactness ratio: within / between (lower is better for compactness)
  # Separation ratio: between / within (higher is better for separation)
  compactness_ratio <- if (between_mean > 0) {
    within_mean / between_mean
  } else {
    NA_real_
  }
  
  list(
    within_cluster_dist = within_mean,
    between_cluster_dist = between_mean,
    compactness_ratio = compactness_ratio,
    separation_ratio = if (within_mean > 0) between_mean / within_mean else NA_real_
  )
}

#' Print method for anchor_assessment
#'
#' @param x Object to print.
#' @param ... Additional arguments (unused).
#' @export
print.anchor_assessment <- function(x, ...) {
  cat("== anchorCorrectR: correction assessment ==\n")
  
  # Display biology type if detected
  if (!is.null(x$biology_type)) {
    cat(sprintf("\nBiology variable type: %s\n", x$biology_type))
  }
  
  if (!is.null(x$n_samples_used)) {
    cat(sprintf("\nUsing %d unique samples for metrics", x$n_samples_used))
    if (!is.null(x$n_samples_excluded) && x$n_samples_excluded > 0) {
      cat(sprintf(" (%d replicated samples excluded)\n", x$n_samples_excluded))
    } else {
      cat("\n")
    }
  }
  
  # === Block 1: Batch Removal Metrics ===
  cat("\n----------------------------------------------------------------------------\n")
  cat("  BATCH REMOVAL METRICS\n")
  cat("----------------------------------------------------------------------------\n")
  
  # RMSE replicates
  if (!is.null(x$rmse_replicates)) {
    cat("\n* RMSE (Technical Replicates)\n")
    cat(sprintf("  Before: %.3f\n", x$rmse_replicates$before))
    cat(sprintf("  After:  %.3f\n", x$rmse_replicates$after))
    cat(sprintf("  Delta:  %.3f (lower is better)\n", x$rmse_replicates$after - x$rmse_replicates$before))
  }
  
  # PVCA batch variance
  if (!is.null(x$pvca)) {
    cat("\n* PVCA Batch Variance\n")
    cat(sprintf("  Before: %.3f\n", x$pvca$before$batch_variance))
    cat(sprintf("  After:  %.3f\n", x$pvca$after$batch_variance))
    cat(sprintf("  Delta:  %.3f (negative delta = improvement)\n", 
                x$pvca$after$batch_variance - x$pvca$before$batch_variance))
  }
  
  # === Block 2: Biological Preservation Metrics ===
  cat("\n----------------------------------------------------------------------------\n")
  cat("  BIOLOGICAL PRESERVATION METRICS\n")
  cat("----------------------------------------------------------------------------\n")
  
  # PVCA biology variance
  if (!is.null(x$pvca) && !is.na(x$pvca$before$biology_variance)) {
    cat("\n* PVCA Biology Variance\n")
    cat(sprintf("  Before: %.3f\n", x$pvca$before$biology_variance))
    cat(sprintf("  After:  %.3f\n", x$pvca$after$biology_variance))
    cat(sprintf("  Delta:  %.3f (positive delta = better preservation)\n", 
                x$pvca$after$biology_variance - x$pvca$before$biology_variance))
  }
  
  # kNN Jaccard
  cat("\n* kNN Jaccard (neighborhood preservation)\n")
  cat(sprintf("  Value: %.3f (higher = better preservation of local structure)\n", x$knn_jaccard))
  
  # HVG Overlap
  cat("\n* HVG Overlap (highly variable genes)\n")
  cat(sprintf("  Value: %.3f (higher = better preservation of variable genes)\n", x$hvg_overlap))
  
  # Biology concordance
  if (!is.null(x$biology_concordance) && x$biology_concordance$metric != "none") {
    bc <- x$biology_concordance
    
    if (bc$metric == "logFC_correlation") {
      cat("\n* Biology Concordance (logFC correlation)\n")
      cat(sprintf("  Value: %.3f (higher = better preservation of differential expression)\n", bc$value))
      
    } else if (bc$metric == "ANOVA_concordance") {
      cat("\n* Biology Concordance (ANOVA)\n")
      cat(sprintf("  F-statistic correlation: %.3f\n", bc$f_stat_correlation))
      cat(sprintf("  Eta-squared correlation: %.3f\n", bc$eta_squared_correlation))
      cat(sprintf("  Mean eta2 before: %.3f, after: %.3f\n", 
                  bc$mean_eta_squared_before, bc$mean_eta_squared_after))
      
    } else if (bc$metric == "continuous_concordance") {
      cat("\n* Biology Concordance (continuous variable)\n")
      cat(sprintf("  Correlation preservation: %.3f\n", bc$correlation_preservation))
      cat(sprintf("  R2 preservation: %.3f\n", bc$r2_preservation))
      cat(sprintf("  Mean |correlation| before: %.3f, after: %.3f\n",
                  bc$mean_abs_cor_before, bc$mean_abs_cor_after))
    }
  }
  
  # === Block 3: Clustering Structure Metrics ===
  cat("\n----------------------------------------------------------------------------\n")
  cat("  CLUSTERING STRUCTURE METRICS\n")
  cat("----------------------------------------------------------------------------\n")
  
  # Adjusted Rand Index
  if (!is.null(x$adjusted_rand_index)) {
    cat("\n* Adjusted Rand Index (ARI)\n")
    cat(sprintf("  Value: %.3f (n_clusters = %d)\n", x$adjusted_rand_index, x$ari_n_clusters))
    cat("  Interpretation: Higher values indicate better preservation of clustering structure\n")
  }
  
  # Cluster compactness/separation
  if (!is.null(x$cluster_compactness_separation)) {
    cs_before <- x$cluster_compactness_separation$before
    cs_after <- x$cluster_compactness_separation$after
    if (!is.na(cs_before$compactness_ratio)) {
      cat("\n* Cluster Compactness\n")
      cat(sprintf("  Before: %.3f (within/between cluster distance ratio)\n", cs_before$compactness_ratio))
      cat(sprintf("  After:  %.3f\n", cs_after$compactness_ratio))
      cat(sprintf("  Delta:  %.3f (small delta = preserved compactness)\n", 
                  cs_after$compactness_ratio - cs_before$compactness_ratio))
      
      cat("\n* Cluster Separation\n")
      cat(sprintf("  Before: %.3f (between/within cluster distance ratio)\n", cs_before$separation_ratio))
      cat(sprintf("  After:  %.3f\n", cs_after$separation_ratio))
      cat(sprintf("  Delta:  %.3f (small delta = preserved separation)\n", 
                  cs_after$separation_ratio - cs_before$separation_ratio))
    }
  }
  
  invisible(x)
}

#' Merge multiple assess_correction results
#'
#' Combines results from multiple runs of \code{assess_correction()} into a single
#' merged result object. This is useful for aggregating assessments across different
#' datasets, correction methods, or parameter settings.
#'
#' @param ... One or more \code{anchor_assessment} objects returned by \code{assess_correction()}.
#'   Can also pass a single list of assessment objects.
#' @param run_names Optional character vector of names for each run. If NULL, runs are
#'   numbered sequentially (run_1, run_2, ...).
#'
#' @return A list with the same structure as \code{assess_correction()} output, but with:
#' \itemize{
#'   \item \code{summary}: data.frame with added \code{run} column identifying the source run.
#'   \item \code{knn_jaccard}: numeric vector with one value per run.
#'   \item \code{hvg_overlap}: numeric vector with one value per run.
#'   \item \code{biology_concordance}: list of biology concordance objects (one per run).
#'   \item \code{biology_consistency_concordance}: data.frame with biology consistency (knn_jaccard, hvg_overlap) 
#'     and concordance metrics (varies by biology type) with run column.
#'   \item \code{biology_type}: character vector with biology types (should be consistent).
#'   \item \code{n_samples_used}: numeric vector with sample counts per run.
#'   \item \code{n_samples_excluded}: numeric vector with excluded sample counts per run.
#'   \item \code{pvca}: list of PVCA results (one per run).
#'   \item \code{adjusted_rand_index}: numeric vector with ARI values per run.
#'   \item \code{ari_n_clusters}: numeric vector with number of clusters used per run.
#'   \item \code{cluster_compactness_separation}: list of cluster metrics (one per run).
#'   \item \code{comprehensive_summary}: data.frame with detailed metrics grouped by category with interpretations.
#'   \item \code{consolidated_table}: data.frame with all metrics in wide format (one row per run).
#'   \item \code{n_runs}: number of runs merged.
#' }
#' The result has class \code{anchor_assessment_merged}.
#'
#' @examples
#' \dontrun{
#' # Example with two assessment runs
#' result1 <- assess_correction(x_before1, x_after1, batch1, biology1)
#' result2 <- assess_correction(x_before2, x_after2, batch2, biology2)
#' 
#' # Merge with default names
#' merged <- merge_assessments(result1, result2)
#' 
#' # Merge with custom names
#' merged <- merge_assessments(result1, result2, 
#'                             run_names = c("dataset1", "dataset2"))
#' 
#' # Merge from a list
#' results_list <- list(result1, result2, result3)
#' merged <- merge_assessments(results_list)
#' }
#' @export
merge_assessments <- function(..., run_names = NULL) {
  # Collect all arguments
  args <- list(...)
  
  # If first argument is a list, use it directly
  if (length(args) == 1 && is.list(args[[1]]) && 
      !inherits(args[[1]], "anchor_assessment")) {
    assessments <- args[[1]]
  } else {
    assessments <- args
  }
  
  # Validate inputs
  if (length(assessments) == 0) {
    stop("No assessment objects provided")
  }
  
  # Check that all are anchor_assessment objects
  for (i in seq_along(assessments)) {
    if (!inherits(assessments[[i]], "anchor_assessment")) {
      stop(sprintf("Element %d is not an anchor_assessment object", i))
    }
  }
  
  n_runs <- length(assessments)
  
  # Generate run names if not provided
  if (is.null(run_names)) {
    run_names <- paste0("run_", seq_len(n_runs))
  } else {
    if (length(run_names) != n_runs) {
      stop(sprintf("Length of run_names (%d) does not match number of assessments (%d)", 
                   length(run_names), n_runs))
    }
  }
  
  # Check biology_type consistency (warn if inconsistent)
  biology_types <- sapply(assessments, function(x) {
    if (is.null(x$biology_type)) return(NA_character_) else x$biology_type
  })
  unique_types <- unique(biology_types[!is.na(biology_types)])
  if (length(unique_types) > 1) {
    warning(sprintf("Inconsistent biology_type across runs: %s. Using first non-NA type.", 
                    paste(unique_types, collapse = ", ")))
  }
  merged_biology_type <- if (length(unique_types) > 0) unique_types[1] else NULL
  
  # Merge data.frames with run identifiers
  merge_df_with_run <- function(df_list, run_col_name = "run") {
    result_list <- list()
    for (i in seq_along(df_list)) {
      df <- df_list[[i]]
      if (!is.null(df) && is.data.frame(df)) {
        df[[run_col_name]] <- run_names[i]
        result_list[[i]] <- df
      }
    }
    do.call(rbind, result_list)
  }
  
  # Collect numeric vectors
  knn_jaccard_vec <- sapply(assessments, function(x) x$knn_jaccard)
  names(knn_jaccard_vec) <- run_names
  
  hvg_overlap_vec <- sapply(assessments, function(x) x$hvg_overlap)
  names(hvg_overlap_vec) <- run_names
  
  # Collect biology concordance (keep as list)
  biology_concordance_list <- lapply(assessments, function(x) x$biology_concordance)
  names(biology_concordance_list) <- run_names
  
  # Collect sample counts
  n_samples_used_vec <- sapply(assessments, function(x) x$n_samples_used)
  names(n_samples_used_vec) <- run_names
  
  n_samples_excluded_vec <- sapply(assessments, function(x) x$n_samples_excluded)
  names(n_samples_excluded_vec) <- run_names
  
  # Collect new metrics
  pvca_list <- lapply(assessments, function(x) x$pvca)
  names(pvca_list) <- run_names
  
  adjusted_rand_index_vec <- sapply(assessments, function(x) {
    if (is.null(x$adjusted_rand_index)) NA_real_ else x$adjusted_rand_index
  })
  names(adjusted_rand_index_vec) <- run_names
  
  ari_n_clusters_vec <- sapply(assessments, function(x) {
    if (is.null(x$ari_n_clusters)) NA_integer_ else x$ari_n_clusters
  })
  names(ari_n_clusters_vec) <- run_names
  
  cluster_compactness_separation_list <- lapply(assessments, function(x) x$cluster_compactness_separation)
  names(cluster_compactness_separation_list) <- run_names
  
  # Create biology consistency and concordance table
  biology_consistency_concordance <- data.frame(
    run = run_names,
    knn_jaccard = knn_jaccard_vec,
    hvg_overlap = hvg_overlap_vec,
    stringsAsFactors = FALSE
  )
  
  # Extract biology concordance metrics based on type
  concordance_metrics <- list()
  for (i in seq_along(biology_concordance_list)) {
    bc <- biology_concordance_list[[i]]
    if (is.null(bc) || bc$metric == "none") {
      concordance_metrics[[i]] <- list(
        biology_concordance_metric = NA_character_,
        biology_concordance_value = NA_real_,
        biology_f_stat_correlation = NA_real_,
        biology_eta_squared_correlation = NA_real_,
        biology_correlation_preservation = NA_real_,
        biology_r2_preservation = NA_real_
      )
    } else if (bc$metric == "logFC_correlation") {
      concordance_metrics[[i]] <- list(
        biology_concordance_metric = "logFC_correlation",
        biology_concordance_value = bc$value,
        biology_f_stat_correlation = NA_real_,
        biology_eta_squared_correlation = NA_real_,
        biology_correlation_preservation = NA_real_,
        biology_r2_preservation = NA_real_
      )
    } else if (bc$metric == "ANOVA_concordance") {
      concordance_metrics[[i]] <- list(
        biology_concordance_metric = "ANOVA_concordance",
        biology_concordance_value = NA_real_,
        biology_f_stat_correlation = bc$f_stat_correlation,
        biology_eta_squared_correlation = bc$eta_squared_correlation,
        biology_correlation_preservation = NA_real_,
        biology_r2_preservation = NA_real_
      )
    } else if (bc$metric == "continuous_concordance") {
      concordance_metrics[[i]] <- list(
        biology_concordance_metric = "continuous_concordance",
        biology_concordance_value = NA_real_,
        biology_f_stat_correlation = NA_real_,
        biology_eta_squared_correlation = NA_real_,
        biology_correlation_preservation = bc$correlation_preservation,
        biology_r2_preservation = bc$r2_preservation
      )
    } else {
      concordance_metrics[[i]] <- list(
        biology_concordance_metric = NA_character_,
        biology_concordance_value = NA_real_,
        biology_f_stat_correlation = NA_real_,
        biology_eta_squared_correlation = NA_real_,
        biology_correlation_preservation = NA_real_,
        biology_r2_preservation = NA_real_
      )
    }
  }
  
  # Convert to data.frame and merge with biology_consistency_concordance
  concordance_df <- do.call(rbind, lapply(concordance_metrics, function(x) {
    data.frame(
      biology_concordance_metric = x$biology_concordance_metric,
      biology_concordance_value = x$biology_concordance_value,
      biology_f_stat_correlation = x$biology_f_stat_correlation,
      biology_eta_squared_correlation = x$biology_eta_squared_correlation,
      biology_correlation_preservation = x$biology_correlation_preservation,
      biology_r2_preservation = x$biology_r2_preservation,
      stringsAsFactors = FALSE
    )
  }))
  
  biology_consistency_concordance <- cbind(biology_consistency_concordance, concordance_df)
  
  # Create comprehensive summary table with interpretations
  comprehensive_summary <- .create_comprehensive_summary(
    assessments, run_names,
    knn_jaccard_vec, hvg_overlap_vec, biology_consistency_concordance,
    pvca_list, adjusted_rand_index_vec, ari_n_clusters_vec,
    cluster_compactness_separation_list, merged_biology_type
  )
  
  # Create consolidated results table (all metrics in one wide table)
  consolidated_table <- .create_consolidated_table(
    assessments, run_names,
    knn_jaccard_vec, hvg_overlap_vec, biology_consistency_concordance,
    pvca_list, adjusted_rand_index_vec, ari_n_clusters_vec,
    cluster_compactness_separation_list
  )
  
  # Build merged result
  res <- list(
    knn_jaccard = knn_jaccard_vec,
    hvg_overlap = hvg_overlap_vec,
    biology_concordance = biology_concordance_list,
    biology_consistency_concordance = biology_consistency_concordance,
    biology_type = merged_biology_type,
    n_samples_used = n_samples_used_vec,
    n_samples_excluded = n_samples_excluded_vec,
    pvca = pvca_list,
    adjusted_rand_index = adjusted_rand_index_vec,
    ari_n_clusters = ari_n_clusters_vec,
    cluster_compactness_separation = cluster_compactness_separation_list,
    comprehensive_summary = comprehensive_summary,
    consolidated_table = consolidated_table,
    n_runs = n_runs,
    run_names = run_names
  )
  
  class(res) <- "anchor_assessment_merged"
  res
}

#' Create comprehensive summary table with interpretations
#' @keywords internal
.create_comprehensive_summary <- function(
    assessments, run_names,
    knn_jaccard_vec, hvg_overlap_vec, biology_consistency_concordance,
    pvca_list, adjusted_rand_index_vec, ari_n_clusters_vec,
    cluster_compactness_separation_list, biology_type
) {
  n_runs <- length(run_names)
  
  # Initialize result list
  result_rows <- list()
  
  # Helper to add metric rows
  add_metric <- function(category, metric_name, values, interpretation, 
                         better_direction = "higher", format_func = function(x) sprintf("%.3f", x)) {
    for (i in seq_len(n_runs)) {
      result_rows[[length(result_rows) + 1]] <- data.frame(
        category = category,
        metric = metric_name,
        run = run_names[i],
        value = format_func(values[i]),
        interpretation = interpretation,
        better_direction = better_direction,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # 1. Batch Removal Metrics
  # RMSE replicates
  rmse_values <- sapply(assessments, function(x) {
    if (!is.null(x$rmse_replicates)) {
      x$rmse_replicates$after
    } else NA_real_
  })
  add_metric("Batch Removal", "RMSE (replicates)", rmse_values,
             "Lower = better batch correction (replicates should be more similar)",
             "lower")
  
  # PVCA batch variance
  pvca_batch_before <- sapply(pvca_list, function(x) {
    if (!is.null(x$before)) x$before$batch_variance else NA_real_
  })
  pvca_batch_after <- sapply(pvca_list, function(x) {
    if (!is.null(x$after)) x$after$batch_variance else NA_real_
  })
  pvca_batch_delta <- pvca_batch_after - pvca_batch_before
  
  for (i in seq_len(n_runs)) {
    result_rows[[length(result_rows) + 1]] <- data.frame(
      category = "Batch Removal",
      metric = "PVCA Batch Variance",
      run = run_names[i],
      value = sprintf("before=%.3f, after=%.3f (delta=%.3f)", 
                      pvca_batch_before[i], pvca_batch_after[i], pvca_batch_delta[i]),
      interpretation = "Lower after = better batch removal. Negative delta = improvement.",
      better_direction = "lower",
      stringsAsFactors = FALSE
    )
  }
  
  # 2. Biological Preservation Metrics
  # PVCA biology variance
  pvca_biol_before <- sapply(pvca_list, function(x) {
    if (!is.null(x$before) && !is.na(x$before$biology_variance)) {
      x$before$biology_variance
    } else NA_real_
  })
  pvca_biol_after <- sapply(pvca_list, function(x) {
    if (!is.null(x$after) && !is.na(x$after$biology_variance)) {
      x$after$biology_variance
    } else NA_real_
  })
  pvca_biol_delta <- pvca_biol_after - pvca_biol_before
  
  if (any(!is.na(pvca_biol_before))) {
    for (i in seq_len(n_runs)) {
      if (!is.na(pvca_biol_before[i])) {
        result_rows[[length(result_rows) + 1]] <- data.frame(
          category = "Biological Preservation",
          metric = "PVCA Biology Variance",
          run = run_names[i],
          value = sprintf("before=%.3f, after=%.3f (delta=%.3f)", 
                          pvca_biol_before[i], pvca_biol_after[i], pvca_biol_delta[i]),
          interpretation = "Higher after = better biology preservation. Positive delta = improvement.",
          better_direction = "higher",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Biology consistency metrics
  add_metric("Biological Preservation", "kNN Jaccard", knn_jaccard_vec,
             "Higher = better preservation of local neighborhood structure (0-1 scale)",
             "higher")
  
  add_metric("Biological Preservation", "HVG Overlap", hvg_overlap_vec,
             "Higher = better preservation of highly variable genes (0-1 scale)",
             "higher")
  
  # Biology concordance (varies by type)
  if (!is.null(biology_consistency_concordance)) {
    for (i in seq_len(n_runs)) {
      bc_row <- biology_consistency_concordance[biology_consistency_concordance$run == run_names[i], ]
      if (nrow(bc_row) > 0) {
        metric_type <- bc_row$biology_concordance_metric[1]
        if (!is.na(metric_type) && metric_type != "none") {
          if (metric_type == "logFC_correlation") {
            val <- bc_row$biology_concordance_value[1]
            result_rows[[length(result_rows) + 1]] <- data.frame(
              category = "Biological Preservation",
              metric = "Biology Concordance (logFC correlation)",
              run = run_names[i],
              value = sprintf("%.3f", val),
              interpretation = "Higher = better preservation of differential expression patterns",
              better_direction = "higher",
              stringsAsFactors = FALSE
            )
          } else if (metric_type == "ANOVA_concordance") {
            f_val <- bc_row$biology_f_stat_correlation[1]
            eta_val <- bc_row$biology_eta_squared_correlation[1]
            result_rows[[length(result_rows) + 1]] <- data.frame(
              category = "Biological Preservation",
              metric = "Biology Concordance (F-stat correlation)",
              run = run_names[i],
              value = sprintf("%.3f", f_val),
              interpretation = "Higher = better preservation of group differences",
              better_direction = "higher",
              stringsAsFactors = FALSE
            )
            result_rows[[length(result_rows) + 1]] <- data.frame(
              category = "Biological Preservation",
              metric = "Biology Concordance (eta2 correlation)",
              run = run_names[i],
              value = sprintf("%.3f", eta_val),
              interpretation = "Higher = better preservation of effect sizes",
              better_direction = "higher",
              stringsAsFactors = FALSE
            )
          } else if (metric_type == "continuous_concordance") {
            cor_val <- bc_row$biology_correlation_preservation[1]
            r2_val <- bc_row$biology_r2_preservation[1]
            result_rows[[length(result_rows) + 1]] <- data.frame(
              category = "Biological Preservation",
              metric = "Biology Concordance (correlation preservation)",
              run = run_names[i],
              value = sprintf("%.3f", cor_val),
              interpretation = "Higher = better preservation of gene-variable associations",
              better_direction = "higher",
              stringsAsFactors = FALSE
            )
            result_rows[[length(result_rows) + 1]] <- data.frame(
              category = "Biological Preservation",
              metric = "Biology Concordance (R2 preservation)",
              run = run_names[i],
              value = sprintf("%.3f", r2_val),
              interpretation = "Higher = better preservation of linear model fits",
              better_direction = "higher",
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }
  
  # 3. Clustering Structure Metrics
  for (i in seq_len(n_runs)) {
    if (!is.na(adjusted_rand_index_vec[i])) {
      result_rows[[length(result_rows) + 1]] <- data.frame(
        category = "Clustering Structure",
        metric = "Adjusted Rand Index (ARI)",
        run = run_names[i],
        value = sprintf("%.3f (n_clusters=%d)", adjusted_rand_index_vec[i], ari_n_clusters_vec[i]),
        interpretation = "Higher = better preservation of clustering structure (range: -1 to 1, 1 = perfect match)",
        better_direction = "higher",
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Cluster compactness/separation
  for (i in seq_len(n_runs)) {
    cs <- cluster_compactness_separation_list[[i]]
    if (!is.null(cs) && !is.na(cs$before$compactness_ratio)) {
      result_rows[[length(result_rows) + 1]] <- data.frame(
        category = "Clustering Structure",
        metric = "Cluster Compactness Ratio",
        run = run_names[i],
        value = sprintf("before=%.3f, after=%.3f (delta=%.3f)", 
                        cs$before$compactness_ratio, cs$after$compactness_ratio,
                        cs$after$compactness_ratio - cs$before$compactness_ratio),
        interpretation = "Lower = tighter clusters. Small delta = preserved compactness.",
        better_direction = "lower",
        stringsAsFactors = FALSE
      )
      result_rows[[length(result_rows) + 1]] <- data.frame(
        category = "Clustering Structure",
        metric = "Cluster Separation Ratio",
        run = run_names[i],
        value = sprintf("before=%.3f, after=%.3f (delta=%.3f)", 
                        cs$before$separation_ratio, cs$after$separation_ratio,
                        cs$after$separation_ratio - cs$before$separation_ratio),
        interpretation = "Higher = better separated clusters. Small delta = preserved separation.",
        better_direction = "higher",
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Combine all rows
  comprehensive_df <- do.call(rbind, result_rows)
  
  # Add sample counts
  n_samples_used <- sapply(assessments, function(x) x$n_samples_used)
  n_samples_excluded <- sapply(assessments, function(x) x$n_samples_excluded)
  
  for (i in seq_len(n_runs)) {
    result_rows[[length(result_rows) + 1]] <- data.frame(
      category = "Sample Information",
      metric = "Samples Used",
      run = run_names[i],
      value = sprintf("%d", n_samples_used[i]),
      interpretation = "Number of unique samples used for metrics",
      better_direction = "N/A",
      stringsAsFactors = FALSE
    )
    if (n_samples_excluded[i] > 0) {
      result_rows[[length(result_rows) + 1]] <- data.frame(
        category = "Sample Information",
        metric = "Samples Excluded",
        run = run_names[i],
        value = sprintf("%d", n_samples_excluded[i]),
        interpretation = "Number of replicated samples excluded",
        better_direction = "N/A",
        stringsAsFactors = FALSE
      )
    }
  }
  
  comprehensive_df <- do.call(rbind, result_rows)
  comprehensive_df
}

#' Create consolidated table with all metrics in wide format
#' Structure: one row per metric, one column per run
#' @keywords internal
.create_consolidated_table <- function(
    assessments, run_names,
    knn_jaccard_vec, hvg_overlap_vec, biology_consistency_concordance,
    pvca_list, adjusted_rand_index_vec, ari_n_clusters_vec,
    cluster_compactness_separation_list
) {
  n_runs <- length(run_names)
  
  # Collect all metric values
  metrics_list <- list()
  
  # 1. RMSE replicates
  rmse_values <- sapply(assessments, function(x) {
    if (!is.null(x$rmse_replicates)) {
      x$rmse_replicates$after
    } else NA_real_
  })
  metrics_list[["rmse_replicates"]] <- rmse_values
  
  # 2. PVCA batch variance (after correction)
  pvca_batch_after <- sapply(pvca_list, function(x) {
    if (!is.null(x$after)) x$after$batch_variance else NA_real_
  })
  metrics_list[["pvca_batch_variance"]] <- pvca_batch_after
  
  # PVCA batch variance delta
  pvca_batch_before <- sapply(pvca_list, function(x) {
    if (!is.null(x$before)) x$before$batch_variance else NA_real_
  })
  metrics_list[["pvca_batch_variance_delta"]] <- pvca_batch_after - pvca_batch_before
  
  # 3. PVCA biology variance (after correction)
  pvca_biol_after <- sapply(pvca_list, function(x) {
    if (!is.null(x$after) && !is.na(x$after$biology_variance)) {
      x$after$biology_variance
    } else NA_real_
  })
  metrics_list[["pvca_biology_variance"]] <- pvca_biol_after
  
  # PVCA biology variance delta
  pvca_biol_before <- sapply(pvca_list, function(x) {
    if (!is.null(x$before) && !is.na(x$before$biology_variance)) {
      x$before$biology_variance
    } else NA_real_
  })
  metrics_list[["pvca_biology_variance_delta"]] <- pvca_biol_after - pvca_biol_before
  
  # 4. kNN Jaccard
  metrics_list[["knn_jaccard"]] <- knn_jaccard_vec
  
  # 5. HVG Overlap
  metrics_list[["hvg_overlap"]] <- hvg_overlap_vec
  
  # 6. Adjusted Rand Index
  metrics_list[["adjusted_rand_index"]] <- adjusted_rand_index_vec
  
  # 7. Cluster compactness (after correction)
  compactness_after <- sapply(cluster_compactness_separation_list, function(x) {
    if (!is.null(x) && !is.na(x$after$compactness_ratio)) {
      x$after$compactness_ratio
    } else NA_real_
  })
  metrics_list[["cluster_compactness"]] <- compactness_after
  
  # Cluster compactness delta
  compactness_before <- sapply(cluster_compactness_separation_list, function(x) {
    if (!is.null(x) && !is.na(x$before$compactness_ratio)) {
      x$before$compactness_ratio
    } else NA_real_
  })
  metrics_list[["cluster_compactness_delta"]] <- compactness_after - compactness_before
  
  # 8. Cluster separation (after correction)
  separation_after <- sapply(cluster_compactness_separation_list, function(x) {
    if (!is.null(x) && !is.na(x$after$separation_ratio)) {
      x$after$separation_ratio
    } else NA_real_
  })
  metrics_list[["cluster_separation"]] <- separation_after
  
  # Cluster separation delta
  separation_before <- sapply(cluster_compactness_separation_list, function(x) {
    if (!is.null(x) && !is.na(x$before$separation_ratio)) {
      x$before$separation_ratio
    } else NA_real_
  })
  metrics_list[["cluster_separation_delta"]] <- separation_after - separation_before
  
  # 9. Biology concordance metrics (varies by type)
  biology_logfc_correlation <- rep(NA_real_, n_runs)
  biology_f_stat_correlation <- rep(NA_real_, n_runs)
  biology_eta_squared_correlation <- rep(NA_real_, n_runs)
  biology_correlation_preservation <- rep(NA_real_, n_runs)
  biology_r2_preservation <- rep(NA_real_, n_runs)
  
  if (!is.null(biology_consistency_concordance)) {
    for (i in seq_len(n_runs)) {
      bc_row <- biology_consistency_concordance[biology_consistency_concordance$run == run_names[i], ]
      if (nrow(bc_row) > 0) {
        metric_type <- bc_row$biology_concordance_metric[1]
        if (!is.na(metric_type) && metric_type != "none") {
          if (metric_type == "logFC_correlation") {
            biology_logfc_correlation[i] <- bc_row$biology_concordance_value[1]
          } else if (metric_type == "ANOVA_concordance") {
            biology_f_stat_correlation[i] <- bc_row$biology_f_stat_correlation[1]
            biology_eta_squared_correlation[i] <- bc_row$biology_eta_squared_correlation[1]
          } else if (metric_type == "continuous_concordance") {
            biology_correlation_preservation[i] <- bc_row$biology_correlation_preservation[1]
            biology_r2_preservation[i] <- bc_row$biology_r2_preservation[1]
          }
        }
      }
    }
  }
  
  # Add biology concordance metrics (only add non-NA ones)
  if (any(!is.na(biology_logfc_correlation))) {
    metrics_list[["biology_logfc_correlation"]] <- biology_logfc_correlation
  }
  if (any(!is.na(biology_f_stat_correlation))) {
    metrics_list[["biology_f_stat_correlation"]] <- biology_f_stat_correlation
  }
  if (any(!is.na(biology_eta_squared_correlation))) {
    metrics_list[["biology_eta_squared_correlation"]] <- biology_eta_squared_correlation
  }
  if (any(!is.na(biology_correlation_preservation))) {
    metrics_list[["biology_correlation_preservation"]] <- biology_correlation_preservation
  }
  if (any(!is.na(biology_r2_preservation))) {
    metrics_list[["biology_r2_preservation"]] <- biology_r2_preservation
  }
  
  # Build result data.frame: one row per metric, one column per run
  result_df <- data.frame(
    Metric = names(metrics_list),
    stringsAsFactors = FALSE
  )
  
  # Add BEFORE column (only for metrics that have before values)
  before_values <- rep(NA_real_, length(metrics_list))
  names(before_values) <- names(metrics_list)
  
  # Set before values for metrics that have them
  before_values[["rmse_replicates"]] <- if (!is.null(assessments[[1]]$rmse_replicates)) {
    assessments[[1]]$rmse_replicates$before
  } else NA_real_
  
  before_values[["pvca_batch_variance"]] <- pvca_batch_before[1]  # Use first assessment's before value
  before_values[["pvca_biology_variance"]] <- pvca_biol_before[1]  # Use first assessment's before value
  before_values[["cluster_compactness"]] <- compactness_before[1]  # Use first assessment's before value
  before_values[["cluster_separation"]] <- separation_before[1]  # Use first assessment's before value
  
  # Delta metrics don't have before values
  before_values[["pvca_batch_variance_delta"]] <- NA_real_
  before_values[["pvca_biology_variance_delta"]] <- NA_real_
  before_values[["cluster_compactness_delta"]] <- NA_real_
  before_values[["cluster_separation_delta"]] <- NA_real_
  
  # Metrics without before values
  before_values[["knn_jaccard"]] <- NA_real_
  before_values[["hvg_overlap"]] <- NA_real_
  before_values[["adjusted_rand_index"]] <- NA_real_
  
  # Biology concordance metrics don't have before values
  if ("biology_logfc_correlation" %in% names(metrics_list)) {
    before_values[["biology_logfc_correlation"]] <- NA_real_
  }
  if ("biology_f_stat_correlation" %in% names(metrics_list)) {
    before_values[["biology_f_stat_correlation"]] <- NA_real_
  }
  if ("biology_eta_squared_correlation" %in% names(metrics_list)) {
    before_values[["biology_eta_squared_correlation"]] <- NA_real_
  }
  if ("biology_correlation_preservation" %in% names(metrics_list)) {
    before_values[["biology_correlation_preservation"]] <- NA_real_
  }
  if ("biology_r2_preservation" %in% names(metrics_list)) {
    before_values[["biology_r2_preservation"]] <- NA_real_
  }
  
  result_df$BEFORE <- before_values
  
  # Add columns for each run
  for (i in seq_len(n_runs)) {
    result_df[[run_names[i]]] <- sapply(metrics_list, function(x) x[i])
  }
  
  result_df
}

#' Print method for merged assessments
#'
#' @param x Object to print.
#' @param ... Additional arguments (unused).
#' @export
print.anchor_assessment_merged <- function(x, ...) {
  cat("===============================================================================\n")
  cat("  anchorCorrectR: Merged Batch Correction Assessment\n")
  cat("===============================================================================\n")
  cat(sprintf("\nNumber of runs merged: %d\n", x$n_runs))
  cat(sprintf("Run names: %s\n", paste(x$run_names, collapse = ", ")))
  
  if (!is.null(x$biology_type)) {
    cat(sprintf("Biology variable type: %s\n", x$biology_type))
  }
  
  # Print comprehensive summary if available
  if (!is.null(x$comprehensive_summary) && nrow(x$comprehensive_summary) > 0) {
    cat("\n")
    cat("===============================================================================\n")
    cat("  COMPREHENSIVE SUMMARY WITH INTERPRETATIONS\n")
    cat("===============================================================================\n")
    
    # Group by category
    categories <- unique(x$comprehensive_summary$category)
    for (cat in categories) {
      cat(sprintf("\n----------------------------------------------------------------------------\n"))
      cat(sprintf("  %s\n", toupper(cat)))
      cat("----------------------------------------------------------------------------\n")
      
      cat_df <- x$comprehensive_summary[x$comprehensive_summary$category == cat, ]
      
      # Print by metric
      metrics <- unique(cat_df$metric)
      for (met in metrics) {
        met_df <- cat_df[cat_df$metric == met, ]
        cat(sprintf("\n  * %s\n", met))
        cat(sprintf("    Interpretation: %s\n", met_df$interpretation[1]))
        cat(sprintf("    Better direction: %s\n", met_df$better_direction[1]))
        cat("    Values:\n")
        
        # Print values for each run
        for (i in seq_len(nrow(met_df))) {
          cat(sprintf("      %-15s: %s\n", met_df$run[i], met_df$value[i]))
        }
      }
    }
    
    cat("\n")
    cat("===============================================================================\n")
    cat("  END OF SUMMARY\n")
    cat("===============================================================================\n")
  }
  
  # Print consolidated table
  if (!is.null(x$consolidated_table) && nrow(x$consolidated_table) > 0) {
    cat("\n")
    cat("===============================================================================\n")
    cat("  CONSOLIDATED RESULTS TABLE (All Metrics)\n")
    cat("===============================================================================\n")
    print(x$consolidated_table, row.names = FALSE)
    cat("\n")
  }
  
  invisible(x)
}