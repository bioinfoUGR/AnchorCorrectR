# Location-shift correction using replicate anchors across batches.
# Supports mean or median centers. Works on log; returns in input format.

#' Row-wise mean or median of a matrix
#' @keywords internal
.row_center <- function(M, center = c("mean", "median")) {
  center <- match.arg(center)
  if (is.null(dim(M))) {
    return(as.numeric(M))
  }
  if (ncol(M) == 1L) {
    return(as.numeric(M[, 1L]))
  }
  if (center == "mean") {
    return(rowMeans(M, na.rm = TRUE))
  }
  apply(M, 1L, stats::median, na.rm = TRUE)
}

#' Location-shift batch correction using replicate anchors
#'
#' Estimates per-batch offsets from the deviation of anchor samples relative to
#' their cross-batch center (mean or median). Offsets are then recentered
#' (sum-to-zero / median-centered, or relative to a reference batch) and applied
#' to all samples in each batch.
#'
#' @param x matrix genes x samples (counts or log)
#' @param batch factor, length ncol(x)
#' @param sample_id factor of replicate IDs (IDs in >=2 batches act as anchors)
#' @param input_type "auto","counts","log"
#' @param center "mean" (default) or "median"; statistic used for (1) each
#'   anchor's cross-batch center, (2) aggregating within-batch replicate
#'   contributions, (3) aggregating offsets across anchors, and (4) recentering
#'   offsets across batches when \code{ref_batch} is \code{NULL}.
#' @param ref_batch optional single batch label to use as uncorrected reference.
#'   When specified, this batch will have zero correction and other batches
#'   will be adjusted relative to it. If NULL, offsets are centered across
#'   batches (\code{center = "mean"}: sum-to-zero; \code{center = "median"}:
#'   subtract the row-wise median of batch offsets).
#' @param clip_counts logical; clip negatives to zero when returning counts
#' @param round_counts logical; round counts on return
#' @param verbose logical; print diagnostic messages
#' @return corrected matrix in the same format (counts or log) as the input
#' @export
#' @examples
#' set.seed(42)
#' counts <- matrix(rpois(50 * 20, 100), nrow = 50)
#' colnames(counts) <- paste0("S", seq_len(20))
#' rownames(counts) <- paste0("Gene", seq_len(50))
#' batch <- factor(rep(c("A", "B"), each = 10))
#' counts[, batch == "B"] <- counts[, batch == "B"] * 1.5
#' sample_id <- factor(rep(seq_len(10), 2))
#' corrected_mean <- correct_shift(counts, batch, sample_id, center = "mean",
#'                                 input_type = "counts", verbose = FALSE)
#' corrected_med <- correct_shift(counts, batch, sample_id, center = "median",
#'                                input_type = "counts", verbose = FALSE)
#' dim(corrected_mean)
correct_shift <- function(
  x, batch, sample_id,
  input_type = c("auto", "counts", "log"),
  center = c("mean", "median"),
  ref_batch = NULL,
  clip_counts = TRUE,
  round_counts = FALSE,
  verbose = TRUE
) {
  input_type <- match.arg(input_type)
  center <- match.arg(center)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  batch <- droplevels(factor(batch))
  sample_id <- droplevels(factor(sample_id))
  B <- levels(batch)
  G <- nrow(x)

  if (input_type == "auto") {
    input_type <- detect_input_type(x)
    if (verbose) {
      message(sprintf("Shift (%s): auto-detected input type: %s", center, input_type))
    }
  }

  tf <- make_transform_cpm_log1p(x, input_type = input_type)
  X <- tf$forward(x)

  ref_batch <- validate_ref_batch(ref_batch, batch, sample_id)

  anchor_info <- identify_anchors(batch, sample_id, min_batches = 2)
  anchor_ids <- anchor_info$anchor_ids

  if (verbose) {
    message(sprintf(
      "Shift (%s): using %d anchor sample IDs across %d observations",
      center,
      length(anchor_ids),
      length(anchor_info$anchor_indices)
    ))
  }

  off_gb <- matrix(0, nrow = G, ncol = length(B), dimnames = list(rownames(X), B))

  for (b in B) {
    cols_b <- which(batch == b)
    if (!length(cols_b)) next
    s_in_b <- intersect(anchor_ids, as.character(sample_id[cols_b]))
    if (!length(s_in_b)) {
      if (verbose) {
        warning(sprintf(
          "Batch '%s' has no anchor samples - will receive zero correction",
          b
        ), call. = FALSE)
      }
      next
    }

    idxs <- lapply(s_in_b, function(sid) which(sample_id == sid))
    mu_sid <- lapply(idxs, function(ix) {
      .row_center(X[, ix, drop = FALSE], center = center)
    })

    contrib_list <- vector("list", length(s_in_b))
    n_keep <- 0L
    for (k in seq_along(s_in_b)) {
      ix_all <- idxs[[k]]
      ix_b <- ix_all[batch[ix_all] == b]
      if (!length(ix_b)) next
      n_keep <- n_keep + 1L
      contrib <- X[, ix_b, drop = FALSE] - mu_sid[[k]]
      contrib_list[[n_keep]] <- .row_center(contrib, center = center)
    }
    if (n_keep == 0L) next
    contrib_mat <- do.call(cbind, contrib_list[seq_len(n_keep)])
    off_gb[, b] <- .row_center(contrib_mat, center = center)
  }

  if (!is.null(ref_batch)) {
    jref <- match(ref_batch, B)
    off_gb <- sweep(off_gb, 1, off_gb[, jref], "-")
    off_gb[, jref] <- 0
    if (verbose) {
      message(sprintf(
        "Shift (%s): reference batch '%s' will remain uncorrected",
        center, ref_batch
      ))
    }
  } else {
    off_center <- .row_center(off_gb, center = center)
    off_gb <- sweep(off_gb, 1, off_center, "-")
    if (verbose) {
      message(sprintf(
        "Shift (%s): centering offsets across batches to preserve global gene levels",
        center
      ))
    }
  }

  Xc <- X
  for (j in seq_along(B)) {
    cols <- which(batch == B[j])
    if (!length(cols)) next
    off <- off_gb[, j]
    off[!is.finite(off)] <- 0
    Xc[, cols] <- sweep(Xc[, cols, drop = FALSE], 1, off, "-")
  }

  if (tf$input_type == "counts") {
    out <- tf$inverse_counts_from_log(Xc)
    out <- .post_counts_guard(out, clip_counts, round_counts)
    return(out)
  }
  Xc
}
