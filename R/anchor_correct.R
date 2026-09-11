# Unified entry point. Default method: "shift".
# Internally uses log processing and returns in the same format as input.

#' Anchor-based batch correction
#'
#' Methods:
#' - "shift" (default): replicate-anchored location-shift removal (mean or median).
#' - "ridge": glmnet ridge on batch dummies.
#' - "combat": anchor-aware ComBat with empirical Bayes shrinkage.
#'
#' All methods process internally on the log scale (log1p-CPM if counts),
#' returning the output in the same format as the input (counts or log).
#'
#' @param x numeric matrix genes x samples (counts or log)
#' @param batch factor of batch labels (length = ncol(x))
#' @param sample_id factor of replicate IDs (IDs in >=2 batches act as anchors)
#' @param method correction method: "shift", "ridge", or "combat"
#' @param input_type "auto","counts","log"
#' @param ref_batch Optional single batch label to treat as the uncorrected reference.
#'   When set, offsets are recentered so the reference batch has zero shift.
#'   The reference batch will remain unchanged after correction.
#'   Supported for "shift", "ridge", and "combat".
#'   WARNING: If ref_batch has no anchor samples, corrections may be poorly estimated.
#' @param clip_counts logical; clip negatives to zero when returning counts
#' @param round_counts logical; round counts on return
#' @param verbose logical; print diagnostic messages
#' @param ... extra args passed to method-specific functions
#'   (e.g. \code{center = "mean"} or \code{center = "median"} for \code{method = "shift"}).
#' @return corrected matrix in the same format (counts or log) as the input
#' @export
#' @examples
#' set.seed(123)
#' n_genes <- 100
#' n_samples_per_batch <- 8
#' counts <- matrix(rpois(n_genes * n_samples_per_batch * 2, lambda = 100),
#'                  nrow = n_genes)
#' colnames(counts) <- paste0("S", seq_len(ncol(counts)))
#' rownames(counts) <- paste0("Gene", seq_len(n_genes))
#' batch <- factor(rep(c("Batch1", "Batch2"), each = n_samples_per_batch))
#' sample_id <- factor(rep(seq_len(n_samples_per_batch), 2))
#' corrected <- anchor_correct(counts, batch, sample_id, method = "shift",
#'                             center = "mean", verbose = FALSE)
#' corrected_med <- anchor_correct(counts, batch, sample_id, method = "shift",
#'                                 center = "median", verbose = FALSE)
#' dim(corrected)
anchor_correct <- function(
    x, batch, sample_id,
    method = c("shift", "ridge", "combat"),
    input_type = c("auto", "counts", "log"),
    ref_batch = NULL,
    clip_counts = TRUE,
    round_counts = TRUE,
    verbose = TRUE,
    ...
) {
  method <- match.arg(method)
  input_type <- match.arg(input_type)

  x <- as.matrix(x)
  storage.mode(x) <- "double"
  stopifnot(
    "Number of columns in x must match length of batch" = ncol(x) == length(batch),
    "Length of batch must match length of sample_id" = length(batch) == length(sample_id)
  )

  batch <- droplevels(factor(batch))
  sample_id <- droplevels(factor(sample_id))

  if (input_type == "auto") {
    input_type <- detect_input_type(x)
    if (verbose) message(sprintf("Auto-detected input type: %s", input_type))
  }

  ref_batch <- validate_ref_batch(ref_batch, batch, sample_id)

  anchor_info <- identify_anchors(batch, sample_id, min_batches = 2)
  if (verbose) {
    message(sprintf(
      "Found %d anchor sample IDs across %d observations",
      anchor_info$coverage$n_anchor_samples,
      anchor_info$coverage$n_anchor_observations
    ))
    if (length(anchor_info$coverage$batches_without_anchors) > 0) {
      message(sprintf(
        "Note: %d batches have no anchors: %s",
        length(anchor_info$coverage$batches_without_anchors),
        paste(anchor_info$coverage$batches_without_anchors, collapse = ", ")
      ))
    }
  }

  if (method == "shift") {
    return(correct_shift(
      x = x, batch = batch, sample_id = sample_id,
      input_type = input_type,
      ref_batch = ref_batch,
      clip_counts = clip_counts,
      round_counts = round_counts,
      verbose = verbose,
      ...
    ))
  }

  if (method == "ridge") {
    tf <- make_transform_cpm_log1p(x, input_type = input_type)
    Xlog <- tf$forward(x)

    fit <- fit_anchor_ridge(
      Xlog = Xlog, batch = batch, sample_id = sample_id,
      verbose = verbose, ...
    )
    Xc <- apply_correction_ridge(
      fit = fit, Xlog = Xlog, batch = batch, ref_batch = ref_batch
    )

    if (tf$input_type == "counts") {
      out <- tf$inverse_counts_from_log(Xc)
      out <- .post_counts_guard(out, clip_counts, round_counts)
      return(out)
    }
    return(Xc)
  }

  if (method == "combat") {
    return(correct_combat_anchor(
      x = x, batch = batch, sample_id = sample_id,
      input_type = input_type,
      ref_batch = ref_batch,
      clip_counts = clip_counts,
      round_counts = round_counts,
      verbose = verbose,
      ...
    ))
  }

  stop("Unknown method: ", method)
}
