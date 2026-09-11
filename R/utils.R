# Shared utility functions for anchorCorrectR
# ============================================

#' Detect input type automatically
#' 
#' Heuristic: non-negative values that are (almost) integers are treated as raw
#' counts — including sparse or low-depth data where both max and high quantiles
#' can be small. Log-scale matrices are usually continuous on a bounded range;
#' those are classified after excluding integer-like count matrices.
#' 
#' If detection is wrong for your data (e.g. expected counts as floats), set
#' \code{input_type} explicitly in \code{anchor_correct()} / \code{make_transform_cpm_log1p()}.
#' 
#' @param x numeric matrix
#' @return "counts" or "log"
#' @keywords internal
detect_input_type <- function(x) {
  x <- as.matrix(x)
  x_finite <- x[is.finite(x)]
  
  if (length(x_finite) == 0) {
    stop("Matrix contains no finite values")
  }
  
  is_non_negative <- all(x_finite >= 0)
  if (!is_non_negative) {
    return("log")
  }
  
  # Proportion of values that round to themselves (float-safe integer test)
  prop_integer <- mean(abs(x_finite - round(x_finite)) < 1e-5)
  is_integer_like <- prop_integer >= 0.97
  
  # Raw counts first: sparse / UMI / shallow sequencing often have max <= 50 and
  # q99 < 30; an early "q99 < 30 => log" rule misclassified those as log.
  if (is_integer_like) {
    return("counts")
  }
  
  # Continuous non-negative (e.g. log1p-CPM, log-expression, or rounded floats)
  return("log")
}

#' Validate reference batch
#' 
#' Checks that ref_batch exists in the batch factor and issues warnings
#' about anchor coverage.
#' 
#' @param ref_batch character or NULL
#' @param batch factor
#' @param sample_id factor
#' @return validated ref_batch (character) or NULL
#' @keywords internal
validate_ref_batch <- function(ref_batch, batch, sample_id) {
  if (is.null(ref_batch)) {
    return(NULL)
  }
  
  ref_batch <- as.character(ref_batch)
  batch_levels <- levels(as.factor(batch))
  
  # Check existence
  if (!ref_batch %in% batch_levels) {
    stop(sprintf(
      "ref_batch '%s' not found in batch levels: %s",
      ref_batch,
      paste(batch_levels, collapse = ", ")
    ))
  }
  
  # Check if ref_batch has any samples
  ref_samples <- sample_id[batch == ref_batch]
  if (length(ref_samples) == 0) {
    stop(sprintf(
      "ref_batch '%s' has no samples in the dataset",
      ref_batch
    ))
  }
  
  # Check anchor coverage
  tab <- table(sample_id, batch) > 0
  anchor_ids <- rownames(tab)[rowSums(tab) >= 2]
  ref_has_anchors <- any(as.character(ref_samples) %in% anchor_ids)
  
  if (!ref_has_anchors) {
    warning(sprintf(
      paste0("ref_batch '%s' has no anchor samples (replicates spanning multiple batches).\n",
             "  Batch effects relative to this batch may not be well-estimated.\n",
             "  Consider using a batch with anchor samples as reference."),
      ref_batch
    ))
  }
  
  return(ref_batch)
}

#' Identify anchor samples and check coverage
#' 
#' @param batch factor
#' @param sample_id factor
#' @param min_batches minimum number of batches a sample must appear in to be an anchor
#' @return list with anchor_ids, anchor_indices, coverage_report
#' @keywords internal
identify_anchors <- function(batch, sample_id, min_batches = 2) {
  # Clean inputs - trim whitespace and handle case sensitivity
  batch_clean <- trimws(as.character(batch))
  sample_id_clean <- trimws(as.character(sample_id))
  
  batch <- factor(batch_clean)
  sample_id <- factor(sample_id_clean)
  
  # Cross-tabulation
  tab <- table(sample_id, batch) > 0
  anchor_ids <- rownames(tab)[rowSums(tab) >= min_batches]
  
  if (length(anchor_ids) == 0) {
    # Diagnostic: show the sample distribution
    cat("\n=== ANCHOR DIAGNOSTIC ===\n")
    cat("Sample ID distribution across batches:\n")
    print(table(sample_id, batch))
    cat("\nUnique sample_ids per batch:\n")
    for (b in levels(batch)) {
      sids <- unique(sample_id[batch == b])
      cat(sprintf("  Batch %s: %s\n", b, paste(sids, collapse=", ")))
    }
    cat("========================\n\n")
    
    stop(sprintf(
      paste0("No anchor samples found. Samples must appear in at least %d batches to serve as anchors.\n",
             "  Check that your sample_id labels are consistent across batches.\n",
             "  Anchor samples (biological/technical replicates) are required for batch correction.\n",
             "  See diagnostic output above for sample distribution."),
      min_batches
    ))
  }
  
  anchor_indices <- which(sample_id %in% anchor_ids)
  
  # Generate coverage report
  batches_with_anchors <- colnames(tab)[colSums(tab[anchor_ids, , drop = FALSE]) > 0]
  batches_without_anchors <- setdiff(levels(batch), batches_with_anchors)
  
  n_anchors_per_batch <- colSums(tab[anchor_ids, , drop = FALSE])
  
  coverage_report <- list(
    n_anchor_samples = length(anchor_ids),
    n_anchor_observations = length(anchor_indices),
    batches_with_anchors = batches_with_anchors,
    batches_without_anchors = batches_without_anchors,
    n_anchors_per_batch = n_anchors_per_batch
  )
  
  # Warnings for poor coverage
  if (length(batches_without_anchors) > 0) {
    # Additional diagnostic for batches without anchors
    cat("\n=== BATCHES WITHOUT ANCHORS DIAGNOSTIC ===\n")
    for (b in batches_without_anchors) {
      sids_in_batch <- unique(sample_id[batch == b])
      cat(sprintf("Batch %s:\n", b))
      cat(sprintf("  Sample IDs in this batch: %s\n", paste(sids_in_batch, collapse=", ")))
      # Check if any of these sample_ids exist in other batches
      for (sid in sids_in_batch) {
        other_batches <- unique(batch[sample_id == sid & batch != b])
        if (length(other_batches) > 0) {
          cat(sprintf("    '%s' also appears in: %s\n", sid, paste(other_batches, collapse=", ")))
        } else {
          cat(sprintf("    '%s' is unique to this batch (not an anchor)\n", sid))
        }
      }
    }
    cat("==========================================\n\n")
    
    warning(sprintf(
      paste0("The following batches have no anchor samples and cannot be directly corrected: %s\n",
             "  These batches will receive zero correction. Consider adding replicate samples.\n",
             "  See diagnostic output above for details."),
      paste(batches_without_anchors, collapse = ", ")
    ))
  }
  
  if (any(n_anchors_per_batch < 3)) {
    sparse_batches <- names(n_anchors_per_batch)[n_anchors_per_batch < 3]
    warning(sprintf(
      paste0("Batches with fewer than 3 anchor observations may have unstable corrections: %s\n",
             "  Consider adding more replicates for better estimation."),
      paste(sparse_batches, collapse = ", ")
    ))
  }
  
  list(
    anchor_ids = anchor_ids,
    anchor_indices = anchor_indices,
    coverage = coverage_report
  )
}

#' Post-process counts (guard against invalid values)
#' 
#' @param x numeric matrix
#' @param clip_counts logical; clip negative values to zero
#' @param round_counts logical; round to integers
#' @return processed matrix
#' @keywords internal
.post_counts_guard <- function(x, clip_counts = TRUE, round_counts = TRUE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x[!is.finite(x)] <- 0
  if (clip_counts) x[x < 0] <- 0
  if (round_counts) x <- round(x)
  x
}