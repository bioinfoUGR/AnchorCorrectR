#!/usr/bin/env Rscript
# Example: Modifying batch correction to preserve variance
# This script shows how to modify correction methods to retain more variability

#' Modified anchor-ComBat with variance preservation options
#' 
#' This is a modified version of correct_combat_anchor that adds options to:
#' 1. Disable scale correction (location-only correction)
#' 2. Preserve original variance after correction
#' 
#' To use this, you can either:
#' - Copy this function and modify correct_combat.R
#' - Or use the wrapper functions below
#' 
#' @param apply_scale_correction logical; if FALSE, only location correction is applied
#' @param preserve_variance logical; if TRUE, scales corrected data to preserve original variance
correct_combat_anchor_variance_preserving <- function(
  x, batch, sample_id,
  input_type = c("auto","counts","log"),
  ref_batch = NULL,
  clip_counts = TRUE,
  round_counts = FALSE,
  apply_scale_correction = FALSE,  # NEW: Control scale correction
  preserve_variance = FALSE,        # NEW: Preserve original variance
  verbose = TRUE
) {
  # This would be the modified version of correct_combat_anchor
  # For now, this is a placeholder showing the key modifications needed
  
  # Key modifications needed in correct_combat.R:
  # 1. Add parameters to function signature (lines 50-57)
  # 2. Modify Step 6 (line 275) to conditionally apply scale correction
  # 3. Add Step 8 after unstandardization to preserve variance if requested
  
  stop("This is a template. See VARIANCE_PRESERVATION_ANALYSIS.md for implementation details.")
}

# ============================================================================
# WRAPPER FUNCTIONS: Use these without modifying the package code
# ============================================================================

#' Wrapper: Location-only ComBat correction (preserves variance better)
#' 
#' Applies only location correction, skipping scale correction.
#' This preserves more variance than full ComBat correction.
#' 
#' @param x matrix genes x samples
#' @param batch factor of batch labels
#' @param sample_id factor of sample IDs
#' @param ... other arguments passed to correct_combat_anchor
#' @return corrected matrix
#' @export
correct_combat_location_only <- function(x, batch, sample_id, ...) {
  # For now, use shift which already does location-only correction
  # This is equivalent to ComBat with scale correction disabled
  message("Using shift method (location-only correction) to preserve variance...")
  return(anchor_correct(x, batch, sample_id, method = "shift", ...))
}

#' Wrapper: Variance-preserving correction
#' 
#' Applies correction and then scales to preserve original variance.
#' WARNING: This may reintroduce some batch effects.
#' 
#' @param x matrix genes x samples
#' @param batch factor of batch labels  
#' @param sample_id factor of sample IDs
#' @param method correction method: "shift", "ridge", or "combat"
#' @param ... other arguments
#' @return corrected matrix with preserved variance
#' @export
correct_with_variance_preservation <- function(x, batch, sample_id, 
                                                 method = "shift", ...) {
  # Apply correction
  Xc <- anchor_correct(x, batch, sample_id, method = method, ...)
  
  # Transform to log scale for variance computation
  input_type <- detect_input_type(x)
  tf <- make_transform_cpm_log1p(x, input_type = input_type)
  X_before <- tf$forward(x)
  X_after <- tf$forward(Xc)
  
  # Compute variance before and after
  var_before <- apply(X_before, 1, var, na.rm = TRUE)
  var_after <- apply(X_after, 1, var, na.rm = TRUE)
  
  # Compute scaling factor to preserve variance
  variance_scale <- sqrt(var_before / pmax(var_after, 1e-8))
  variance_scale[!is.finite(variance_scale)] <- 1
  
  # Scale to preserve variance
  X_after_scaled <- sweep(X_after, 1, variance_scale, "*")
  
  # Recenter to preserve means
  mean_before <- rowMeans(X_before, na.rm = TRUE)
  mean_after <- rowMeans(X_after_scaled, na.rm = TRUE)
  X_after_scaled <- sweep(X_after_scaled, 1, (mean_after - mean_before), "-")
  
  # Transform back to original format
  if (input_type == "counts") {
    out <- tf$inverse_counts_from_log(X_after_scaled)
    return(out)
  } else {
    return(X_after_scaled)
  }
}

#' Wrapper: Partial correction (less aggressive)
#' 
#' Applies correction with a damping factor to preserve more variance.
#' 
#' @param x matrix genes x samples
#' @param batch factor of batch labels
#' @param sample_id factor of sample IDs
#' @param method correction method
#' @param correction_strength numeric between 0 and 1; 1 = full correction, 0.5 = half correction
#' @param ... other arguments
#' @return partially corrected matrix
#' @export
correct_partial <- function(x, batch, sample_id, method = "shift",
                            correction_strength = 0.8, ...) {
  # Apply full correction
  Xc_full <- anchor_correct(x, batch, sample_id, method = method, ...)
  
  # Transform to log scale
  input_type <- detect_input_type(x)
  tf <- make_transform_cpm_log1p(x, input_type = input_type)
  X_before <- tf$forward(x)
  X_after <- tf$forward(Xc_full)
  
  # Interpolate between original and fully corrected
  X_partial <- X_before * (1 - correction_strength) + X_after * correction_strength
  
  # Transform back
  if (input_type == "counts") {
    out <- tf$inverse_counts_from_log(X_partial)
    return(out)
  } else {
    return(X_partial)
  }
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if (FALSE) {
  # Example 1: Use shift (already preserves variance well)
  corrected <- anchor_correct(x, batch, sample_id, method = "shift")
  
  # Example 2: Use location-only wrapper
  corrected <- correct_combat_location_only(x, batch, sample_id)
  
  # Example 3: Use variance-preserving wrapper
  corrected <- correct_with_variance_preservation(x, batch, sample_id, method = "shift")
  
  # Example 4: Use partial correction (80% correction strength)
  corrected <- correct_partial(x, batch, sample_id, method = "shift", 
                                correction_strength = 0.8)
  
  # Example 5: Compare variance retention
  library(anchorCorrectR)
  
  # Before correction
  X_before <- make_transform_cpm_log1p(x, input_type = "auto")$forward(x)
  var_before <- apply(X_before, 1, var)
  
  # After correction (shift)
  corrected_ms <- anchor_correct(x, batch, sample_id, method = "shift")
  X_after_ms <- make_transform_cpm_log1p(corrected_ms, input_type = "auto")$forward(corrected_ms)
  var_after_ms <- apply(X_after_ms, 1, var)
  
  # After correction (ComBat - if available)
  corrected_cb <- anchor_correct(x, batch, sample_id, method = "combat")
  X_after_cb <- make_transform_cpm_log1p(corrected_cb, input_type = "auto")$forward(corrected_cb)
  var_after_cb <- apply(X_after_cb, 1, var)
  
  # Compare variance retention
  cat("Shift variance retention:", mean(var_after_ms / var_before, na.rm = TRUE), "\n")
  cat("ComBat variance retention:", mean(var_after_cb / var_before, na.rm = TRUE), "\n")
  
  # MDS plots
  library(stats)
  dist_before <- dist(t(X_before))
  mds_before <- cmdscale(dist_before, k = 2)
  
  dist_after_ms <- dist(t(X_after_ms))
  mds_after_ms <- cmdscale(dist_after_ms, k = 2)
  
  dist_after_cb <- dist(t(X_after_cb))
  mds_after_cb <- cmdscale(dist_after_cb, k = 2)
  
  par(mfrow = c(1, 3))
  plot(mds_before, col = as.numeric(batch), main = "Before", pch = 19)
  plot(mds_after_ms, col = as.numeric(batch), main = "After: Mean-Shift", pch = 19)
  plot(mds_after_cb, col = as.numeric(batch), main = "After: ComBat", pch = 19)
}

