# Helper: CPM + log1p transform (counts) and identity (log).
# The inverse returns CPM-like values (not integer counts).

#' Create forward/inverse transforms for count or log data
#' 
#' This function creates transformation functions for converting between
#' count and log-scale data. For counts, it applies CPM normalization
#' followed by log1p transformation. For log data, it acts as identity.
#' 
#' @param x numeric matrix (genes x samples)
#' @param input_type "auto", "counts", or "log"
#' @param scale_factor CPM scaling factor (default 1e6)
#' @param pseudocount additional pseudocount before log (default 0)
#' @return list with:
#'   - input_type: detected or specified type
#'   - forward: function to transform to log scale
#'   - inverse_counts_from_log: function to transform back to count-like values
#' @keywords internal
#' @examples
#' \dontrun{
#' # For count data
#' counts <- matrix(rpois(1000, 50), nrow = 100)
#' tf <- make_transform_cpm_log1p(counts, input_type = "counts")
#' log_data <- tf$forward(counts)
#' counts_back <- tf$inverse_counts_from_log(log_data)
#' 
#' # For log data
#' log_data <- matrix(rnorm(1000, 5, 2), nrow = 100)
#' tf <- make_transform_cpm_log1p(log_data, input_type = "log")
#' log_out <- tf$forward(log_data)  # Identity transform
#' }
make_transform_cpm_log1p <- function(x, input_type = c("auto","counts","log"),
                                     scale_factor = 1e6, pseudocount = 0) {
  input_type <- match.arg(input_type)
  x <- as.matrix(x)
  storage.mode(x) <- "double"

  # Use centralized auto-detection
  if (input_type == "auto") {
    input_type <- detect_input_type(x)
  }

  # Forward transform: counts to log
  forward_counts <- function(mat) {
    mat <- as.matrix(mat)
    storage.mode(mat) <- "double"
    lib <- colSums(mat, na.rm = TRUE)
    lib[lib <= 0 | !is.finite(lib)] <- 1
    cpm <- t(t(mat) / lib) * scale_factor
    log1p(pmax(cpm + pseudocount, 0))
  }

  # Inverse transform: log back to count-like values
  inverse_counts_from_log <- function(logmat) {
    # Back to CPM-like values, then scale by library sizes
    tmp <- expm1(as.matrix(logmat))
    tmp <- pmax(tmp - pseudocount, 0)  # remove pseudocount
    # Scale back using original library sizes
    lib <- colSums(x, na.rm = TRUE)
    lib[lib <= 0] <- 1
    counts <- t(t(tmp) * lib / scale_factor)
    counts
  }

  # Identity transforms for log data
  identity_transform <- function(m) {
    m <- as.matrix(m)
    storage.mode(m) <- "double"
    m
  }

  list(
    input_type = input_type,
    forward = if (input_type == "counts") forward_counts else identity_transform,
    inverse_counts_from_log = inverse_counts_from_log
  )
}