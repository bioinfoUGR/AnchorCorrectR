# Helper file: MDS utilities for anchorCorrectR
# ------------------------------------------------

#' Compute sample-level MDS on log-scale data
#'
#' Transforms input to log scale (CPM + log1p if counts), computes a sample-sample
#' distance (Euclidean or 1 - Pearson correlation), then runs classical MDS
#' (\code{cmdscale}) to obtain low-dimensional coordinates.
#'
#' @param x Matrix (genes x samples), counts or log-scale.
#' @param input_type One of \code{c("auto","counts","log")}. If "auto",
#'   the transform is chosen based on the range of values.
#' @param distance One of \code{c("euclidean","1-corr")}. If "1-corr", uses
#'   \code{1 - cor(X)} as a distance on samples.
#' @param k Integer, number of MDS dimensions to return (default 2).
#'
#' @return A \code{data.frame} with columns \code{sample}, \code{Dim1}, \code{Dim2}, ... up to \code{Dim k}.
#' @export
mds_samples <- function(x, input_type = c("auto","counts","log"),
                        distance = c("euclidean","1-corr"), k = 2) {
  input_type <- match.arg(input_type)
  distance <- match.arg(distance)

  tr <- make_transform_cpm_log1p(x, input_type = input_type)
  Xlog <- tr$forward(x)
  Xlog[!is.finite(Xlog)] <- 0

  # distance on samples
  if (distance == "euclidean") {
    d <- stats::dist(t(Xlog))
  } else {
    cmat <- stats::cor(Xlog, use = "pairwise.complete.obs")
    cmat[!is.finite(cmat)] <- 0
    d <- stats::as.dist(1 - pmax(pmin(cmat, 1), -1))
  }

  # try MDS; fall back to 1D and pad to k columns if needed
  fit <- try(stats::cmdscale(d, k = k, eig = FALSE), silent = TRUE)
  if (inherits(fit, "try-error")) fit <- stats::cmdscale(d, k = 1, eig = FALSE)

  df <- as.data.frame(fit)
  if (ncol(df) < k) for (j in (ncol(df)+1):k) df[[j]] <- 0
  colnames(df) <- paste0("Dim", seq_len(ncol(df)))
  df$sample <- if (!is.null(colnames(Xlog))) colnames(Xlog) else paste0("S", seq_len(nrow(df)))
  df[, c("sample", paste0("Dim", seq_len(k)))]
}

#' Plot MDS before vs after, side by side, colored by batch (optionally shaped by biology)
#'
#' Produces two panels (before/after) with identical axis limits so the layouts are comparable.
#' If ggplot2 is available, returns a ggplot object with facet panels; otherwise, draws base R plots.
#'
#' @param x_before Matrix (genes x samples) before correction.
#' @param x_after  Matrix (genes x samples) after correction.
#' @param batch    Factor/character of length \code{ncol(x_before)} used for point color.
#' @param biology  Optional factor/character for point shape.
#' @param input_type One of \code{c("auto","counts","log")}.
#' @param distance One of \code{c("euclidean","1-corr")}.
#' @param k Integer, MDS dimensions to compute (2 recommended for plotting).
#' @param same_limits Logical; enforce the same x/y limits across both panels (default TRUE).
#'
#' @return A ggplot object (if ggplot2 is installed) or \code{NULL} (after drawing base plots).
#' @export
plot_mds_before_after <- function(x_before, x_after, batch, biology = NULL,
                                  input_type = c("auto","counts","log"),
                                  distance = c("euclidean","1-corr"), k = 2,
                                  same_limits = TRUE) {
  input_type <- match.arg(input_type)
  distance <- match.arg(distance)
  stopifnot(ncol(x_before) == ncol(x_after))
  stopifnot(length(batch) == ncol(x_before))
  if (!is.null(biology)) stopifnot(length(biology) == ncol(x_before))

  batch <- as.factor(batch)
  biology <- if (is.null(biology)) NULL else as.factor(biology)

  # Make names unique
  colnames(x_before) <- make.unique(colnames(x_before))
  colnames(x_after) <- make.unique(colnames(x_after))

  # Compute MDS for both states
  mds_b <- mds_samples(x_before, input_type = input_type, distance = distance, k = k)
  mds_a <- mds_samples(x_after,  input_type = input_type, distance = distance, k = k)
  mds_b$state <- "before"
  mds_a$state <- "after"
  mds <- rbind(mds_b, mds_a)

  # Reorder states
  mds$state <- factor(mds$state, levels = c("before", "after"))

  # Attach aesthetics (map by column name)
  colmap <- colnames(x_before)
  mds$batch <- batch[match(mds$sample, colmap)]
  if (!is.null(biology)) mds$biology <- biology[match(mds$sample, colmap)]

  # Unified limits
  if (same_limits) {
    xr <- range(mds$Dim1, na.rm = TRUE)
    yr <- range(mds$Dim2, na.rm = TRUE)
  } else {
    xr <- yr <- NULL
  }

  # ggplot2 path
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2 <- asNamespace("ggplot2")
    aes_map <- ggplot2$aes(x = .data$Dim1, y = .data$Dim2, color = .data$batch)
    if (!is.null(biology)) {
      aes_map <- ggplot2$aes(x = .data$Dim1, y = .data$Dim2, color = .data$batch, shape = .data$biology)
    }
    gg <- ggplot2$ggplot(mds, aes_map) +
      ggplot2$geom_point(size = 2, alpha = 0.9) +
      ggplot2$facet_grid(rows = NULL, cols = ggplot2$vars(state)) +  # side by side
      ggplot2$labs(title = "MDS before vs after",
                   x = "MDS 1", y = "MDS 2", color = "Batch",
                   shape = if (!is.null(biology)) "Biology" else NULL) +
      ggplot2$theme_bw() +
      ggplot2$theme(panel.grid = ggplot2$element_blank())
    if (same_limits) {
      gg <- gg + ggplot2$xlim(xr) + ggplot2$ylim(yr)
    }
    return(gg)
  }

  # Base R fallback: side-by-side panels with shared limits
  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))
  graphics::par(mfrow = c(1,2))

  if (is.null(xr)) {
    xr_b <- range(mds_b$Dim1, na.rm = TRUE)
    yr_b <- range(mds_b$Dim2, na.rm = TRUE)
    xr_a <- range(mds_a$Dim1, na.rm = TRUE)
    yr_a <- range(mds_a$Dim2, na.rm = TRUE)
  } else {
    xr_b <- xr_a <- xr
    yr_b <- yr_a <- yr
  }

  # before
  cols_b <- as.integer(as.factor(mds_b$batch))
  pch_b <- if (!is.null(biology)) as.integer(as.factor(mds_b$biology)) else 16
  graphics::plot(mds_b$Dim1, mds_b$Dim2, xlim = xr_b, ylim = yr_b,
                 col = cols_b, pch = pch_b, main = "MDS - before",
                 xlab = "MDS 1", ylab = "MDS 2")
  graphics::legend("topright", legend = levels(batch),
                   col = seq_along(levels(batch)), pch = 16, title = "Batch", cex = 0.8)

  # after
  cols_a <- as.integer(as.factor(mds_a$batch))
  pch_a <- if (!is.null(biology)) as.integer(as.factor(mds_a$biology)) else 16
  graphics::plot(mds_a$Dim1, mds_a$Dim2, xlim = xr_a, ylim = yr_a,
                 col = cols_a, pch = pch_a, main = "MDS - after",
                 xlab = "MDS 1", ylab = "MDS 2")
  graphics::legend("topright", legend = levels(batch),
                   col = seq_along(levels(batch)), pch = 16, title = "Batch", cex = 0.8)

  invisible(NULL)
}
