#' Detect Outliers in Gene Expression Data using PCA and MDS
#'
#' Iteratively identifies outlier samples that deviate beyond a threshold
#' in either PCA or MDS space. Stops when no new outliers are detected.
#'
#' @param data_matrix Numeric matrix of expression values (genes in rows, samples in columns).
#'   Row names should be gene IDs, column names sample IDs.
#' @param n_sd Numeric. Number of standard deviations for outlier threshold (default 5).
#' @param min_gene_count Numeric. Minimum count value for gene filtering (default 10).
#' @param min_gene_samples Integer. Minimum samples with min_gene_count for gene inclusion (default 5).
#' @param min_total_count Numeric. Minimum total count per gene for inclusion (default 0).
#' @param max_iterations Integer. Maximum outlier removal iterations (default 20).
#' @param verbose Logical. Print iteration progress (default TRUE).
#'
#' @return List containing:
#'   \item{pc_plot}{ggplot object for PCA results}
#'   \item{mds_plot}{ggplot object for MDS results}
#'   \item{filter_table}{data.frame with sample IDs and outlier status}
#'   \item{outliers}{character vector of outlier sample IDs}
#'   \item{n_iterations}{number of iterations performed}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- detect_outliers_pca_mds(
#'   data_matrix = counts,
#'   n_sd = 4,
#'   verbose = TRUE
#' )
#' print(result$pc_plot)
#' table(result$filter_table$filter_outliers)
#' }
detect_outliers_pca_mds <- function(
    data_matrix,
    n_sd = 5,
    min_gene_count = 10,
    min_gene_samples = 5,
    min_total_count = 0,
    max_iterations = 20,
    verbose = TRUE
) {
  
  # Input validation
  .validate_outlier_inputs(data_matrix, n_sd, min_gene_count, 
                           min_gene_samples, max_iterations)
  
  # Initialize
  sample_names <- colnames(data_matrix)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample_", seq_len(ncol(data_matrix)))
    colnames(data_matrix) <- sample_names
  }
  
  # Transpose and filter genes
  mat <- t(data_matrix)
  mat <- .filter_low_genes(mat, min_total_count, min_gene_count, min_gene_samples)
  
  if (ncol(mat) == 0) {
    stop("No genes passed filtering criteria. Consider relaxing thresholds.")
  }
  
  # Store initial data for plotting
  initial_pc <- .calculate_pca(mat)
  initial_mds <- .calculate_mds(mat)
  
  # Iterative outlier detection
  outlier_result <- .iterative_outlier_detection(
    mat, n_sd, max_iterations, verbose
  )
  
  # Create filter table
  filter_table <- data.frame(
    sample_id = sample_names,
    filter_outliers = ifelse(
      sample_names %in% outlier_result$outliers,
      "Outlier",
      "Passed"
    ),
    stringsAsFactors = FALSE
  )
  
  # Mark outliers in initial coordinates
  initial_pc$outliers <- rownames(initial_pc) %in% outlier_result$outliers
  initial_mds$outliers <- rownames(initial_mds) %in% outlier_result$outliers
  
  # Create plots
  pc_plot <- .create_outlier_plot(
    initial_pc, "PC1", "PC2",
    title = sprintf("PCA Outlier Detection (threshold = %g SD)", n_sd)
  )
  
  mds_plot <- .create_outlier_plot(
    initial_mds, "MDS1", "MDS2",
    title = sprintf("MDS Outlier Detection (threshold = %g SD)", n_sd)
  )
  
  if (verbose) {
    message(sprintf(
      "Outlier detection complete: %d outliers found in %d iterations",
      length(outlier_result$outliers),
      outlier_result$n_iterations
    ))
  }
  
  list(
    pc_plot = pc_plot,
    mds_plot = mds_plot,
    filter_table = filter_table,
    outliers = outlier_result$outliers,
    n_iterations = outlier_result$n_iterations
  )
}

# ---- Helper Functions ----

#' Validate inputs for outlier detection
#' @keywords internal
.validate_outlier_inputs <- function(data_matrix, n_sd, min_gene_count, 
                                     min_gene_samples, max_iterations) {
  if (!is.matrix(data_matrix) && !is.data.frame(data_matrix)) {
    stop("data_matrix must be a matrix or data frame")
  }
  if (ncol(data_matrix) < 3) {
    stop("Need at least 3 samples for outlier detection")
  }
  if (!is.numeric(n_sd) || n_sd <= 0) {
    stop("n_sd must be a positive number")
  }
  if (!is.numeric(min_gene_count) || min_gene_count < 0) {
    stop("min_gene_count must be non-negative")
  }
  if (!is.numeric(min_gene_samples) || min_gene_samples < 0) {
    stop("min_gene_samples must be non-negative")
  }
  if (!is.numeric(max_iterations) || max_iterations < 1) {
    stop("max_iterations must be at least 1")
  }
}

#' Filter genes with low counts
#' @keywords internal
.filter_low_genes <- function(mat, min_total, min_count, min_samples) {
  # Remove genes with zero total counts
  keep_total <- colSums(mat) > min_total
  
  # Remove genes not expressed above threshold in enough samples
  keep_samples <- colSums(mat >= min_count) >= min_samples
  
  mat[, keep_total & keep_samples, drop = FALSE]
}

#' Calculate PCA coordinates
#' @keywords internal
.calculate_pca <- function(mat) {
  pca_result <- stats::prcomp(mat, scale. = TRUE)
  pc_values <- as.data.frame(pca_result$x[, 1:2])
  colnames(pc_values) <- c("PC1", "PC2")
  pc_values
}

#' Calculate MDS coordinates
#' @keywords internal
.calculate_mds <- function(mat) {
  dist_mat <- stats::dist(mat)
  mds_result <- stats::cmdscale(dist_mat, k = 2)
  mds_values <- as.data.frame(mds_result)
  colnames(mds_values) <- c("MDS1", "MDS2")
  mds_values
}

#' Identify outliers in 2D coordinates
#' @keywords internal
.identify_outliers_2d <- function(coords, n_sd) {
  dim1_mean <- mean(coords[, 1])
  dim1_sd <- stats::sd(coords[, 1])
  dim2_mean <- mean(coords[, 2])
  dim2_sd <- stats::sd(coords[, 2])
  
  outliers_dim1 <- abs(coords[, 1] - dim1_mean) > n_sd * dim1_sd
  outliers_dim2 <- abs(coords[, 2] - dim2_mean) > n_sd * dim2_sd
  
  rownames(coords)[outliers_dim1 | outliers_dim2]
}

#' Iteratively detect and remove outliers
#' @keywords internal
.iterative_outlier_detection <- function(mat, n_sd, max_iterations, verbose) {
  all_outliers <- character(0)
  iteration <- 0
  
  while (iteration < max_iterations) {
    iteration <- iteration + 1
    
    # Calculate coordinates
    pc_coords <- .calculate_pca(mat)
    mds_coords <- .calculate_mds(mat)
    
    # Find outliers
    pc_outliers <- .identify_outliers_2d(pc_coords, n_sd)
    mds_outliers <- .identify_outliers_2d(mds_coords, n_sd)
    new_outliers <- unique(c(pc_outliers, mds_outliers))
    
    if (length(new_outliers) == 0) {
      if (verbose) message(sprintf("  Iteration %d: no new outliers found", iteration))
      break
    }
    
    if (verbose) {
      message(sprintf(
        "  Iteration %d: found %d outliers (%d from PCA, %d from MDS)",
        iteration, length(new_outliers), length(pc_outliers), length(mds_outliers)
      ))
    }
    
    # Add to cumulative list
    all_outliers <- c(all_outliers, new_outliers)
    
    # Remove outliers and re-filter genes
    mat <- mat[!rownames(mat) %in% new_outliers, , drop = FALSE]
    mat <- .filter_low_genes(mat, 0, 10, 5)
    
    if (nrow(mat) < 3) {
      warning("Fewer than 3 samples remain after outlier removal")
      break
    }
  }
  
  if (iteration >= max_iterations) {
    warning(sprintf(
      "Reached maximum iterations (%d). Consider increasing max_iterations or adjusting n_sd.",
      max_iterations
    ))
  }
  
  list(
    outliers = unique(all_outliers),
    n_iterations = iteration
  )
}

#' Create outlier plot
#' @keywords internal
.create_outlier_plot <- function(data, x_var, y_var, title) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; returning NULL for plot")
    return(NULL)
  }
  
  ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_var]], y = .data[[y_var]], 
                                     color = .data$outliers)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.7) +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "skyblue4", "TRUE" = "tomato3"),
      labels = c("FALSE" = "Passed", "TRUE" = "Outlier"),
      name = "Status"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = title, x = x_var, y = y_var) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}