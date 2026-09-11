#' Detect Sex Mismatch and Contamination in RNA-seq Samples
#'
#' Evaluates potential sex mismatches and contamination based on XIST and
#' Y chromosome gene expression. Uses a contamination zone defined by an
#' angular range around a diagonal line.
#'
#' @param counts Numeric matrix with genes in rows, samples in columns.
#'   Row names should be gene IDs (ENSEMBL format supported).
#' @param metadata Data frame with columns:
#'   \itemize{
#'     \item \code{sample_id}: Sample identifiers matching colnames(counts)
#'     \item \code{sex}: Biological sex ("Male" or "Female", case-insensitive)
#'   }
#' @param xist_gene Character. Gene ID for XIST (default: "ENSG00000229807").
#' @param chry_genes Character vector of Y chromosome gene IDs.
#'   If NULL, must be provided.
#' @param contamination_slope Either "estimated" (default) or numeric angle in degrees.
#'   Defines the center of the contamination zone.
#' @param contamination_fraction Numeric (0-1). Fraction of 90 degrees defining
#'   the contamination zone width (default 0.3 = 27 degrees).
#' @param shared_score_threshold Numeric (0-1). Minimum score for classification
#'   (default 0.6).
#' @param plot_title Character. Title for the plot.
#' @param id_col Character. Name of sample ID column in metadata (default "sample_id").
#' @param sex_col Character. Name of sex column in metadata (default "sex").
#'
#' @return List with:
#' \describe{
#'   \item{plot}{ggplot object showing contamination/mismatch analysis}
#'   \item{status_table}{data.frame with sample status classifications}
#'   \item{cpm_matrix}{CPM-normalized expression matrix}
#'   \item{summary}{Summary statistics of classifications}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- detect_sex_mismatch(
#'   counts = counts_matrix,
#'   metadata = metadata_df,
#'   chry_genes = chry_gene_ids,
#'   plot_title = "Sample QC"
#' )
#' print(result$plot)
#' table(result$status_table$status)
#' }
detect_sex_mismatch <- function(
    counts,
    metadata,
    xist_gene = "ENSG00000229807",
    chry_genes = default_chrY_genes,
    contamination_slope = "estimated",
    contamination_fraction = 0.3,
    shared_score_threshold = 0.6,
    plot_title = "Sex Mismatch and Contamination Analysis",
    id_col = "sample_id",
    sex_col = "sex"
) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required for detect_sex_mismatch().", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required for detect_sex_mismatch().", call. = FALSE)
  }
  `%>%` <- dplyr::`%>%`

  # Input validation
  .validate_sex_mismatch_inputs(
    counts, metadata, xist_gene, chry_genes, 
    contamination_fraction, id_col, sex_col
  )
  
  # Prepare data
  counts <- .strip_gene_versions(counts)
  counts_cpm <- edgeR::cpm(round(counts))
  rownames(counts_cpm) <- rownames(counts)
  
  # Extract expression values
  expr_data <- .extract_sex_chromosome_expression(
    counts_cpm, metadata, xist_gene, chry_genes, id_col, sex_col
  )
  
  # Calculate contamination zone parameters
  zone_params <- .calculate_contamination_zone(
    expr_data, contamination_slope, contamination_fraction
  )
  
  # Classify samples (now returns anchors too)
  classification_result <- .classify_samples(
    expr_data, zone_params, shared_score_threshold
  )
  
  # Create visualization with anchor points
  plot_obj <- .create_sex_mismatch_plot(
    expr_data, classification_result$classifications, 
    zone_params, classification_result$anchors, plot_title
  )
  
  # Summary statistics
  summary_stats <- classification_result$classifications %>%
    dplyr::count(.data$status) %>%
    dplyr::mutate(
      percent = round(100 * .data$n / sum(.data$n), 1),
      label = paste0(.data$status, " (n=", .data$n, ", ", .data$percent, "%)")
    )
  
  list(
    plot = plot_obj,
    status_table = classification_result$classifications,
    cpm_matrix = counts_cpm,
    summary = summary_stats
  )
}

# ---- Helper Functions ----

#' Validate inputs for sex mismatch detection
#' @keywords internal
.validate_sex_mismatch_inputs <- function(counts, metadata, xist_gene, 
                                          chry_genes, contamination_fraction,
                                          id_col, sex_col) {
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("counts must be a matrix or data frame")
  }
  
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }
  
  if (!id_col %in% colnames(metadata)) {
    stop(sprintf("Column '%s' not found in metadata", id_col))
  }
  
  if (!sex_col %in% colnames(metadata)) {
    stop(sprintf("Column '%s' not found in metadata", sex_col))
  }
  
  if (is.null(chry_genes) || length(chry_genes) == 0) {
    stop("chry_genes must be provided as a character vector of gene IDs")
  }
  
  if (!is.numeric(contamination_fraction) || 
      contamination_fraction <= 0 || contamination_fraction > 1) {
    stop("contamination_fraction must be between 0 and 1")
  }
  
  if (!xist_gene %in% rownames(counts)) {
    warning(sprintf(
      "XIST gene '%s' not found in count matrix. Check gene ID format.",
      xist_gene
    ))
  }
  
  missing_chry <- sum(!chry_genes %in% rownames(counts))
  if (missing_chry > 0) {
    warning(sprintf(
      "%d of %d Y chromosome genes not found in count matrix",
      missing_chry, length(chry_genes)
    ))
  }
}

#' Strip gene version numbers from Ensembl IDs
#' @keywords internal
.strip_gene_versions <- function(counts) {
  if (any(grepl("\\.", rownames(counts)))) {
    rownames(counts) <- sub("\\..*$", "", rownames(counts))
  }
  counts
}

#' Extract and normalize XIST and Y chromosome expression
#' @keywords internal
.extract_sex_chromosome_expression <- function(counts_cpm, metadata, 
                                               xist_gene, chry_genes,
                                               id_col, sex_col) {
  `%>%` <- dplyr::`%>%`

  # Extract expression
  xist_expr <- if (xist_gene %in% rownames(counts_cpm)) {
    counts_cpm[xist_gene, ]
  } else {
    rep(0, ncol(counts_cpm))
  }
  
  chry_present <- chry_genes[chry_genes %in% rownames(counts_cpm)]
  chry_expr <- if (length(chry_present) > 0) {
    colSums(counts_cpm[chry_present, , drop = FALSE])
  } else {
    rep(0, ncol(counts_cpm))
  }
  
  # Prepare metadata with expression values
  result <- metadata %>%
    dplyr::mutate(
      xist_raw = xist_expr[match(.data[[id_col]], names(xist_expr))],
      chry_raw = chry_expr[match(.data[[id_col]], names(chry_expr))],
      sex = tolower(as.character(.data[[sex_col]]))
    ) %>%
    dplyr::filter(!is.na(.data$xist_raw), !is.na(.data$chry_raw))
  
  # Normalize to [0, 1]
  max_xist <- max(result$xist_raw, na.rm = TRUE)
  max_chry <- max(result$chry_raw, na.rm = TRUE)
  
  if (max_xist == 0) max_xist <- 1
  if (max_chry == 0) max_chry <- 1
  
  result %>%
    dplyr::mutate(
      xist_norm = .data$xist_raw / max_xist,
      chry_norm = .data$chry_raw / max_chry,
      sex_numeric = dplyr::case_when(
        .data$sex == "male" ~ 1L,
        .data$sex == "female" ~ 2L,
        TRUE ~ 0L
      )
    )
}

#' Calculate contamination zone parameters
#' @keywords internal
.calculate_contamination_zone <- function(data, contamination_slope, 
                                         contamination_fraction) {
  `%>%` <- dplyr::`%>%`

  contamination_width <- contamination_fraction * 90
  
  if (contamination_slope == "estimated") {
    # Estimate from median expression in each sex
    male_data <- data %>% dplyr::filter(.data$sex_numeric == 1)
    female_data <- data %>% dplyr::filter(.data$sex_numeric == 2)
    
    x_median <- median(male_data$xist_norm, na.rm = TRUE)
    y_median <- median(female_data$chry_norm, na.rm = TRUE)
    
    # Calculate angle
    hypotenuse <- sqrt(x_median^2 + y_median^2)
    angle_rad <- asin(y_median / hypotenuse)
    center_angle <- atan(angle_rad) * (180 / pi)
  } else {
    center_angle <- as.numeric(contamination_slope)
    if (is.na(center_angle)) {
      stop("contamination_slope must be 'estimated' or a numeric value in degrees")
    }
  }
  
  # Calculate slope values
  lower_angle <- center_angle - contamination_width / 2
  upper_angle <- center_angle + contamination_width / 2
  
  list(
    center_slope = tan(center_angle * pi / 180),
    lower_slope = tan(lower_angle * pi / 180),
    upper_slope = tan(upper_angle * pi / 180),
    center_angle = center_angle,
    width = contamination_width
  )
}

#' Classify samples based on expression and contamination zone
#' @keywords internal
.classify_samples <- function(data, zone_params, threshold) {
  `%>%` <- dplyr::`%>%`

  # Naive sex classification
  data <- data %>%
    dplyr::mutate(
      expr_sex_naive = dplyr::if_else(
        .data$chry_norm > .data$xist_norm * zone_params$center_slope,
        1L, 2L
      )
    )
  
  # Calculate anchor points for contamination zone
  male_classified <- data %>% 
    dplyr::filter(.data$sex_numeric == 1, .data$expr_sex_naive == 1)
  female_classified <- data %>% 
    dplyr::filter(.data$sex_numeric == 2, .data$expr_sex_naive == 2)
  
  x_anchor <- min(
    median(male_classified$xist_norm, na.rm = TRUE),
    0
    #min(female_classified$xist_norm, na.rm = TRUE) - 0.01
  )
  y_anchor <- min(
    median(male_classified$chry_norm, na.rm = TRUE),
    0
    #min(male_classified$chry_norm, na.rm = TRUE) - 0.01
  )
  
  # Classify contamination and sex
  classifications <- data %>%
    dplyr::mutate(
      xist_corrected = .data$xist_norm - x_anchor,
      lower_bound = .data$xist_corrected * zone_params$lower_slope + y_anchor,
      upper_bound = .data$xist_corrected * zone_params$upper_slope + y_anchor,
      middle_line = .data$xist_corrected * zone_params$center_slope + y_anchor,
      contaminated = dplyr::if_else(
        .data$chry_norm > .data$lower_bound & 
          .data$chry_norm < .data$upper_bound,
        "yes", "no"
      ),
      expr_sex = dplyr::if_else(
        .data$chry_norm < .data$middle_line,
        2L, 1L
      ),
      mismatch = dplyr::case_when(
        .data$sex_numeric == 0 ~ "unknown",
        .data$expr_sex == .data$sex_numeric ~ "no",
        TRUE ~ "yes"
      ),
      status = dplyr::case_when(
        .data$contaminated == "yes" & .data$mismatch == "yes" ~ 
          "Contaminated and sex mismatch",
        .data$contaminated == "yes" ~ "Likely contaminated",
        .data$mismatch == "yes" ~ "Sex mismatch",
        TRUE ~ "Passed"
      )
    )
  
  # Return both classifications and anchors
  list(
    classifications = classifications,
    anchors = list(x_anchor = x_anchor, y_anchor = y_anchor)
  )
}

#' Create sex mismatch visualization
#' @keywords internal
.create_sex_mismatch_plot <- function(data, classifications, zone_params, 
                                      anchors, title) {
  `%>%` <- dplyr::`%>%`

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 not available; returning NULL for plot")
    return(NULL)
  }
  
  # Prepare exclusion zone data with anchor corrections
  max_norm <- max(c(data$xist_norm, data$chry_norm), na.rm = TRUE)
  x_seq <- seq(0, max_norm, length.out = 100)
  
  # Calculate the intercepts for the anchor-corrected lines
  # Formula: y = (x - x_anchor) * slope + y_anchor
  # Expanded: y = x * slope + (y_anchor - x_anchor * slope)
  lower_intercept <- anchors$y_anchor - anchors$x_anchor * zone_params$lower_slope
  upper_intercept <- anchors$y_anchor - anchors$x_anchor * zone_params$upper_slope
  middle_intercept <- anchors$y_anchor - anchors$x_anchor * zone_params$center_slope
  
  zone_data <- data.frame(x = x_seq) %>%
    dplyr::mutate(
      lower_bound = .data$x * zone_params$lower_slope + lower_intercept,
      upper_bound = .data$x * zone_params$upper_slope + upper_intercept,
      middle_line = .data$x * zone_params$center_slope + middle_intercept
    )
  
  # Color palette
  status_colors <- c(
    "Passed" = "#1f77b4",
    "Likely contaminated" = "#ff7f0e",
    "Sex mismatch" = "#d62728",
    "Contaminated and sex mismatch" = "#9467bd"
  )
  
  # Create plot
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = zone_data,
      ggplot2::aes(x = .data$x, ymin = .data$lower_bound, ymax = .data$upper_bound),
      alpha = 0.2, fill = "orange"
    ) +
    ggplot2::geom_line(
      data = zone_data,
      ggplot2::aes(x = .data$x, y = .data$middle_line),
      linetype = "dashed", color = "blue", linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = classifications,
      ggplot2::aes(
        x = .data$xist_norm, 
        y = .data$chry_norm,
        color = .data$status,
        shape = factor(.data$sex_numeric)
      ),
      size = 2.5, alpha = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = status_colors,
      name = "Status"
    ) +
    ggplot2::scale_shape_manual(
      values = c("1" = 17, "2" = 16, "0" = 4),
      labels = c("1" = "Male", "2" = "Female", "0" = "Unknown"),
      name = "Reported Sex"
    ) +
    ggplot2::coord_cartesian(xlim = c(0, max_norm), ylim = c(0, max_norm)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = title,
      x = "XIST Normalized Expression",
      y = "Chr Y Genes Normalized Expression"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}