#' Compare RNA and DNA VCF Files for Sample Agreement
#'
#' Compares variants between RNA-seq and DNA/genotype VCF files to assess
#' sample concordance. Computes shared and unique variant counts for each
#' RNA-DNA pair and classifies samples based on agreement scores.
#'
#' @param rna_dir Path to folder containing RNA VCF files.
#' @param dna_dir Path to folder containing DNA/genotype VCF files.
#' @param find.endpattern Logical, whether to filter files by regex
#' @param endpattern_file Regex pattern for VCF filenames (default ".vcf\\.gz$")
#' @param n.cores Number of cores for parallel computation (default 1)
#' @param id_pattern_substring Logical, use substring of filename as ID
#' @param id_length Length of substring for sample ID
#' @param correct_refalt Logical, correct REF/ALT swaps if detected
#' @param shared_score_threshold Numeric threshold used for classification (default 0.6)
#' @param verbose Logical, print progress messages
#'

#' @return List with:
#' \describe{
#'   \item{comparisons}{data.frame with all pairwise comparison statistics}
#'   \item{summary}{Summary table of sample classifications}
#'   \item{filter_table}{Per-sample classification results}
#'   \item{plot}{Combined ggplot visualization (bar + boxplot + pie)}
#'   \item{rna_matrix}{normalized genotype matrices from RNA samples}
#'   \item{dna_matrix}{normalized genotype matrices from DNA samples}
#' }
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' result <- vcf_compare_fast(
#'   rna_dir = rna_dir,
#'   dna_dir = dna_dir,
#'   find.endpattern = TRUE,
#'   endpattern_file = "\\.vcf\\.gz$",
#'   n.cores = 4,
#'   id_pattern_substring = TRUE,
#'   id_length = 7,
#'   correct_refalt = TRUE,
#'   shared_score_threshold = 0.6,
#'   verbose = TRUE
#' )
#' }
vcf_compare_fast <- function(
    rna_dir, dna_dir, find.endpattern = TRUE, endpattern_file = ".vcf\\.gz$",
    n.cores = 1, id_pattern_substring = TRUE, id_length = 7,
    correct_refalt = TRUE, shared_score_threshold = 0.6, verbose = TRUE
) {
  for (pkg in c("vcfR", "dplyr", "foreach", "doParallel", "parallel", "tidyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required for vcf_compare_fast().", pkg), call. = FALSE)
    }
  }
  `%dopar%` <- foreach::`%dopar%`
  `%:%` <- foreach::`%:%`
  `%>%` <- dplyr::`%>%`

  if (verbose) {
    message("Detecting RNA and DNA VCF files...")
    message("Classification shared_score_threshold = ", shared_score_threshold)
  }
  
  # --- Detect files ---
  rna_files <- list.files(rna_dir, pattern = if (find.endpattern) endpattern_file else NULL, full.names = TRUE)
  dna_files <- list.files(dna_dir, pattern = if (find.endpattern) endpattern_file else NULL, full.names = TRUE)
  if (length(rna_files) == 0 || length(dna_files) == 0)
    stop("No RNA or DNA VCF files found")
  
  if (verbose) message("Reading and filtering VCF files...")
  
  # --- Read VCFs ---
  rna_vcfs <- lapply(rna_files, .read_and_filter_vcf, filter_pass = TRUE)
  dna_vcfs <- lapply(dna_files, .read_and_filter_vcf, filter_pass = FALSE)
  
  # --- Assign IDs ---
  rna_names <- basename(rna_files)
  dna_names <- basename(dna_files)
  if (id_pattern_substring) {
    rna_ids <- substr(rna_names, 1, id_length)
    dna_ids <- substr(dna_names, 1, id_length)
  } else {
    rna_ids <- rna_names
    dna_ids <- dna_names
  }
  names(rna_vcfs) <- rna_ids
  names(dna_vcfs) <- dna_ids
  
  # --- Correct REF/ALT swaps if requested ---
  if (correct_refalt) {
    if (verbose) message("Correcting REF/ALT swaps for shared samples...")
    shared <- intersect(names(rna_vcfs), names(dna_vcfs))
    for (sid in shared)
      rna_vcfs[[sid]] <- .correct_refalt_pairs(rna_vcfs[[sid]], dna_vcfs[[sid]])
  }
  
  if (verbose) message("Creating genotype matrices...")
  
  rna_mat <- .normalize_genotypes(.make_geno_matrix(rna_vcfs, rna_names))
  dna_mat <- .normalize_genotypes(.make_geno_matrix(dna_vcfs, dna_names))
  
  # --- Restrict to shared variant sites ---
  shared_sites <- intersect(rownames(rna_mat), rownames(dna_mat))
  rna_mat <- rna_mat[shared_sites, , drop = FALSE]
  dna_mat <- dna_mat[shared_sites, , drop = FALSE]
  
  if (verbose) message("Running RNA-DNA sample comparisons...")
  
  # --- Parallel setup ---
  if (n.cores > 1) {
    if (.Platform$OS.type == "unix") {
      doParallel::registerDoParallel(cores = n.cores)
    } else {
      cl <- parallel::makeCluster(n.cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }
  }
  
  # --- Compute all pairwise RNA-DNA similarities ---
  comparisons <- foreach::foreach(i = seq_len(ncol(rna_mat)), .combine = rbind) %:%
    foreach::foreach(j = seq_len(ncol(dna_mat)), .combine = rbind) %dopar% {
      rna_gt <- rna_mat[, i]
      dna_gt <- dna_mat[, j]
      both <- !is.na(rna_gt) & !is.na(dna_gt)
      match_gt <- rna_gt[both] == dna_gt[both]
      
      data.frame(
        rna_sample_full = colnames(rna_mat)[i],
        dna_sample_full = colnames(dna_mat)[j],
        rna_id = names(rna_vcfs)[i],
        dna_id = names(dna_vcfs)[j],
        n_shared_sites = sum(both),
        n_equal = sum(match_gt),
        similarity = ifelse(sum(both) > 0, sum(match_gt) / sum(both), NA_real_)
      )
    }
  
  comparisons <- as.data.frame(comparisons, stringsAsFactors = FALSE) %>%
    dplyr::mutate(match_id = rna_id == dna_id)
  
  if (verbose) message("Classifying RNA samples by similarity...")
  
  classifications <- .classify_rna_samples(comparisons, shared_score_threshold, verbose)
  
  if (verbose) message("Creating visualization plots...")
  
  plot_obj <- .create_vcf_comparison_plot(comparisons, classifications, shared_score_threshold)
  
  if (verbose) message("RNA-DNA comparison complete!")
  
  list(
    comparisons = comparisons,
    summary = classifications$summary,
    filter_table = classifications$filter_table,
    plot = plot_obj,
    rna_matrix = rna_mat,
    dna_matrix = dna_mat
  )
}

# ------------------------------------------------------------------------------
# --- Helper functions ---
# ------------------------------------------------------------------------------

# --- Read and filter VCF ---
.read_and_filter_vcf <- function(path, filter_pass = FALSE) {
  v <- vcfR::read.vcfR(path, verbose = FALSE)
  fix <- v@fix
  gt <- vcfR::extract.gt(v, element = "GT")
  rm(v)
  
  # Keep only SNPs with single alleles (no indels, no multiallelic ALT)
  keep <- nchar(fix[, "REF"]) == 1 & nchar(fix[, "ALT"]) == 1 & !grepl(",", fix[, "ALT"])
  if (filter_pass) keep <- keep & (fix[, "FILTER"] == "PASS")
  
  fix <- fix[keep, , drop = FALSE]
  gt <- gt[keep, , drop = FALSE]
  cbind(fix, GT = gt)
}

# --- Correct REF/ALT swaps between RNA and DNA ---
.correct_refalt_pairs <- function(rna_vcf, dna_vcf) {
  rna_vcf <- as.data.frame(rna_vcf, stringsAsFactors = FALSE)
  dna_vcf <- as.data.frame(dna_vcf, stringsAsFactors = FALSE)
  
  merged <- merge(
    rna_vcf[, c("CHROM", "POS", "REF", "ALT")],
    dna_vcf[, c("CHROM", "POS", "REF", "ALT")],
    by = c("CHROM", "POS"),
    suffixes = c(".RNA", ".DNA")
  )
  
  inverted <- merged$REF.RNA == merged$ALT.DNA & merged$ALT.RNA == merged$REF.DNA
  if (any(inverted)) {
    pos_fix <- merged$POS[inverted]
    rows <- which(rna_vcf$POS %in% pos_fix)
    
    # Swap REF/ALT for inverted sites
    tmp <- rna_vcf[rows, "REF"]
    rna_vcf[rows, "REF"] <- rna_vcf[rows, "ALT"]
    rna_vcf[rows, "ALT"] <- tmp
    
    # Flip genotype encoding (0<->1)
    gt_cols <- setdiff(colnames(rna_vcf), c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"))
    if (length(gt_cols) > 0) {
      for (col in gt_cols) {
        gt <- rna_vcf[rows, col]
        gt <- gsub("\\|", "/", gt)
        gt <- gsub("0", "x", gt)
        gt <- gsub("1", "0", gt)
        gt <- gsub("x", "1", gt)
        rna_vcf[rows, col] <- gt
      }
    }
  }
  rna_vcf
}

# --- Build genotype matrix from list of VCFs ---
.make_geno_matrix <- function(vcf_list, sample_names = names(vcf_list)) {
  all_sites <- unique(do.call(c, lapply(vcf_list, function(x)
    paste(x[, "CHROM"], x[, "POS"], x[, "REF"], x[, "ALT"], sep = "_")
  )))
  
  mat <- matrix(NA, nrow = length(all_sites), ncol = length(vcf_list),
                dimnames = list(all_sites, sample_names))
  
  for (i in seq_along(vcf_list)) {
    x <- as.data.frame(vcf_list[[i]], stringsAsFactors = FALSE)
    ids <- paste(x[, "CHROM"], x[, "POS"], x[, "REF"], x[, "ALT"], sep = "_")
    gt_cols <- setdiff(colnames(x), c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO", "ID"))
    if (length(gt_cols) > 1) gt_cols <- gt_cols[1]
    gts <- x[, gt_cols]
    idx <- match(ids, all_sites)
    mat[idx[!is.na(idx)], i] <- gts[!is.na(idx)]
  }
  mat
}

# --- Normalize genotypes: remove phasing, sort alleles, unify format ---
.normalize_genotypes <- function(gt_matrix) {
  mat <- as.data.frame(gt_matrix, stringsAsFactors = FALSE)
  mat[] <- lapply(mat, function(x) gsub("\\|", "/", x)) # remove phasing
  mat[] <- lapply(mat, function(x)
    sapply(x, function(g) {
      if (is.na(g)) return(NA)
      alleles <- strsplit(g, "/")[[1]]
      if (length(alleles) != 2) return(g)
      paste0(sort(alleles), collapse = "/")
    })
  )
  mat <- as.matrix(mat)
  if (!is.null(rownames(gt_matrix))) rownames(mat) <- rownames(gt_matrix)
  if (!is.null(colnames(gt_matrix))) colnames(mat) <- colnames(gt_matrix)
  mat
}

.classify_rna_samples <- function(comparisons, threshold, verbose) {
  # Find best match for each RNA sample
  top_matches <- comparisons %>%
    dplyr::group_by(rna_sample_full) %>%
    dplyr::slice_max(order_by = similarity, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  
  # Identify categories
  perfect_high <- top_matches %>% dplyr::filter(rna_id == dna_id, similarity > threshold) %>% dplyr::pull(rna_sample_full)
  perfect_low  <- top_matches %>% dplyr::filter(rna_id == dna_id, similarity <= threshold) %>% dplyr::pull(rna_sample_full)
  
  multi_hits <- comparisons %>%
    dplyr::filter(similarity > threshold) %>%
    dplyr::group_by(rna_sample_full) %>%
    dplyr::summarise(n_hits = dplyr::n(), .groups = "drop") %>%
    dplyr::filter(n_hits > 1) %>%
    dplyr::pull(rna_sample_full)
  
  multi_hits <- setdiff(multi_hits, perfect_high)
  
  multi_with_match <- comparisons %>%
    dplyr::filter(rna_sample_full %in% multi_hits, similarity > threshold, rna_id == dna_id) %>%
    dplyr::pull(rna_sample_full) %>%
    unique()
  multi_no_match <- setdiff(multi_hits, multi_with_match)
  
  available_dna <- unique(comparisons$dna_id)
  no_dna_info <- setdiff(unique(comparisons$rna_id), available_dna)
  
  low_score <- comparisons %>%
    dplyr::group_by(rna_sample_full) %>%
    dplyr::summarise(max_score = max(similarity, na.rm = TRUE), .groups = "drop") %>%
    dplyr::filter(max_score < threshold) %>%
    dplyr::pull(rna_sample_full)
  low_score <- setdiff(low_score, no_dna_info)
  
  all_rna <- unique(comparisons$rna_sample_full)
  filter_table <- data.frame(rna_sample_full = all_rna, stringsAsFactors = FALSE)
  
  filter_table <- filter_table %>%
    dplyr::mutate(
      classification = dplyr::case_when(
        rna_sample_full %in% perfect_high ~ "Pass: High concordance",
        rna_sample_full %in% perfect_low ~ "Pass: Low concordance",
        rna_sample_full %in% multi_with_match ~ "Multiple matches (includes own ID)",
        rna_sample_full %in% multi_no_match ~ "Multiple matches (no own ID)",
        rna_sample_full %in% low_score ~ "Low concordance (no good match)",
        rna_sample_full %in% no_dna_info ~ "No matching DNA sample",
        TRUE ~ "Unclassified"
      )
    )
  
  summary_table <- filter_table %>%
    dplyr::count(classification) %>%
    dplyr::mutate(
      percent = round(100 * n / sum(n), 1),
      label = paste0(classification, " (n=", n, ", ", percent, "%)")
    )
  
  if (verbose) {
    message("\nClassification summary:")
    print(summary_table)
  }
  
  list(filter_table = filter_table, summary = summary_table)
}


# --- Visualization plot combining barplot, boxplot and pie ---
.create_vcf_comparison_plot <- function(comparisons, classifications, threshold) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE) ||
      !requireNamespace("dplyr", quietly = TRUE)) {
    warning("ggplot2, patchwork, and dplyr required for plotting; returning NULL")
    return(NULL)
  }

  `%>%` <- dplyr::`%>%`

  plot_data <- comparisons %>%
    dplyr::left_join(classifications$filter_table, by = "rna_sample_full") %>%
    dplyr::group_by(rna_sample_full) %>%
    dplyr::mutate(
      is_max = similarity == max(similarity, na.rm = TRUE),
      same_id = rna_id == dna_id
    ) %>%
    dplyr::ungroup()

  class_colors <- c(
    "Pass: High concordance" = "darkseagreen",
    "Pass: Low concordance" = "#8dd3c7",
    "Multiple matches (includes own ID)" = "darkseagreen1",
    "Multiple matches (no own ID)" = "lightsalmon3",
    "Low concordance (no good match)" = "mistyrose2",
    "No matching DNA sample" = "lightsalmon",
    "Unclassified" = "#cccccc"
  )

  bar_plot <- plot_data %>%
    dplyr::distinct(rna_sample_full, n_shared_sites, classification) %>%
    ggplot2::ggplot(ggplot2::aes(
      x = rna_sample_full, y = n_shared_sites, fill = classification
    )) +
    ggplot2::geom_col() +
    ggplot2::facet_grid(~ classification, scales = "free_x", space = "free") +
    ggplot2::scale_fill_manual(values = class_colors) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::labs(
      x = NULL, y = "Total RNA SNPs",
      title = "RNA-DNA Concordance Analysis"
    )

  box_plot <- plot_data %>%
    ggplot2::ggplot(ggplot2::aes(x = rna_sample_full, y = similarity)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_point(
      data = plot_data %>% dplyr::filter(is_max | same_id),
      ggplot2::aes(color = interaction(is_max, same_id)),
      size = 2, alpha = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = threshold, linetype = "dashed",
      color = "tomato3", linewidth = 0.5
    ) +
    ggplot2::facet_grid(~ classification, scales = "free_x", space = "free") +
    ggplot2::scale_color_manual(
      values = c(
        "FALSE.TRUE" = "seagreen3",
        "TRUE.FALSE" = "pink2",
        "TRUE.TRUE" = "seagreen3"
      ),
      labels = c(
        "FALSE.TRUE" = "Matching ID",
        "TRUE.FALSE" = "Best match (different ID)",
        "TRUE.TRUE" = "Best match (same ID)"
      ),
      name = "Point Type"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "RNA Sample", y = "Concordance Score")

  pie_plot <- classifications$summary %>%
    ggplot2::ggplot(ggplot2::aes(x = "", y = n, fill = classification)) +
    ggplot2::geom_col(width = 1, color = "white") +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::scale_fill_manual(values = class_colors) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::labs(title = "Sample Distribution", fill = "Classification")

  design <- "
    AAC
    BBC
    BBC
  "
  bar_plot + box_plot + pie_plot + patchwork::plot_layout(design = design)
}
