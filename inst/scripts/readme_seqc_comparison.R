#!/usr/bin/env Rscript
# Rebuild README SEQC metrics + MDS including sva::ComBat_seq.

suppressPackageStartupMessages({
  library(anchorCorrectR)
  library(sva)
  library(ggplot2)
})

out_dir <- file.path("man", "figures", "readme")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- as.matrix(read.delim(
  system.file("extdata", "counts.tsv.gz", package = "anchorCorrectR"),
  row.names = 1, check.names = FALSE
))
meta <- read.delim(
  system.file("extdata", "metadata.tsv.gz", package = "anchorCorrectR"),
  check.names = FALSE
)
stopifnot(identical(colnames(counts), meta$ID))

pick_ab <- function(meta_ab, sample_label, n_ilm, n_solid) {
  ilm <- meta_ab$ID[meta_ab$Sample == sample_label & meta_ab$Platform == "ILLUMINA"]
  sol <- meta_ab$ID[meta_ab$Sample == sample_label & meta_ab$Platform == "ABI_SOLID"]
  c(sample(ilm, n_ilm), sample(sol, n_solid))
}

# Each scenario re-runs from set.seed(123) like the README (same anchors; A/B draws
# start from the RNG state after sampling those anchors).
design_for <- function(scenario) {
  set.seed(123)
  cd <- meta[meta$Sample %in% c("C", "D"), ]
  cd_anchor_ids <- names(which(tapply(cd$Platform, cd$IDR, function(z) {
    all(c("ILLUMINA", "ABI_SOLID") %in% z)
  })))
  keep_idr <- sample(cd_anchor_ids, 5)
  meta_cd <- cd[cd$IDR %in% keep_idr, ]

  ab <- meta[meta$Sample %in% c("A", "B"), ]
  ids_ab <- switch(
    scenario,
    balanced   = c(pick_ab(ab, "A", 8, 7),  pick_ab(ab, "B", 8, 7)),
    mild       = c(pick_ab(ab, "A", 11, 4), pick_ab(ab, "B", 11, 4)),
    strong     = c(pick_ab(ab, "A", 14, 1), pick_ab(ab, "B", 14, 1)),
    confounded = c(pick_ab(ab, "A", 15, 0), pick_ab(ab, "B", 0, 15)),
    stop("Unknown scenario")
  )
  meta_ab <- ab[ab$ID %in% ids_ab, ]
  stopifnot(nrow(meta_ab) == 30)
  meta_use <- rbind(meta_cd, meta_ab)
  meta_use <- meta_use[match(intersect(meta_use$ID, colnames(counts)), meta_use$ID), ]
  list(keep_idr = keep_idr, meta_use = meta_use)
}

scenarios <- c("balanced", "mild", "strong", "confounded")

extract_row <- function(a, scenario, method) {
  data.frame(
    batch_var_before = a$pvca$before$batch_variance,
    batch_var_after  = a$pvca$after$batch_variance,
    bio_var_before   = a$pvca$before$biology_variance,
    bio_var_after    = a$pvca$after$biology_variance,
    knn_jaccard      = a$knn_jaccard,
    hvg_overlap      = a$hvg_overlap,
    rmse_before      = a$rmse_replicates$before,
    rmse_after       = a$rmse_replicates$after,
    ari              = a$adjusted_rand_index,
    scenario         = scenario,
    method           = method,
    stringsAsFactors = FALSE
  )
}

run_one <- function(scenario) {
  message("=== scenario: ", scenario, " ===")
  des <- design_for(scenario)
  keep_idr <<- des$keep_idr
  message("  Anchors: ", paste(des$keep_idr, collapse = ", "))
  meta_use <- des$meta_use
  counts_use <- counts[, meta_use$ID, drop = FALSE]
  batch     <- factor(meta_use$Platform)
  biology   <- factor(meta_use$Sample)
  sample_id <- factor(meta_use$IDR)

  outs <- list()
  outs$shift <- anchor_correct(
    counts_use, batch, sample_id, method = "shift", center = "mean",
    input_type = "counts", ref_batch = "ILLUMINA", verbose = FALSE
  )
  outs$ridge <- anchor_correct(
    counts_use, batch, sample_id, method = "ridge",
    input_type = "counts", ref_batch = "ILLUMINA", verbose = FALSE
  )
  outs$combat <- anchor_correct(
    counts_use, batch, sample_id, method = "combat",
    input_type = "counts", ref_batch = "ILLUMINA", verbose = FALSE
  )

  message("  running ComBat_seq (sva) with group = biology ...")
  outs$combatseq <- sva::ComBat_seq(
    counts = counts_use,
    batch = as.character(batch),
    group = as.character(biology)
  )
  storage.mode(outs$combatseq) <- "double"
  dimnames(outs$combatseq) <- dimnames(counts_use)

  assessments <- lapply(names(outs), function(m) {
    message("  assess: ", m)
    assess_correction(
      counts_use, outs[[m]], batch, biology, sample_id,
      input_type = "counts", verbose = FALSE
    )
  })
  names(assessments) <- names(outs)

  merged <- merge_assessments(
    assessments$shift, assessments$ridge, assessments$combat, assessments$combatseq,
    run_names = c("shift", "ridge", "combat", "combatseq")
  )
  write.csv(
    merged$consolidated_table,
    file.path(out_dir, paste0("consolidated_", scenario, ".csv")),
    row.names = FALSE
  )

  # MDS for all methods in balanced; shift+combatseq in imbalance; confounded key methods
  plot_methods <- switch(
    scenario,
    balanced = c("shift", "ridge", "combat", "combatseq"),
    mild = c("shift", "combatseq"),
    strong = c("shift", "combatseq"),
    confounded = c("shift", "ridge", "combatseq"),
    names(outs)
  )
  for (m in plot_methods) {
    p <- plot_mds_before_after(
      counts_use, outs[[m]], batch, biology, input_type = "counts"
    )
    ggsave(
      filename = file.path(out_dir, sprintf("mds_%s_%s.png", scenario, m)),
      plot = p, width = 9, height = 4.2, dpi = 150
    )
  }

  rows <- do.call(rbind, lapply(names(assessments), function(m) {
    extract_row(assessments[[m]], scenario, m)
  }))
  list(rows = rows, assessments = assessments, outs = outs, merged = merged)
}

results <- lapply(scenarios, run_one)
names(results) <- scenarios

summary_df <- do.call(rbind, lapply(results, `[[`, "rows"))
summary_df$batch_delta <- summary_df$batch_var_after - summary_df$batch_var_before
summary_df$bio_delta   <- summary_df$bio_var_after - summary_df$bio_var_before
summary_df$rmse_delta  <- summary_df$rmse_after - summary_df$rmse_before
summary_df$scenario <- factor(summary_df$scenario, levels = scenarios)
summary_df$method <- factor(
  summary_df$method,
  levels = c("shift", "ridge", "combat", "combatseq")
)

write.csv(summary_df, file.path(out_dir, "summary_metrics.csv"), row.names = FALSE)
saveRDS(
  list(anchors = keep_idr, summary = summary_df),
  file.path(out_dir, "analysis_results.rds")
)

# Print markdown-ready table rows
fmt <- function(x) sprintf("%.3f", x)
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf(
    "| %s | %s | %s → **%s** | %s → %s | %s → **%s** | %s | %s |\n",
    r$scenario, r$method,
    fmt(r$batch_var_before), fmt(r$batch_var_after),
    fmt(r$bio_var_before), fmt(r$bio_var_after),
    fmt(r$rmse_before), fmt(r$rmse_after),
    fmt(r$knn_jaccard), fmt(r$hvg_overlap)
  ))
}

message("Done. Wrote metrics and MDS under ", out_dir)
