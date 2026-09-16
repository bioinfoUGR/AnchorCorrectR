#!/usr/bin/env Rscript
# Strong-imbalance SEQC: effect of number of C/D replicate IDRs.
#
# n = number of cross-platform IDRs (biological anchors), NOT libraries.
#   2  -> 1 C + 1 D   (4 libraries: each IDR on both platforms)
#   3  -> 2 C + 1 D
#   5  -> 3 C + 2 D
#   10 -> 5 C + 5 D
# Nested: larger sets contain smaller ones.

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

# Balanced nested C/D IDR pools
set.seed(123)
cd <- meta[meta$Sample %in% c("C", "D"), ]
cd_anchor_ids <- names(which(tapply(cd$Platform, cd$IDR, function(z) {
  all(c("ILLUMINA", "ABI_SOLID") %in% z)
})))
c_ids <- sample(cd_anchor_ids[grepl("^C", cd_anchor_ids)], 5)
d_ids <- sample(cd_anchor_ids[grepl("^D", cd_anchor_ids)], 5)

# Nested compositions: n -> (n_C, n_D)
comp <- list(
  `2`  = list(c = 1L, d = 1L),
  `3`  = list(c = 2L, d = 1L),
  `5`  = list(c = 3L, d = 2L),
  `10` = list(c = 5L, d = 5L)
)

anchor_sets <- lapply(comp, function(cd_n) {
  c(c_ids[seq_len(cd_n$c)], d_ids[seq_len(cd_n$d)])
})

message("C pool: ", paste(c_ids, collapse = ", "))
message("D pool: ", paste(d_ids, collapse = ", "))
for (nm in names(anchor_sets)) {
  a <- anchor_sets[[nm]]
  message(sprintf(
    "n=%s -> %s  (C=%d D=%d, libraries=%d)",
    nm, paste(a, collapse = ", "),
    sum(grepl("^C", a)), sum(grepl("^D", a)), 2L * length(a)
  ))
}

# Strong A/B imbalance fixed across all n
ab <- meta[meta$Sample %in% c("A", "B"), ]
ids_ab <- c(pick_ab(ab, "A", 14, 1), pick_ab(ab, "B", 14, 1))
meta_ab <- ab[ab$ID %in% ids_ab, ]
stopifnot(nrow(meta_ab) == 30)

methods <- c("shift", "ridge", "combat", "combatseq")
n_vec <- as.integer(names(anchor_sets))

extract_row <- function(a, n_rep, method, anchors) {
  data.frame(
    n_replicates = n_rep,
    n_C = sum(grepl("^C", anchors)),
    n_D = sum(grepl("^D", anchors)),
    n_libraries = 2L * length(anchors),
    method = method,
    anchors = paste(anchors, collapse = ","),
    batch_var_before = a$pvca$before$batch_variance,
    batch_var_after  = a$pvca$after$batch_variance,
    bio_var_before   = a$pvca$before$biology_variance,
    bio_var_after    = a$pvca$after$biology_variance,
    knn_jaccard      = a$knn_jaccard,
    hvg_overlap      = a$hvg_overlap,
    rmse_before      = a$rmse_replicates$before,
    rmse_after       = a$rmse_replicates$after,
    ari              = a$adjusted_rand_index,
    stringsAsFactors = FALSE
  )
}

# Helper: MDS long data for custom multi-panel titles
mds_long <- function(x_before, x_after, batch, biology, method_label) {
  mb <- mds_samples(x_before, input_type = "counts")
  ma <- mds_samples(x_after, input_type = "counts")
  mb$state <- "before"
  ma$state <- "after"
  m <- rbind(mb, ma)
  m$batch <- batch[match(m$sample, colnames(x_before))]
  m$biology <- biology[match(m$sample, colnames(x_before))]
  m$method <- method_label
  m$panel <- paste(method_label, m$state, sep = " — ")
  m
}

rows <- list()
outs_n2 <- list()
counts_n2 <- NULL
batch_n2 <- biology_n2 <- NULL

for (n_rep in n_vec) {
  keep_idr <- anchor_sets[[as.character(n_rep)]]
  message("=== n_replicates = ", n_rep,
          " (1 unit = 1 cross-platform IDR; composition ",
          sum(grepl("^C", keep_idr)), "C+",
          sum(grepl("^D", keep_idr)), "D) ===")

  meta_cd <- cd[cd$IDR %in% keep_idr, ]
  stopifnot(nrow(meta_cd) == 2L * n_rep)
  meta_use <- rbind(meta_cd, meta_ab)
  meta_use <- meta_use[match(intersect(meta_use$ID, colnames(counts)), meta_use$ID), ]
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
  message("  ComBat_seq ...")
  outs$combatseq <- sva::ComBat_seq(
    counts = counts_use,
    batch = as.character(batch),
    group = as.character(biology)
  )
  storage.mode(outs$combatseq) <- "double"
  dimnames(outs$combatseq) <- dimnames(counts_use)

  for (m in methods) {
    message("  assess: ", m)
    a <- assess_correction(
      counts_use, outs[[m]], batch, biology, sample_id,
      input_type = "counts", verbose = FALSE
    )
    rows[[length(rows) + 1L]] <- extract_row(a, n_rep, m, keep_idr)
  }

  # Individual method MDS with method-specific overall title
  title_map <- c(
    shift = "MDS before vs after — shift (mean)",
    ridge = "MDS before vs after — ridge",
    combat = "MDS before vs after — anchor-aware ComBat",
    combatseq = "MDS before vs after — ComBat_seq (sva)"
  )
  for (m in methods) {
    p <- plot_mds_before_after(
      counts_use, outs[[m]], batch, biology,
      input_type = "counts",
      title = sprintf("%s · %d replicates (%dC+%dD)",
                      title_map[[m]], n_rep,
                      sum(grepl("^C", keep_idr)), sum(grepl("^D", keep_idr)))
    )
    ggsave(
      filename = file.path(out_dir, sprintf("mds_strong_nrep%d_%s.png", n_rep, m)),
      plot = p, width = 9, height = 4.4, dpi = 150
    )
  }

  if (n_rep == 2L) {
    outs_n2 <- outs
    counts_n2 <- counts_use
    batch_n2 <- batch
    biology_n2 <- biology
  }
}

summary_df <- do.call(rbind, rows)
summary_df$method <- factor(summary_df$method, levels = methods)
summary_df$batch_delta <- summary_df$batch_var_after - summary_df$batch_var_before
summary_df$bio_delta   <- summary_df$bio_var_after - summary_df$bio_var_before
summary_df$rmse_delta  <- summary_df$rmse_after - summary_df$rmse_before

write.csv(
  summary_df,
  file.path(out_dir, "summary_metrics_nreplicates_strong.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    c_pool = c_ids, d_pool = d_ids, anchor_sets = anchor_sets,
    summary = summary_df
  ),
  file.path(out_dir, "analysis_nreplicates_strong.rds")
)

# Trend plot
long <- rbind(
  transform(summary_df, metric = "Batch variance (after)", value = batch_var_after),
  transform(summary_df, metric = "Biology variance (after)", value = bio_var_after),
  transform(summary_df, metric = "Anchor RMSE (after)", value = rmse_after)
)
long$metric <- factor(
  long$metric,
  levels = c("Batch variance (after)", "Biology variance (after)", "Anchor RMSE (after)")
)

p_trend <- ggplot(long, aes(x = n_replicates, y = value, color = method)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = n_vec) +
  labs(
    title = "Strong imbalance: effect of number of replicate IDRs",
    subtitle = "n = 2 is 1C+1D (not 2C+2D); each IDR is still on both platforms",
    x = "Number of replicate IDRs (C+D)",
    y = NULL,
    color = "Method"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "nreplicates_strong_trends.png"),
  plot = p_trend, width = 10, height = 4.0, dpi = 150
)

# Combined ComBat vs ComBat_seq at n=2 with overall + panel titles
stopifnot(length(outs_n2) > 0)
mds_cmp <- rbind(
  mds_long(counts_n2, outs_n2$combat, batch_n2, biology_n2, "Anchor-aware ComBat"),
  mds_long(counts_n2, outs_n2$combatseq, batch_n2, biology_n2, "ComBat_seq (sva)")
)
mds_cmp$method <- factor(
  mds_cmp$method,
  levels = c("Anchor-aware ComBat", "ComBat_seq (sva)")
)
mds_cmp$state <- factor(mds_cmp$state, levels = c("before", "after"))
mds_cmp$panel <- factor(
  paste(mds_cmp$method, mds_cmp$state, sep = " — "),
  levels = c(
    "Anchor-aware ComBat — before",
    "Anchor-aware ComBat — after",
    "ComBat_seq (sva) — before",
    "ComBat_seq (sva) — after"
  )
)

p_cmp <- ggplot(mds_cmp, aes(x = Dim1, y = Dim2, color = batch, shape = biology)) +
  geom_point(size = 2.2, alpha = 0.9) +
  facet_wrap(~ panel, nrow = 2) +
  labs(
    title = "Strong imbalance · 2 replicates (1C + 1D) · ComBat vs ComBat_seq",
    x = "MDS 1", y = "MDS 2", color = "Batch", shape = "Biology"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave(
  filename = file.path(out_dir, "mds_strong_nrep2_combat_vs_combatseq.png"),
  plot = p_cmp, width = 9.5, height = 7.2, dpi = 150
)

fmt <- function(x) sprintf("%.3f", x)
cat("\nMarkdown table rows:\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf(
    "| %d (%dC+%dD) | %s | %s → %s | %s → %s | %s → %s | %s | %s |\n",
    r$n_replicates, r$n_C, r$n_D, r$method,
    fmt(r$batch_var_before), fmt(r$batch_var_after),
    fmt(r$bio_var_before), fmt(r$bio_var_after),
    fmt(r$rmse_before), fmt(r$rmse_after),
    fmt(r$knn_jaccard), fmt(r$hvg_overlap)
  ))
}

message("Done.")
