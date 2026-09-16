#!/usr/bin/env Rscript
# Regenerate README scenario MDS plots with overall + panel titles.

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

pick_ab <- function(meta_ab, sample_label, n_ilm, n_solid) {
  ilm <- meta_ab$ID[meta_ab$Sample == sample_label & meta_ab$Platform == "ILLUMINA"]
  sol <- meta_ab$ID[meta_ab$Sample == sample_label & meta_ab$Platform == "ABI_SOLID"]
  c(sample(ilm, n_ilm), sample(sol, n_solid))
}

method_label <- c(
  shift = "shift (mean)",
  ridge = "ridge",
  combat = "anchor-aware ComBat",
  combatseq = "ComBat_seq (sva)"
)

# Plots referenced in the cleaned README
wanted <- list(
  balanced   = c("shift", "ridge", "combat", "combatseq"),
  strong     = c("shift", "combatseq"),
  confounded = c("shift", "ridge", "combatseq")
)

for (scenario in names(wanted)) {
  message("=== ", scenario, " ===")
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
    strong     = c(pick_ab(ab, "A", 14, 1), pick_ab(ab, "B", 14, 1)),
    confounded = c(pick_ab(ab, "A", 15, 0), pick_ab(ab, "B", 0, 15)),
    stop(scenario)
  )
  meta_ab <- ab[ab$ID %in% ids_ab, ]
  meta_use <- rbind(meta_cd, meta_ab)
  meta_use <- meta_use[match(intersect(meta_use$ID, colnames(counts)), meta_use$ID), ]
  counts_use <- counts[, meta_use$ID, drop = FALSE]
  batch     <- factor(meta_use$Platform)
  biology   <- factor(meta_use$Sample)
  sample_id <- factor(meta_use$IDR)

  outs <- list()
  for (m in setdiff(wanted[[scenario]], "combatseq")) {
    args <- list(
      x = counts_use, batch = batch, sample_id = sample_id,
      method = m, input_type = "counts",
      ref_batch = "ILLUMINA", verbose = FALSE
    )
    if (m == "shift") args$center <- "mean"
    outs[[m]] <- do.call(anchor_correct, args)
  }
  if ("combatseq" %in% wanted[[scenario]]) {
    message("  ComBat_seq ...")
    outs$combatseq <- sva::ComBat_seq(
      counts = counts_use,
      batch = as.character(batch),
      group = as.character(biology)
    )
    storage.mode(outs$combatseq) <- "double"
    dimnames(outs$combatseq) <- dimnames(counts_use)
  }

  for (m in wanted[[scenario]]) {
    ttl <- sprintf("MDS before vs after — %s · %s", method_label[[m]], scenario)
    p <- plot_mds_before_after(
      counts_use, outs[[m]], batch, biology,
      input_type = "counts", title = ttl
    )
    f <- file.path(out_dir, sprintf("mds_%s_%s.png", scenario, m))
    ggsave(f, p, width = 9, height = 4.4, dpi = 150)
    message("  wrote ", basename(f))
  }
}

message("Done regenerating titled scenario MDS plots.")
