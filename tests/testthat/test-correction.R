helper_simulate_batches <- function(n_genes = 80, n_per_batch = 6, n_batches = 3,
                                    seed = 1) {
  set.seed(seed)
  sample_id <- factor(rep(seq_len(n_per_batch), times = n_batches))
  batch <- factor(rep(paste0("B", seq_len(n_batches)), each = n_per_batch))
  counts <- matrix(
    rpois(n_genes * length(sample_id), lambda = 60),
    nrow = n_genes,
    dimnames = list(
      paste0("Gene", seq_len(n_genes)),
      paste0("S", seq_along(sample_id))
    )
  )
  mult <- setNames(seq(1, 2, length.out = n_batches), levels(batch))
  counts <- sweep(counts, 2, mult[as.character(batch)], `*`)
  list(counts = counts, batch = batch, sample_id = sample_id)
}

test_that("correct_shift mean and median return same dimensions", {
  sim <- helper_simulate_batches()
  out_mean <- correct_shift(
    sim$counts, sim$batch, sim$sample_id,
    center = "mean", input_type = "counts", verbose = FALSE
  )
  out_med <- correct_shift(
    sim$counts, sim$batch, sim$sample_id,
    center = "median", input_type = "counts", verbose = FALSE
  )
  expect_equal(dim(out_mean), dim(sim$counts))
  expect_equal(dim(out_med), dim(sim$counts))
  expect_equal(rownames(out_mean), rownames(sim$counts))
  expect_equal(colnames(out_mean), colnames(sim$counts))
})

test_that("anchor_correct shift matches correct_shift", {
  sim <- helper_simulate_batches()
  a <- anchor_correct(
    sim$counts, sim$batch, sim$sample_id,
    method = "shift", center = "mean",
    input_type = "counts", verbose = FALSE, round_counts = FALSE
  )
  b <- correct_shift(
    sim$counts, sim$batch, sim$sample_id,
    center = "mean", input_type = "counts",
    verbose = FALSE, round_counts = FALSE
  )
  expect_equal(a, b)
})

test_that("detect_input_type classifies integer counts", {
  m <- matrix(c(0, 1, 2, 10, 3, 4), nrow = 2)
  expect_equal(anchorCorrectR:::detect_input_type(m), "counts")
})

test_that("ref_batch leaves reference columns unchanged on log scale", {
  sim <- helper_simulate_batches()
  tf <- anchorCorrectR:::make_transform_cpm_log1p(sim$counts, input_type = "counts")
  X <- tf$forward(sim$counts)
  # work on log input for exact equality
  out <- correct_shift(
    X, sim$batch, sim$sample_id,
    center = "mean", input_type = "log",
    ref_batch = "B1", verbose = FALSE
  )
  cols_ref <- which(sim$batch == "B1")
  expect_equal(out[, cols_ref], X[, cols_ref])
})

test_that("combat correction returns finite matrix", {
  sim <- helper_simulate_batches(n_genes = 40)
  out <- correct_combat_anchor(
    sim$counts, sim$batch, sim$sample_id,
    input_type = "counts", verbose = FALSE, round_counts = FALSE
  )
  expect_equal(dim(out), dim(sim$counts))
  expect_true(all(is.finite(out)))
})
