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

test_that("ebayes location keeps more signal with more anchors", {
  set.seed(1)
  gamma <- cbind(rnorm(200, 0.8, 0.2), rnorm(200, -0.8, 0.2))
  n_few <- matrix(2, nrow = 200, ncol = 2)
  n_many <- matrix(30, nrow = 200, ncol = 2)
  g_few <- anchorCorrectR:::.ebayes_shrink_location(gamma, n_few, verbose = FALSE)
  g_many <- anchorCorrectR:::.ebayes_shrink_location(gamma, n_many, verbose = FALSE)
  # More anchors => less shrinkage toward 0
  expect_gt(mean(abs(g_many)), mean(abs(g_few)))
  # With many anchors, posterior should stay close to the raw estimate (~ n/(n+1))
  expect_gt(mean(abs(g_many) / pmax(abs(gamma), 1e-8)), 0.9)
  expect_lt(mean(abs(g_few) / pmax(abs(gamma), 1e-8)), 0.8)
})

test_that("combat reduces batch mean differences with many anchors", {
  set.seed(2)
  sim <- helper_simulate_batches(n_genes = 120, n_per_batch = 20, n_batches = 2)
  tf <- anchorCorrectR:::make_transform_cpm_log1p(sim$counts, input_type = "counts")
  before <- tf$forward(sim$counts)
  out <- correct_combat_anchor(
    sim$counts, sim$batch, sim$sample_id,
    input_type = "counts", ref_batch = "B1",
    verbose = FALSE, round_counts = FALSE
  )
  after <- tf$forward(out)
  batch_gap <- function(X, batch) {
    lev <- levels(batch)
    mean(abs(
      rowMeans(X[, batch == lev[1], drop = FALSE]) -
        rowMeans(X[, batch == lev[2], drop = FALSE])
    ))
  }
  expect_lt(batch_gap(after, sim$batch), 0.25 * batch_gap(before, sim$batch))
})
