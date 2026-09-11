# anchorCorrectR

<p align="center">
  <img src="man/figures/anchorCorrectR-logo.png" alt="anchorCorrectR logo" width="420"/>
</p>

**Anchor-based batch correction for gene expression data.**

`anchorCorrectR` removes technical batch effects when the same biological samples
(or known replicates) are measured in more than one batch. Those cross-batch
replicates act as **anchors**: the package estimates per-batch shifts from them
and applies the correction to all samples.

It supports three correction methods (location shift with mean or median,
ridge, and anchor-aware ComBat), works on counts or log-scale matrices, and
includes tools to assess correction quality and run common RNA-seq QC checks.

## Why anchors?

Standard batch-correction methods estimate batch effects from all samples in a
batch. When biology and batch are confounded, that can remove real biological
signal.

With anchors, batch effects are estimated by comparing **the same biological
sample** across batches. That isolates technical shift more cleanly and is
especially useful when batches are unbalanced or partially confounded with
biology.

## Installation

Install from GitHub:

```r
# install.packages("remotes")
remotes::install_github("gbarturen/anchorCorrectR")
```

Or from a local clone:

```r
install.packages("path/to/anchorCorrectR", repos = NULL, type = "source")
```

After Bioconductor acceptance:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("anchorCorrectR")
```

### Dependencies

**Imports:** `stats`, `utils`, `methods`

**Suggested (methods / QC / docs):** `glmnet` (ridge), `ggplot2`, `FNN`,
`cluster`, `edgeR`, `dplyr`, `tibble`, `tidyr`, `patchwork`, `scales`, `vcfR`,
`foreach`, `doParallel`, `parallel`, `ShortRead`, `Matrix`, `BiocStyle`,
`knitr`, `rmarkdown`, `testthat`

Install selected Bioconductor suggestions if needed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("ShortRead", "edgeR"))
```

## Quick start

```r
library(anchorCorrectR)

# Expression matrix: genes x samples (counts or log-scale)
# batch:      factor of technical batch labels
# sample_id:  factor of biological sample IDs
#             IDs present in >= 2 batches are used as anchors

corrected <- anchor_correct(
  x         = counts,
  batch     = batch,
  sample_id = sample_id,
  method    = "shift",       # default; also "ridge" or "combat"
  center    = "mean"         # or "median" for robust location shift
)
```

Output is returned in the **same scale** as the input (counts or log).

## Simulated example

```r
library(anchorCorrectR)

set.seed(123)
n_genes  <- 500
n_per_batch <- 12
n_batches   <- 3

# Shared biological sample IDs across batches (all are anchors here)
sample_id <- factor(rep(seq_len(n_per_batch), times = n_batches))
batch     <- factor(rep(paste0("Batch", seq_len(n_batches)), each = n_per_batch))

# Baseline counts + artificial batch multipliers
counts <- matrix(
  rpois(n_genes * length(sample_id), lambda = 80),
  nrow = n_genes,
  dimnames = list(paste0("Gene", seq_len(n_genes)),
                  paste0("S", seq_along(sample_id)))
)
batch_mult <- c(Batch1 = 1.0, Batch2 = 1.8, Batch3 = 0.6)
counts <- sweep(counts, 2, batch_mult[as.character(batch)], `*`)

# Mean-centered location shift (good default; preserves variance well)
corrected <- anchor_correct(
  counts, batch, sample_id,
  method = "shift",
  center = "mean",
  input_type = "counts"
)

# Median-centered location shift (more robust to outlier anchors)
corrected_med <- anchor_correct(
  counts, batch, sample_id,
  method = "shift",
  center = "median",
  input_type = "counts"
)

# Optional: leave one batch unchanged and correct relative to it
corrected_ref <- anchor_correct(
  counts, batch, sample_id,
  method = "shift",
  ref_batch = "Batch1",
  input_type = "counts"
)
```

## Correction methods

All methods process internally on the log scale (log1p-CPM if the input is
counts) and return data in the input scale.

| Method | Function | Description |
|--------|----------|-------------|
| `"shift"` (default) | `correct_shift()` | Per-gene, per-batch offsets from anchor deviations vs their cross-batch center (`mean` or `median`). Simple, fast, variance-preserving. |
| `"ridge"` | `fit_anchor_ridge()` + `apply_correction_ridge()` | Ridge regression (`glmnet`, α = 0) of expression on batch dummies, fitted on anchors. Requires the **glmnet** package. |
| `"combat"` | `correct_combat_anchor()` | Anchor-aware ComBat: location (+ scale) parameters from anchors with empirical Bayes shrinkage across genes. |

```r
# Location shift (recommended starting point)
out_mean <- anchor_correct(counts, batch, sample_id, method = "shift", center = "mean")
out_med  <- anchor_correct(counts, batch, sample_id, method = "shift", center = "median")

# Ridge (needs glmnet)
out_ridge <- anchor_correct(counts, batch, sample_id, method = "ridge")

# Anchor-aware ComBat
out_combat <- anchor_correct(counts, batch, sample_id, method = "combat")
```

You can also call the method-specific functions directly:

```r
out_mean <- correct_shift(counts, batch, sample_id, center = "mean", input_type = "counts")
out_med  <- correct_shift(counts, batch, sample_id, center = "median", input_type = "counts")
out_cb   <- correct_combat_anchor(counts, batch, sample_id, input_type = "counts")
```

### Mean vs median (`center`)

For `method = "shift"` (or `correct_shift()`), `center` controls the location
statistic used to:

1. Summarize each anchor across batches  
2. Aggregate within-batch replicate contributions  
3. Aggregate offsets across anchors in a batch  
4. Recenter offsets across batches when `ref_batch` is `NULL`

| `center` | Behavior |
|----------|----------|
| `"mean"` (default) | Classical mean-shift; sum-to-zero recentering |
| `"median"` | Robust to outlier anchors; median-centered recentering |

### Choosing a method

- **shift** — default for most use cases; location-only correction that
  tends to preserve biological variance. Prefer `center = "median"` if a few
  anchors may be unreliable.
- **ridge** — useful when you want regularized batch estimates and can install
  `glmnet`.
- **combat** — closer to classical ComBat (location and scale, EB shrinkage),
  but parameters are estimated from anchors rather than all samples.

### Reference batch

Set `ref_batch` so that one batch is left uncorrected and others are adjusted
relative to it:

```r
corrected <- anchor_correct(
  counts, batch, sample_id,
  method = "shift",
  ref_batch = "Batch1"
)
```

If `ref_batch` has no anchors, corrections may be poorly estimated — a warning
is issued.

## Input types (counts vs log)

By default, `input_type = "auto"` tries to detect whether `x` is raw counts or
already on a log scale.

| Value | When to use |
|-------|-------------|
| `"auto"` | Default heuristic (integer-like non-negative → counts). |
| `"counts"` | Raw or count-like matrices (including decimal expected counts). |
| `"log"` | Already on the scale you want to correct (e.g. log1p-CPM). |

**Important:** Many RNA-seq pipelines produce **non-integer** count-like values
(RSEM/Salmon expected counts, aggregated replicates, etc.). Auto-detection may
label those as `"log"` and skip CPM + log1p. For count-like data with decimals,
set the type explicitly:

```r
corrected <- anchor_correct(
  x, batch, sample_id,
  input_type = "counts"
)
```

## Assessing correction quality

Compare before/after matrices with biology labels (optional) and sample IDs:

```r
# biology: factor/character/numeric biological variable (optional)
assessment <- assess_correction(
  x_before  = counts,
  x_after   = corrected,
  batch     = batch,
  biology   = biology,      # optional
  sample_id = sample_id,    # enables replicate RMSE
  input_type = "counts"
)

print(assessment)
```

Metrics include:

- **PVCA** — variance attributed to batch vs biology (before/after)
- **kNN Jaccard** — neighborhood stability before vs after
- **HVG overlap** — highly variable gene concordance
- **Biology concordance** — how well biological structure is retained
- **Adjusted Rand Index** — clustering structure before vs after
- **Replicate RMSE** — cross-batch distance of anchors (should decrease)

Compare several methods:

```r
a_mean   <- assess_correction(counts, out_mean,   batch, biology, sample_id)
a_med    <- assess_correction(counts, out_med,    batch, biology, sample_id)
a_ridge  <- assess_correction(counts, out_ridge,  batch, biology, sample_id)
a_combat <- assess_correction(counts, out_combat, batch, biology, sample_id)

merged <- merge_assessments(
  a_mean, a_med, a_ridge, a_combat,
  run_names = c("shift_mean", "shift_median", "ridge", "combat")
)
print(merged)
merged$consolidated_table
```

### MDS before / after

```r
# Coordinates only
mds <- mds_samples(corrected, input_type = "counts", distance = "euclidean")

# Side-by-side MDS colored by batch (shape by biology if provided)
# Requires ggplot2
p <- plot_mds_before_after(
  x_before = counts,
  x_after  = corrected,
  batch    = batch,
  biology  = biology,
  input_type = "counts"
)
print(p)
```

## QC utilities

### Outlier detection (PCA + MDS)

```r
outliers <- detect_outliers_pca_mds(
  data_matrix = counts,
  n_sd = 5
)
outliers$outliers
outliers$pc_plot
outliers$mds_plot
```

### Sex mismatch / contamination

Uses XIST and chrY gene expression vs reported sex:

```r
sex_qc <- detect_sex_mismatch(
  counts   = counts,
  metadata = data.frame(
    sample_id = colnames(counts),
    sex       = reported_sex   # "Male" / "Female"
  )
)
sex_qc$filter_table
sex_qc$plot
```

### RNA vs DNA SNP agreement

Compare RNA-derived VCFs to DNA/genotype VCFs (requires **vcfR**):

```r
snp_qc <- vcf_compare_fast(
  rna_dir = "path/to/rna_vcfs",
  dna_dir = "path/to/dna_vcfs",
  n.cores = 4
)
snp_qc$summary
snp_qc$plot
```

### Bootstrap FASTQ replicates

Generate in silico technical replicates by resampling reads (requires
**ShortRead**):

```r
generate_bootstrap_fastq(
  input_dir  = "path/to/fastq",
  file_names = c("sampleA", "sampleB"),  # prefixes; supports paired-end
  output_dir = "path/to/bootstrap_out",
  reps       = 3,
  percent    = 0.5,
  seed       = 1
)
```

## API overview

| Function | Role |
|----------|------|
| `anchor_correct()` | Unified entry point for all correction methods |
| `correct_shift()` | Location-shift correction (`center = "mean"` or `"median"`) |
| `correct_combat_anchor()` | Anchor-aware ComBat |
| `fit_anchor_ridge()` / `apply_correction_ridge()` | Ridge fit and apply |
| `assess_correction()` | Before/after quality metrics |
| `merge_assessments()` | Combine assessments across methods/runs |
| `mds_samples()` / `plot_mds_before_after()` | MDS visualization |
| `detect_outliers_pca_mds()` | Iterative PCA/MDS outlier detection |
| `detect_sex_mismatch()` | Sex label vs expression QC |
| `default_chrY_genes()` | Default chrY gene set for sex QC |
| `vcf_compare_fast()` | RNA–DNA VCF concordance |
| `generate_bootstrap_fastq()` | Bootstrap FASTQ replicates |

## Requirements for correction

1. A numeric matrix with **genes in rows** and **samples in columns**.
2. `batch` and `sample_id` of length `ncol(x)`.
3. At least some `sample_id` values present in **two or more batches**
   (anchors). Batches without anchors receive zero correction for those
   methods that rely on anchors in that batch.

## Citation / license

This package is licensed under the **GNU General Public License v3 (GPL-3)**.
See the [`LICENSE`](LICENSE) file and `DESCRIPTION` for details.

```
Package: anchorCorrectR
Authors:
  María Rivas-Torrubia <maria.rivas@genyo.es>
  Ana Brescia-Zapata <ana.brescia@genyo.es>
  Guillermo Barturen <gbarturen@ugr.es>
License: GPL-3
```
