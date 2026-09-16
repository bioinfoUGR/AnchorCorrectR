# anchorCorrectR News

## anchorCorrectR 0.99.0

* Initial public development version prepared for GitHub and Bioconductor submission.
* Location-shift correction via `correct_shift()` / `anchor_correct(method = "shift")`
  with `center = "mean"` or `center = "median"`.
* Anchor-aware ComBat and ridge regression methods.
* **Ridge fits batch coefficients on anchor observations only** (same principle as
  shift/ComBat), so non-anchor biology cannot leak into platform estimates.
* Assessment utilities (`assess_correction()`, MDS helpers) and QC helpers
  (outliers, sex mismatch, VCF concordance, FASTQ bootstrap).
* Package documentation, vignette, and unit tests added.
