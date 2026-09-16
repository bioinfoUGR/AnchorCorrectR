# anchorCorrectR News

## anchorCorrectR 0.99.0

* Initial public development version prepared for GitHub and Bioconductor submission.
* Location-shift correction via `correct_shift()` / `anchor_correct(method = "shift")`
  with `center = "mean"` or `center = "median"`.
* Anchor-aware ComBat and ridge regression methods.
* Assessment utilities (`assess_correction()`, MDS helpers) and QC helpers
  (outliers, sex mismatch, VCF concordance, FASTQ bootstrap).
* Package documentation, vignette, and unit tests added.
