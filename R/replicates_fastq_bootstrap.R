# FASTQ bootstrap utilities for anchorCorrectR
# ===========================================

#' Generate bootstrap replicates of FASTQ files
#'
#' Randomly resamples reads with replacement from specified FASTQ/FASTQ.gz files
#' and writes bootstrap replicates to `output_dir`.
#'
#' @param input_dir Directory containing input FASTQ or FASTQ.gz files.
#' @param file_names Character vector of sample prefixes. For each prefix the
#'   function searches `input_dir` for files named `prefix.fastq[.gz]` (single
#'   end) or `prefix_1.fastq[.gz]`/`prefix_2.fastq[.gz]` (paired-end). Presence
#'   of both `_1` and `_2` files triggers paired processing automatically.
#' @param output_dir Directory where bootstrap FASTQ files will be written. Will
#'   be created recursively if it does not exist.
#' @param reps Integer; number of bootstrap replicates to generate per file.
#' @param percent Numeric in (0, 1]; proportion of reads sampled (with
#'   replacement) for each bootstrap replicate.
#' @param merge_inputs Logical; when `TRUE` (default) reads from all inputs are
#'   merged prior to bootstrapping and one merged replicate (per requested
#'   `reps`) is written out. Set to `FALSE` to generate replicates separately for
#'   each input file or pair.
#' @param merge_label Character scalar used as the base file name when
#'   `merge_inputs = TRUE`.
#' @param use_parallel Logical; when `TRUE` (default) bootstrap replicates are
#'   generated in parallel via `parallel::mclapply()`. Set to `FALSE` to force
#'   sequential execution.
#' @param workers Integer; number of parallel workers to use when
#'   `use_parallel = TRUE`. Defaults to `parallel::detectCores() - 1`.
#' @param chunk_size Integer; number of reads pulled per chunk when streaming
#'   FASTQ data. Larger values reduce disk I/O but increase peak memory while
#'   generating each replicate.
#' @param seed Optional integer to set a reproducible RNG seed.
#'
#' @return Invisibly returns `TRUE` if all replicates are generated.
#' @examples
#' \dontrun{
#' generate_bootstrap_fastq(
#'   input_dir = "fastq_inputs",
#'   file_names = c("sample1", "sample2"),
#'   output_dir = "replicates_fastq_bootstrap",
#'   reps = 5,
#'   percent = 0.25,
#'   seed = 123
#' )
#' }
#' @export
generate_bootstrap_fastq <- function(
    input_dir,
    file_names,
    output_dir,
    reps = 3,
    percent = 0.5,
    merge_inputs = TRUE,
    merge_label = "merged",
    use_parallel = TRUE,
    workers = NULL,
    chunk_size = 100000,
    seed = NULL) {

  if (!requireNamespace("ShortRead", quietly = TRUE)) {
    stop("Package 'ShortRead' is required for generate_bootstrap_fastq().", call. = FALSE)
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required for generate_bootstrap_fastq().", call. = FALSE)
  }

  stopifnot(dir.exists(input_dir))
  stopifnot(is.numeric(percent), percent > 0, percent <= 1)
  stopifnot(is.numeric(reps), reps >= 1)
  stopifnot(is.logical(merge_inputs), length(merge_inputs) == 1)
  stopifnot(is.character(merge_label), length(merge_label) == 1)
  stopifnot(is.logical(use_parallel), length(use_parallel) == 1)
  stopifnot(is.numeric(chunk_size), length(chunk_size) == 1)

  if (use_parallel) {
    if (is.null(workers)) {
      workers <- max(1L, parallel::detectCores() - 1L)
    }
    workers <- as.integer(workers)
    if (is.na(workers) || workers < 1L) {
      stop("`workers` must be a positive integer when `use_parallel = TRUE`.")
    }
    message(sprintf("Using %d parallel workers for bootstrap replicates", workers))
  } else {
    workers <- 1L
  }

  chunk_size <- as.integer(chunk_size)
  if (is.na(chunk_size) || chunk_size <= 0L) {
    stop("`chunk_size` must be a positive integer.")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  specs <- normalize_fastq_inputs(file_names, input_dir)
  specs <- lapply(
    specs,
    annotate_sample_metadata,
    input_dir = input_dir,
    chunk_size = chunk_size
  )

  total_steps <- if (merge_inputs) reps else reps * length(specs)
  progress <- new_progress_bar(total_steps)
  on.exit(progress$close(), add = TRUE)

  runner <- make_replicate_runner(use_parallel, workers, progress)

  if (merge_inputs) {
    run_merged_bootstrap(
      specs = specs,
      percent = percent,
      reps = reps,
      input_dir = input_dir,
      output_dir = output_dir,
      merge_label = merge_label,
      runner = runner,
      chunk_size = chunk_size
    )
  } else {
    for (spec in specs) {
      process_sample_replicates(
        spec = spec,
        percent = percent,
        reps = reps,
        input_dir = input_dir,
        output_dir = output_dir,
        runner = runner,
        chunk_size = chunk_size
      )
    }
  }

  invisible(TRUE)
}

run_merged_bootstrap <- function(specs,
                                 percent,
                                 reps,
                                 input_dir,
                                 output_dir,
                                 merge_label,
                                 runner,
                                 chunk_size) {
  if (length(specs) == 0) {
    stop("No samples supplied for merging.")
  }

  pairing <- vapply(specs, function(s) s$paired, logical(1))
  if (length(unique(pairing)) != 1) {
    stop("Cannot merge single-end and paired-end inputs together.")
  }

  read_counts <- vapply(specs, function(s) s$read_count, numeric(1))
  total_reads <- sum(read_counts)
  if (total_reads == 0) {
    warning("No reads found after merging; nothing to write.")
    runner(reps, function(...) NULL)
    return(invisible(TRUE))
  }

  compress_flag <- any(unlist(lapply(specs, `[[`, "compress")))
  is_paired <- pairing[1]
  message(
    "Processing merged inputs as ",
    if (is_paired) "paired-end" else "single-end"
  )

  runner(
    reps,
    function(r) {
      indices <- simulate_bootstrap_indices(total_reads, percent)
      split_idx <- split_indices_by_sample(indices, read_counts)
      out_info <- build_merged_output_paths(
        output_dir = output_dir,
        merge_label = merge_label,
        rep_id = r,
        is_paired = is_paired,
        compress_flag = compress_flag
      )

      mode <- "w"
      for (i in seq_along(specs)) {
        local_idx <- split_idx[[i]]
        if (length(local_idx) == 0L) {
          next
        }
        write_indices_to_fastq(
          input_paths = file.path(input_dir, specs[[i]]$files),
          indices = local_idx,
          output_paths = out_info$paths,
          compress_flags = out_info$compress,
          chunk_size = chunk_size,
          mode = mode
        )
        mode <- "a"
      }

      invisible(NULL)
    }
  )

  invisible(TRUE)
}

process_sample_replicates <- function(spec,
                                      percent,
                                      reps,
                                      input_dir,
                                      output_dir,
                                      runner,
                                      chunk_size) {
  message(
    "Processing ",
    spec$label,
    if (spec$paired) " (paired-end)" else " (single-end)"
  )

  if (spec$read_count == 0) {
    warning(sprintf("No reads found in '%s'; skipping.", spec$label))
    runner(reps, function(...) NULL)
    return(invisible(NULL))
  }

  runner(
    reps,
    function(r) {
      indices <- simulate_bootstrap_indices(spec$read_count, percent)
      out_info <- build_sample_output_paths(spec$files, output_dir, r)
      write_indices_to_fastq(
        input_paths = file.path(input_dir, spec$files),
        indices = indices,
        output_paths = out_info$paths,
        compress_flags = out_info$compress,
        chunk_size = chunk_size,
        mode = "w"
      )
      invisible(NULL)
    }
  )

  invisible(NULL)
}

annotate_sample_metadata <- function(spec, input_dir, chunk_size) {
  fq_paths <- file.path(input_dir, spec$files)
  missing <- spec$files[!file.exists(fq_paths)]
  if (length(missing) > 0) {
    stop(sprintf(
      "The following files were not found in %s: %s",
      input_dir,
      paste(missing, collapse = ", ")
    ))
  }

  read_counts <- vapply(
    fq_paths,
    count_fastq_reads,
    numeric(1),
    chunk_size = chunk_size
  )

  if (spec$paired && diff(range(read_counts)) != 0) {
    stop(sprintf(
      "Paired files for '%s' do not have the same number of reads.",
      spec$label
    ))
  }

  spec$read_count <- read_counts[1]
  spec$compress <- grepl("\\.gz$", spec$files)
  spec
}

count_fastq_reads <- function(path, chunk_size) {
  streamer <- ShortRead::FastqStreamer(path, n = chunk_size)
  on.exit(close(streamer), add = TRUE)
  total <- 0L
  repeat {
    chunk <- ShortRead::yield(streamer)
    len <- length(chunk)
    if (len == 0L) {
      break
    }
    total <- total + len
  }
  total
}

simulate_bootstrap_indices <- function(read_count, percent) {
  n_draw <- ceiling(percent * read_count)
  if (n_draw <= 0L) {
    stop("Computed draw size must be positive. Check `percent` and read counts.")
  }
  sort(sample.int(read_count, n_draw, replace = TRUE))
}

build_sample_output_paths <- function(files, output_dir, rep_id) {
  base_names <- basename(files)
  paths <- file.path(
    output_dir,
    sprintf("%s.rep%02d.fastq", base_names, rep_id)
  )
  list(
    paths = paths,
    compress = grepl("\\.gz$", base_names)
  )
}

build_merged_output_paths <- function(output_dir,
                                      merge_label,
                                      rep_id,
                                      is_paired,
                                      compress_flag) {
  suffixes <- if (is_paired) c("_R1", "_R2") else ""
  file_ext <- if (compress_flag) ".fastq.gz" else ".fastq"
  paths <- file.path(
    output_dir,
    sprintf("%s%s.rep%02d%s", merge_label, suffixes, rep_id, file_ext)
  )
  list(
    paths = paths,
    compress = rep(compress_flag, length(paths))
  )
}

split_indices_by_sample <- function(indices, read_counts) {
  bounds <- c(0L, cumsum(as.integer(read_counts)))
  sample_ids <- findInterval(indices, bounds, rightmost.closed = TRUE)
  local_idx <- indices - bounds[sample_ids]
  unname(split(
    local_idx,
    factor(sample_ids, levels = seq_along(read_counts)),
    drop = FALSE
  ))
}

write_indices_to_fastq <- function(input_paths,
                                   indices,
                                   output_paths,
                                   compress_flags,
                                   chunk_size,
                                   mode = "w") {
  if (!length(indices)) {
    return(invisible(NULL))
  }

  streamers <- lapply(
    input_paths,
    function(p) ShortRead::FastqStreamer(p, n = chunk_size)
  )
  on.exit(lapply(streamers, close), add = TRUE)

  idx_ptr <- 1L
  total_idx <- length(indices)
  chunk_start <- 1L
  mode_current <- mode

  repeat {
    chunks <- lapply(streamers, ShortRead::yield)
    chunk_len <- length(chunks[[1]])
    if (chunk_len == 0L) {
      break
    }

    chunk_end <- chunk_start + chunk_len - 1L

    if (idx_ptr <= total_idx && indices[idx_ptr] <= chunk_end) {
      end_ptr <- idx_ptr
      while (end_ptr <= total_idx && indices[end_ptr] <= chunk_end) {
        end_ptr <- end_ptr + 1L
      }
      rel_idx <- indices[idx_ptr:(end_ptr - 1L)] - chunk_start + 1L
      for (i in seq_along(chunks)) {
        ShortRead::writeFastq(
          chunks[[i]][rel_idx],
          output_paths[[i]],
          mode = mode_current,
          compress = compress_flags[[i]]
        )
      }
      mode_current <- "a"
      idx_ptr <- end_ptr
      if (idx_ptr > total_idx) {
        break
      }
    }

    chunk_start <- chunk_end + 1L
  }

  invisible(NULL)
}

normalize_fastq_inputs <- function(file_names, input_dir) {
  if (!is.character(file_names) || length(file_names) == 0) {
    stop("`file_names` must be a non-empty character vector of prefixes.")
  }

  specs <- vector("list", length(file_names))
  for (i in seq_along(file_names)) {
    prefix <- file_names[[i]]
    if (!nzchar(prefix)) {
      stop("Prefixes in `file_names` must be non-empty strings.")
    }
    specs[[i]] <- resolve_fastq_prefix(input_dir, prefix)
  }

  specs
}

make_replicate_runner <- function(use_parallel, workers, progress) {
  function(reps, task) {
    reps <- as.integer(reps)
    if (is.na(reps) || reps <= 0L) {
      return(invisible(NULL))
    }

    if (!use_parallel || reps == 1L) {
      lapply(
        seq_len(reps),
        function(r) {
          task(r)
          progress$tick()
        }
      )
    } else {
      # Use mc.set.seed = TRUE to ensure each worker gets proper RNG
      # Use mc.preschedule = FALSE for I/O-bound tasks (better for file operations)
      # Don't update progress inside workers (not fork-safe)
      actual_workers <- min(workers, reps)
      res <- parallel::mclapply(
        seq_len(reps),
        task,
        mc.cores = actual_workers,
        mc.set.seed = TRUE,
        mc.preschedule = FALSE
      )
      # Update progress after all workers complete
      progress$tick(length(res))
      res
    }
  }
}

new_progress_bar <- function(total_steps) {
  total_steps <- as.integer(total_steps)
  if (is.na(total_steps) || total_steps <= 0L) {
    return(list(
      tick = function(...) NULL,
      close = function() NULL
    ))
  }

  pb <- utils::txtProgressBar(min = 0, max = total_steps, style = 3)
  env <- new.env(parent = emptyenv())
  env$val <- 0L

  list(
    tick = function(n = 1L) {
      n <- as.integer(n)
      if (is.na(n) || n <= 0L) {
        return(invisible(NULL))
      }
      env$val <- min(total_steps, env$val + n)
      utils::setTxtProgressBar(pb, env$val)
      invisible(NULL)
    },
    close = function() {
      utils::setTxtProgressBar(pb, total_steps)
      close(pb)
    }
  )
}

resolve_fastq_prefix <- function(input_dir, prefix) {
  exts <- c(".fastq.gz", ".fastq")
  single_candidate <- NULL

  for (ext in exts) {
    r1 <- paste0(prefix, "_1", ext)
    r2 <- paste0(prefix, "_2", ext)
    path1 <- file.path(input_dir, r1)
    path2 <- file.path(input_dir, r2)
    if (file.exists(path1) && file.exists(path2)) {
      return(list(files = c(r1, r2), paired = TRUE, label = prefix))
    }
    if (file.exists(path1) && is.null(single_candidate)) {
      single_candidate <- r1
    }
  }

  for (ext in exts) {
    candidate <- paste0(prefix, ext)
    if (file.exists(file.path(input_dir, candidate))) {
      return(list(files = candidate, paired = FALSE, label = prefix))
    }
  }

  if (!is.null(single_candidate)) {
    warning(sprintf(
      "Found %s but not matching '_2' file; treating as single-end.",
      single_candidate
    ))
    return(list(files = single_candidate, paired = FALSE, label = prefix))
  }

  stop(sprintf(
    "Could not find FASTQ files for prefix '%s' in %s.",
    prefix,
    input_dir
  ))
}


