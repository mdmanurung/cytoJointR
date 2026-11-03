# Input layer: every constructor ends in cj_from_counts(), which owns validation.
# Vocabulary follows GLOSSARY.md; rules follow the bundle's 01_DATA_CONTRACT.md.

#' Build a `cj_data` object from exact pattern counts
#'
#' Validates a long table of raw integer pattern counts, fills absent patterns
#' with 0 so all 2^J patterns exist, and pairs every stimulated sample with the
#' control sample of the same donor, visit and parent population.
#'
#' @param counts Long data frame: `sample_id`, one logical or 0/1 column per
#'   marker, and `count` (exact raw integers). Percentages, overlapping marginal
#'   tables and background-subtracted counts are rejected.
#' @param samples One row per sample: `sample_id`, `donor_id`, `visit` (a factor
#'   with explicit levels), `parent_population`, `condition` (`"control"` or the
#'   stimulation name) and optionally `N`, the parent count. When `N` is given
#'   and a sample has no all-negative row, that count is `N` minus the other
#'   patterns, and it must be a nonnegative integer.
#' @param markers Character vector of functional marker names (1 to 8). The
#'   first marker is the low bit of the pattern index.
#'
#' @details
#' A stimulated sample without a control sample gets `control_sample_id = NA`.
#' More than one control sample in a donor x visit x parent population cell is
#' an error: v0.1 has no replicate pooling policy. One control may serve several
#' stimulations and is stored once. Samples with `N = 0` are kept and listed in
#' `provenance$unobserved`; all samples at `N = 0` is an error. Samples whose
#' all-negative count was absent and not reconstructed (no `N`) are listed in
#' `provenance$all_negative_assumed_zero`.
#'
#' @return A `cj_data` object.
#' @examples
#' counts <- expand.grid(sample_id = c("D1_ctrl", "D1_cmv"), IFNG = 0:1, TNF = 0:1,
#'                       stringsAsFactors = FALSE)
#' counts$count <- c(900, 800, 50, 90, 30, 60, 20, 50)
#' samples <- data.frame(
#'   sample_id = c("D1_ctrl", "D1_cmv"), donor_id = "D1",
#'   visit = factor("v1"), parent_population = "CD4",
#'   condition = c("control", "CMV"))
#' x <- cj_from_counts(counts, samples, c("IFNG", "TNF"))
#' x
#' @export
cj_from_counts <- function(counts, samples, markers) {
  bits <- cj_bits(markers)
  stopifnot(is.data.frame(counts), is.data.frame(samples))
  key_cols <- c("sample_id", "donor_id", "visit", "parent_population", "condition")

  # Samples table ------------------------------------------------------------
  if (length(miss <- setdiff(key_cols, names(samples))))
    stop("samples lacks: ", paste(miss, collapse = ", "), call. = FALSE)
  if (anyNA(samples[key_cols])) stop("samples has NA in: ", paste(key_cols, collapse = ", "), call. = FALSE)
  if (anyDuplicated(samples$sample_id)) {
    if (nrow(unique(samples)) > length(unique(samples$sample_id)))
      stop("sample metadata varies within a sample_id", call. = FALSE)
    stop("duplicate sample_id in samples", call. = FALSE)
  }
  if (!is.factor(samples$visit))
    stop("visit must be a factor with explicit levels, e.g. factor(visit, levels = c('v1', 'v2'))",
         call. = FALSE)
  samples$condition <- as.character(samples$condition)
  bad <- samples$condition != CONTROL_CONDITION & tolower(trimws(samples$condition)) == CONTROL_CONDITION
  if (any(bad) || any(tolower(samples$condition) %in% c("unstim", "unstimulated")))
    stop("condition must be exactly '", CONTROL_CONDITION, "' for control samples; ",
         "a stimulation may not be named like a control", call. = FALSE)
  for (k in c("sample_id", "donor_id", "parent_population")) samples[[k]] <- as.character(samples[[k]])
  S <- nrow(samples)
  K <- nrow(bits)

  # Counts table -------------------------------------------------------------
  if (length(miss <- setdiff(c("sample_id", markers, "count"), names(counts))))
    stop("counts lacks: ", paste(miss, collapse = ", "), call. = FALSE)
  if (length(pct <- grep("^(percent|pct|freq|prop)", names(counts), ignore.case = TRUE, value = TRUE)))
    stop("counts has percentage-like column(s) ", paste(pct, collapse = ", "),
         "; supply raw integer counts", call. = FALSE)
  meta_in_counts <- intersect(setdiff(names(samples), "sample_id"), names(counts))
  if (length(extra <- setdiff(names(counts), c("sample_id", markers, "count", names(samples)))))
    stop("counts has unknown column(s): ", paste(extra, collapse = ", "), call. = FALSE)
  counts$sample_id <- as.character(counts$sample_id)
  if (length(unk <- setdiff(counts$sample_id, samples$sample_id)))
    stop("counts has sample_id not in samples: ", paste(utils::head(unk, 5), collapse = ", "), call. = FALSE)
  si <- match(counts$sample_id, samples$sample_id)
  for (k in meta_in_counts)
    if (any(as.character(counts[[k]]) != as.character(samples[[k]][si])))
      stop("sample metadata varies within a sample_id: ", k, call. = FALSE)

  cnt <- counts$count
  if (!is.numeric(cnt) || anyNA(cnt) || any(!is.finite(cnt)) || any(cnt < 0) || any(cnt != round(cnt)))
    stop("count must be nonnegative integers without NA (no fractions or percentages)", call. = FALSE)
  if (any(cnt > .Machine$integer.max)) stop("count exceeds the integer range", call. = FALSE)
  Z <- as.matrix(counts[markers])
  if (anyNA(Z) || !all(Z %in% c(0, 1))) stop("marker columns must be logical or 0/1 without NA", call. = FALSE)
  col <- drop(Z %*% 2^(seq_along(markers) - 1L)) + 1
  if (anyDuplicated(data.frame(si, col))) stop("duplicate pattern keys within a sample", call. = FALSE)

  # Fill absent patterns with 0 ----------------------------------------------
  M <- matrix(0, S, K)
  seen <- matrix(FALSE, S, K)
  M[cbind(si, col)] <- cnt
  seen[cbind(si, col)] <- TRUE

  # Parent count N: reconstruct the all-negative pattern or check the sum -----
  N <- if ("N" %in% names(samples)) as.numeric(samples$N) else rep(NA_real_, S)
  if (any(N < 0 | N != round(N), na.rm = TRUE)) stop("N must be nonnegative integers", call. = FALSE)
  rest <- rowSums(M[, -1L, drop = FALSE])
  recon <- !is.na(N) & !seen[, 1L]
  if (any(recon & !rowSums(seen) & N > 0))
    stop("samples with N > 0 but no count rows: ",
         paste(samples$sample_id[recon & !rowSums(seen) & N > 0], collapse = ", "), call. = FALSE)
  if (any(resid_bad <- recon & N - rest < 0))
    stop("N is smaller than the summed patterns (negative all-negative count) in: ",
         paste(samples$sample_id[resid_bad], collapse = ", "), call. = FALSE)
  M[recon, 1L] <- (N - rest)[recon]
  if (any(mismatch <- !is.na(N) & rowSums(M) != N))
    stop("pattern counts do not sum to N in: ", paste(samples$sample_id[mismatch], collapse = ", "),
         call. = FALSE)
  samples$N <- rowSums(M)
  if (all(samples$N == 0)) stop("all samples have N = 0; nothing to infer from", call. = FALSE)

  counts_mat <- M
  storage.mode(counts_mat) <- "integer"
  dimnames(counts_mat) <- list(samples$sample_id, pattern_labels(bits))
  samples <- samples[c(SAMPLE_COLUMNS, setdiff(names(samples), SAMPLE_COLUMNS))]
  rownames(samples) <- NULL

  contrasts <- derive_contrasts(samples)

  new_cj_data(counts_mat, samples, contrasts, bits, provenance = list(
    source = "counts",
    unobserved = samples$sample_id[samples$N == 0],
    all_negative_reconstructed = samples$sample_id[recon],
    all_negative_assumed_zero = samples$sample_id[!seen[, 1L] & !recon]))
}

#' Keep one parent population of a `cj_data` object
#'
#' v0.1 fits one parent population at a time ([cj_design()] and [cj_fit()] stop
#' on more than one). `cj_subset()` keeps the samples, counts and contrasts of
#' one population and re-runs the `cj_data` validators.
#'
#' @param x A `cj_data` object.
#' @param parent_population One population name found in `x$samples`.
#' @return A `cj_data` object.
#' @export
#' @examples
#' sim <- cj_simulate_counts(n_donors = 3, seed = 1)
#' cj_subset(sim$data, "CD4")$samples$parent_population |> unique()
cj_subset <- function(x, parent_population) {
  stopifnot(inherits(x, "cj_data"), length(parent_population) == 1L)
  keep <- as.character(x$samples$parent_population) == parent_population
  if (!any(keep)) stop("no samples in parent population: ", parent_population, call. = FALSE)
  samples <- x$samples[keep, , drop = FALSE]
  rownames(samples) <- NULL
  contrasts <- x$contrasts[as.character(x$contrasts$parent_population) == parent_population, ,
                           drop = FALSE]
  rownames(contrasts) <- NULL
  carry_intensity(new_cj_data(x$counts[keep, , drop = FALSE], samples, contrasts, x$bits, x$provenance),
                  x, keep)
}

# Long count table from a K-pattern code vector per sample; shared by the adapters.
counts_long <- function(sample_id, code, n, bits) {
  data.frame(sample_id = sample_id, bits[code + 1L, , drop = FALSE], count = n,
             check.names = FALSE, row.names = NULL)
}

#' Build a `cj_data` object from a per-cell table
#'
#' Compresses cells to per-sample pattern counts with [cj_compress()], then
#' calls [cj_from_counts()].
#'
#' @param cells One row per cell: `sample_id`, one fluorescence column per
#'   marker and, unless `thresholds` is given, `<marker>_pos` binary calls.
#'   Donor and condition come from `samples`; if `cells` also has `donor_id` or
#'   `condition` they must agree.
#' @param markers Marker names, as in [cj_from_counts()].
#' @param samples Sample table, as in [cj_from_counts()].
#' @param thresholds Optional named numeric vector (names equal `markers`);
#'   positive means fluorescence above the threshold.
#' @return A `cj_data` object.
#' @examples
#' cells <- data.frame(sample_id = rep(c("c", "s"), each = 4),
#'                     A = c(0, 0, 2, 2, 2, 2, 0, 2), B = c(0, 2, 0, 2, 2, 2, 0, 0))
#' samples <- data.frame(sample_id = c("c", "s"), donor_id = "D1", visit = factor("v1"),
#'                       parent_population = "CD4", condition = c("control", "X"))
#' cj_from_events(cells, c("A", "B"), samples, thresholds = c(A = 1, B = 1))
#' @export
cj_from_events <- function(cells, markers, samples, thresholds = NULL) {
  stopifnot(is.data.frame(cells), is.data.frame(samples), "sample_id" %in% names(cells))
  if (length(unk <- setdiff(cells$sample_id, samples$sample_id)))
    stop("cells has sample_id not in samples: ", paste(utils::head(unk, 5), collapse = ", "), call. = FALSE)
  m <- match(cells$sample_id, samples$sample_id)
  for (k in c("donor_id", "condition")) {
    if (k %in% names(cells) && any(as.character(cells[[k]]) != as.character(samples[[k]][m])))
      stop(k, " in cells disagrees with samples", call. = FALSE)
    cells[[k]] <- samples[[k]][m]
  }
  cmp <- cj_compress(cells, markers, thresholds = thresholds, metadata = character())
  s <- cmp$summary
  cj_from_counts(counts_long(s$sample_id, s$pattern_id, s$n, cmp$bits), samples, markers)
}

#' Build a `cj_data` object from a flowWorkspace GatingSet
#'
#' Counts, per sample, the parent-population events that are positive in each
#' marker's gate, and tabulates the exact patterns. Requires the flowWorkspace
#' package.
#'
#' @param gs A `GatingSet`.
#' @param node Path of the parent population node.
#' @param markers Named character vector. Names are marker names, values are the
#'   node paths of each marker's positive gate (the "boolean-gate route").
#' @param samples Sample table, as in [cj_from_counts()]; `sample_id` must equal
#'   the sample names in `gs`. `N` is set to the parent count when absent and
#'   must equal it when given.
#' @return A `cj_data` object.
#' @examples
#' \dontrun{
#' cj_from_gatingset(gs, "/live/CD4", c(IFNG = "/live/CD4/IFNG+", TNF = "/live/CD4/TNF+"),
#'                   samples)
#' }
#' @export
cj_from_gatingset <- function(gs, node, markers, samples) {
  if (!requireNamespace("flowWorkspace", quietly = TRUE))
    stop("cj_from_gatingset() needs the flowWorkspace package", call. = FALSE)
  if (is.null(names(markers)) || anyNA(names(markers)))
    stop("markers must be a named vector: names = markers, values = positive-gate nodes", call. = FALSE)
  if (length(unk <- setdiff(samples$sample_id, flowWorkspace::sampleNames(gs))))
    stop("samples not in the GatingSet: ", paste(unk, collapse = ", "), call. = FALSE)
  if (length(extra <- setdiff(flowWorkspace::sampleNames(gs), samples$sample_id)))
    message("GatingSet samples not in `samples` and ignored: ", paste(extra, collapse = ", "))
  per_sample <- lapply(as.character(samples$sample_id), function(id) {
    gh <- gs[[id]]
    parent <- flowWorkspace::gh_pop_get_indices(gh, node)
    pos <- lapply(markers, function(g) flowWorkspace::gh_pop_get_indices(gh, g))
    tabulate_patterns(parent, pos)
  })
  K <- 2L^length(markers)
  counts <- counts_long(rep(as.character(samples$sample_id), each = K), rep(seq_len(K) - 1L, nrow(samples)),
                        unlist(per_sample), cj_bits(names(markers)))
  cj_from_counts(counts, samples, names(markers))
}

# Pattern counts (length 2^J, pattern k at k + 1) among parent events; pos is a
# list of logical vectors over all events, one per marker, first marker low bit.
tabulate_patterns <- function(parent, pos) {
  code <- Reduce(`+`, Map(function(p, j) 2^(j - 1L) * p[parent], pos, seq_along(pos)))
  tabulate(code + 1L, nbins = 2L^length(pos))
}

#' Build a `cj_data` object from polyICSFlow cell counts
#'
#' Accepts the table returned by polyICSFlow `extractCellCounts()`: a factor or
#' character `MarkerComb` (marker names each followed by `+` or `-`, joined with
#' no separator, e.g. `"IFNg+TNF-"`) and the cell count `n`. Other columns
#' (`MarkerCount`, `Functionality`) are ignored. The table covers one sample, so
#' for several samples add a `sample_id` column (for example by row-binding the
#' per-sample tables).
#'
#' @param x The polyICSFlow count table.
#' @param samples Sample table, as in [cj_from_counts()]. One row if `x` has no
#'   `sample_id` column.
#' @param markers Marker names in the order used in `MarkerComb`.
#' @return A `cj_data` object.
#' @examples
#' x <- data.frame(MarkerComb = c("IFNg-TNF-", "IFNg+TNF-", "IFNg+TNF+"), n = c(90L, 6L, 4L))
#' samples <- data.frame(sample_id = "s1", donor_id = "D1", visit = factor("v1"),
#'                       parent_population = "CD4", condition = "control", N = 100)
#' cj_from_polyicsflow(x, samples, c("IFNg", "TNF"))
#' @export
cj_from_polyicsflow <- function(x, samples, markers) {
  stopifnot(is.data.frame(x), all(c("MarkerComb", "n") %in% names(x)))
  if (!"sample_id" %in% names(x)) {
    if (nrow(samples) != 1L) stop("x has no sample_id column; samples must have exactly one row", call. = FALSE)
    x$sample_id <- samples$sample_id
  }
  bits <- cj_bits(markers)
  expected <- apply(bits, 1L, function(b) paste0(markers, ifelse(b == 1L, "+", "-"), collapse = ""))
  code <- match(as.character(x$MarkerComb), expected) - 1L
  if (anyNA(code))
    stop("MarkerComb not matching markers (order matters): ",
         paste(utils::head(unique(as.character(x$MarkerComb)[is.na(code)]), 5), collapse = ", "), call. = FALSE)
  cj_from_counts(counts_long(as.character(x$sample_id), code, x$n, bits), samples, markers)
}
