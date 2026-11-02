# Reference implementation for a cell-level data.frame.
# Prototype: independent of RTMB. R execution tests are mandatory.
# Fluorescence is assumed ALREADY transformed and quality controlled.
#' Reference utilities for cytoJointR development
#'
#' Provides deterministic data-frame compression, conditional Gaussian log
#' likelihoods and marginal-matched polyfunctionality scores. Full joint
#' hierarchical inference remains unimplemented.
#'
#' This reference version does not contain the planned RTMB joint hierarchical
#' inference engine. Results from the reference functions are not empirical
#' Bayesian posterior estimates.
#'
#' @param markers Character vector of marker names.
#' @param cells One event per row in a data.frame, with sample metadata.
#' @param sample Sample identifier column.
#' @param donor Donor identifier column.
#' @param condition Stimulation-condition column.
#' @param positive Named or ordered columns containing binary positivity calls.
#'   For the likelihood functions, whether the gate is the positive population.
#' @param thresholds Named transformed-fluorescence thresholds, ordered like
#'   markers.
#' @param metadata Additional sample-level metadata fields.
#' @param x Vector of transformed fluorescence values.
#' @param mu Gaussian mean.
#' @param sigma Gaussian standard deviation.
#' @param threshold Fixed transformed positivity cutoff.
#' @param n Number of cells or counts for one category.
#' @param sx Sum of transformed fluorescence.
#' @param ss Sum of squared transformed fluorescence.
#' @param p Probability vector of mutually exclusive patterns.
#' @param p_u Unstimulated pattern-probability vector.
#' @param p_s Stimulated pattern-probability vector.
#' @param bits Matrix of binary marker patterns.
#' @param tol Convergence tolerance for IPF.
#' @param max_iter Maximum iterations for IPF.
#' @param donors Number of donors in the synthetic experiment.
#' @param cells_per_sample Acquired cells per observed donor-condition sample.
#' @param regime Antigen-specific or polyclonal simulation regime.
#' @param missing_controls Fraction of deliberately missing control samples.
#' @param seed Random seed.
#' @param min_functional Deprecated compatibility argument, ignored with a
#'   warning when supplied; inspect positivity denominators and support
#'   diagnostics.
#' @return Depending on the function: binary pattern matrix, compressed object,
#'   numeric log likelihood, pattern-probability vector or named functional
#'   scores.
#' @name cytoJointR-reference
#' @aliases cytoJointR-reference
#' @export
cj_bits <- function(markers) {
  stopifnot(is.character(markers), length(markers) >= 1L,
            length(markers) <= 8L, !anyDuplicated(markers))
  out <- as.matrix(expand.grid(rep(list(0:1), length(markers))))
  colnames(out) <- markers
  storage.mode(out) <- "integer"
  out
}

#' @rdname cytoJointR-reference
#' @export
cj_compress <- function(cells, markers, sample = "sample_id",
                        donor = "donor_id", condition = "condition",
                        positive = paste0(markers, "_pos"),
                        thresholds = NULL,
                        metadata = c("cell_subset", "stimulus", "visit", "batch", "group")) {
  stopifnot(is.data.frame(cells), length(markers) >= 1L, length(markers) <= 8L,
            all(c(sample, donor, condition, markers) %in% names(cells)))
  metadata <- intersect(metadata, names(cells))
  meta_cols <- unique(c(sample, donor, condition, metadata))
  if (anyNA(cells[, meta_cols, drop = FALSE])) stop("Missing sample metadata")
  X <- as.matrix(cells[, markers, drop = FALSE]); storage.mode(X) <- "double"
  if (any(!is.finite(X))) stop("Fluorescence must be finite and transformed")
  if (all(positive %in% names(cells))) {
    Z <- as.matrix(cells[, positive, drop = FALSE])
    if (anyNA(Z) || !all(Z %in% c(0, 1, FALSE, TRUE))) stop("Positivity calls must be binary")
    Z <- matrix(as.integer(Z), ncol = length(markers))
    if (!is.null(thresholds)) {
      if (!identical(names(thresholds), markers))
        stop("Threshold names/order must match markers")
      implied <- sweep(X, 2L, thresholds, `>`)
      if (!identical(unname(Z), unname(matrix(as.integer(implied), ncol=length(markers)))))
        stop("Positivity columns contradict fluorescence thresholds")
    }
  } else {
    if (is.null(thresholds) || !identical(names(thresholds), markers))
      stop("Supply all positivity columns or named thresholds ordered as markers")
    Z <- sweep(X, 2L, thresholds, `>`)
    storage.mode(Z) <- "integer"
  }
  names(cells)[names(cells) == sample] <- ".cj_sample"
  if (sample != ".cj_sample") sample <- ".cj_sample"
  sample_levels <- unique(as.character(cells[[sample]]))
  ids <- match(as.character(cells[[sample]]), sample_levels)
  S <- length(sample_levels)
  K <- as.integer(2^length(markers))
  # Lowest marker has lowest bit; matches cj_bits() expand.grid ordering.
  code <- as.integer(drop(Z %*% (2^(seq_along(markers)-1L))))
  bin <- (ids - 1L) * K + code + 1L
  counts <- tabulate(bin, nbins = S*K)
  sums <- matrix(0, S*K, length(markers), dimnames = list(NULL, markers))
  squares <- sums
  for (j in seq_along(markers)) {
    rs <- rowsum(X[,j], bin, reorder = TRUE)
    rq <- rowsum(X[,j]^2, bin, reorder = TRUE)
    sums[as.integer(rownames(rs)), j] <- rs[,1]
    squares[as.integer(rownames(rq)), j] <- rq[,1]
  }
  idx_first <- match(sample_levels, as.character(cells[[sample]]))
  M <- cells[idx_first, unique(c(sample, donor, condition, metadata)), drop = FALSE]
  # Refuse ambiguous sample-level fields.
  for (field in setdiff(names(M), sample)) {
    u <- tapply(as.character(cells[[field]]), ids, function(x) length(unique(x)))
    if (any(u != 1L)) stop("Inconsistent sample metadata: ", field)
  }
  bits <- cj_bits(markers)
  # Sample-major, then pattern-minor row ordering.
  compressed <- data.frame(
    sample_id = rep(sample_levels, each = K),
    pattern_id = rep(seq_len(K) - 1L, times = S),
    n = counts,
    sums, setNames(as.data.frame(squares), paste0("ss_", markers)),
    check.names = FALSE
  )
  names(compressed)[4L:(3L+length(markers))] <- paste0("sx_", markers)
  # Duplicate the sample metadata; preserve original sample key for joins.
  for (field in setdiff(names(M), sample))
    compressed[[field]] <- rep(M[[field]], each = K)
  stopifnot(sum(compressed$n) == nrow(cells),
            all(colSums(matrix(compressed$n, nrow=K)) == tabulate(ids, nbins=S)))
  structure(list(summary = compressed, bits = bits, markers = markers,
                 sample_metadata = M, thresholds = thresholds,
                 n_cells = nrow(cells)), class = "cj_compressed")
}
