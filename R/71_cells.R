# From a cell data frame to a cj_data that also keeps per-pattern intensity
# sums. Reuses cj_compress() and cj_from_counts(); adds no second compressor.

#' Build a `cj_data` object from a cell data frame, with fluorescence summaries
#'
#' Transforms the marker columns, applies the gate's thresholds, tabulates exact
#' patterns per sample (through [cj_compress()] and [cj_from_counts()]) and keeps
#' what the fluorescence scores need.
#'
#' @details
#' The returned `cj_data` has three extra fields, ignored by the count model:
#' * `intensity`: `thresholds`, `transform` (label), and two arrays
#'   `[sample, pattern, marker]`, `sx` (sum of `X - t` over the events of that
#'   sample and pattern) and `ss` (sum of `(X - t)^2`). The pattern count `n` is
#'   `x$counts`. `X` is the transformed fluorescence and `t` the threshold, so
#'   `sx / n` is the mean excess over the threshold. For a marker-negative
#'   pattern it is negative.
#' * `gate`: the `cj_gate`.
#' * `events`: sample index and transformed `X` of every event, kept only when
#'   `keep_events = TRUE`. [cj_fluorescence_scores()] needs it for the
#'   cell-level bootstrap. It adds one numeric matrix of events by markers.
#'
#' @param cells Data frame, one row per event: `sample_id` and one fluorescence
#'   column per marker. Events must already be in the parent population.
#' @param markers Marker names mapped to columns, as in [cj_gate()]. `NULL`
#'   takes them from `gate`.
#' @param samples Sample table, as in [cj_from_counts()]. Its `condition` may use
#'   another label for controls, named in `control`.
#' @param gate A [cj_gate()]. `NULL` builds one with `cj_gate(cells, markers,
#'   samples, transform = transform, control = control, ...)`.
#' @param transform Required when `gate` is `NULL` (see [cj_gate()]). When a gate is
#'   given, its transform is used and a `transform` here must be `NULL`.
#' @param control Control labels in `samples$condition`; with a `gate`, the
#'   gate's own.
#' @param keep_events Keep the events for the bootstrap.
#' @param ... Passed to [cj_gate()] when `gate` is `NULL`.
#' @return A `cj_data` object with the extra fields above.
#' @export
#' @examples
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 300, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' samples$visit <- factor("v1")
#' samples$parent_population <- "CD4"
#' x <- cj_from_cells(sim$cells, c("IFNG", "TNF", "IL2"), samples,
#'                    transform = "none", control = "unstim", method = "quantile")
#' x
cj_from_cells <- function(cells, markers = NULL, samples, gate = NULL, transform = NULL,
                          control = "control", keep_events = TRUE, ...) {
  stopifnot(is.data.frame(cells), is.data.frame(samples), "sample_id" %in% names(cells))
  if (is.null(gate)) {
    if (is.null(markers)) stop("give markers or a gate", call. = FALSE)
    gate <- cj_gate(cells, markers, samples, transform = transform, control = control, ...)
  } else {
    stopifnot(inherits(gate, "cj_gate"))
    if (!is.null(transform))
      stop("the gate fixes the transform (", gate$transform$label, "); leave transform = NULL",
           call. = FALSE)
    if (!is.null(markers) && !identical(names(gate_markers(markers)), names(gate$markers)))
      stop("markers differ from the gate's markers", call. = FALSE)
    control <- gate$control
  }
  mk <- gate$markers
  th <- gate$thresholds[names(mk)]
  samples <- recode_control(samples, control)
  X <- gate_matrix(cells, mk, gate$transform)
  ev <- data.frame(sample_id = as.character(cells$sample_id), X, check.names = FALSE)
  m <- match(ev$sample_id, samples$sample_id)
  if (anyNA(m)) stop("cells has sample_id not in samples: ",
                     paste(utils::head(unique(ev$sample_id[is.na(m)]), 5), collapse = ", "), call. = FALSE)
  ev$donor_id <- as.character(samples$donor_id)[m]
  ev$condition <- samples$condition[m]
  x <- cj_from_events(ev, names(mk), samples, thresholds = th)
  s <- cj_compress(ev, names(mk), thresholds = th, metadata = character())$summary

  # Shift the sums to X - t: sum(X - t) = sx - n t; sum((X - t)^2) = ss - 2 t sx + n t^2.
  K <- nrow(x$bits)
  lv <- unique(s$sample_id)
  stopifnot("compressed counts match the cj_data counts" =
              identical(unname(x$counts[lv, ]), unname(t(matrix(as.integer(s$n), nrow = K)))))
  dn <- list(x$samples$sample_id, colnames(x$counts), names(mk))
  sx <- ss <- array(0, vapply(dn, length, 1L), dimnames = dn)
  for (j in names(mk)) {
    a <- matrix(s[[paste0("sx_", j)]], nrow = K)
    q <- matrix(s[[paste0("ss_", j)]], nrow = K)
    n <- matrix(s$n, nrow = K)
    sx[lv, , j] <- t(a - n * th[[j]])
    ss[lv, , j] <- t(q - 2 * th[[j]] * a + n * th[[j]]^2)
  }
  x$intensity <- list(thresholds = th, transform = gate$transform$label, sx = sx, ss = ss)
  x$gate <- gate
  if (keep_events) x$events <- list(sample = match(ev$sample_id, x$samples$sample_id), X = X)
  x$provenance$source <- "cells"
  x$provenance$gate <- list(method = gate$method, method_used = gate$method_used,
                            transform = gate$transform$label, n_control = gate$n_control)
  x
}

# Carry the intensity fields of x into out = cj_subset(x), keeping the rows in `keep`.
carry_intensity <- function(out, x, keep) {
  if (is.null(x$intensity)) return(out)
  out$intensity <- x$intensity
  out$intensity$sx <- x$intensity$sx[keep, , , drop = FALSE]
  out$intensity$ss <- x$intensity$ss[keep, , , drop = FALSE]
  out$gate <- x$gate
  if (!is.null(x$events)) {
    new_row <- cumsum(keep)
    ok <- keep[x$events$sample]
    out$events <- list(sample = new_row[x$events$sample[ok]], X = x$events$X[ok, , drop = FALSE])
  }
  out
}
