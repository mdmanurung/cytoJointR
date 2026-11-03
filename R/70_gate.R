# Pooled-control gating for the adoption API. One threshold per marker, found on
# all control (unstimulated) events pooled, on an explicitly stated transform.
# The result feeds cj_from_cells(). Vocabulary follows GLOSSARY.md.

# Every optional-package check goes through here, so tests can mock it. (Loading
# openCyto aborts R on some hosts because it pulls in flowWorkspace.)
has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

GATE_METHODS <- c("tail", "mindensity", "quantile", "manual")

# Named character vector marker name -> data-frame column. Unnamed means the
# column name is the marker name.
gate_markers <- function(markers) {
  if (!is.character(markers) || !length(markers)) stop("markers must be a character vector", call. = FALSE)
  if (is.null(names(markers))) names(markers) <- markers
  if (anyNA(names(markers)) || any(!nzchar(names(markers))) || anyDuplicated(names(markers)))
    stop("marker names must be unique and non-empty", call. = FALSE)
  markers
}

# The transform is never assumed. "none" declares that cells are already on the
# analysis scale; "asinh" needs its cofactor; a function is used as given.
gate_transform <- function(transform, cofactor = NULL) {
  if (is.null(transform))
    stop("transform is required. Use transform = \"none\" if the cells are already transformed, ",
         "transform = \"asinh\" with cofactor = <number>, or a function of one argument.",
         call. = FALSE)
  if (is.function(transform))
    return(list(fn = transform, label = attr(transform, "label") %||% "custom function"))
  if (identical(transform, "none")) return(list(fn = identity, label = "none"))
  if (identical(transform, "asinh")) {
    if (!is.numeric(cofactor) || length(cofactor) != 1L || !is.finite(cofactor) || cofactor <= 0)
      stop("transform = \"asinh\" needs a positive cofactor, e.g. cofactor = 5", call. = FALSE)
    return(list(fn = function(x) asinh(x / cofactor), label = sprintf("asinh(x / %g)", cofactor)))
  }
  stop("transform must be \"none\", \"asinh\" or a function", call. = FALSE)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

# Samples whose condition is one of the control labels become "control".
recode_control <- function(samples, control) {
  samples$condition <- as.character(samples$condition)
  samples$condition[samples$condition %in% control] <- CONTROL_CONDITION
  samples
}

# Transformed marker matrix (events x markers), columns named by marker.
gate_matrix <- function(cells, markers, transform) {
  if (length(miss <- setdiff(markers, names(cells))))
    stop("cells lacks column(s): ", paste(miss, collapse = ", "), call. = FALSE)
  X <- matrix(NA_real_, nrow(cells), length(markers), dimnames = list(NULL, names(markers)))
  for (m in names(markers)) X[, m] <- as.numeric(transform$fn(as.numeric(cells[[markers[[m]]]])))
  bad <- colSums(!is.finite(X))
  if (any(bad > 0))
    stop("transformed values are not finite (", paste0(names(bad)[bad > 0], ": ", bad[bad > 0],
         collapse = ", "), "). Fix the transform or remove those events; nothing is dropped silently.",
         call. = FALSE)
  X
}

# Tail gate on a pooled control sample: the main peak is the negative
# population; its spread is estimated from the half below the mode (mirrored,
# so positive events in the upper tail do not inflate it); the cutoff sits
# qnorm(1 - tail_p) spreads above the mode.
gate_cut_tail <- function(x, tail_p) {
  d <- stats::density(x)
  m <- d$x[which.max(d$y)]
  low <- x[x <= m]
  if (length(low) < 10L) stop("too few control events below the density mode to estimate the spread", call. = FALSE)
  s <- sqrt(mean((low - m)^2))
  if (!is.finite(s) || s <= 0) stop("the control events have no spread below the mode", call. = FALSE)
  m + stats::qnorm(1 - tail_p) * s
}

# openCyto::gate_mindensity() on one pooled channel; returns the cut point.
gate_cut_mindensity <- function(x, marker, ...) {
  fn <- getExportedValue("openCyto", "gate_mindensity")
  fr <- flowCore::flowFrame(matrix(x, ncol = 1L, dimnames = list(NULL, marker)))
  g <- fn(fr, channel = marker, ...)
  as.numeric(unname(g@min)[1L])
}

#' Gate cytokines with one pooled threshold per marker
#'
#' Finds one positivity threshold per marker on all control (unstimulated)
#' events pooled across donors, visits and samples. The threshold is on the
#' transformed scale, and the transform must be stated.
#'
#' @details
#' Methods:
#' * `"tail"` (default): the main peak of the pooled control density is the
#'   negative population. Its spread is the root mean square of the events below
#'   the mode about the mode, and the threshold is `mode + qnorm(1 - tail_p) *
#'   spread`. This is the idea of the openCyto/flowDensity "tailgate" for
#'   cytokines. It runs in base R. (The `gate_tail()` function of older openCyto
#'   releases is not in current openCyto, so nothing is called there.)
#' * `"mindensity"`: `openCyto::gate_mindensity()` on the pooled control events,
#'   with `...` passed on. It needs openCyto and flowCore. Without openCyto the
#'   call warns and falls back to `"quantile"`; `method_used` records it. A
#'   valley gate only works when the pooled control has a visible positive mode.
#' * `"quantile"`: the `probs` quantile of the pooled control events (default
#'   0.999).
#' * `"manual"`: the named threshold vector in `thresholds`, on the transformed
#'   scale. It is also used when only `thresholds` is given.
#'
#' All control events are pooled; one gate does not follow donors or visits. The
#' gate depends on the control events only, never on stimulated ones.
#' Events are assumed to be already inside the parent population (for example
#' live CD4 T cells) and quality controlled.
#'
#' @param cells Data frame, one row per event, with `sample_id` and one
#'   fluorescence column per marker.
#' @param markers Marker names mapped to data-frame columns: a named character
#'   vector (`c(IFNG = "BV786-A")`) or, if unnamed, the column names.
#' @param samples Sample table with `sample_id` and `condition`, as in
#'   [cj_from_counts()].
#' @param method One of `"tail"`, `"mindensity"`, `"quantile"`, `"manual"`.
#' @param thresholds For `"manual"`: named numeric vector (one value per marker)
#'   on the transformed scale.
#' @param transform Required. `"none"` (cells already transformed), `"asinh"`
#'   (with `cofactor`), or a one-argument function. There is no default.
#' @param cofactor Positive number for `transform = "asinh"`: `asinh(x / cofactor)`.
#' @param control Values of `samples$condition` that mean unstimulated. They
#'   are recoded to `"control"` in [cj_from_cells()].
#' @param probs Quantile for `"quantile"`.
#' @param tail_p One-sided tail probability for `"tail"`; 0.001 puts the cut
#'   3.09 spreads above the mode.
#' @param ... Passed to `openCyto::gate_mindensity()` for `"mindensity"`.
#' @return A `cj_gate` list: `thresholds` (named, transformed scale), `markers`,
#'   `method`, `method_used` (per marker), `transform` (`fn` and `label`),
#'   `control`, `n_control` (events used), `control_values` (the pooled control
#'   events, transformed) and `notes`. Draw it with [cj_plot_gate()].
#' @export
#' @examples
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 400, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' g <- cj_gate(sim$cells, c("IFNG", "TNF", "IL2"), samples, method = "quantile",
#'              transform = "none", control = "unstim")
#' g
cj_gate <- function(cells, markers, samples, method = c("tail", "mindensity", "quantile", "manual"),
                    thresholds = NULL, transform = NULL, cofactor = NULL, control = "control",
                    probs = 0.999, tail_p = 0.001, ...) {
  stopifnot(is.data.frame(cells), is.data.frame(samples), "sample_id" %in% names(cells),
            all(c("sample_id", "condition") %in% names(samples)))
  if (missing(method) && !is.null(thresholds)) method <- "manual"
  method <- match.arg(method)
  if (method != "manual" && !is.null(thresholds))
    stop("thresholds are used only with method = \"manual\"", call. = FALSE)
  markers <- gate_markers(markers)
  tr <- gate_transform(transform, cofactor)
  samples <- recode_control(samples, control)
  is_ctrl <- as.character(cells$sample_id) %in% as.character(samples$sample_id[samples$condition == CONTROL_CONDITION])
  X <- gate_matrix(cells[is_ctrl, , drop = FALSE], markers, tr)
  notes <- character()

  if (method == "manual") {
    if (!is.numeric(thresholds) || is.null(names(thresholds)) ||
        !setequal(names(thresholds), names(markers)) || anyDuplicated(names(thresholds)))
      stop("thresholds must be a named numeric vector with one value per marker: ",
           paste(names(markers), collapse = ", "), call. = FALSE)
    th <- thresholds[names(markers)]
    used <- setNames(rep("manual", length(th)), names(th))
  } else {
    if (nrow(X) < 100L)
      stop("only ", nrow(X), " control events; a pooled gate needs at least 100", call. = FALSE)
    if (method == "mindensity" && !(has_pkg("openCyto") && has_pkg("flowCore"))) {
      warning("openCyto or flowCore is not available; using the ", probs, " control quantile instead",
              call. = FALSE)
      notes <- "mindensity requested but openCyto/flowCore missing: quantile used"
      eff <- "quantile"
    } else eff <- method
    th <- vapply(names(markers), function(m) switch(eff,
      tail = gate_cut_tail(X[, m], tail_p),
      quantile = stats::quantile(X[, m], probs, names = FALSE),
      mindensity = gate_cut_mindensity(X[, m], m, ...)), numeric(1))
    used <- setNames(rep(eff, length(th)), names(th))
  }
  structure(list(thresholds = th, markers = markers, method = method, method_used = used,
                 transform = tr, control = control, n_control = nrow(X), control_values = X,
                 notes = notes), class = "cj_gate")
}

#' @rdname cj_gate
#' @param x A `cj_gate`.
#' @export
print.cj_gate <- function(x, ...) {
  cat("cj_gate: ", length(x$thresholds), " markers, ", x$n_control, " pooled control events\n", sep = "")
  cat("Transform: ", x$transform$label, "\n", sep = "")
  print(data.frame(marker = names(x$thresholds), column = unname(x$markers[names(x$thresholds)]),
                   threshold = unname(x$thresholds), method = unname(x$method_used),
                   row.names = NULL))
  if (length(x$notes)) cat("Note: ", x$notes, "\n", sep = "")
  invisible(x)
}

#' Plot the pooled control density with its threshold
#'
#' One panel per marker: the density of the pooled control events on the
#' transformed scale and a vertical line at the threshold. Check every panel
#' before using the gate; a threshold inside the negative peak or far from the
#' tail is a gating error.
#'
#' @param gate A `cj_gate`.
#' @return A ggplot2 object.
#' @export
#' @examples
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 400, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' g <- cj_gate(sim$cells, c("IFNG", "TNF"), samples, method = "quantile",
#'              transform = "none", control = "unstim")
#' cj_plot_gate(g)
cj_plot_gate <- function(gate) {
  stopifnot(inherits(gate, "cj_gate"))
  X <- gate$control_values
  if (!NROW(X)) stop("this gate holds no control events to plot", call. = FALSE)
  m <- names(gate$thresholds)
  lv <- paste0(m, " (", gate$method_used[m], ")")
  d <- data.frame(marker = factor(rep(lv, each = nrow(X)), levels = lv), value = as.vector(X[, m]))
  th <- data.frame(marker = factor(lv, levels = lv), threshold = unname(gate$thresholds[m]))
  ggplot2::ggplot(d, ggplot2::aes(.data$value)) +
    ggplot2::geom_density(fill = "grey80", colour = "grey30") +
    ggplot2::geom_vline(data = th, ggplot2::aes(xintercept = .data$threshold), colour = "firebrick") +
    ggplot2::facet_wrap(~marker, scales = "free") +
    ggplot2::labs(x = paste0("Fluorescence (", gate$transform$label, ")"), y = "Density",
                  title = paste0("Pooled control events (n = ", gate$n_control, ")"))
}
