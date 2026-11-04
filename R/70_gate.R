# Pooled-control gating for the adoption API. One threshold per marker, found on
# all control (unstimulated) events pooled, on an explicitly stated transform.
# The result feeds cj_from_cells(). Vocabulary follows GLOSSARY.md.

# Every optional-package check goes through here, so tests can mock it. (Loading
# openCyto aborts R on some hosts because it pulls in flowWorkspace.)
has_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

GATE_METHODS <- c("tail", "mindensity", "quantile", "manual", "auto")

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
    ok <- is.numeric(cofactor) && length(cofactor) >= 1L && all(is.finite(cofactor)) && all(cofactor > 0)
    if (!ok) stop("transform = \"asinh\" needs a positive cofactor, e.g. cofactor = 5", call. = FALSE)
    if (length(cofactor) > 1L) {
      if (is.null(names(cofactor))) stop("a cofactor vector must be named by marker", call. = FALSE)
      return(list(fn = NULL, by_marker = lapply(cofactor, function(cf) function(x) asinh(x / cf)),
                  label = paste0("asinh(x / cofactor), cofactor per marker: ",
                                 paste0(names(cofactor), "=", signif(cofactor, 3), collapse = ", "))))
    }
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
  for (m in names(markers)) {
    fn <- if (is.null(transform$by_marker)) transform$fn else
      transform$by_marker[[m]] %||% stop("no cofactor for marker ", m, call. = FALSE)
    X[, m] <- as.numeric(fn(as.numeric(cells[[markers[[m]]]])))
  }
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
  pk <- gate_peak(x)
  pk[["mode"]] + stats::qnorm(1 - tail_p) * pk[["spread"]]
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
#' events pooled across donors, visits and samples, or one per batch with
#' `by`. The threshold is on the transformed scale, and the transform must be
#' stated.
#'
#' @details
#' Methods:
#' * `"tail"` (default): the main peak of the pooled control density is the
#'   negative population. Its spread is the root mean square of the events below
#'   the mode about the mode, and the threshold is `mode + qnorm(1 - tail_p) *
#'   spread`. This is the idea of the openCyto/flowDensity "tailgate" for
#'   cytokines. It runs in base R. (The `gate_tail()` function of older openCyto
#'   releases is not in current openCyto, so nothing is called there.) It assumes
#'   a symmetric negative peak; with a heavy upper tail the cut falls inside the tail.
#' * `"mindensity"`: `openCyto::gate_mindensity()` on the pooled control events,
#'   with `...` passed on. It needs openCyto and flowCore. Without openCyto the
#'   call warns and falls back to `"quantile"`; `method_used` records it. A
#'   valley gate only works when the pooled control has a visible positive mode.
#' * `"quantile"`: the `probs` quantile of the pooled control events (default
#'   0.999).
#' * `"manual"`: the named threshold vector in `thresholds`, on the transformed
#'   scale. It is also used when only `thresholds` is given.
#' * `"auto"`: per marker, `"tail"` or `"quantile"`, chosen from the shape of the
#'   pooled control distribution alone (see below). `method_used` and `$auto`
#'   record the choice. **Experimental.** It helps in simulation when the negative
#'   peak has a heavy upper tail, but the held-out comparison with expert gates did
#'   not support it (see the `from-cells` vignette): it lowered agreement on the
#'   markers where it switched. Use the `$auto` table as a diagnostic, and prefer
#'   `"tail"` with `by` and `cofactor = "auto"`.
#'
#' **The automatic rule.** The tail rule extrapolates a Gaussian from the lower
#' half of the negative peak. The statistic is the *tail ratio*: the distance
#' from the mode to the 99th percentile of the control events, divided by
#' `qnorm(0.99)` times the one-sided spread. It is about 1 for a symmetric peak and
#' grows when the upper tail is heavier than the lower half, whether from a skewed
#' negative population or from weakly positive events in the control. A marker
#' with a tail ratio above 1.5 gets the `probs` quantile (the tail extrapolation
#' is contradicted by the data), unless more than 3% of the control events lie
#' above the tail cut. That much is more likely a positive population than a
#' negative tail, the quantile would sit inside it, and the tail cut is kept, with
#' the reason in `$auto`. (3% is where, in simulation, the quantile stops beating
#' the tail cut.) Any other marker gets the tail cut. The 1.5 is
#' above the 99th percentile of the Gaussian null at 1000 or more events (1.46 at
#' 1000, 1.33 at 2000, 1.13 at 1e5 events, simulated), so a clean marker keeps the
#' tail cut. With fewer than 1000 control events for a marker the tail cut is
#' kept, with the reason in `$auto`. The constants were fixed from this null and
#' from simulated mixtures with clean and heavy tails, before any comparison with
#' expert gates; the simulator is [cj_simulate_gate_cells()]. The rule reads control events only. It does not handle a
#' control with a distinct positive mode (use `"mindensity"` there; the
#' `high_control_positive` flag warns).
#'
#' **Batches.** With `by`, the cut of every batch is found on that batch's
#' control events, with the same method, and applied to that batch's samples.
#' `thresholds` stays the pooled cut, the reference level; [cj_from_cells()]
#' shifts each event by its batch's cut minus the pooled cut, so downstream
#' code still compares one value with one threshold. A batch with fewer than
#' `min_batch_events` control events, or none, takes the pooled cut and the call
#' warns. For `"auto"`, the shape is judged on control events centred on each
#' batch's mode, so batch offsets do not look like a heavy tail.
#'
#' **QC flags** (`$qc`, shown by `print()` and [cj_plot_gate()]), from control
#' events only. Limits are in `qc_thresholds`.
#' * `cut_in_broad_tail`: the control density at the cut is above 5% of the peak
#'   density (a clean cut 3.09 spreads above the mode sits near 1%).
#' * `high_control_positive`: more than 1% of the control events lie above the
#'   tail-rule cut. A clean control has 0.1% by construction. Unstimulated
#'   cytokine-positive events, or a heavy negative tail, exceed this. (The
#'   percentage above the cut actually used is also in `$qc`; for a quantile cut
#'   it is fixed near `100 * (1 - probs)`.)
#' * `batch_cut_spread`: with `by`, the range of the batch cuts is more than one
#'   negative-peak spread. The batches then disagree about where the negative
#'   population ends, so a pooled cut would be wrong for some of them.
#'
#' All control events are pooled within a batch (or entirely, without `by`). The
#' gate depends on the control events only, never on stimulated ones.
#' Events are assumed to be already inside the parent population (for example
#' live CD4 T cells) and quality controlled.
#'
#' @param cells Data frame, one row per event, with `sample_id` and one
#'   fluorescence column per marker.
#' @param markers Marker names mapped to data-frame columns: a named character
#'   vector (`c(IFNG = "BV786-A")`) or, if unnamed, the column names.
#' @param samples Sample table with `sample_id` and `condition`, as in
#'   [cj_from_counts()]; plus the `by` column if used.
#' @param method One of `"tail"`, `"mindensity"`, `"quantile"`, `"manual"`, `"auto"`.
#' @param thresholds For `"manual"`: named numeric vector (one value per marker)
#'   on the transformed scale.
#' @param transform Required. `"none"` (cells already transformed), `"asinh"`
#'   (with `cofactor`), or a one-argument function. There is no default.
#' @param cofactor Positive number for `transform = "asinh"`: `asinh(x / cofactor)`.
#'   A vector named by marker gives one cofactor per marker. `"auto"` takes it
#'   from the width of the negative peak of the raw control events ([cj_cofactor()]).
#' @param control Values of `samples$condition` that mean unstimulated. They
#'   are recoded to `"control"` in [cj_from_cells()].
#' @param probs Quantile for `"quantile"`.
#' @param tail_p One-sided tail probability for `"tail"`; 0.001 puts the cut
#'   3.09 spreads above the mode.
#' @param by `NULL` (default, one pooled gate) or the name of a `samples` column
#'   (for example a staining batch or screen). Not used with `"manual"`.
#' @param min_batch_events Minimum control events per batch (default 1000).
#' @param qc_thresholds Named numeric vector overriding any of
#'   `max_control_pos` (0.01), `tail_density` (0.05) and `batch_spread` (1).
#' @param ... Passed to `openCyto::gate_mindensity()` for `"mindensity"`.
#' @return A `cj_gate` list: `thresholds` (named, transformed scale), `markers`,
#'   `method`, `method_used` (per marker), `transform` (`fn` and `label`),
#'   `control`, `n_control` (events used), `control_values` (the pooled control
#'   events, transformed), `notes`, `qc` (one row per marker: flags and the
#'   numbers behind them), `auto` (tail ratio and choice per marker, for
#'   `"auto"`), and for `by`: `by`, `batch_thresholds` (batch by marker),
#'   `batch_of_sample`, `batch_n_control`, `control_batch`. Draw it with [cj_plot_gate()].
#' @export
#' @examples
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 400, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' g <- cj_gate(sim$cells, c("IFNG", "TNF", "IL2"), samples, method = "quantile",
#'              transform = "none", control = "unstim")
#' g
#'
#' # Heavy-tailed marker: "auto" switches that marker to the control quantile.
#' s <- cj_simulate_gate_cells(c("clean", "skewed"), c("A", "B"), seed = 1)
#' cj_gate(s$cells, c("A", "B"), s$samples, method = "auto", transform = "none",
#'         control = "unstim")
cj_gate <- function(cells, markers, samples, method = c("tail", "mindensity", "quantile", "manual", "auto"),
                    thresholds = NULL, transform = NULL, cofactor = NULL, control = "control",
                    probs = 0.999, tail_p = 0.001, by = NULL, min_batch_events = 1000L,
                    qc_thresholds = NULL, ...) {
  stopifnot(is.data.frame(cells), is.data.frame(samples), "sample_id" %in% names(cells),
            all(c("sample_id", "condition") %in% names(samples)))
  if (missing(method) && !is.null(thresholds)) method <- "manual"
  method <- match.arg(method)
  if (method != "manual" && !is.null(thresholds))
    stop("thresholds are used only with method = \"manual\"", call. = FALSE)
  if (method == "manual" && !is.null(by))
    stop("by is not used with method = \"manual\": the thresholds are given", call. = FALSE)
  markers <- gate_markers(markers)
  qc <- GATE_QC
  if (!is.null(qc_thresholds)) {
    if (!is.numeric(qc_thresholds) || !all(names(qc_thresholds) %in% names(GATE_QC)))
      stop("qc_thresholds must be a named numeric vector with names in: ", paste(names(GATE_QC), collapse = ", "),
           call. = FALSE)
    qc[names(qc_thresholds)] <- qc_thresholds
  }
  samples <- recode_control(samples, control)
  is_ctrl <- as.character(cells$sample_id) %in% as.character(samples$sample_id[samples$condition == CONTROL_CONDITION])
  if (identical(cofactor, "auto")) {
    if (!identical(transform, "asinh")) stop("cofactor = \"auto\" needs transform = \"asinh\"", call. = FALSE)
    if (sum(is_ctrl) < 100L) stop("only ", sum(is_ctrl), " control events; a pooled gate needs at least 100", call. = FALSE)
    cofactor <- vapply(names(markers), function(m)
      gate_cofactor(as.numeric(cells[[markers[[m]]]][is_ctrl])), numeric(1))
  }
  tr <- gate_transform(transform, cofactor)
  X <- gate_matrix(cells[is_ctrl, , drop = FALSE], markers, tr)
  notes <- character()
  bmap <- bev <- NULL
  if (!is.null(by)) {
    bmap <- gate_batch_map(samples, by)
    bev <- unname(bmap[as.character(cells$sample_id[is_ctrl])])
    if (anyNA(bev)) stop("cells has sample_id not in samples", call. = FALSE)
  }
  auto <- batch <- NULL

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
      eff <- setNames(rep("quantile", ncol(X)), colnames(X))
    } else if (method == "auto") {
      xc <- X
      if (!is.null(bev)) {
        enough <- bev %in% names(which(table(bev) >= min_batch_events))
        if (any(enough)) {
          xc <- X[enough, , drop = FALSE]
          for (b in unique(bev[enough])) {
            i <- bev[enough] == b
            xc[i, ] <- sweep(xc[i, , drop = FALSE], 2L, apply(xc[i, , drop = FALSE], 2L, function(v)
              gate_peak(v)[["mode"]]))
          }
        }
      }
      auto <- gate_auto_choose(xc, tail_p)
      eff <- setNames(auto$choice, auto$marker)
    } else eff <- setNames(rep(method, ncol(X)), colnames(X))
    th <- vapply(names(markers), function(m) gate_cut(eff[[m]], X[, m], m, probs, tail_p, ...), numeric(1))
    used <- eff[names(th)]
    if (!is.null(bev)) {
      all_b <- sort(unique(unname(bmap)))
      batch <- gate_batch_cuts(X, bev, all_b, th, eff, min_batch_events, probs, tail_p, ...)
      few <- all_b[!all_b %in% batch$ok]
      if (length(few))
        warning(length(few), " batch(es) have fewer than ", min_batch_events, " control events (",
                paste0(few, ": ", batch$n[few], collapse = ", "), "); they use the pooled cut",
                call. = FALSE)
      if (length(few)) notes <- c(notes, paste0("batches on the pooled cut (too few control events): ",
                                                 paste(few, collapse = ", ")))
    }
  }

  # QC flags: cut and tail-rule cut that apply to each control event.
  cut_ev <- matrix(rep(th, each = nrow(X)), nrow(X), dimnames = list(NULL, names(th)))
  tail_pooled <- vapply(names(markers), function(m)
    tryCatch(gate_cut_tail(X[, m], tail_p), error = function(e) NA_real_), numeric(1))
  tailcut_ev <- matrix(rep(tail_pooled, each = nrow(X)), nrow(X), dimnames = list(NULL, names(th)))
  if (!is.null(batch)) {
    cut_ev[] <- batch$cut[bev, names(th)]
    tc <- batch$tailcut
    tc[is.na(tc)] <- matrix(rep(tail_pooled, each = nrow(tc)), nrow(tc))[is.na(tc)]
    tailcut_ev[] <- tc[bev, names(th)]
  }
  qcdf <- gate_qc(X, cut_ev, tailcut_ev, batch, qc)
  qcdf$flags <- gate_flag_text(qcdf)
  out <- list(thresholds = th, markers = markers, method = method, method_used = used,
              transform = tr, control = control, n_control = nrow(X), control_values = X,
              notes = notes, qc = qcdf, auto = auto, qc_thresholds = qc)
  if (!is.null(batch)) {
    out$by <- by
    out$batch_thresholds <- batch$cut[, names(th), drop = FALSE]
    out$batch_of_sample <- bmap
    out$batch_n_control <- batch$n
    out$control_batch <- bev
  }
  structure(out, class = "cj_gate")
}

#' @rdname cj_gate
#' @param x A `cj_gate`.
#' @export
print.cj_gate <- function(x, ...) {
  cat("cj_gate: ", length(x$thresholds), " markers, ", x$n_control, " pooled control events",
      if (!is.null(x$by)) paste0(", cuts per ", x$by, " (", nrow(x$batch_thresholds), " batches)"), "\n", sep = "")
  cat("Transform: ", x$transform$label, "\n", sep = "")
  d <- data.frame(marker = names(x$thresholds), column = unname(x$markers[names(x$thresholds)]),
                  threshold = unname(x$thresholds), method = unname(x$method_used), row.names = NULL)
  if (!is.null(x$qc)) {
    q <- x$qc[match(d$marker, x$qc$marker), ]
    d$pct_ctrl_pos <- signif(q$pct_control_pos_tailrule, 3)
    if (!is.null(x$by)) d$batch_spread <- signif(q$batch_cut_spread, 3)
    d$flags <- q$flags
  }
  print(d)
  if (!is.null(x$qc) && any(nzchar(x$qc$flags)))
    cat("QC: pct_ctrl_pos is the % of control events above the tail-rule cut; flags are listed per marker.\n")
  if (length(x$notes)) cat("Note: ", x$notes, "\n", sep = "")
  invisible(x)
}

#' Plot the pooled control density with its threshold
#'
#' One panel per marker: the density of the pooled control events on the
#' transformed scale and a vertical line at the threshold. Dashed lines are the
#' per-batch cuts of a gate made with `by`. The panel title lists the method and
#' any QC flag (`cut_in_broad_tail`, `high_control_positive`,
#' `batch_cut_spread`). Check every panel before using the gate; a threshold
#' inside the negative peak or far from the tail is a gating error.
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
  fl <- if (is.null(gate$qc)) rep("", length(m)) else gate$qc$flags[match(m, gate$qc$marker)]
  lv <- paste0(m, " (", gate$method_used[m], ")", ifelse(nzchar(fl), paste0("\n", fl), ""))
  d <- data.frame(marker = factor(rep(lv, each = nrow(X)), levels = lv), value = as.vector(X[, m]))
  th <- data.frame(marker = factor(lv, levels = lv), threshold = unname(gate$thresholds[m]))
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$value)) +
    ggplot2::geom_density(fill = "grey80", colour = "grey30") +
    ggplot2::geom_vline(data = th, ggplot2::aes(xintercept = .data$threshold), colour = "firebrick")
  if (!is.null(gate$batch_thresholds)) {
    bt <- gate$batch_thresholds[, m, drop = FALSE]
    bd <- data.frame(marker = factor(rep(lv, each = nrow(bt)), levels = lv), threshold = as.vector(bt))
    p <- p + ggplot2::geom_vline(data = bd, ggplot2::aes(xintercept = .data$threshold),
                                 colour = "firebrick", linetype = "dashed", alpha = 0.6)
  }
  p + ggplot2::facet_wrap(~marker, scales = "free") +
    ggplot2::labs(x = paste0("Fluorescence (", gate$transform$label, ")"), y = "Density",
                  title = paste0("Pooled control events (n = ", gate$n_control, ")"))
}
