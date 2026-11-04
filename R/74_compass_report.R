# COMPASS calls with cytoJointR magnitudes beside them: one table per contrast.
# Docs: vignettes/from-cells.Rmd ("Magnitudes beside COMPASS calls").

# Score thresholds for the epsilon probabilities in the report. They change only
# `*_p_positive`; the estimates and intervals do not depend on them. RMS and
# marker magnitude use 1e-4 (one cell in 10,000, the benchmark's epsilon), CPS 0.01.
REPORT_EPSILON <- c(RMS = 1e-4, poly_absolute = 1e-4, poly_marginal = 1e-4, CPS = 0.01,
                    CCS = 1e-4, breadth = 0.01, marker_magnitude = 1e-4)

CONCORDANCE_LEVELS <- c("response_with_magnitude", "no_response_no_magnitude",
                        "response_without_magnitude", "magnitude_without_response",
                        "not_assessed")

# Wide RMS/CPS columns from a cj_scores() table.
report_global <- function(sc, which) {
  out <- NULL
  for (s in which) {
    r <- sc[sc$score == s & is.na(sc$marker), ]
    r <- data.frame(contrast_id = r$contrast_id, estimate = r$estimate, lower = r$lower,
                    upper = r$upper, p_positive = r$p_positive, status = r$status)
    names(r)[-1] <- paste0(tolower(s), c("", "_lower", "_upper", "_p_positive", "_status"))
    out <- if (is.null(out)) r else merge(out, r, by = "contrast_id", all = TRUE, sort = FALSE)
  }
  out
}

#' COMPASS calls with cytoJointR magnitudes and intervals
#'
#' COMPASS says who responds. This table adds how large the response is, with
#' an interval, so one row answers both questions. Per donor, stimulation,
#' visit and parent population it joins:
#'
#' 1. COMPASS `compass_FS`, `compass_PFS`, the response probabilities and
#'    `control_source` (from [cj_run_compass()]);
#' 2. cytoJointR RMS and CPS with credible intervals, and per-marker
#'    magnitudes in `attr(, "markers")`;
#' 3. when `x` is a cell data frame (so fluorescence is available), NICR and
#'    Poly-NICR per marker with bootstrap intervals, also in `attr(, "markers")`;
#' 4. a `concordance` column that compares the COMPASS call with the RMS
#'    interval.
#'
#' @section Magnitude model:
#' `magnitude_model = "c0"` (default) is the conjugate model
#' C0: an independent Dirichlet(`a0` + counts) per
#' sample, no Stan, no pooling. It is fast and, in the benchmark simulations,
#' ranks donor RMS far better than COMPASS FS. Its intervals describe the
#' sampled donor-visit only and are narrow when the response is rare: they are
#' not a donor-population interval. `"m1dm"` calls `cj_fit(model = "m1dm",
#' level = "sample")`, which pools donors and needs CmdStan; give extra
#' arguments in `fit_args`. A `cj_fit` given as `x` is used as it is (for the
#' magnitudes) and `magnitude_model` is ignored.
#'
#' @section Concordance:
#' The COMPASS call is `compass_response_degree1 > response_threshold`
#' (`compass_call`). The magnitude is positive when the lower end of the RMS
#' interval exceeds `magnitude_floor` (`magnitude_positive`). `concordance` is
#' one of:
#'
#' * `"response_with_magnitude"`: call and positive magnitude;
#' * `"no_response_no_magnitude"`: neither;
#' * `"response_without_magnitude"`: COMPASS calls a response, the RMS interval
#'   includes `magnitude_floor`;
#' * `"magnitude_without_response"`: no COMPASS call, the RMS interval lies
#'   above `magnitude_floor`;
#' * `"not_assessed"`: the call or the magnitude is missing (for example a
#'   contrast without a control sample, where no magnitude is defined, or a
#'   failed COMPASS run).
#'
#' `discordant` is `TRUE` for the two mismatch classes. A discordant row asks
#' for a look, not a verdict: COMPASS's call and a posterior interval answer
#' different questions.
#'
#' @param x A `cj_data`, a `cj_fit`, or a cell data frame (then `samples` is
#'   required, as in [cj_run_compass()]).
#' @param markers,samples,gate,transform,control Used only when `x` is a cell
#'   data frame; see [cj_from_cells()].
#' @param magnitude_model `"c0"` (default, no Stan) or `"m1dm"`.
#' @param missing_control Passed to [cj_run_compass()]. Rows borrowed from
#'   another control get COMPASS numbers but no magnitude (`"not_assessed"`).
#' @param level Interval level for the magnitudes (and the NICR bootstrap).
#' @param epsilon Named score thresholds for the `*_p_positive` columns, as in
#'   [cj_scores()]. Default: RMS and marker magnitude 1e-4, CPS 0.01.
#' @param response_threshold COMPASS response probability above which a donor
#'   is called a responder.
#' @param magnitude_floor The RMS interval must lie above this to count as a
#'   positive magnitude.
#' @param n_boot Cell bootstrap replicates for NICR and Poly-NICR. `0` gives
#'   point estimates only. Needs a cell data frame (or a `cj_data` built with
#'   `keep_events = TRUE`).
#' @param a0 Dirichlet prior pseudo-count of C0.
#' @param n_draws Posterior draws of C0.
#' @param seed Seed for C0 and the bootstrap.
#' @param fit_args List of further arguments to [cj_fit()] when
#'   `magnitude_model = "m1dm"`.
#' @param compass Optional result of [cj_run_compass()] on the same data, to
#'   skip running COMPASS again.
#' @param compass_args List of further arguments to [cj_run_compass()]
#'   (for example `iterations`, `replications`).
#' @param gate_args List of extra arguments to [cj_gate()] for a cell data frame.
#' @return A data frame of class `cj_compass_report`, one row per contrast:
#'   the contrast keys; `compass_FS`, `compass_PFS`, `compass_response_degree1`,
#'   `compass_response_degree2`, `control_source` and the other COMPASS columns;
#'   `compass_call`; `rms`, `rms_lower`, `rms_upper`, `rms_p_positive`,
#'   `rms_status` and the same for `cps`; `magnitude_model`;
#'   `magnitude_positive`; `concordance`; `discordant`. Attributes: `markers`
#'   (one row per contrast and marker: `marker_magnitude` with interval, and,
#'   from cells, `nicr`, `poly_nicr` and their `_lo`/`_hi` bootstrap bounds),
#'   `settings`, and `compass` (the [cj_run_compass()] table).
#' @export
#' @examplesIf requireNamespace("COMPASS", quietly = TRUE)
#' \donttest{
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 400, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' samples$visit <- factor("v1")
#' samples$parent_population <- "CD4"
#' rep <- cj_compass_report(sim$cells, c("IFNG", "TNF", "IL2"), samples,
#'                          transform = "none", control = "unstim",
#'                          gate_args = list(method = "quantile"), n_boot = 20,
#'                          compass_args = list(iterations = 200, replications = 2))
#' rep[c("donor_id", "compass_response_degree1", "rms", "rms_lower", "rms_upper", "concordance")]
#' }
cj_compass_report <- function(x, markers = NULL, samples = NULL, gate = NULL, transform = NULL,
                              control = "control", magnitude_model = c("c0", "m1dm"),
                              missing_control = c("drop", "donor", "population"),
                              level = 0.9, epsilon = REPORT_EPSILON,
                              response_threshold = 0.5, magnitude_floor = 0,
                              n_boot = 0L, a0 = 0.5, n_draws = 4000L, seed = 1L,
                              fit_args = list(), compass = NULL, compass_args = list(),
                              gate_args = list()) {
  magnitude_model <- match.arg(magnitude_model)
  missing_control <- match.arg(missing_control)
  is_fit <- inherits(x, "cj_fit")
  from_cells <- is.data.frame(x)
  if (from_cells && n_boot > 0) gate_args <- c(list(keep_events = TRUE), gate_args)
  d <- do.call(as_cj_data, c(list(x, markers, samples, gate, transform, control), gate_args))
  if (is.null(compass)) {
    compass <- do.call(cj_run_compass, c(list(d, missing_control = missing_control), compass_args))
  }
  comp <- as.data.frame(compass)

  # Magnitudes ----------------------------------------------------------------
  draws <- if (is_fit) {
    x
  } else if (magnitude_model == "c0") {
    cj_fit_dirichlet(d, a0 = a0, n_draws = n_draws, seed = seed)
  } else {
    do.call(cj_fit, c(list(d, model = "m1dm", level = "sample", seed = seed), fit_args))
  }
  used_model <- if (is_fit) "cj_fit" else magnitude_model
  sc <- suppressMessages(cj_scores(draws, epsilon, level = level, ccs_gate = NULL))
  mk <- cj_markers(draws, epsilon[["marker_magnitude"]], level = level)
  mag <- report_global(sc, c("RMS", "CPS"))

  # Join to the COMPASS table, which keeps one row per contrast it ran -------
  keys <- d$contrasts[, c("contrast_id", "donor_id", "visit", "parent_population", "stimulation")]
  comp_cols <- setdiff(names(comp), c("donor_id", "visit", "parent_population", "stimulation"))
  out <- merge(keys, comp[comp_cols], by = "contrast_id", all.x = TRUE, sort = FALSE)
  out <- merge(out, mag, by = "contrast_id", all.x = TRUE, sort = FALSE)
  out <- out[match(keys$contrast_id, out$contrast_id), ]
  stopifnot("report rows match contrasts" = identical(out$contrast_id, keys$contrast_id))
  rownames(out) <- NULL

  out$compass_call <- out$compass_response_degree1 > response_threshold
  out$magnitude_model <- used_model
  out$magnitude_positive <- out$rms_lower > magnitude_floor
  out$concordance <- factor(
    ifelse(is.na(out$compass_call) | is.na(out$magnitude_positive), "not_assessed",
      ifelse(out$compass_call & out$magnitude_positive, "response_with_magnitude",
        ifelse(!out$compass_call & !out$magnitude_positive, "no_response_no_magnitude",
          ifelse(out$compass_call, "response_without_magnitude", "magnitude_without_response")))),
    levels = CONCORDANCE_LEVELS)
  out$discordant <- ifelse(out$concordance == "not_assessed", NA,
                           out$concordance %in% c("response_without_magnitude",
                                                  "magnitude_without_response"))

  # Per-marker table: marker magnitude, plus NICR from cells -------------------
  mtab <- data.frame(contrast_id = mk$contrast_id, marker = as.character(mk$marker),
                     marker_magnitude = mk$estimate, marker_magnitude_lower = mk$lower,
                     marker_magnitude_upper = mk$upper)
  if (!is.null(d$intensity)) {
    fl <- suppressMessages(cj_fluorescence_scores(d, n_boot = n_boot, level = level, seed = seed))
    keep <- intersect(c("contrast_id", "marker", "nicr", "nicr_lo", "nicr_hi",
                        "poly_nicr", "poly_nicr_lo", "poly_nicr_hi"), names(fl))
    mtab <- merge(mtab, fl[keep], by = c("contrast_id", "marker"), all.x = TRUE, sort = FALSE)
  }
  mtab <- mtab[order(match(mtab$contrast_id, keys$contrast_id), match(mtab$marker, d$markers)), ]
  rownames(mtab) <- NULL

  attr(out, "markers") <- mtab
  attr(out, "settings") <- list(magnitude_model = used_model, level = level, epsilon = epsilon,
                                response_threshold = response_threshold,
                                magnitude_floor = magnitude_floor, n_boot = n_boot, seed = seed)
  attr(out, "compass") <- comp
  class(out) <- c("cj_compass_report", "data.frame")
  out
}

#' @export
print.cj_compass_report <- function(x, ...) {
  s <- attr(x, "settings")
  cat("<cj_compass_report> ", nrow(x), " contrast(s); magnitude model ", s$magnitude_model,
      ", ", 100 * s$level, "% intervals; COMPASS call = response > ", s$response_threshold, "\n", sep = "")
  print(table(concordance = x$concordance), ...)
  invisible(x)
}

#' Forest plot of COMPASS calls and magnitudes
#'
#' One row per contrast. The position is the cytoJointR magnitude with its
#' credible interval; the point colour is the COMPASS response probability
#' (the midpoint is the call threshold). Discordant rows (see
#' [cj_compass_report()]) are drawn as triangles.
#'
#' @param report Output of [cj_compass_report()].
#' @param score `"rms"` (default) or `"cps"`.
#' @param order Order rows by `"magnitude"` (default), or keep the table order.
#' @param max_rows Show at most this many contrasts (the largest magnitudes
#'   first); `NULL` shows all.
#' @return A ggplot object.
#' @export
#' @examples
#' rep <- data.frame(contrast_id = paste0("d", 1:6), donor_id = paste0("d", 1:6),
#'                   compass_response_degree1 = c(0.9, 0.8, 0.6, 0.4, 0.2, 0.05),
#'                   rms = c(0.02, 0.015, 0.001, 0.004, 0.0002, 0),
#'                   rms_lower = c(0.015, 0.01, -0.0005, 0.002, -0.0001, -0.0002),
#'                   rms_upper = c(0.025, 0.02, 0.003, 0.006, 0.0006, 0.0002),
#'                   discordant = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE))
#' attr(rep, "settings") <- list(response_threshold = 0.5, level = 0.9,
#'                               magnitude_model = "c0")
#' cj_plot_compass_report(rep)
cj_plot_compass_report <- function(report, score = c("rms", "cps"), order = c("magnitude", "table"),
                                   max_rows = NULL) {
  score <- match.arg(score)
  order <- match.arg(order)
  s <- attr(report, "settings")
  d <- data.frame(label = report$contrast_id, response = report$compass_response_degree1,
                  estimate = report[[score]], lower = report[[paste0(score, "_lower")]],
                  upper = report[[paste0(score, "_upper")]],
                  discordant = report$discordant %in% TRUE)
  d <- d[!is.na(d$estimate), ]
  if (order == "magnitude") d <- d[order(d$estimate, decreasing = TRUE), ]
  if (!is.null(max_rows)) d <- utils::head(d, max_rows)
  d$label <- factor(d$label, levels = rev(unique(d$label)))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$estimate, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey40") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
                           orientation = "y", width = 0.25, colour = "grey55") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$response, shape = .data$discordant), size = 2.4) +
    ggplot2::scale_colour_gradient2(low = "#2c7bb6", mid = "grey85", high = "#d7191c",
                                    midpoint = s$response_threshold %||% 0.5, limits = c(0, 1),
                                    name = "COMPASS\nresponse") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = "discordant") +
    ggplot2::labs(x = paste0(toupper(score), " (", 100 * (s$level %||% 0.9), "% interval; ",
                             s$magnitude_model %||% "model", ")"),
                  y = NULL)
}
