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
#' `magnitude_model = "c0eb"` (default) is an empirical-Bayes pooled conjugate
#' model, no Stan, 3 to 5 seconds for 24 donors. Per donor and visit it models the
#' any-positive fraction of the control, `q`, and the fold change `rho` of each
#' stimulation, both log-normal across donors. The across-donor means and
#' standard deviations are fitted from the data by marginal likelihood on a
#' grid, per parent population and stimulation, with a floor of 0.2 on the log
#' scale. The stimulated counts also inform `q`, so the baseline, which drives
#' the true RMS when responses are rare, is estimated from both samples. Within
#' the positive patterns each sample keeps C0's Dirichlet(`a0` + counts). Extra
#' arguments (`sd_min`, `pool_patterns`, `null_mass`) go in `fit_args`. A parent
#' population with fewer than 8 donor-visit blocks gets a warning: the
#' across-donor hyperparameters are then poorly determined (use `"c0"`).
#' Every control sample of a donor and visit measures the same baseline and
#' shares one `q` draw.
#'
#' Because the fold change is pooled, a donor with no positive cell in either
#' sample gets an interval that reflects the cohort, not the donor. When the
#' cohort responds strongly the interval can lie wholly above 0 (for example
#' 3e-5 to 3e-4 for a 10,000-cell sample), because the model says "if this donor
#' responds like the others, the effect is about this small". Its lower end stays
#' below the default `magnitude_floor` of 1e-4, so the donor is not
#' `magnitude_positive` and no `magnitude_without_response` flag is raised; with
#' `magnitude_floor = 0` it would be. A mixture with a point mass at no response
#' (`fit_args = list(null_mass = TRUE)`) did not change this when every donor
#' responds (the fitted mass was 0.001), and in dev simulations it cut false
#' flags at the cost of power (rare, weak: true positives 0.64 to 0.53, false
#' flags among non-responders 0.19 to 0.10, pooled Spearman -0.02), so it is off.
#'
#' `"c0"` is the plain conjugate model: an independent Dirichlet(`a0` + counts)
#' per sample, no pooling. `"m1dm"` calls `cj_fit(model = "m1dm", level =
#' "sample")`, which pools donors in Stan and needs CmdStan (14 to 24 seconds for
#' 24 donors; arguments in `fit_args`, for example `level = "structured"`). A
#' `cj_fit` given as `x` is used as it is for the magnitudes, and
#' `magnitude_model` is ignored.
#'
#' In the benchmark simulations (24 donors, 10,000 parent cells, 3 or 4 markers,
#' seeds 11 to 20, held out), Spearman with the true donor RMS, pooled over
#' scenarios and seeds, was:
#'
#' | scenario | COMPASS FS | C0 | C0-EB (default) | m1dm |
#' |---|---|---|---|---|
#' | rare, weak | 0.50 | 0.69 | 0.82 | 0.82 |
#' | moderate | 0.72 | 0.98 | 0.99 | 0.99 |
#'
#' Within one replicate at the rare, weak effect, the gain over FS was +0.06
#' (95% interval -0.00 to 0.12 over seeds) for C0, +0.21 (0.15 to 0.27) for
#' C0-EB and +0.20 (0.15 to 0.26) for m1dm. The 90% intervals covered the true
#' RMS in 0.85 to 0.94 of donor contrasts for C0-EB (0.89 to 0.94 at the rare,
#' weak effect), 0.92 to 0.96 for C0 and 0.91 to 0.98 for m1dm. The intervals
#' describe the sampled donor-visit: against a structured truth with
#' sample-level noise (generator B) they covered 0.43. **CPS from C0 or C0-EB is
#' not reliable** (dev Spearman -0.6 to 0.7, coverage 0.2 to 1.00 across
#' scenarios): take CPS from `"m1dm"`. NICR and Poly-NICR bootstrap intervals
#' cover cell sampling within a sample only.
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
#' @param magnitude_model `"c0eb"` (default, no Stan), `"c0"` (no Stan) or `"m1dm"`.
#' @param missing_control Passed to [cj_run_compass()]. Rows borrowed from
#'   another control get COMPASS numbers but no magnitude (`"not_assessed"`).
#' @param level Interval level for the magnitudes (and the NICR bootstrap).
#' @param epsilon Named score thresholds for the `*_p_positive` columns, as in
#'   [cj_scores()]. Default: RMS and marker magnitude 1e-4, CPS 0.01.
#' @param response_threshold COMPASS response probability above which a donor
#'   is called a responder.
#' @param magnitude_floor The RMS interval must lie above this to count as a
#'   positive magnitude. Default `NULL`: the RMS entry of `epsilon` (1e-4), so
#'   "positive magnitude" is P(RMS > epsilon) > 0.95 for a 90% interval, the
#'   benchmark's detection rule. A floor of 0 flags more donors, including
#'   non-responders whose sampled stimulated and control wells differ by chance.
#' @param n_boot Cell bootstrap replicates for NICR and Poly-NICR. `0` gives
#'   point estimates only. Needs a cell data frame (or a `cj_data` built with
#'   `keep_events = TRUE`).
#' @param a0 Dirichlet prior pseudo-count of C0 and C0-EB.
#' @param n_draws Posterior draws of C0 and C0-EB.
#' @param seed Seed for C0, C0-EB and the bootstrap.
#' @param fit_args List of further arguments to [cj_fit()] when
#'   `magnitude_model = "m1dm"`, or to the C0-EB fit (`sd_min`, `pool_patterns`).
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
                              control = "control", magnitude_model = c("c0eb", "c0", "m1dm"),
                              missing_control = c("drop", "donor", "population"),
                              level = 0.9, epsilon = REPORT_EPSILON,
                              response_threshold = 0.5, magnitude_floor = NULL,
                              n_boot = 0L, a0 = 0.5, n_draws = 4000L, seed = 1L,
                              fit_args = list(), compass = NULL, compass_args = list(),
                              gate_args = list()) {
  magnitude_model <- match.arg(magnitude_model)
  epsilon <- unlist(utils::modifyList(as.list(REPORT_EPSILON), as.list(epsilon)))
  if (is.null(magnitude_floor)) magnitude_floor <- epsilon[["RMS"]]
  if (inherits(x, "cj_fit") && (!missing(magnitude_model) || !missing(a0) || !missing(n_draws) || length(fit_args)))
    warning("x is a cj_fit: magnitude_model, a0, n_draws and fit_args are ignored", call. = FALSE)
  else if (magnitude_model == "m1dm" && (!missing(a0) || !missing(n_draws)))
    warning("a0 and n_draws are ignored by magnitude_model = \"m1dm\"", call. = FALSE)
  else if (magnitude_model == "c0" && length(fit_args))
    warning("fit_args is ignored by magnitude_model = \"c0\"", call. = FALSE)
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
  } else if (magnitude_model == "c0eb") {
    do.call(cj_fit_dirichlet_eb, c(list(d, a0 = a0, n_draws = n_draws, seed = seed), fit_args))
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
