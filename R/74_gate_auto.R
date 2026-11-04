# Helpers for cj_gate(): negative-peak statistics, the automatic tail/quantile
# choice, the asinh cofactor rule, per-batch cuts, QC flags, and a simulator of
# control distributions with known truth. Everything here is label-free: it reads
# control events only, never stimulated events and never an expert gate.

# A priori constants (fixed before any real-data comparison; see ?cj_gate).
# ratio:      the 1.5 sits above the 99th percentile of the Gaussian null of the tail
#             ratio for n >= 1000 events (1.46 at n = 1000, 1.33 at 2000, 1.13 at 1e5).
# max_excess: the share of control events above the tail-rule cut. The 0.999 quantile
#             tolerates about 0.1% positives in the control; beyond 3% the surplus is a
#             population, not a tail, and the quantile sits inside it. Simulation (dev,
#             scripts/98_gating-dev-grid.R) puts the crossover where the quantile stops
#             beating the tail cut at 2.5% to 6% depending on how far the positives lie.
# min_events: below this the 99th percentile is too noisy to judge a shape.
# COFACTOR_K: the cut (mode + 3.09 spreads) then lies near asinh argument 1.
GATE_AUTO <- c(ratio = 1.5, max_excess = 0.03, min_events = 1000, prob = 0.99)
GATE_COFACTOR_K <- 3
GATE_QC <- c(max_control_pos = 0.01, tail_density = 0.05, batch_spread = 1)

# Mode and one-sided spread of the main peak. The spread is the root mean square
# of the events below the mode about the mode, so positive events in the upper
# tail do not enter it.
gate_peak <- function(x) {
  d <- stats::density(x)
  m <- d$x[which.max(d$y)]
  low <- x[x <= m]
  if (length(low) < 10L) stop("too few control events below the density mode to estimate the spread", call. = FALSE)
  s <- sqrt(mean((low - m)^2))
  if (!is.finite(s) || s <= 0) stop("the control events have no spread below the mode", call. = FALSE)
  c(mode = m, spread = s)
}

# Label-free shape statistic: how far the 99th percentile lies above the mode,
# in units of what the lower-half spread predicts for a symmetric (Gaussian) peak.
# About 1 for a clean peak; above 1 when the upper tail is heavier than the lower half.
gate_tail_ratio <- function(x, prob = GATE_AUTO[["prob"]]) {
  pk <- gate_peak(x)
  (stats::quantile(x, prob, names = FALSE) - pk[["mode"]]) / (stats::qnorm(prob) * pk[["spread"]])
}

# Cut by an effective method on one vector of control events.
gate_cut <- function(eff, x, marker, probs, tail_p, ...) {
  switch(eff,
    tail = gate_cut_tail(x, tail_p),
    quantile = stats::quantile(x, probs, names = FALSE),
    mindensity = gate_cut_mindensity(x, marker, ...))
}

# Per-marker choice between "tail" and "quantile" (rule in ?cj_gate). `xc` are control events
# (batch-centred when gating by batch, so offsets do not look like a heavy tail).
gate_auto_choose <- function(xc, tail_p) {
  out <- lapply(colnames(xc), function(m) {
    x <- xc[, m]
    n <- sum(is.finite(x))
    row <- function(r, ex, choice, reason)
      data.frame(marker = m, n = n, tail_ratio = r, excess = ex, choice = choice, reason = reason,
                 stringsAsFactors = FALSE)
    if (n < GATE_AUTO[["min_events"]])
      return(row(NA_real_, NA_real_, "tail", "too few control events to judge the shape"))
    r <- tryCatch(gate_tail_ratio(x), error = function(e) NA_real_)
    ex <- tryCatch(mean(x > gate_cut_tail(x, tail_p)), error = function(e) NA_real_)
    if (is.na(r) || is.na(ex)) return(row(r, ex, "tail", "tail ratio not computable"))
    if (r <= GATE_AUTO[["ratio"]]) return(row(r, ex, "tail", "clean negative peak"))
    if (ex > GATE_AUTO[["max_excess"]])
      return(row(r, ex, "tail", "too much control above the tail cut: probable positive population, tail cut kept"))
    row(r, ex, "quantile", "heavy upper tail")
  })
  do.call(rbind, out)
}

# Asinh cofactor from the width of the negative peak on the raw scale: k times the
# one-sided spread of the main peak. The peak is located on the lower 90% of the
# events, so a long positive tail does not coarsen the density estimate.
gate_cofactor <- function(x, k = GATE_COFACTOR_K) {
  x <- x[is.finite(x)]
  if (length(x) > 2e5) x <- x[round(seq(1, length(x), length.out = 2e5))]
  core <- x[x <= stats::quantile(x, 0.9, names = FALSE)]
  k * gate_peak(core)[["spread"]]
}

# gate_cofactor() with the marker named in any error.
gate_cofactor_marker <- function(marker, x, k = GATE_COFACTOR_K) {
  tryCatch(gate_cofactor(as.numeric(x), k),
           error = function(e) stop("cofactor for marker ", marker, ": ", conditionMessage(e), call. = FALSE))
}

# Per-event batch labels of the control events, and the sample -> batch map.
gate_batch_map <- function(samples, by) {
  if (!is.character(by) || length(by) != 1L || !by %in% names(samples))
    stop("by must be the name of one column of samples", call. = FALSE)
  b <- as.character(samples[[by]])
  stats::setNames(b, as.character(samples$sample_id))
}

# Per-batch cuts for every marker. Batches with fewer than `min_events` control
# events (or none) take the pooled cut. Returns the cut matrix (batch x marker)
# and the tail-rule cut and negative-peak spread per batch (NA for fallback batches).
gate_batch_cuts <- function(X, bev, all_batches, pooled, eff, min_events, probs, tail_p, ...) {
  mk <- colnames(X)
  dn <- list(all_batches, mk)
  cut <- tailcut <- spread <- matrix(NA_real_, length(all_batches), length(mk), dimnames = dn)
  n <- as.integer(table(factor(bev, all_batches)))
  names(n) <- all_batches
  ok <- all_batches[n >= min_events]
  for (b in all_batches) {
    use <- b %in% ok
    for (m in mk) {
      cut[b, m] <- if (use) gate_cut(eff[[m]], X[bev == b, m], m, probs, tail_p, ...) else pooled[[m]]
      if (use) {
        tailcut[b, m] <- tryCatch(gate_cut_tail(X[bev == b, m], tail_p), error = function(e) NA_real_)
        spread[b, m] <- tryCatch(gate_peak(X[bev == b, m])[["spread"]], error = function(e) NA_real_)
      }
    }
  }
  list(cut = cut, tailcut = tailcut, spread = spread, n = n, ok = ok)
}

# QC flags per marker, from control events only. `cut_ev` and `tailcut_ev` are
# the cut and the tail-rule cut that apply to each event (its batch's, if any).
#  * cut_in_broad_tail:   density at the cut is above qc[["tail_density"]] of the peak density.
#  * high_control_positive: control events above the tail-rule cut exceed qc[["max_control_pos"]].
#  * batch_cut_spread:    range of the batch cuts over the typical negative-peak spread
#                         exceeds qc[["batch_spread"]] (only with `by`).
gate_qc <- function(X, cut_ev, tailcut_ev, batch, qc, fallback = character()) {
  na_row <- function(m) data.frame(marker = m, cut_density_ratio = NA_real_, pct_control_pos = NA_real_,
    pct_control_pos_tailrule = NA_real_, batch_cut_spread = NA_real_, cut_in_broad_tail = NA,
    high_control_positive = NA, batch_cut_spread_flag = NA, batch_without_control = length(fallback) > 0L,
    batch_fallback = paste(fallback, collapse = ","), stringsAsFactors = FALSE)
  rows <- lapply(colnames(X), function(m) {
    # Too few control events (a manual gate needs none): no flags, NA numbers.
    if (nrow(X) < 2L) return(na_row(m))
    tryCatch(gate_qc_marker(X, m, cut_ev, tailcut_ev, batch, qc, fallback), error = function(e) na_row(m))
  })
  do.call(rbind, rows)
}

gate_qc_marker <- function(X, m, cut_ev, tailcut_ev, batch, qc, fallback) {
  {
    xc <- X[, m] - cut_ev[, m]
    if (length(xc) > 2e5) xc <- xc[round(seq(1, length(xc), length.out = 2e5))]
    d <- stats::density(xc)
    dens <- stats::approx(d$x, d$y, xout = 0, rule = 2)$y / max(d$y)
    if (0 < min(d$x) || 0 > max(d$x)) dens <- 0
    pos_tail <- if (anyNA(tailcut_ev[, m])) NA_real_ else mean(X[, m] > tailcut_ev[, m])
    spread <- NA_real_
    if (!is.null(batch)) {
      cuts <- batch$cut[batch$ok, m]
      sp <- stats::median(batch$spread[batch$ok, m], na.rm = TRUE)
      spread <- if (length(cuts) > 1L && is.finite(sp)) (max(cuts) - min(cuts)) / sp else NA_real_
    }
    data.frame(marker = m, cut_density_ratio = dens, pct_control_pos = 100 * mean(X[, m] > cut_ev[, m]),
               pct_control_pos_tailrule = 100 * pos_tail, batch_cut_spread = spread,
               cut_in_broad_tail = dens > qc[["tail_density"]],
               high_control_positive = isTRUE(pos_tail > qc[["max_control_pos"]]),
               batch_cut_spread_flag = isTRUE(spread > qc[["batch_spread"]]),
               batch_without_control = length(fallback) > 0L, batch_fallback = paste(fallback, collapse = ","),
               stringsAsFactors = FALSE)
  }
}

# Flag names set for one marker, as a short string.
gate_flag_text <- function(qc) {
  f <- cbind(cut_in_broad_tail = qc$cut_in_broad_tail, high_control_positive = qc$high_control_positive,
             batch_cut_spread = qc$batch_cut_spread_flag, batch_without_control = qc$batch_without_control)
  f[is.na(f)] <- FALSE
  vapply(seq_len(nrow(f)), function(i) paste(colnames(f)[f[i, ]], collapse = ", "), "")
}

#' Suggest an asinh cofactor from the width of the negative peak
#'
#' A label-free rule for the cofactor of `asinh(x / cofactor)`. It reads the
#' control events only. For each marker, the main peak of the raw control
#' density is the negative population, its one-sided spread `s` is the root mean
#' square of the events below the mode about the mode (as in the `"tail"` gate),
#' and the cofactor is `k * s`. The default `k = 3` puts the negative peak and the
#' tail cut (mode plus 3.09 spreads) in the near-linear zone of the transform, at
#' an asinh argument of about 1, while positive events spread out on the log-like
#' zone. The peak is located on the lower 90% of the events so that a long
#' positive tail does not coarsen the density estimate.
#'
#' The rule needs a negative peak near zero on the raw scale (compensated
#' linear values). It is a starting point, not a verdict: check the gate plot.
#' Pass `cofactor = "auto"` to [cj_gate()] to use it there.
#'
#' @param cells,markers,samples,control As in [cj_gate()].
#' @param k Multiplier of the negative-peak spread.
#' @return A named numeric vector, one cofactor per marker.
#' @export
#' @examples
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 400, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' cj_cofactor(sim$cells, c("IFNG", "TNF"), samples, control = "unstim")
cj_cofactor <- function(cells, markers, samples, control = "control", k = GATE_COFACTOR_K) {
  stopifnot(is.data.frame(cells), is.data.frame(samples), "sample_id" %in% names(cells),
            is.numeric(k), length(k) == 1L, k > 0)
  markers <- gate_markers(markers)
  samples <- recode_control(samples, control)
  is_ctrl <- as.character(cells$sample_id) %in% as.character(samples$sample_id[samples$condition == CONTROL_CONDITION])
  if (sum(is_ctrl) < 100L) stop("only ", sum(is_ctrl), " control events", call. = FALSE)
  vapply(names(markers), function(m) gate_cofactor_marker(m, cells[[markers[[m]]]][is_ctrl], k), numeric(1))
}

#' Simulate cells with controlled control-distribution shapes
#'
#' A simulator for testing gates against a known truth. Each marker has a
#' negative population and a positive population; `truth_<marker>` says which an
#' event came from. Controls hold `control_pos` positives and stimulated samples
#' `stim_pos`. Samples are spread over `batches`, each shifting the channel
#' offset (and optionally the width) of all its events.
#'
#' Profiles (`profile`, one per marker, recycled):
#' * `"clean"`: negative `N(0, 1)`.
#' * `"skewed"`: negative `N(0, 1)` plus an exponential upper tail of mean 0.9 on
#'   a random 45% of events (a heavy, right-skewed negative peak).
#' * `"very_skewed"`: as `"skewed"`, but 80% of events and a mean of 2.5, so more than 1%
#'   of the control lies beyond 7 spreads.
#' * `"t"`: negative scaled Student `t` with 4 degrees of freedom (heavy both sides).
#' * `"wide_background"`: clean negative, but 7% of the control events are weak
#'   positives (`N(3.5, 1.2)`), as for markers with real unstimulated expression.
#'
#' Strong positives are `N(pos_mu, 3)`. In `"wide_background"` the weak positives
#' also exist in every stimulated sample.
#'
#' @param profile Character vector of profiles, one per marker (recycled).
#' @param markers Marker names.
#' @param samples_per_batch Number of control and of stimulated samples per batch.
#' @param cells_per_sample Events per sample.
#' @param batches Number of batches.
#' @param batch_shift Standard deviation of the batch offset in negative-peak
#'   units (offsets are drawn once per batch and marker, seeded).
#' @param batch_scale Standard deviation of the log width factor per batch.
#' @param control_pos,stim_pos Fraction of positives in controls and in stimulated
#'   samples (a vector with one value per marker, recycled). `control_pos = NULL`
#'   uses 0.0005, or 0.07 for `"wide_background"`.
#' @param pos_mu Mean of the strong positives, in negative-peak units (default 14).
#' @param seed Integer seed.
#' @return A list: `cells` (data frame with `sample_id`, `donor_id`, `condition`
#'   (`"unstim"`/`"stim"`), `batch`, marker columns and `truth_<marker>`) and
#'   `samples` (the sample table).
#' @export
#' @examples
#' s <- cj_simulate_gate_cells(c("clean", "skewed"), c("A", "B"), seed = 1)
#' table(s$cells$condition)
cj_simulate_gate_cells <- function(profile = "clean", markers = c("A", "B", "C", "D"),
                                   samples_per_batch = 4L, cells_per_sample = 3000L, batches = 1L,
                                   batch_shift = 0, batch_scale = 0, control_pos = NULL,
                                   stim_pos = 0.2, pos_mu = 14, seed = 1L) {
  profile <- match.arg(profile, c("clean", "skewed", "very_skewed", "t", "wide_background"), several.ok = TRUE)
  profile <- rep_len(profile, length(markers))
  if (is.null(control_pos)) control_pos <- ifelse(profile == "wide_background", 0.07, 0.0005)
  control_pos <- rep_len(control_pos, length(markers))
  stim_pos <- rep_len(stim_pos, length(markers))
  set.seed(seed)
  bnames <- sprintf("B%d", seq_len(batches))
  off <- matrix(stats::rnorm(batches * length(markers), 0, batch_shift), batches, dimnames = list(bnames, markers))
  wid <- exp(matrix(stats::rnorm(batches * length(markers), 0, batch_scale), batches, dimnames = list(bnames, markers)))
  grid <- expand.grid(rep = seq_len(samples_per_batch), condition = c("unstim", "stim"), batch = bnames,
                      stringsAsFactors = FALSE)
  grid$sample_id <- sprintf("S%03d_%s", seq_len(nrow(grid)), grid$condition)
  grid$donor_id <- sprintf("%s_D%02d", grid$batch, grid$rep)
  one <- function(g) {
    n <- cells_per_sample
    d <- data.frame(sample_id = rep(g$sample_id, n), donor_id = rep(g$donor_id, n),
                    condition = rep(g$condition, n), batch = rep(g$batch, n))
    for (j in seq_along(markers)) {
      m <- markers[[j]]
      p <- if (g$condition == "stim") stim_pos[[j]] else control_pos[[j]]
      pos <- stats::runif(n) < p
      neg <- switch(profile[[j]],
        t = stats::rt(n, 4),
        skewed = stats::rnorm(n) + (stats::runif(n) < 0.45) * stats::rexp(n, 1 / 0.9),
        very_skewed = stats::rnorm(n) + (stats::runif(n) < 0.8) * stats::rexp(n, 1 / 2.5),
        stats::rnorm(n))
      strong <- pos
      if (profile[[j]] == "wide_background") {
        # weak positives (N(3.5, 1.2)) at the control rate in every sample; stimulation adds strong ones
        strong <- if (g$condition == "stim") stats::runif(n) < stim_pos[[j]] else rep(FALSE, n)
        weak <- stats::runif(n) < control_pos[[j]]
        pos <- strong | weak
        x <- ifelse(strong, stats::rnorm(n, pos_mu, 3), ifelse(weak, stats::rnorm(n, 3.5, 1.2), neg))
      } else x <- ifelse(pos, stats::rnorm(n, pos_mu, 3), neg)
      d[[m]] <- off[g$batch, m] + wid[g$batch, m] * x
      d[[paste0("truth_", m)]] <- pos
    }
    d
  }
  cells <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) one(grid[i, ])))
  list(cells = cells, samples = grid[c("sample_id", "donor_id", "condition", "batch")])
}
