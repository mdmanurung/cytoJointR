# Adapter to the real COMPASS package (docs/adr/0001). Native COMPASSResult
# objects are returned untouched; COMPASS FS/PFS/response probabilities are
# COMPASS quantities, not cytoJointR sign probabilities.

# COMPASS default category filter (COMPASS::COMPASS), applied to stimulated counts.
COMPASS_DEFAULT_FILTER <- function(x) colSums(x > 5) > 2

# "IFNG&!TNF&IL2": COMPASS's own pattern-name format, markers in bits order.
compass_pattern_names <- function(bits) {
  apply(bits, 1L, function(b) paste0(ifelse(b == 1L, "", "!"), colnames(bits), collapse = "&"))
}

compass_require <- function() {
  if (!requireNamespace("COMPASS", quietly = TRUE))
    stop("The COMPASS package is required; install it to use cj_compass()", call. = FALSE)
}

# A priori choices, not tuned on any data: two controls of one donor "differ" when a
# two-sample proportion test of the any-positive fraction gives p < 0.01 and the
# fractions differ by at least 25% (fold). The fold keeps a huge-N test from flagging
# a trivial difference; the test keeps a small-N fold from flagging noise.
COMPASS_SHIFT_ALPHA <- 0.01
COMPASS_SHIFT_FOLD <- 1.25

# "ok", "differs" or "untestable" for the controls of one donor (rows of counts).
compass_shift_check <- function(counts, bits) {
  if (nrow(counts) < 2L) return("untestable")
  n <- rowSums(counts)
  pos <- n - counts[, which(rowSums(bits) == 0)]
  differs <- vapply(utils::combn(nrow(counts), 2L, simplify = FALSE), function(i) {
    fold <- max(pos[i] / n[i]) / min(pos[i] / n[i])
    p <- suppressWarnings(stats::prop.test(pos[i], n[i])$p.value)
    isTRUE(p < COMPASS_SHIFT_ALPHA && fold >= COMPASS_SHIFT_FOLD)
  }, NA)
  if (any(differs)) "differs" else "ok"
}

# Decide the control sample of every stimulated sample. Returns the contrasts to run
# (with control_source, control_fallback and a control_sample_id that is NA for
# "population") and the contrasts left out with the reason.
# `force_missing` lists contrast_ids whose control is treated as absent (pseudo-missing fits
# used by the calibration).
compass_resolve_controls <- function(x, missing_control, force_missing = character()) {
  ct <- x$contrasts
  ct$control_sample_id[ct$contrast_id %in% force_missing] <- NA
  none <- is.na(ct$control_sample_id)
  ct$control_source <- ifelse(none, NA_character_, "paired")
  ct$control_fallback <- FALSE
  ct$control_shift_check <- NA_character_
  if (missing_control != "drop") {
    ctrl <- x$samples[x$samples$condition == CONTROL_CONDITION, ]
    lev <- levels(x$samples$visit)
    for (i in which(none)) {
      cand <- ctrl[ctrl$donor_id == ct$donor_id[i] & ctrl$parent_population == ct$parent_population[i] &
                    ctrl$visit != ct$visit[i], ]
      # A control with no cells cannot stand in for one; leave it out.
      cand <- cand[rowSums(x$counts[cand$sample_id, , drop = FALSE]) > 0, ]
      if (missing_control == "donor" && nrow(cand)) {
        pos <- match(as.character(cand$visit), lev)
        # Nearest visit by level order; a tie goes to the earlier visit.
        pick <- cand[order(abs(pos - match(as.character(ct$visit[i]), lev)), pos)[1], ]
        ct$control_sample_id[i] <- pick$sample_id
        ct$control_source[i] <- paste0("donor_visit:", pick$visit)
        ct$control_shift_check[i] <- compass_shift_check(x$counts[cand$sample_id, , drop = FALSE], x$bits)
      } else {
        ct$control_source[i] <- "population"
        ct$control_fallback[i] <- missing_control == "donor"
      }
    }
  }
  drop <- is.na(ct$control_source)
  # A run in which every donor lacks a control has no population control to borrow.
  gid <- paste(ct$stimulation, ct$visit, ct$parent_population, sep = " | ")
  nocontrol <- !drop & unname(tapply(is.na(ct$control_sample_id), gid, all)[gid])
  keep_cols <- c("contrast_id", "donor_id", "visit", "parent_population", "stimulation")
  excluded <- ct[drop | nocontrol, keep_cols]
  if (nrow(excluded))
    excluded$reason <- ifelse(nocontrol[drop | nocontrol], "no donor in this run has a control sample",
                              "no control sample")
  list(contrasts = ct[!drop & !nocontrol, ], excluded = excluded)
}

# Counts of one COMPASS run (one stimulation, visit and parent population) after COMPASS's
# category filter. A sample without a control gets an all-zero row: COMPASS then uses
# the population control.
compass_prepare_group <- function(cg, x, ord, nm, null_col, category_filter) {
  n_s <- x$counts[cg$stim_sample_id, ord, drop = FALSE]
  n_u <- x$counts[match(cg$control_sample_id, rownames(x$counts)), ord, drop = FALSE]
  n_u[is.na(cg$control_sample_id), ] <- 0L
  n_u_raw <- n_u
  dimnames(n_s) <- dimnames(n_u) <- list(cg$donor_id, nm)
  # Same filter, same null-category bookkeeping as COMPASS::COMPASS().
  del <- integer()
  if (!is.null(category_filter)) {
    del <- setdiff(which(!match.fun(category_filter)(n_s)), null_col)
    if (length(del)) {
      n_s[, null_col] <- n_s[, null_col] + rowSums(n_s[, del, drop = FALSE])
      n_u[, null_col] <- n_u[, null_col] + rowSums(n_u[, del, drop = FALSE])
      n_s <- n_s[, -del, drop = FALSE]; n_u <- n_u[, -del, drop = FALSE]
    }
  }
  list(n_s = n_s, n_u = n_u, n_u_raw = n_u_raw, del = del)
}

# SimpleCOMPASS() on a prepared group; a COMPASSResult or the error message.
compass_run_group <- function(prep, cg, verbose, ...) {
  meta <- data.frame(donor_id = cg$donor_id, stimulation = cg$stimulation,
                     visit = cg$visit, parent_population = cg$parent_population,
                     control_source = cg$control_source)
  # SimpleCOMPASS() skips the posterior that PosteriorDiff() and friends read;
  # COMPASS() computes it with the same unexported helper.
  run <- function() {
    f <- COMPASS::SimpleCOMPASS(prep$n_s, prep$n_u, meta, "donor_id", ...)
    f$fit$posterior <- utils::getFromNamespace("compute_posterior", "COMPASS")(f)
    f
  }
  tryCatch({
    res <- NULL
    if (verbose) res <- run() else utils::capture.output(res <- suppressMessages(run()))
    res
  }, error = function(e) conditionMessage(e))
}

#' Fit COMPASS to a cytoJointR data set
#'
#' Runs one `COMPASS::SimpleCOMPASS()` per stimulation, visit and parent
#' population, using each donor's stimulated and control sample as the
#' `n_s` and `n_u` rows. The native `COMPASSResult` objects are returned in
#' `fits`, so every COMPASS function (`COMPASS::scores()`,
#' `COMPASS::PosteriorDiff()`, `plot()`, ...) works on them.
#'
#' By default (`missing_control = "drop"`) contrasts without a control sample
#' cannot enter COMPASS; they are listed in `excluded`. The other two settings
#' keep them (see below). `SimpleCOMPASS()` has no category filter, so the adapter applies
#' COMPASS's own default from `COMPASS::COMPASS()`,
#' `function(x) colSums(x > 5) > 2` on the stimulated counts, and moves removed
#' patterns' cells to the all-negative pattern as COMPASS does. What was removed
#' is reported in `filters`. Pass `category_filter = NULL` to keep every pattern.
#'
#' COMPASS FS, PFS and response probabilities have COMPASS's definitions and
#' range. `compass_response_degree*` is `COMPASS::Response()`: the mean, over
#' the cell subsets with at least that many markers, of the posterior
#' probability that the subset responds. They are not cytoJointR sign probabilities and should not be
#' compared with `p_positive` as if they were.
#'
#' @section Missing control samples:
#' COMPASS integrates out each donor's control proportions, so an all-zero
#' `n_u` row adds nothing to the likelihood and that donor's stimulated counts
#' are compared with the population control distribution that the other donors
#' define. `missing_control` chooses what a stimulated sample without a control
#' of its own visit gets:
#' * `"drop"`: nothing; the contrast is excluded (default, unchanged).
#' * `"donor"`: the same donor's control from the nearest other visit with the
#'   same parent population. Nearest is by order of the visit levels; a tie goes
#'   to the earlier visit. Controls with no cells are skipped. A donor with no
#'   usable control at another visit of the same parent population falls back
#'   to `"population"`, flagged by `control_fallback = TRUE`.
#' * `"population"`: an all-zero control row.
#'
#' Every `scores` row has `control_source`: `"paired"`, `"donor_visit:<visit>"`
#' or `"population"`. A `"population"` row is scored against the other donors'
#' controls, so it is only as good as the assumption that this donor's
#' unstimulated pattern distribution resembles theirs (`docs/adr/0003`). A
#' borrowed control carries any real visit effect on the background. A run in
#' which no donor has a control is skipped and listed in `excluded`.
#' `COMPASS::PosteriorDiff()` on a `"population"` row is relative to that
#' population control. The category filter reads only stimulated counts, so
#' zero rows do not change which patterns it removes.
#'
#' @section Calibration (experimental):
#' `"donor"` and `"population"` print a message saying so. In simulation (40
#' donors, 4 markers, 50% responders, 20 seeds, controls missing at random,
#' `runs/2026-11-07_compass-missing-control/`) the donor-level call "max mean
#' gamma > 0.5" kept its ranking but not its error rate for non-paired rows.
#' Among non-responders:
#'
#' | Option | Setting | FPR, non-paired | FPR, paired |
#' |---|---|---|---|
#' | `"population"` | moderate effect | 0.19 to 0.39 | 0.05 to 0.10 |
#' | `"population"` | rare, weak effect | 0.01 to 0.03, TPR 0.20 to 0.70 (paired 0.64) | 0.01 |
#' | `"donor"` | stable background | 0.06 to 0.17 (moderate) | 0.05 to 0.08 |
#' | `"donor"` | background x1.5 at visit 2 | 0.4 to 0.55 pooled; 0.81 to 0.94 for the later visit | 0.06 to 0.19 |
#'
#' Use them to rank donors, not to call responders. For `"donor"`, the donor's
#' available controls at other visits are compared first (any-positive fraction,
#' two-sample proportion test and a fold difference, both fixed a priori and not
#' tuned on benchmark data). `control_shift_check` is `"ok"`, `"differs"` or
#' `"untestable"` (only one control to borrow), and `"differs"` or `"untestable"`
#' gives a warning. See `docs/adr/0003-compass-missing-control.md`.
#'
#' @section Calibrated calls (experimental):
#' `missing_control_calibrate = TRUE` maps the max mean gamma of non-paired rows onto the
#' paired scale without responder labels. In each COMPASS run, paired donors are
#' pseudo-removed in folds of 5 (at most `max_pseudo` per partition, chosen at random), the
#' control is resolved by the same `missing_control` rule, and the run is refitted. A pseudo
#' donor then has two gammas: one with its control missing and one with its own control.
#' `pseudo_reps` random partitions are drawn and their points pooled. An isotonic regression
#' of the paired gamma on the pseudo-missing gamma gives the map. It is applied to the run's
#' non-paired rows of the same kind (`"donor_visit"` or `"population"`), and a calibrated
#' gamma never exceeds the raw gamma. Use `compass_call` (calibrated gamma above 0.5), not
#' raw max gamma above 0.5. `scores` gets `compass_max_gamma`, `compass_gamma_cal`,
#' `compass_call`, `calibration` (`"paired"`, `"isotonic"`, `"quantile"`, `"fallback"` or
#' `"switched"`) and `compass_max_gamma_switch`; `$calibration$table` has one line per run and
#' kind (with `pseudo_auc`, `excess`, `switched`), `$calibration$pseudo` has the pseudo pairs,
#' `$calibration$n_failed` counts failed refits and `$calibration$switched` lists switched runs.
#'
#' A map needs at least 15 distinct pseudo donors, so at least 15 paired donors in the run
#' (with a usable borrowed control for `"donor"`), and at least 3 pseudo donors with a paired
#' gamma above 0.5 and 3 not above it. Otherwise the run uses a fallback cut: raw gamma 0.85 is
#' sent to 0.5. That cut comes from simulation only (the cut at which non-paired FPR met paired
#' FPR on the dev seeds of `runs/2026-11-07_compass-missing-control/`); it depends on that
#' design, and it costs power at a rare, weak effect. The assumption is that the pseudo-missing
#' donors resemble the truly missing ones (controls missing at random) and that the pseudo
#' gamma ranks donors like the paired gamma. The second is checked: when the AUC of
#' pseudo-missing gamma for the paired call is below 0.8, a warning says the calibrated calls
#' there are unreliable.
#'
#' **Visit effect with `"donor"`.** If the background differs between visits, a control
#' borrowed from another visit is biased and no map repairs it. `control_shift_check` cannot
#' tell: with two visits it is always `"untestable"`, for every borrowed row. The detector is
#' the pseudo AUC. With `donor_switch = TRUE` (default), a `"donor"` run whose pseudo AUC is
#' below 0.7, with at least 5 pseudo donors with a paired call, is refitted with the
#' same-visit `"population"` control. That control is an all-zero row, which COMPASS compares
#' with the other donors' controls of the same stimulation, visit and parent population, so it
#' shares the visit effect. Those rows are calibrated with population-rule pseudo fits and
#' marked `calibration = "switched"`; `compass_max_gamma_switch` holds their gamma. In
#' simulation (clean test, seeds 31 to 40) this took the later-visit FPR from 0.78 to 0.97 down to
#' 0.05 to 0.07 at 10 and 25% missing and 0.21 at 50% (moderate effect), and to 0.04 to 0.47 at the rare, weak
#' effect, where power falls to the population level (0.45 to 0.8). A run is switched without a shift in
#' 0 to 6% of the borrowed rows. Switching is partial: a shift whose pseudo AUC stays above 0.7, or with
#' fewer than 5 pseudo donors with a paired call, is not repaired (the pseudo-AUC warning still applies). Cost: about `pseudo_reps * max_pseudo / 5 + 1` COMPASS fits per run with
#' non-paired rows, 12 to 14 times the plain call with the defaults, plus a second set for switched
#' runs. Numbers: `docs/adr/0003-compass-missing-control.md`.
#'
#' @param x A `cj_data` (or a `cj_fit`, whose `data` is used).
#' @param category_filter COMPASS category filter; the default is COMPASS's own.
#' @param verbose Show COMPASS's console output.
#' @param missing_control What to do with a stimulated sample that has no control
#'   sample at its visit: `"drop"`, `"donor"` or `"population"`. See *Missing
#'   control samples*.
#' @param missing_control_calibrate `FALSE` (default), `TRUE` or `"isotonic"`, or
#'   `"quantile"`. With `missing_control = "donor"` or `"population"`, map the
#'   max mean gamma of the non-paired rows onto the paired scale with a label-free
#'   empirical null from the same run. See *Calibrated calls*.
#' @param max_pseudo Most paired donors per COMPASS run that are pseudo-removed for the
#'   calibration in one partition (default 20). Each 5 donors cost one extra COMPASS fit per run, so the
#'   calibrated call costs about `1 + max_pseudo / 5` fits per run with non-paired rows. A
#'   map needs at least 15 pseudo donors, so a cap below 15 always uses the fallback cut.
#' @param donor_switch With `"donor"` and calibration, switch a run whose borrowed controls cannot rank
#'   donors (pseudo AUC below 0.7, with at least 5 pseudo donors with a paired call) to the same-visit
#'   `"population"` control. See *Calibrated calls*.
#' @param pseudo_reps Number of random partitions of the paired donors into pseudo-removal folds
#'   (default 1). The pseudo points of all repeats are pooled for the map.
#' @param calibrate_seed Seed for choosing the pseudo-removed donors. The caller's RNG
#'   state is left as it was.
#' @param ... Passed to `COMPASS::SimpleCOMPASS()`, e.g. `iterations`,
#'   `replications`, `seed`. COMPASS's defaults apply when omitted.
#' @return A `cj_compass` list with `fits` (named list of `COMPASSResult`),
#'   `scores` (tidy table with `compass_FS`, `compass_PFS`,
#'   `compass_response_degree1`, `compass_response_degree2`, `control_source`,
#'   `control_fallback`, `control_shift_check`), `groups` (one
#'   row per COMPASS run with its status), `filters` (what the category filter
#'   removed per run) and `excluded` (contrasts left out, with the reason).
#' @export
#' @examplesIf requireNamespace("COMPASS", quietly = TRUE)
#' \donttest{
#' d <- cj_example_draws(n_donor = 4, markers = c("IFNG", "TNF"), visits = "v1",
#'                       stimulations = "CMV", n_draw = 5)
#' counts <- t(vapply(seq_len(nrow(d$samples)), function(i)
#'   rmultinom(1, d$samples$N[i], colMeans(d$p[, i, ]))[, 1], numeric(4)))
#' storage.mode(counts) <- "integer"
#' dimnames(counts) <- list(d$samples$sample_id, dimnames(d$p)[[3]])
#' x <- cytoJointR:::new_cj_data(counts, d$samples, d$contrasts, d$bits)
#' cc <- cj_compass(x, iterations = 200, replications = 2)
#' cc$scores
#' }
cj_compass <- function(x, category_filter = COMPASS_DEFAULT_FILTER, verbose = FALSE,
                       missing_control = c("drop", "donor", "population"),
                       missing_control_calibrate = FALSE, max_pseudo = 20L, calibrate_seed = 1L, pseudo_reps = 3L, donor_switch = TRUE, ...) {
  compass_require()
  # SimpleCOMPASS() calls set.seed(): give the caller's RNG state back on exit.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else
    suppressWarnings(rm(".Random.seed", envir = globalenv())), add = TRUE)
  missing_control <- match.arg(missing_control)
  if (missing_control == "drop" && !isFALSE(missing_control_calibrate))
    warning("missing_control_calibrate is ignored with missing_control = \"drop\": no row runs ",
            "without a paired control", call. = FALSE)
  cal_method <- if (isFALSE(missing_control_calibrate)) NA_character_ else
    if (isTRUE(missing_control_calibrate)) "isotonic" else
      match.arg(missing_control_calibrate, c("isotonic", "quantile"))
  if (missing_control == "drop") cal_method <- NA_character_
  if (inherits(x, "cj_fit")) x <- x$data
  if (!inherits(x, "cj_data")) stop("x must be a cj_data or cj_fit", call. = FALSE)
  if (grepl("[&!+*?()\\[\\]{}|^$\\\\]", paste(x$markers, collapse = ""), perl = TRUE))
    stop("COMPASS parses marker names with '&', '!' and a regular expression; ",
         "rename markers to use only letters, digits, '_', '.' or '-'", call. = FALSE)
  ct <- compass_resolve_controls(x, missing_control)
  excluded <- ct$excluded
  ct <- ct$contrasts
  n_alt <- sum(ct$control_source != "paired")
  if (missing_control != "drop" && is.na(cal_method))
    message("missing_control = \"", missing_control, "\" is experimental: use it to rank donors. ",
            "In simulation (40 donors, 4 markers) the gamma > 0.5 call is not calibrated for non-paired ",
            "rows. False-positive rate among non-responders: ",
            if (missing_control == "population")
              "0.19 to 0.39 at a moderate effect (paired 0.05 to 0.10); power falls at a rare, weak effect."
            else
              "0.06 to 0.17 when the background is stable across visits, 0.4 to 0.55 pooled and 0.81 to 0.94 for the later visit when it rises by 1.5-fold (paired 0.05 to 0.10).",
            " Pass missing_control_calibrate = TRUE (needs 15 or more paired donors per run) for a calibrated call. ",
            "See docs/adr/0003-compass-missing-control.md.")
  bad <- ct$control_shift_check %in% c("differs", "untestable")
  if (any(bad))
    warning(sum(ct$control_shift_check %in% "differs"), " borrowed control(s) come from a visit whose background ",
            "differs from the donor's other controls and ", sum(ct$control_shift_check %in% "untestable"),
            " could not be checked (the donor has one control); see control_shift_check", call. = FALSE)
  if (n_alt)
    message(n_alt, " of ", nrow(ct), " contrast(s) ran without a paired control: ",
            sum(startsWith(ct$control_source, "donor_visit")), " with the donor's control from another visit, ",
            sum(ct$control_source == "population"), " against the population control (",
            sum(ct$control_fallback), " of them a fallback from \"donor\")")

  # COMPASS wants the all-negative pattern last; the rest by degree.
  degree <- rowSums(x$bits)
  ord <- c(which(degree > 0)[order(degree[degree > 0])], which(degree == 0))
  nm <- compass_pattern_names(x$bits)[ord]
  null_col <- length(nm)
  gid <- paste(ct$stimulation, ct$visit, ct$parent_population, sep = " | ")
  fits <- list(); groups <- list(); filters <- list(); tidy <- list()
  for (g in unique(gid)) {
    cg <- ct[gid == g, ]
    cg <- cg[order(cg$donor_id), ]
    prep <- compass_prepare_group(cg, x, ord, nm, null_col, category_filter)
    del <- prep$del
    n_u_raw <- prep$n_u_raw
    key <- cg[1, c("stimulation", "visit", "parent_population")]
    filters[[g]] <- cbind(group = g, key, n_patterns = length(nm), n_removed = length(del),
                          removed = paste(nm[del], collapse = "; "),
                          cells_moved_stimulated = sum(x$counts[cg$stim_sample_id, ord[del]]),
                          cells_moved_control = sum(n_u_raw[, del]),
                          n_zero_control = sum(is.na(cg$control_sample_id)),
                          row.names = NULL)
    fit <- compass_run_group(prep, cg, verbose, ...)
    ok <- !is.character(fit)
    groups[[g]] <- cbind(group = g, key, n_donors = nrow(cg),
                         status = if (ok) "fitted" else paste("failed:", fit), row.names = NULL)
    if (!ok) next
    fits[[g]] <- fit
    # Response() returns a one-column matrix (rows = donors), so index its column
    # by donor; indexing the matrix itself by name gives NA.
    resp <- function(deg) if (length(x$markers) >= deg)
      suppressMessages(COMPASS::Response(fit, degree = deg))[cg$donor_id, 1] else
        rep(NA_real_, nrow(cg))
    tidy[[g]] <- data.frame(
      contrast_id = cg$contrast_id, donor_id = cg$donor_id, stimulation = cg$stimulation,
      visit = cg$visit, parent_population = cg$parent_population,
      compass_FS = unname(COMPASS::FunctionalityScore(fit)[cg$donor_id]),
      compass_PFS = unname(COMPASS::PolyfunctionalityScore(fit)[cg$donor_id]),
      compass_response_degree1 = unname(resp(1)),
      compass_response_degree2 = unname(resp(2)),
      control_source = cg$control_source, control_fallback = cg$control_fallback,
      control_shift_check = cg$control_shift_check,
      compass_max_gamma = unname(compass_max_gamma(fit)[cg$donor_id]))
  }
  bind <- function(l) if (length(l)) do.call(rbind, unname(l))
  scores <- bind(tidy)
  calibration <- NULL
  if (!is.na(cal_method) && !is.null(scores)) {
    pseudo <- if (n_alt) compass_pseudo_fits(x, ct, fits, gid, ord, nm, null_col, category_filter, verbose,
                                  missing_control = missing_control, max_pseudo = max_pseudo,
                                  calibrate_seed = calibrate_seed, reps = pseudo_reps, ...)
    gs <- paste(scores$stimulation, scores$visit, scores$parent_population, sep = " | ")
    if (is.null(pseudo)) pseudo <- data.frame(group = character(), donor_id = character(), cls = character(),
                                              fold = integer(), rep = integer(), gamma_pseudo = numeric(), gamma_paired = numeric())
    n_failed <- if (is.null(attr(pseudo, "n_failed"))) 0L else attr(pseudo, "n_failed")
    if (n_failed) message(n_failed, " pseudo-missing refit(s) failed or could not run; the map uses the folds ",
                          "that did (see $calibration$n_failed)")
    cal <- compass_calibrate_scores(scores, gs, pseudo, cal_method)
    scores <- cal$scores
    if (!is.null(cal$info)) cal$info$switched <- FALSE
    calibration <- list(table = cal$info, pseudo = pseudo, n_failed = n_failed, switched = NULL)
    # "donor" groups whose borrowed control cannot rank donors take the same-visit population control.
    if (donor_switch && missing_control == "donor" && !is.null(cal$info)) {
      bad <- cal$info$group[cal$info$control_class == "donor_visit" & !is.na(cal$info$pseudo_auc) &
                              cal$info$pseudo_auc < COMPASS_CAL_SWITCH_AUC &
                              cal$info$n_pseudo_above >= COMPASS_CAL_SWITCH_POS]
      if (length(bad)) {
        sw <- compass_switch_groups(bad, scores, gs, x, ord, nm, null_col, category_filter, verbose,
                                    cal_method, max_pseudo, calibrate_seed, pseudo_reps, ...)
        scores <- sw$scores
        calibration$switched <- sw$table
        calibration$n_failed <- calibration$n_failed + sw$n_failed
        cal$info$switched <- cal$info$group %in% sw$table$group & cal$info$control_class == "donor_visit"
        calibration$table <- cal$info
        message(sum(scores$calibration %in% "switched"), " borrowed-control row(s) in ", nrow(sw$table),
                " group(s) switched to the same-visit population control: the borrowed control did not rank ",
                "donors like the paired control (pseudo AUC < ", COMPASS_CAL_SWITCH_AUC, ")")
      }
    }
    n_fb <- sum(scores$calibration %in% "fallback")
    low <- if (is.null(cal$info)) FALSE else
      cal$info$used_map & !cal$info$switched & !is.na(cal$info$pseudo_auc) & cal$info$pseudo_auc < COMPASS_CAL_MIN_AUC
    if (n_alt) message("missing_control_calibrate: ", sum(scores$calibration %in% cal_method), " row(s) mapped onto the ",
            "paired scale (", cal_method, "), ", n_fb, " on the fallback cut (too few paired donors); ",
            "use compass_call, not max gamma > 0.5. Calibrated in simulation only (40 donors, 4 markers): ",
            "non-paired FPR was within 0.05 of the paired FPR in 12 of 12 cells on a clean held-out test (largest gap 0.05), but \"donor\" at a visit with a ",
            "shifted background is only partly repaired (switch); see docs/adr/0003-compass-missing-control.md.")
    if (any(low))
      warning("In ", sum(low), " group(s) a pseudo-missing control does not rank donors like the paired ",
              "control (AUC < ", COMPASS_CAL_MIN_AUC, "); calibrated calls there are unreliable, see ",
              "$calibration$table", call. = FALSE)
  }
  structure(list(fits = fits, scores = scores, groups = bind(groups), calibration = calibration,
                 filters = bind(filters), excluded = excluded,
                 note = paste("COMPASS FS, PFS and response probabilities are COMPASS",
                              "quantities, not cytoJointR sign probabilities.")),
            class = "cj_compass")
}

#' Draw COMPASS's own heatmap for one run
#'
#' Calls `plot()` on a native `COMPASSResult`, which draws COMPASS's mean-gamma
#' heatmap.
#'
#' @param x Output of [cj_compass()].
#' @param group Name of a run in `names(x$fits)`; default the first.
#' @param ... Passed to COMPASS's `plot()` method.
#' @return What COMPASS's `plot()` returns, invisibly.
#' @export
cj_plot_compass <- function(x, group = names(x$fits)[1], ...) {
  compass_require()
  if (!group %in% names(x$fits))
    stop("group must be one of: ", paste(names(x$fits), collapse = ", "), call. = FALSE)
  invisible(plot(x$fits[[group]], ...))
}
