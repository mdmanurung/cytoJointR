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
compass_resolve_controls <- function(x, missing_control) {
  ct <- x$contrasts
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
#' @section Convergence diagnostics:
#' COMPASS's `replications` are not independent chains. They are extra passes of
#' one Gibbs chain, thinned by `replications`, and the model runs one chain per
#' call. `n_runs > 1` therefore refits with different seeds and treats the fits
#' as chains. For each donor, from the per-iteration subset indicators `gamma`:
#' * `compass_gamma_rhat`, `compass_gamma_ess`: rank-normalised R-hat and bulk
#'   ESS ([posterior::rhat()], [posterior::ess_bulk()]) of the mean of `gamma`
#'   over the non-null subsets;
#' * `compass_FS_rhat`, `compass_FS_ess`, `compass_PFS_rhat`, `compass_PFS_ess`:
#'   the same for COMPASS's FS and PFS formulas applied to each iteration;
#' * `compass_runs_called`: runs in which the donor's call (max over subsets of
#'   mean `gamma` > 0.5) is positive, and `compass_max_gamma`, the mean over runs
#'   of that maximum;
#' * `compass_call` (the reported call, that of the first run) and
#'   `compass_call_unsupported`: `TRUE` when at most half of the runs make the
#'   reported call, so the reported call is probably the odd one out;
#' * `unstable_call`: `TRUE` when the runs disagree on the call (some, not all,
#'   are positive; `compass_call_disagree`), or when the largest of the three
#'   R-hats exceeds `rhat_max` (`compass_rhat_high`).
#'
#' Which flag to read (simulation: 40 donors, 4 markers, 12 data sets with weak
#' effects, 10x-iterations reference from 4 seeds; 18 of 480 donors had a
#' reported call that differed from the reference at COMPASS's default
#' iterations):
#' * `unstable_call` finds the donors whose call is not reproducible. It caught
#'   18 of 18 wrong reported calls, but 52 donors were flagged: in most of them
#'   the reported call happened to be the right one (precision 0.35).
#' * `compass_call_unsupported` with 5 runs flagged 19 donors, 18 of them wrong
#'   (recall 1.0, precision 0.95). Use it to decide which calls to recheck. With
#'   4 runs a 2-2 tie counts as unsupported and the precision falls to 0.50, so
#'   `diagnostics = TRUE` uses 5.
#' * With only 1000 iterations (0.1x the default) `compass_call_unsupported`
#'   found 62% of the wrong calls (precision 0.68) and `unstable_call` 95%
#'   (precision 0.47): run more iterations first.
#'
#' `rhat_max` is `Inf` by default, so R-hat is reported but does not flag. At the
#' default iterations an R-hat above 1.05 flagged 27% of donors (44% of the
#' RV144 dev donors) and found no wrong call that run disagreement had missed.
#' R-hat of 0/1 draws stays high for donors whose call is stable. Use
#' `rhat_max = 1.05` for the stricter check.
#'
#' `diagnostics$summary` counts the flagged donors per run group and gives the
#' worst R-hat, the lowest ESS and the seconds the group took. The cost is about
#' `n_runs` times that of one fit. An R-hat of 1 with an ESS of `NA` means every
#' run sat on one value. The thresholds are fixed a priori (1.05, 0.5), not tuned.
#' A flagged donor needs more `iterations` (or a closer look); an unflagged donor
#' only shows that these runs agree.
#'
#' @param x A `cj_data` (or a `cj_fit`, whose `data` is used).
#' @param category_filter COMPASS category filter; the default is COMPASS's own.
#' @param verbose Show COMPASS's console output.
#' @param missing_control What to do with a stimulated sample that has no control
#'   sample at its visit: `"drop"`, `"donor"` or `"population"`. See *Missing
#'   control samples*.
#' @param n_runs Number of COMPASS fits per run group, with seeds `seed`,
#'   `seed + 1`, ... (`seed` defaults to COMPASS's 100). The first fit is the
#'   one reported in `fits` and in the score columns, so `n_runs = 1` is the old
#'   behaviour. With `n_runs > 1` the extra fits feed the convergence columns
#'   described in *Convergence diagnostics* and are then dropped.
#' @param diagnostics `TRUE` is `n_runs = max(n_runs, 5)`. Five, because an odd number
#'   of runs never ties on the call.
#' @param rhat_max A donor is also `unstable_call` when its largest R-hat exceeds
#'   this. `Inf` (default) never flags on R-hat; the usual limit is 1.05.
#' @param ... Passed to `COMPASS::SimpleCOMPASS()`, e.g. `iterations`,
#'   `replications`, `seed`. COMPASS's defaults apply when omitted.
#' @return A `cj_compass` list with `fits` (named list of `COMPASSResult`),
#'   `scores` (tidy table with `compass_FS`, `compass_PFS`,
#'   `compass_response_degree1`, `compass_response_degree2`, `control_source`,
#'   `control_fallback`, `control_shift_check`, and with `n_runs > 1` the
#'   convergence columns), `groups` (one
#'   row per COMPASS run with its status), `filters` (what the category filter
#'   removed per run), `excluded` (contrasts left out, with the reason) and
#'   `diagnostics` (`NULL`, or with `n_runs > 1` a list: `n_runs`, `seeds`,
#'   `summary`). A run group whose extra runs fail keeps its fit and scores, has
#'   `NA` convergence columns, and says so in `groups$status`.
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
                       n_runs = 1L, diagnostics = FALSE, rhat_max = Inf, ...) {
  compass_require()
  missing_control <- match.arg(missing_control)
  n_runs <- compass_n_runs(n_runs, diagnostics)
  dots <- list(...)
  seed0 <- if (is.null(dots[["seed"]])) 100L else dots[["seed"]] # COMPASS's default seed
  dots[["seed"]] <- NULL
  if (inherits(x, "cj_fit")) x <- x$data
  if (!inherits(x, "cj_data")) stop("x must be a cj_data or cj_fit", call. = FALSE)
  if (grepl("[&!+*?()\\[\\]{}|^$\\\\]", paste(x$markers, collapse = ""), perl = TRUE))
    stop("COMPASS parses marker names with '&', '!' and a regular expression; ",
         "rename markers to use only letters, digits, '_', '.' or '-'", call. = FALSE)
  ct <- compass_resolve_controls(x, missing_control)
  excluded <- ct$excluded
  ct <- ct$contrasts
  n_alt <- sum(ct$control_source != "paired")
  if (missing_control != "drop")
    message("missing_control = \"", missing_control, "\" is experimental: use it to rank donors. ",
            "In simulation (40 donors, 4 markers) the gamma > 0.5 call is not calibrated for non-paired ",
            "rows. False-positive rate among non-responders: ",
            if (missing_control == "population")
              "0.19 to 0.39 at a moderate effect (paired 0.05 to 0.10); power falls at a rare, weak effect."
            else
              "0.06 to 0.17 when the background is stable across visits, 0.4 to 0.55 pooled and 0.81 to 0.94 for the later visit when it rises by 1.5-fold (paired 0.05 to 0.10).",
            " See docs/adr/0003-compass-missing-control.md.")
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
  fits <- list(); groups <- list(); filters <- list(); tidy <- list(); diag <- list()
  seconds <- c()
  for (g in unique(gid)) {
    cg <- ct[gid == g, ]
    cg <- cg[order(cg$donor_id), ]
    n_s <- x$counts[cg$stim_sample_id, ord, drop = FALSE]
    # A sample without a control gets an all-zero row: COMPASS then uses the population control.
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
    key <- cg[1, c("stimulation", "visit", "parent_population")]
    filters[[g]] <- cbind(group = g, key, n_patterns = length(nm), n_removed = length(del),
                          removed = paste(nm[del], collapse = "; "),
                          cells_moved_stimulated = sum(x$counts[cg$stim_sample_id, ord[del]]),
                          cells_moved_control = sum(n_u_raw[, del]),
                          n_zero_control = sum(is.na(cg$control_sample_id)),
                          row.names = NULL)
    meta <- data.frame(donor_id = cg$donor_id, stimulation = cg$stimulation,
                       visit = cg$visit, parent_population = cg$parent_population,
                       control_source = cg$control_source)
    # SimpleCOMPASS() skips the posterior that PosteriorDiff() and friends read;
    # COMPASS() computes it with the same unexported helper.
    run <- function(seed = seed0, posterior = TRUE) {
      f <- do.call(COMPASS::SimpleCOMPASS, c(list(n_s, n_u, meta, "donor_id"), dots, list(seed = seed)))
      if (posterior)
        f$fit$posterior <- utils::getFromNamespace("compute_posterior", "COMPASS")(f)
      f
    }
    quiet <- function(expr) if (verbose) expr else {
      utils::capture.output(res <- suppressMessages(expr)); res
    }
    t0 <- Sys.time()
    fit <- tryCatch(quiet(run()), error = function(e) conditionMessage(e))
    ok <- !is.character(fit)
    groups[[g]] <- cbind(group = g, key, n_donors = nrow(cg),
                         status = if (ok) "fitted" else paste("failed:", fit), row.names = NULL)
    if (!ok) next
    fits[[g]] <- fit
    if (n_runs > 1L) {
      # Same data and settings, different seeds; the first run is the one reported above.
      more <- tryCatch(lapply(seed0 + seq_len(n_runs - 1L), function(sd) quiet(run(sd, posterior = FALSE))),
                       error = function(e) conditionMessage(e))
      if (is.character(more)) {
        # The fit and its scores stand; only the convergence columns are missing for this group.
        groups[[g]]$status <- paste0("fitted (extra runs failed: ", more, ")")
      } else {
        diag[[g]] <- cbind(group = g, compass_run_diagnostics(c(list(fit), more), rhat_max))
        seconds[g] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      }
    }
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
      control_shift_check = cg$control_shift_check)
    if (!is.null(diag[[g]]))
      tidy[[g]] <- cbind(tidy[[g]], diag[[g]][match(cg$donor_id, diag[[g]]$donor_id),
                                              setdiff(names(diag[[g]]), c("group", "donor_id"))],
                         row.names = NULL)
  }
  bind <- function(l) if (length(l)) do.call(rbind, unname(l))
  # A group whose extra runs failed has no convergence columns: pad them with NA.
  cols <- unique(unlist(lapply(tidy, names)))
  tidy <- lapply(tidy, function(t) {
    t[setdiff(cols, names(t))] <- NA
    t[cols]
  })
  scores <- bind(tidy)
  diagnostics_out <- NULL
  if (length(diag))
    diagnostics_out <- list(n_runs = n_runs, seeds = seed0 + seq_len(n_runs) - 1L,
                            summary = compass_diagnostics_summary(bind(diag), seconds, n_runs))
  structure(list(fits = fits, scores = scores, groups = bind(groups),
                 filters = bind(filters), excluded = excluded, diagnostics = diagnostics_out,
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
