# Convergence diagnostics for cj_compass(): several COMPASS fits with different seeds,
# compared per donor. COMPASS's `replications` are not independent chains. They are
# `replications` extra passes of one Gibbs chain with one `thin_by = replications` thinning
# (COMPASS-discrete.R): each draw continues the previous one. Independent runs need
# independent seeds, so the runs here differ in `seed`.

# A priori choices, not tuned on any data. 1.05 is the usual loose limit for the
# rank-normalised R-hat of Vehtari et al. (2021); 0.5 is the donor-level call of the
# cytoJointR benchmarks ("max over subsets of mean gamma > 0.5").
COMPASS_GAMMA_CALL <- 0.5
# 5, not 4: with an odd number of runs a call never ties, so "the reported call is in the minority" is defined
# for every donor (4 runs: precision of that flag 0.50 at COMPASS's default iterations, 5 runs: 0.95).
COMPASS_DIAG_MIN_RUNS <- 5L

# R-hat and bulk ESS of one donor's draws (rows = iterations, columns = runs). Runs that all
# sit on one value agree exactly (R-hat 1, ESS NA). Constant runs on different values never
# mix (R-hat Inf).
compass_rhat_ess <- function(draws) {
  if (all(draws == draws[1L])) return(c(rhat = 1, ess = NA_real_))
  within_constant <- all(apply(draws, 2L, function(v) all(v == v[1L])))
  if (within_constant) return(c(rhat = Inf, ess = NA_real_))
  c(rhat = posterior::rhat(draws), ess = posterior::ess_bulk(draws))
}

# Per-iteration draws of the three COMPASS quantities for every donor of one fit. Returns a
# list of matrices (iterations x donors): gamma (mean over the non-null subsets of the 0/1
# subset indicator), FS and PFS (COMPASS's own formulas applied to each draw instead of to
# the mean).
compass_draws <- function(fit) {
  g <- fit$fit$gamma # donors x subsets x iterations, null subset last
  k <- dim(g)[[2L]]
  n_marker <- ncol(fit$fit$categories) - 1L
  degree <- fit$fit$categories[, "Counts"]
  # Sum of gamma over the non-null subsets, donors x iterations.
  s <- apply(g[, -k, , drop = FALSE], c(1L, 3L), sum)
  w <- degree[-k] / choose(n_marker, degree[-k]) / n_marker * (2 / (n_marker + 1))
  pfs <- apply(g[, -k, , drop = FALSE], c(1L, 3L), function(v) sum(v * w))
  list(gamma = t(s) / (k - 1L), FS = t(s) / (2^n_marker - 1), PFS = t(pfs))
}

# Per-donor diagnostics from a list of COMPASSResult objects fitted with different seeds.
compass_run_diagnostics <- function(fits, rhat_max = Inf) {
  draws <- lapply(fits, compass_draws)
  donors <- rownames(fits[[1L]]$fit$gamma)
  stopifnot("all runs have the same donors" = all(vapply(fits, function(f)
    identical(rownames(f$fit$gamma), donors), NA)))
  # Donor-level call of each run: max over non-null subsets of mean gamma.
  max_gamma <- vapply(fits, function(f) {
    mg <- f$fit$mean_gamma
    apply(mg[, -ncol(mg), drop = FALSE], 1L, max)[donors]
  }, numeric(length(donors)))
  if (anyNA(max_gamma)) stop("mean gamma has missing values; the calls cannot be compared", call. = FALSE)
  max_gamma <- matrix(max_gamma, nrow = length(donors), dimnames = list(donors, NULL))
  stat <- function(q) t(vapply(seq_along(donors), function(i)
    compass_rhat_ess(vapply(draws, function(d) d[[q]][, i], numeric(nrow(draws[[1L]][[q]])))),
    numeric(2L)))
  gm <- stat("gamma"); fs <- stat("FS"); pf <- stat("PFS")
  called <- rowSums(max_gamma > COMPASS_GAMMA_CALL)
  out <- data.frame(
    donor_id = donors,
    compass_runs = length(fits),
    compass_runs_called = called,
    compass_max_gamma_mean = rowMeans(max_gamma),
    compass_gamma_rhat = gm[, "rhat"], compass_gamma_ess = gm[, "ess"],
    compass_FS_rhat = fs[, "rhat"], compass_FS_ess = fs[, "ess"],
    compass_PFS_rhat = pf[, "rhat"], compass_PFS_ess = pf[, "ess"],
    row.names = NULL)
  # The raw call of the first run: the one cj_compass() reports. With missing_control_calibrate,
  # cj_compass() replaces compass_call by the calibrated call and keeps this one as compass_call_raw.
  out$compass_call <- max_gamma[, 1L] > COMPASS_GAMMA_CALL
  out$compass_call_raw <- out$compass_call
  out$compass_call_disagree <- out$compass_runs_called > 0L & out$compass_runs_called < length(fits)
  # Runs that make the same call as the reported one; a tie or less means the reported call is not backed.
  support <- ifelse(out$compass_call, out$compass_runs_called, length(fits) - out$compass_runs_called)
  out$compass_call_unsupported <- support <= length(fits) / 2
  out$compass_rhat_high <- pmax(out$compass_gamma_rhat, out$compass_FS_rhat,
                                out$compass_PFS_rhat, na.rm = TRUE) > rhat_max
  out$unstable_call <- out$compass_call_disagree | out$compass_rhat_high
  out
}

# One row per COMPASS run: how many donors are flagged, and why.
compass_diagnostics_summary <- function(diag, seconds, n_runs) {
  rows <- lapply(split(diag, diag$group), function(d) data.frame(
    group = d$group[1L], n_runs = n_runs, n_donors = nrow(d),
    n_unstable = sum(d$unstable_call), n_call_disagree = sum(d$compass_call_disagree),
    n_call_unsupported = sum(d$compass_call_unsupported), n_rhat_high = sum(d$compass_rhat_high),
    max_rhat = suppressWarnings(max(d$compass_gamma_rhat, d$compass_FS_rhat, d$compass_PFS_rhat,
                                    na.rm = TRUE)),
    min_ess = suppressWarnings(min(d$compass_gamma_ess, d$compass_FS_ess, d$compass_PFS_ess,
                                   na.rm = TRUE)),
    seconds = unname(seconds[d$group[1L]]), row.names = NULL))
  do.call(rbind, unname(rows))
}

# Validate n_runs / diagnostics. diagnostics = TRUE means at least 5 runs.
compass_n_runs <- function(n_runs, diagnostics) {
  if (!is.numeric(n_runs) || length(n_runs) != 1L || !is.finite(n_runs) || n_runs < 1 || n_runs != round(n_runs))
    stop("n_runs must be a finite whole number >= 1", call. = FALSE)
  n_runs <- as.integer(n_runs)
  if (isTRUE(diagnostics)) n_runs <- max(n_runs, COMPASS_DIAG_MIN_RUNS)
  if (n_runs %in% 2:3)
    warning("n_runs = ", n_runs, " gives weak diagnostics; use 5 or more", call. = FALSE)
  if (n_runs > 3L && n_runs %% 2L == 0L)
    warning("n_runs = ", n_runs, " is even: a 50/50 split of the runs counts the reported call as ",
            "unsupported (compass_call_unsupported); use an odd number", call. = FALSE)
  if (n_runs > 1L && !requireNamespace("posterior", quietly = TRUE))
    stop("The posterior package is required for n_runs > 1", call. = FALSE)
  n_runs
}
