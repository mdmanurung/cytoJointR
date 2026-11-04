# Label-free calibration of COMPASS calls for rows without a paired control
# (docs/adr/0003). Each paired donor of a run is pseudo-removed (its control set
# to missing, resolved by the same `missing_control` rule), the group is refitted,
# and the pseudo-missing gamma is paired with the donor's paired-control gamma
# from the real fit. No responder labels are used: the target is COMPASS's own
# paired call.

# A priori constants (not tuned on any labels).
COMPASS_CAL_GAMMA <- 0.5      # the call: max mean gamma above this
COMPASS_CAL_FOLD_SIZE <- 5L   # paired donors pseudo-removed per refit
COMPASS_CAL_MIN_PSEUDO <- 15L # fewest pseudo-missing donors per map (dev N-sweep: 10-14 paired donors gave FPR +0.11 over paired)
COMPASS_CAL_MIN_SIDE <- 3L    # fewest pseudo donors with paired gamma above, and not above, 0.5
COMPASS_CAL_MIN_AUC <- 0.8    # below this the pseudo gamma cannot rank donors; warn
COMPASS_CAL_ALPHA <- 0.05     # "null_quantile": target false-call rate among paired non-calls
# Fallback cut on raw max gamma (Option 2). From simulation only: the cut at which
# the non-paired false-positive rate meets the paired rate on the dev seeds of
# runs/2026-11-07_compass-missing-control (see ADR 0003). It depends on that design.
COMPASS_CAL_FALLBACK_CUT <- 0.85

# Max over patterns of mean gamma per donor; the last column is the null pattern.
compass_max_gamma <- function(fit) {
  mg <- fit$fit$mean_gamma
  apply(mg[, -ncol(mg), drop = FALSE], 1L, max)
}

# Run `expr` with a fixed RNG seed and give the caller's RNG state back.
compass_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
    suppressWarnings(rm(".Random.seed", envir = globalenv())))
  set.seed(seed)
  expr
}

# AUC by rank: probability that a pseudo donor above 0.5 at the paired fit has a
# higher pseudo-missing gamma than one that is not.
compass_cal_auc <- function(pop, pos) {
  if (!any(pos) || all(pos)) return(NA_real_)
  r <- rank(pop)
  n1 <- sum(pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * sum(!pos))
}

# Map raw max gamma `new` onto the paired scale. `pop` and `paired` are the pseudo pairs.
# "isotonic": monotone regression of paired gamma on pseudo-missing gamma.
# "quantile": the paired quantile at the rank of `new` among the pseudo-missing gammas.
# "isotonic_call": the same on the paired call (paired gamma above 0.5) instead of paired gamma.
# "null_quantile": the empirical null. Pseudo donors whose paired gamma is not above 0.5 are the
# null; the cut is their (1 - alpha) quantile of pseudo-missing gamma, never below 0.5, and a
# piecewise-linear map sends that cut to 0.5 (0 and 1 stay).
compass_cal_map <- function(new, pop, paired, method = c("isotonic", "quantile", "isotonic_call", "null_quantile"),
                            alpha = COMPASS_CAL_ALPHA) {
  method <- match.arg(method)
  if (method == "null_quantile") {
    cut <- min(max(COMPASS_CAL_GAMMA, unname(stats::quantile(pop[paired <= COMPASS_CAL_GAMMA], 1 - alpha, type = 1))), 0.999)
    return(stats::approx(c(0, cut, 1), c(0, COMPASS_CAL_GAMMA, 1), xout = new)$y)
  }
  if (method == "isotonic_call") paired <- as.numeric(paired > COMPASS_CAL_GAMMA)
  if (method == "quantile")
    return(unname(stats::quantile(paired, stats::ecdf(pop)(new), type = 1)))
  # One distinct pseudo-missing gamma carries no ranking: every row gets the paired mean.
  if (length(unique(pop)) < 2L) return(rep(mean(paired), length(new)))
  o <- order(pop)
  iso <- stats::isoreg(pop[o], paired[o])
  stats::approx(pop[o], iso$yf, xout = new, ties = mean, rule = 2)$y
}

# Option 2: piecewise linear map that sends the fallback cut to 0.5 and keeps 0 and 1.
compass_cal_fallback <- function(g) {
  stats::approx(c(0, COMPASS_CAL_FALLBACK_CUT, 1), c(0, COMPASS_CAL_GAMMA, 1), xout = g)$y
}

# Pseudo-missing fits. `fits` and `ct` are from cj_compass(); returns one row per
# pseudo-missing donor: group, donor, class of its pseudo control, gamma under the
# pseudo-missing fit and gamma under the paired fit.
compass_pseudo_fits <- function(x, ct, fits, gid, ord, nm, null_col, category_filter, verbose,
                                missing_control, max_pseudo, seed, ...) {
  todo <- unique(gid[ct$control_source != "paired" & gid %in% names(fits)])
  pick <- compass_with_seed(seed, lapply(todo, function(g) {
    p <- ct$contrast_id[gid == g & ct$control_source == "paired"]
    p <- p[sample.int(length(p), min(length(p), max_pseudo))]
    split(p, rep_len(seq_len(ceiling(length(p) / COMPASS_CAL_FOLD_SIZE)), length(p)))
  }))
  names(pick) <- todo
  out <- list()
  for (k in seq_len(max(lengths(pick), 0L))) {
    forced <- unlist(lapply(pick, \(f) f[k][[1L]]), use.names = FALSE)
    if (is.null(forced)) next
    ck <- compass_resolve_controls(x, missing_control, force_missing = forced)$contrasts
    gk <- paste(ck$stimulation, ck$visit, ck$parent_population, sep = " | ")
    for (g in todo) {
      ids <- pick[[g]][k][[1L]]
      cg <- ck[gk == g, ]
      if (is.null(ids) || !all(ids %in% cg$contrast_id)) next
      cg <- cg[order(cg$donor_id), ]
      prep <- compass_prepare_group(cg, x, ord, nm, null_col, category_filter)
      fit <- compass_run_group(prep, cg, verbose, ...)
      if (is.character(fit)) next
      d <- cg$donor_id[cg$contrast_id %in% ids]
      out[[length(out) + 1L]] <- data.frame(
        group = g, donor_id = d,
        cls = ifelse(startsWith(cg$control_source[match(d, cg$donor_id)], "donor_visit"),
                     "donor_visit", "population"),
        fold = k, gamma_pseudo = unname(compass_max_gamma(fit)[d]),
        gamma_paired = unname(compass_max_gamma(fits[[g]])[d]))
    }
  }
  if (length(out)) do.call(rbind, out) else NULL
}

# Add the calibrated gamma to the tidy scores. Paired rows keep their raw gamma.
# Returns the scores with compass_gamma_cal, compass_call and calibration, and the
# per group and class table of what was done.
compass_calibrate_scores <- function(scores, gid_scores, pseudo, method, alpha = COMPASS_CAL_ALPHA,
                                     min_pseudo = COMPASS_CAL_MIN_PSEUDO, null_when_sparse = FALSE, fallback = TRUE) {
  scores$compass_gamma_cal <- scores$compass_max_gamma
  scores$calibration <- ifelse(scores$control_source == "paired", "paired", NA_character_)
  info <- list()
  cls <- ifelse(startsWith(scores$control_source, "donor_visit"), "donor_visit", "population")
  for (key in unique(paste(gid_scores, cls)[scores$control_source != "paired"])) {
    rows <- which(paste(gid_scores, cls) == key & scores$control_source != "paired")
    ps <- pseudo[paste(pseudo$group, pseudo$cls) == key, ]
    pos <- ps$gamma_paired > COMPASS_CAL_GAMMA
    enough <- nrow(ps) >= min_pseudo
    ok <- enough && sum(pos) >= COMPASS_CAL_MIN_SIDE && sum(!pos) >= COMPASS_CAL_MIN_SIDE
    # Too few paired calls among the pseudo donors for a map, but the null (paired non-calls) is there:
    # the cut stays at 0.5 unless the null itself produces false calls.
    use_null <- null_when_sparse && !ok && enough && sum(!pos) >= COMPASS_CAL_MIN_SIDE
    auc <- if (nrow(ps)) compass_cal_auc(ps$gamma_pseudo, pos) else NA_real_
    if (ok || use_null) {
      m <- if (ok) method else "null_quantile"
      scores$compass_gamma_cal[rows] <- compass_cal_map(scores$compass_max_gamma[rows],
                                                        ps$gamma_pseudo, ps$gamma_paired, m, alpha)
      scores$calibration[rows] <- m
    } else {
      scores$compass_gamma_cal[rows] <- if (fallback) compass_cal_fallback(scores$compass_max_gamma[rows]) else
        scores$compass_max_gamma[rows]
      scores$calibration[rows] <- if (fallback) "fallback" else "raw"
    }
    info[[key]] <- data.frame(group = sub(" (donor_visit|population)$", "", key),
                              control_class = sub("^.* ", "", key), n_rows = length(rows),
                              n_pseudo = nrow(ps), n_pseudo_above = sum(pos),
                              pseudo_auc = auc, method = scores$calibration[rows[1L]])
  }
  # Calibration only tightens: a non-paired row is never called when the raw rule would not
  # call it. Dev: this removed the extra false calls that the map added at a rare, weak effect.
  np <- scores$control_source != "paired"
  scores$compass_gamma_cal[np] <- pmin(scores$compass_gamma_cal[np], scores$compass_max_gamma[np])
  scores$compass_call <- scores$compass_gamma_cal > COMPASS_CAL_GAMMA
  list(scores = scores, info = if (length(info)) do.call(rbind, unname(info)))
}
