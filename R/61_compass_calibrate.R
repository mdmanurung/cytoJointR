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
COMPASS_CAL_GATE_DONORS <- 6L  # fewest pseudo donors for the small-run gate below (dev: 8 and 10 lose the TPR gain)
COMPASS_CAL_SWITCH_AUC <- 0.7  # donor runs with pseudo AUC below this switch to the same-visit population control
COMPASS_CAL_SWITCH_POS <- 5L   # ... if at least this many pseudo donors have a paired call (AUC is noisy below that)
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
compass_cal_map <- function(new, pop, paired, method = c("isotonic", "quantile")) {
  method <- match.arg(method)
  # One distinct pseudo-missing gamma carries no ranking: every row gets the paired mean.
  if (length(unique(pop)) < 2L) return(rep(mean(paired), length(new)))
  if (method == "quantile")
    return(unname(stats::quantile(paired, stats::ecdf(pop)(new), type = 1)))
  o <- order(pop)
  stats::approx(pop[o], stats::isoreg(pop[o], paired[o])$yf, xout = new, ties = mean, rule = 2)$y
}

# Option 2: piecewise linear map that sends the fallback cut to 0.5 and keeps 0 and 1.
compass_cal_fallback <- function(g) {
  stats::approx(c(0, COMPASS_CAL_FALLBACK_CUT, 1), c(0, COMPASS_CAL_GAMMA, 1), xout = g)$y
}

# Pseudo-missing fits. `fits` and `ct` are from cj_compass(); returns one row per
# pseudo-missing donor: group, donor, class of its pseudo control, gamma under the
# pseudo-missing fit and gamma under the paired fit.
compass_pseudo_fits <- function(x, ct, fits, gid, ord, nm, null_col, category_filter, verbose,
                                missing_control, max_pseudo, calibrate_seed, reps = 1L, ...) {
  todo <- unique(gid[ct$control_source != "paired" & gid %in% names(fits)])
  out <- list()
  n_failed <- 0L # pseudo fits that did not run or failed
  for (r in seq_len(reps)) {
    # A new random partition of the paired donors for each repeat.
    # One seed per group, so a group's partition does not depend on which other groups are in the run.
    pick <- lapply(todo, function(g) compass_with_seed(calibrate_seed + r - 1L, {
      i <- which(gid == g & ct$control_source == "paired")
      # By donor ID, so the draw does not depend on the row order of the data.
      p <- ct$contrast_id[i][order(ct$donor_id[i])]
      p <- p[sample.int(length(p), min(length(p), max_pseudo))]
      split(p, rep_len(seq_len(ceiling(length(p) / COMPASS_CAL_FOLD_SIZE)), length(p)))
    }))
    names(pick) <- todo
    for (k in seq_len(max(lengths(pick), 0L))) {
      forced <- unlist(lapply(pick, \(f) f[k][[1L]]), use.names = FALSE)
      if (is.null(forced)) next
      ck <- compass_resolve_controls(x, missing_control, force_missing = forced)$contrasts
      gk <- paste(ck$stimulation, ck$visit, ck$parent_population, sep = " | ")
      for (g in todo) {
        ids <- pick[[g]][k][[1L]]
        cg <- ck[gk == g, ]
        if (is.null(ids)) next
        if (!all(ids %in% cg$contrast_id)) { n_failed <- n_failed + 1L; next }
        cg <- cg[order(cg$donor_id), ]
        prep <- compass_prepare_group(cg, x, ord, nm, null_col, category_filter)
        fit <- compass_run_group(prep, cg, verbose, ...)
        if (is.character(fit)) { n_failed <- n_failed + 1L; next }
        d <- cg$donor_id[cg$contrast_id %in% ids]
        out[[length(out) + 1L]] <- data.frame(
          group = g, donor_id = d,
          cls = ifelse(startsWith(cg$control_source[match(d, cg$donor_id)], "donor_visit"),
                       "donor_visit", "population"),
          fold = k, rep = r, gamma_pseudo = unname(compass_max_gamma(fit)[d]),
          gamma_paired = unname(compass_max_gamma(fits[[g]])[d]))
      }
    }
  }
  res <- if (length(out)) do.call(rbind, out) else
    data.frame(group = character(), donor_id = character(), cls = character(), fold = integer(),
               rep = integer(), gamma_pseudo = numeric(), gamma_paired = numeric())
  attr(res, "n_failed") <- n_failed
  res
}

# Add the calibrated gamma to the tidy scores. Paired rows keep their raw gamma.
# Returns the scores with compass_gamma_cal, compass_call and calibration, and the
# per group and class table of what was done.
compass_calibrate_scores <- function(scores, gid_scores, pseudo, method, min_pseudo = COMPASS_CAL_MIN_PSEUDO,
                                     fallback = TRUE) {
  if (is.null(pseudo$donor_id)) pseudo$donor_id <- seq_len(nrow(pseudo)) # one donor per point
  scores$compass_gamma_cal <- scores$compass_max_gamma
  scores$compass_max_gamma_switch <- NA_real_
  scores$calibration <- ifelse(scores$control_source == "paired", "paired", NA_character_)
  info <- list()
  cls <- ifelse(startsWith(scores$control_source, "donor_visit"), "donor_visit", "population")
  for (key in unique(paste(gid_scores, cls)[scores$control_source != "paired"])) {
    rows <- which(paste(gid_scores, cls) == key & scores$control_source != "paired")
    ps <- pseudo[paste(pseudo$group, pseudo$cls) == key, ]
    ps <- ps[is.finite(ps$gamma_pseudo) & is.finite(ps$gamma_paired), ]
    pos <- ps$gamma_paired > COMPASS_CAL_GAMMA
    # Donors, not points: a repeated partition gives a donor several points but no new donor.
    n_donor <- length(unique(ps$donor_id))
    n_pos <- length(unique(ps$donor_id[pos]))
    ok <- n_donor >= min_pseudo && n_pos >= COMPASS_CAL_MIN_SIDE && n_donor - n_pos >= COMPASS_CAL_MIN_SIDE
    auc <- if (nrow(ps)) compass_cal_auc(ps$gamma_pseudo, pos) else NA_real_
    if (ok) {
      scores$compass_gamma_cal[rows] <- compass_cal_map(scores$compass_max_gamma[rows],
                                                        ps$gamma_pseudo, ps$gamma_paired, method)
      scores$calibration[rows] <- method
    } else if (n_donor >= COMPASS_CAL_GATE_DONORS &&
               mean(ps$gamma_pseudo > COMPASS_CAL_GAMMA) - mean(ps$gamma_paired > COMPASS_CAL_GAMMA) <= 0) {
      # Too few pseudo donors for a map, but the pseudo-missing control calls no more donors than the paired
      # control does: there is no excess to remove, and the fallback cut would only cost power.
      scores$calibration[rows] <- "raw"
    } else {
      scores$compass_gamma_cal[rows] <- if (fallback) compass_cal_fallback(scores$compass_max_gamma[rows]) else
        scores$compass_max_gamma[rows]
      scores$calibration[rows] <- if (fallback) "fallback" else "raw"
    }
    info[[key]] <- data.frame(group = sub(" (donor_visit|population)$", "", key),
                              control_class = sub("^.* ", "", key), n_rows = length(rows),
                              n_pseudo = n_donor, n_pseudo_above = n_pos,
                              pseudo_auc = auc, method = scores$calibration[rows[1L]], used_map = ok)
  }
  # Calibration only tightens: a non-paired row is never called when the raw rule would not
  # call it. Dev: this removed the extra false calls that the map added at a rare, weak effect.
  np <- scores$control_source != "paired"
  scores$compass_gamma_cal[np] <- pmin(scores$compass_gamma_cal[np], scores$compass_max_gamma[np])
  scores$compass_call <- scores$compass_gamma_cal > COMPASS_CAL_GAMMA
  list(scores = scores, info = if (length(info)) do.call(rbind, unname(info)))
}

# "donor" runs whose borrowed control cannot rank donors (pseudo AUC low, e.g. a visit effect on the
# background): refit the run with the same-visit population control (an all-zero control row is
# compared with the other donors' controls of the same stimulation, visit and parent population, so it
# shares the visit effect), calibrate those rows with the population-rule pseudo fits, and put them in
# place of the borrowed rows. Returns the changed scores and a table of what was switched.
compass_switch_groups <- function(bad, scores, gs, x, ord, nm, null_col, category_filter, verbose,
                                  cal_method, max_pseudo, calibrate_seed, reps, ...) {
  ct_pop <- compass_resolve_controls(x, "population")$contrasts
  gid_pop <- paste(ct_pop$stimulation, ct_pop$visit, ct_pop$parent_population, sep = " | ")
  done <- list()
  for (g in bad) {
    cg <- ct_pop[gid_pop == g, ]
    cg <- cg[order(cg$donor_id), ]
    fit <- compass_run_group(compass_prepare_group(cg, x, ord, nm, null_col, category_filter), cg, verbose, ...)
    if (is.character(fit)) {
      done[[g]] <- data.frame(group = g, n_rows = 0L, fit_status = paste("failed:", fit), n_failed = NA_integer_,
                              n_pseudo = NA_integer_, pseudo_auc = NA_real_, method = NA_character_)
      next
    }
    ps <- compass_pseudo_fits(x, cg, stats::setNames(list(fit), g), rep(g, nrow(cg)), ord, nm, null_col,
                              category_filter, verbose, missing_control = "population",
                              max_pseudo = max_pseudo, calibrate_seed = calibrate_seed, reps = reps, ...)
    in_g <- gs == g
    sc <- scores[in_g, ]
    sc$compass_max_gamma <- unname(compass_max_gamma(fit)[sc$donor_id])
    nonp <- sc$control_source != "paired"
    sc$control_source[nonp] <- "population"
    r <- compass_calibrate_scores(sc, rep(g, nrow(sc)), ps, cal_method)
    take <- which(in_g)[startsWith(scores$control_source[in_g], "donor_visit")]
    from <- match(scores$donor_id[take], sc$donor_id)
    scores$compass_max_gamma_switch[take] <- sc$compass_max_gamma[from]
    scores$compass_gamma_cal[take] <- r$scores$compass_gamma_cal[from]
    scores$calibration[take] <- "switched"
    # The population calibration of this group: all its non-paired rows share one population class.
    i <- r$info[r$info$control_class == "population", ][1L, ]
    done[[g]] <- data.frame(group = g, n_rows = length(take), fit_status = "fitted",
                            n_failed = attr(ps, "n_failed"), n_pseudo = i$n_pseudo, pseudo_auc = i$pseudo_auc,
                            method = i$method)
  }
  scores$compass_call <- scores$compass_gamma_cal > COMPASS_CAL_GAMMA
  list(scores = scores, table = if (length(done)) do.call(rbind, unname(done)))
}
