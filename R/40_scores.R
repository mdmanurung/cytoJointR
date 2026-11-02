# Score layer: per-draw net scores for every contrast, and their summaries.
# The score algebra of cj_score_reference() is vectorised over draws in
# score_contrast() (tested equal to it); this file also summarises.
# Docs: docs/04_SCORE_CONTRACT.md, 06_ATTRIBUTION_AND_OUTPUTS.md, 07_UNCERTAINTY_AND_GROUPS.md.

GLOBAL_SCORES <- c("RMS", "poly_absolute", "poly_marginal", "CPS", "CCS", "breadth")
CONTRIBUTION_SCORES <- c("RMS", "CPS", "CCS")
# CPS and CCS are defined on polyfunctional patterns only; RMS on any positive pattern.
CONTRIBUTION_MIN_DEGREE <- c(RMS = 1L, CPS = 2L, CCS = 2L)
SCORE_UNITS <- c(RMS = "fraction", poly_absolute = "fraction", poly_marginal = "fraction",
                 CPS = "conditional fraction", CCS = "fraction",
                 breadth = "functions per positive cell", marker_magnitude = "fraction")

# Accept a cj_draws or a cj_fit; return the cj_draws.
as_draws <- function(x) {
  if (inherits(x, "cj_fit")) return(x$draws)
  if (inherits(x, "cj_draws")) return(x)
  stop("x must be a cj_draws or a cj_fit", call. = FALSE)
}

# "marker_magnitude:IFNG" -> list(score, marker); other names have no marker.
split_score_col <- function(cn) {
  has_marker <- grepl(":", cn, fixed = TRUE)
  list(score = sub(":.*$", "", cn),
       marker = ifelse(has_marker, sub("^[^:]*:", "", cn), NA_character_))
}

# Epsilon must be named per score; a missing or negative entry is an error.
check_epsilon <- function(epsilon, needed) {
  if (!is.numeric(epsilon) || is.null(names(epsilon)))
    stop("epsilon must be a named numeric vector, one entry per score: ",
         paste(needed, collapse = ", "), call. = FALSE)
  absent <- setdiff(needed, names(epsilon))
  if (length(absent))
    stop("epsilon lacks: ", paste(absent, collapse = ", "), call. = FALSE)
  if (any(!is.finite(epsilon[needed]) | epsilon[needed] < 0))
    stop("epsilon must be finite and nonnegative", call. = FALSE)
  epsilon
}

# Per-draw engine -----------------------------------------------------------

# IPF on all draws at once: Pu, Ps are [draw x pattern], rows sum to 1.
# Same algorithm and tolerance as cj_counterfactual(); a draw stops updating
# at the sweep where it converges, so results match the scalar path. Returns
# Q [draw x pattern] and ok [draw]; ok is FALSE where the scalar path errors
# (invalid input, non-interior margins, no support, no convergence).
ipf_draws <- function(Pu, Ps, bits, tol = 1e-10, max_iter = 10000L) {
  nd <- nrow(Pu); J <- ncol(bits)
  target <- Ps %*% bits
  ok <- rowSums(!(is.finite(Pu) & Pu > 0) | !(is.finite(Ps) & Ps > 0)) == 0L
  ok <- ok & rowSums(!(target > 0 & target < 1)) == 0L
  Q <- Pu
  act <- which(ok)
  for (iter in seq_len(max_iter)) {
    if (!length(act)) break
    Qa <- Q[act, , drop = FALSE]; ta <- target[act, , drop = FALSE]
    bad <- rep(FALSE, length(act))
    for (j in seq_len(J)) {
      plus <- bits[, j] == 1L
      a <- rowSums(Qa[, plus, drop = FALSE])
      bad <- bad | !(a > 0 & a < 1)
      a[bad] <- 0.5
      Qa[, plus] <- Qa[, plus] * ta[, j] / a
      Qa[, !plus] <- Qa[, !plus] * (1 - ta[, j]) / (1 - a)
      Qa <- Qa / rowSums(Qa)
    }
    bad <- bad | rowSums(!is.finite(Qa)) > 0L
    done <- !bad & rowSums(abs(Qa %*% bits - ta) >= tol) == 0L
    Q[act, ] <- Qa
    ok[act[bad]] <- FALSE
    act <- act[!bad & !done]
  }
  ok[act] <- FALSE  # still iterating at max_iter
  list(Q = Q, ok = ok)
}

# cj_score_reference() algebra for all draws of one contrast. Returns G
# [draw x score column], status per draw and (keep_patterns) contrib
# [draw x pattern x RMS/CPS/CCS]. Failed-IPF draws stay as NA rows.
score_contrast <- function(Pu, Ps, bits, keep_patterns = FALSE) {
  Pu <- Pu / rowSums(Pu); Ps <- Ps / rowSums(Ps)
  ipf <- ipf_draws(Pu, Ps, bits)
  deg <- rowSums(bits); pos <- deg >= 1L; poly <- deg >= 2L
  dp <- Ps - Pu
  q_u <- rowSums(Pu[, pos, drop = FALSE]); q_s <- rowSums(Ps[, pos, drop = FALSE])
  cps <- Ps / q_s - Pu / q_u
  ccs <- Ps - ipf$Q
  sub <- function(M, keep) rowSums(M[, keep, drop = FALSE])
  G <- cbind(RMS = sub(dp, pos), poly_absolute = sub(dp, poly),
             poly_marginal = sub(ipf$Q - Pu, poly), CPS = sub(cps, poly),
             CCS = sub(ccs, poly), breadth = rowSums(sweep(cps, 2L, deg, `*`)),
             dp %*% bits)
  colnames(G)[-(1:6)] <- paste0("marker_magnitude:", colnames(bits))
  G[!ipf$ok, ] <- NA_real_
  status <- ifelse(ipf$ok, "ok", "ipf_failed")
  status[ipf$ok & rowSums(!is.finite(G)) > 0L] <- "score_failed"
  contrib <- NULL
  if (keep_patterns) {
    mask <- function(M, keep) { M[, !keep] <- 0; M }
    contrib <- array(NA_real_, c(nrow(Pu), nrow(bits), 3L),
                     dimnames = list(NULL, pattern_labels(bits), CONTRIBUTION_SCORES))
    contrib[, , "RMS"] <- mask(dp, pos)
    contrib[, , "CPS"] <- mask(cps, poly)
    contrib[, , "CCS"] <- mask(ccs, poly)
    contrib[!ipf$ok, , ] <- NA_real_
  }
  list(G = G, status = status, contrib = contrib)
}

# For each contrast: G [draw x score column] and per-draw failure status.
# With keep_patterns, also contrib [draw x pattern x RMS/CPS/CCS].
# Memoised on the last input: cj_scores(), cj_patterns(), cj_markers() and
# cj_visit_change() on the same draws share one pass.
# ponytail: single-slot cache; use a keyed cache if alternating objects matter.
.engine_cache <- new.env(parent = emptyenv())

score_engine <- function(d, keep_patterns = FALSE) {
  hit <- .engine_cache$val
  same <- !is.null(hit) && identical(hit$p, d$p) && identical(hit$bits, d$bits) &&
    identical(hit$contrasts, d$contrasts) && identical(hit$samples, d$samples)
  if (same && (hit$keep || !keep_patterns)) return(hit$eng)
  eng <- compute_engine(d, keep_patterns)
  assign("val", list(p = d$p, bits = d$bits, contrasts = d$contrasts,
                     samples = d$samples, keep = keep_patterns, eng = eng),
         envir = .engine_cache)
  eng
}

compute_engine <- function(d, keep_patterns) {
  p <- d$p
  ct <- d$contrasts
  s_row <- match(ct$stim_sample_id, d$samples$sample_id)
  u_row <- match(ct$control_sample_id, d$samples$sample_id)
  cols <- c(GLOBAL_SCORES, paste0("marker_magnitude:", d$markers))
  out <- lapply(seq_len(nrow(ct)), function(i) {
    if (is.na(u_row[i])) return(list(available = FALSE, cols = cols))
    # The control row is read from the same draw as the stimulated row.
    r <- score_contrast(matrix(p[, u_row[i], ], dim(p)[1L]),
                        matrix(p[, s_row[i], ], dim(p)[1L]), d$bits, keep_patterns)
    colnames(r$G) <- cols
    c(list(available = TRUE, cols = cols), r)
  })
  names(out) <- ct$contrast_id
  out
}

contrast_keys <- function(ct) {
  ct[, c("contrast_id", "donor_id", "visit", "parent_population", "stimulation")]
}

#' Per-draw net scores for every contrast
#'
#' Evaluates every score within each posterior draw, using the stimulated
#' sample and its control sample from the same draw. A control shared by
#' several stimulations therefore gives identical draws to each of its
#' contrasts. The counterfactual (IPF) is solved in every draw.
#'
#' A draw whose computation fails is kept as a row with `status` `"ipf_failed"`
#' (or `"score_failed"`) and `value` `NA`; it is never dropped. A contrast with
#' no control sample gets rows with `status` `"unavailable"`.
#'
#' @param x A `cj_draws` or `cj_fit`.
#' @return A data frame with columns `draw_id`, `contrast_id`, `score`,
#'   `marker` (only for `marker_magnitude`), `value`, `status`. Scores are
#'   `RMS`, `poly_absolute`, `poly_marginal`, `CPS`, `CCS`, `breadth` and
#'   `marker_magnitude`. Values are signed and never clipped.
#' @export
#' @examples
#' d <- cj_example_draws(n_draw = 20)
#' head(cj_score_draws(d))
cj_score_draws <- function(x) {
  d <- as_draws(x)
  eng <- score_engine(d)
  nd <- dim(d$p)[1]
  do.call(rbind, lapply(names(eng), function(id) {
    e <- eng[[id]]
    sc <- split_score_col(e$cols)
    n <- length(e$cols)
    data.frame(
      draw_id = rep(seq_len(nd), times = n), contrast_id = id,
      score = rep(sc$score, each = nd), marker = rep(sc$marker, each = nd),
      value = if (e$available) as.vector(e$G) else NA_real_,
      status = if (e$available) rep(e$status, times = n) else "unavailable")
  }))
}

# Summarise the columns of G (draws x scores) into one row per column.
# Failed draws are counted and reported, never silently removed: above
# `max_failure` the summary is withheld.
summarise_columns <- function(cols, G, bad, eps, level, max_failure) {
  # eps: one epsilon per column (named by nothing, aligned to cols).
  sc <- split_score_col(cols)
  na_row <- c(estimate = NA, lower = NA, upper = NA, epsilon = NA, p_positive = NA,
              p_negative = NA, p_negligible = NA)
  rows <- lapply(seq_along(cols), function(j) {
    if (is.null(G))
      return(c(na_row, n_draws = NA, n_failed = NA, failure_fraction = NA))
    nd <- length(bad); nf <- sum(bad)
    s <- if (nf / nd > max_failure) na_row else
      cj_summarize_signed(G[!bad, j], eps[[j]], level)[names(na_row)]
    c(s, n_draws = nd, n_failed = nf, failure_fraction = nf / nd)
  })
  m <- do.call(rbind, rows)
  status <- if (is.null(G)) "unavailable" else
    if (mean(bad) > max_failure) "failure_tolerance_exceeded" else
    if (any(bad)) "ok_with_failed_draws" else "ok"
  data.frame(score = sc$score, marker = sc$marker, as.data.frame(m, row.names = NULL),
             level = level, status = status, units = unname(SCORE_UNITS[sc$score]))
}

#' Signed score summaries per contrast
#'
#' Posterior mean, equal-tailed interval and the three epsilon probabilities
#' for every score of every contrast, from [cj_summarize_signed()].
#' `p_positive` is P(score > epsilon), `p_negative` is P(score < -epsilon) and
#' `p_negligible` is P(|score| <= epsilon); they sum to one. Intervals come
#' from the per-draw values of the score itself.
#'
#' @param x A `cj_draws` or `cj_fit`.
#' @param epsilon Named numeric vector, one nonnegative entry per score in
#'   `RMS`, `poly_absolute`, `poly_marginal`, `CPS`, `CCS`, `breadth`,
#'   `marker_magnitude`, on that score's units. There is no default.
#' @param level Probability mass of the credible interval.
#' @param max_failure Largest tolerated fraction of failed draws per contrast.
#'   Above it the summary is withheld and `status` says so. Failed draws are
#'   counted in `n_failed` and `failure_fraction`, never hidden.
#' @return A data frame, one row per contrast, score and marker, with
#'   `estimate`, `lower`, `upper`, `epsilon`, `p_positive`, `p_negative`,
#'   `p_negligible`, `n_draws`, `n_failed`, `failure_fraction`, `status`
#'   (`ok`, `ok_with_failed_draws`, `failure_tolerance_exceeded`,
#'   `unavailable`) and `units`.
#' @export
#' @examples
#' d <- cj_example_draws(n_draw = 100)
#' eps <- c(RMS = 0.005, poly_absolute = 0.005, poly_marginal = 0.005, CPS = 0.02,
#'          CCS = 0.005, breadth = 0.05, marker_magnitude = 0.005)
#' head(cj_scores(d, eps))
cj_scores <- function(x, epsilon, level = 0.9, max_failure = 0.01) {
  d <- as_draws(x)
  check_epsilon(epsilon, c(GLOBAL_SCORES, "marker_magnitude"))
  eng <- score_engine(d)
  keys <- contrast_keys(d$contrasts)
  do.call(rbind, lapply(seq_along(eng), function(i) {
    e <- eng[[i]]
    bad <- if (e$available) e$status != "ok"
    cbind(keys[i, ], summarise_columns(e$cols, e$G, bad, epsilon[split_score_col(e$cols)$score], level,
                      max_failure),
          row.names = NULL)
  }))
}

#' Marker-wise magnitude views (overlapping)
#'
#' Net marginal frequency of each functional marker per contrast. A cell
#' positive for several markers counts in each of them, so these views
#' overlap and do not sum to RMS or any other global score.
#'
#' @inheritParams cj_scores
#' @param epsilon One nonnegative number on the fraction scale.
#' @return A data frame like [cj_scores()] for `marker_magnitude`, plus
#'   `overlapping = TRUE`.
#' @export
#' @examples
#' cj_markers(cj_example_draws(n_draw = 50), epsilon = 0.005)
cj_markers <- function(x, epsilon, level = 0.9, max_failure = 0.01) {
  if (length(epsilon) != 1L) stop("epsilon must be one number", call. = FALSE)
  d <- as_draws(x)
  eng <- score_engine(d)
  keys <- contrast_keys(d$contrasts)
  eps <- c(marker_magnitude = unname(epsilon))
  check_epsilon(eps, "marker_magnitude")
  out <- do.call(rbind, lapply(seq_along(eng), function(i) {
    e <- eng[[i]]
    j <- startsWith(e$cols, "marker_magnitude:")
    bad <- if (e$available) e$status != "ok"
    cbind(keys[i, ], summarise_columns(e$cols[j], if (e$available) e$G[, j, drop = FALSE],
                                       bad, rep(eps[[1]], sum(j)), level, max_failure), row.names = NULL)
  }))
  out$overlapping <- TRUE
  out
}

#' Exact-pattern contributions to a global score
#'
#' Each exact pattern's signed contribution to RMS, CPS or CCS, per contrast.
#' Contributions are computed per draw and sum to the global score in every
#' draw, so posterior means add up to the global mean. Interval bounds and
#' medians do not add up; never sum them.
#'
#' Patterns that contribute zero by definition are left out: degree 0 for
#' all scores, and degree 1 for CPS and CCS. CPS contributions are changes in
#' positive-cell shares, so single-positive patterns still shape the
#' denominators.
#'
#' @param x A `cj_draws` or `cj_fit`.
#' @param score One or more of `"RMS"`, `"CPS"`, `"CCS"`.
#' @param epsilon Named numeric vector with one nonnegative entry per `score`.
#' @param level Probability mass of the credible interval.
#' @param max_failure As in [cj_scores()].
#' @return A data frame, one row per contrast, score and pattern, with the
#'   full state label `pattern`, `degree`, posterior means `p_U_mean` and
#'   `p_S_mean`, `estimate`, `lower`, `upper`, the epsilon probabilities and
#'   `status`. For a `cj_fit` it also has the raw counts `n_U`, `n_S` and the
#'   parent totals `N_U`, `N_S`.
#' @export
#' @examples
#' pt <- cj_patterns(cj_example_draws(n_draw = 50), "RMS", c(RMS = 0.005))
#' head(pt)
cj_patterns <- function(x, score = c("RMS", "CPS", "CCS"), epsilon, level = 0.9,
                        max_failure = 0.01) {
  score <- match.arg(score, several.ok = TRUE)
  check_epsilon(epsilon, score)
  d <- as_draws(x)
  eng <- score_engine(d, keep_patterns = TRUE)
  keys <- contrast_keys(d$contrasts)
  labels <- pattern_labels(d$bits)
  degree <- rowSums(d$bits)
  s_row <- match(d$contrasts$stim_sample_id, d$samples$sample_id)
  u_row <- match(d$contrasts$control_sample_id, d$samples$sample_id)
  out <- list()
  for (i in seq_along(eng)) for (sc in score) {
    e <- eng[[i]]
    k <- which(degree >= CONTRIBUTION_MIN_DEGREE[[sc]])
    cols <- labels[k]
    G <- if (e$available) matrix(e$contrib[, k, sc], dim(e$contrib)[1L])  # keep a matrix when one pattern remains
    cell <- summarise_columns(cols, G, if (e$available) e$status != "ok",
                             rep(epsilon[[sc]], length(cols)), level, max_failure)
    cell$score <- sc
    cell$units <- SCORE_UNITS[[sc]]
    cell$marker <- NULL
    cell$pattern <- cols
    cell$degree <- degree[k]
    if (e$available) {
      cell$p_U_mean <- unname(colMeans(matrix(d$p[, u_row[i], k], dim(d$p)[1])))
      cell$p_S_mean <- unname(colMeans(matrix(d$p[, s_row[i], k], dim(d$p)[1])))
      if (inherits(x, "cj_fit")) {
        cn <- x$data$counts
        cell$n_U <- unname(cn[u_row[i], k]); cell$n_S <- unname(cn[s_row[i], k])
        cell$N_U <- d$samples$N[u_row[i]]; cell$N_S <- d$samples$N[s_row[i]]
      }
    }
    out[[length(out) + 1L]] <- cbind(keys[i, ], cell, row.names = NULL)
  }
  out <- do.call(rbind, lapply(out, function(o) {
    for (m in c("p_U_mean", "p_S_mean", "n_U", "n_S", "N_U", "N_S"))
      if (!m %in% names(o)) o[[m]] <- NA_real_
    o
  }))
  if (!inherits(x, "cj_fit")) out[c("n_U", "n_S", "N_U", "N_S")] <- NULL
  out[order(out$contrast_id, out$score, out$degree, out$pattern), ]
}

#' Synthetic posterior draws for examples and tests
#'
#' Builds a small `cj_draws` with log-normal perturbation of known pattern
#' probabilities. Not a model fit.
#'
#' @param n_donor,n_draw Number of donors and draws.
#' @param markers Marker names.
#' @param effect Multiplier on the stimulated shift in marker log-odds. Use a
#'   negative value for a suppressive stimulation.
#' @param visits Visit labels; the second visit has twice the shift.
#' @param stimulations Stimulation names; each shares the donor's control.
#' @param missing_control Contrast ids whose `control_sample_id` is set to `NA`.
#' @param seed Random seed.
#' @return A `cj_draws` object.
#' @export
#' @examples
#' d <- cj_example_draws(n_draw = 10)
#' dim(d$p)
cj_example_draws <- function(n_donor = 3, n_draw = 200, markers = c("IFNG", "TNF", "IL2"),
                             effect = 1, visits = c("v1", "v2"),
                             stimulations = c("CMV", "SEB"), missing_control = character(),
                             seed = 1) {
  set.seed(seed)
  bits <- cj_bits(markers)
  grid <- expand.grid(donor_id = paste0("D", seq_len(n_donor)), visit = visits,
                      condition = c(CONTROL_CONDITION, stimulations),
                      stringsAsFactors = FALSE)
  grid$visit <- factor(grid$visit, levels = visits)
  grid$sample_id <- paste(grid$donor_id, grid$visit, grid$condition, sep = "_")
  samples <- data.frame(sample_id = grid$sample_id, donor_id = grid$donor_id,
                        visit = grid$visit, parent_population = "CD4",
                        condition = grid$condition, N = 1000, stringsAsFactors = FALSE)
  stim <- samples[samples$condition != CONTROL_CONDITION, ]
  contrasts <- data.frame(
    contrast_id = paste(stim$donor_id, stim$visit, stim$condition, sep = "_"),
    donor_id = stim$donor_id, visit = stim$visit, parent_population = "CD4",
    stimulation = stim$condition, stim_sample_id = stim$sample_id,
    control_sample_id = paste(stim$donor_id, stim$visit, CONTROL_CONDITION, sep = "_"),
    stringsAsFactors = FALSE)
  contrasts$control_sample_id[contrasts$contrast_id %in% missing_control] <- NA_character_
  # Log-linear base: negative main effects, one positive pairwise interaction.
  shift <- ifelse(samples$condition == CONTROL_CONDITION, 0,
                  effect * as.integer(samples$visit))
  inter <- bits[, 1] * bits[, 2]
  base <- t(vapply(shift, function(s) drop(bits %*% rep(-2 + s, length(markers))) + 0.5 * inter,
                   numeric(nrow(bits))))
  logp <- array(rep(base, each = n_draw), c(n_draw, nrow(samples), nrow(bits)))
  logp <- logp + array(rnorm(length(logp), sd = 0.1), dim(logp))
  p <- exp(logp)
  p <- p / array(apply(p, c(1, 2), sum), dim(p))
  dimnames(p) <- list(NULL, samples$sample_id, pattern_labels(bits))
  new_cj_draws(p, samples, contrasts, bits,
               list(engine_id = "synthetic", model_id = "example", seed = as.integer(seed),
                    n_draws = as.integer(n_draw), draw_kind = "in_sample_donor_conditional"))
}
