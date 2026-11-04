# Native MIMOSA2 model (ADR 0004): baseline-adjusted positivity calls for one cell subset.
# Written from the design in docs/adr/0004-mimosa2-native.md; no upstream code is used.

MIMOSA2_SAMPLES <- c("s1", "u1", "s0", "u0")
MIMOSA2_LABELS <- c("R1", "R2", "R3", "R4", "NR1", "NR2", "NR3", "NSR1", "NSR2", "NSR3", "NSR4")

# Each component: a list of terms, each term = c(level, samples pooled into one beta-binomial term).
MIMOSA2_COMPONENTS <- list(
  R1 = list(list(1, "s1"), list(2, "s0"), list(3, "u1"), list(4, "u0")),
  R2 = list(list(1, "s1"), list(3, "u1"), list(4, c("s0", "u0"))),
  R3 = list(list(1, c("s1", "s0")), list(3, "u0"), list(4, "u1")),
  R4 = list(list(1, "s1"), list(2, "s0"), list(4, c("u1", "u0"))),
  NR1 = list(list(3, c("s1", "u1")), list(4, c("s0", "u0"))),
  NR2 = list(list(3, c("s1", "u1")), list(2, "s0"), list(4, "u0")),
  NR3 = list(list(4, c("s1", "u1", "s0", "u0"))),
  NSR1 = list(list(1, c("u1", "u0")), list(4, c("s1", "s0"))),
  NSR2 = list(list(1, "u1"), list(4, c("s1", "s0", "u0"))),
  NSR3 = list(list(1, "u0"), list(4, c("s1", "u1", "s0"))),
  NSR4 = list(list(2, c("s1", "s0")), list(4, c("u1", "u0")))
)

# Beta-binomial term for a pooled count n of N cells, beta(mu * phi, (1 - mu) * phi) proportion.
# The binomial coefficient is left out: it is the product over the pooled samples, the same in every
# component, so it cancels in the component probabilities. A sample with N = 0 gives exactly 0.
mimosa2_bb <- function(n, N, mu, phi) {
  a <- mu * phi
  b <- (1 - mu) * phi
  out <- lgamma(a + n) + lgamma(b + N - n) - lgamma(a + b + N) - (lgamma(a) + lgamma(b) - lgamma(a + b))
  ifelse(N == 0, 0, out)
}

# theta: list(mu, log_precision), four values each (mu decreasing), or the 8 unconstrained
# parameters of mimosa2_unpack(). n, N: donors x 4 matrices, columns s1, u1, s0, u0.
# Returns the donors x 11 matrix of component log-likelihoods (no one-sided rule applied).
mimosa2_component_loglik <- function(theta, n, N) {
  if (!is.list(theta)) theta <- mimosa2_unpack(theta)
  n <- as.matrix(n); N <- as.matrix(N)
  out <- matrix(0, nrow(n), length(MIMOSA2_COMPONENTS),
                dimnames = list(rownames(n), MIMOSA2_LABELS))
  for (k in seq_along(MIMOSA2_COMPONENTS)) {
    for (term in MIMOSA2_COMPONENTS[[k]]) {
      s <- term[[2]]; lev <- term[[1]]
      out[, k] <- out[, k] + mimosa2_bb(rowSums(n[, s, drop = FALSE]), rowSums(N[, s, drop = FALSE]),
                                        theta$mu[lev], exp(theta$log_precision[lev]))
    }
  }
  out
}

# par = (a, b1, b2, b3, c1..c4): logit(mu4) = a, logit(mu_k) = logit(mu_{k+1}) + exp(b_k) so the means
# are strictly ordered; log precision nu_k = exp(c_k) > 0 (the Gamma prior lives on nu).
mimosa2_unpack <- function(par) {
  lm <- numeric(4)
  lm[4] <- par[1]
  for (k in 3:1) lm[k] <- lm[k + 1] + exp(par[1 + k])
  list(mu = stats::plogis(lm), log_precision = pmin(exp(par[5:8]), 20))
}

mimosa2_pack <- function(mu, nu) {
  lm <- stats::qlogis(mu)
  c(lm[4], log(pmax(lm[1:3] - lm[2:4], 1e-3)), log(nu))
}

# Proportion-scale contrasts that carry the Normal penalties (mean 1e-4, SD 5e-4).
mimosa2_log_prior <- function(theta) {
  mu <- theta$mu
  contrasts <- c(mu[1] - mu[3] - mu[2] + mu[4], mu[1] - mu[3], mu[2] - mu[4], mu[1] - mu[2], mu[3] - mu[4])
  sum(stats::dgamma(theta$log_precision, shape = 11 / 4, rate = 1 / 2, log = TRUE)) +
    sum(stats::dnorm(contrasts, mean = 1e-4, sd = 5e-4, log = TRUE))
}

# Control proportions with a missing control (total 0) replaced by the donor's other control.
mimosa2_controls <- function(n, N) {
  p <- n / pmax(N, 1)
  miss1 <- N[, "u1"] == 0; miss0 <- N[, "u0"] == 0
  list(u1 = ifelse(miss1 & !miss0, p[, "u0"], p[, "u1"]), u0 = ifelse(miss0 & !miss1, p[, "u1"], p[, "u0"]))
}

mimosa2_delta <- function(n, N) {
  p <- n / pmax(N, 1)
  u <- mimosa2_controls(n, N)
  unname(p[, "s1"] - u$u1 - p[, "s0"] + u$u0)
}

# Newton et al. (2004) q-value: running mean of 1 - P(response) over donors ordered by decreasing P.
# Tied donors share the value at the end of their tie group.
mimosa2_qvalue <- function(p) {
  if (!all(is.finite(p))) return(rep(NA_real_, length(p)))
  o <- order(p, decreasing = TRUE)
  ps <- p[o]
  cm <- cumsum(1 - ps) / seq_along(ps)
  q <- cm[length(ps) - match(ps, rev(ps)) + 1L]
  out <- numeric(length(p)); out[o] <- q
  names(out) <- names(p)
  out
}

mimosa2_init <- function(n, N) {
  p <- n / pmax(N, 1)
  u <- mimosa2_controls(n, N)
  resp <- mimosa2_delta(n, N) > 0 & p[, "s1"] > u$u1
  z <- matrix(0, nrow(n), 11L, dimnames = list(rownames(n), MIMOSA2_LABELS))
  z[resp, 1:4] <- 1 / 4
  z[!resp, 5:11] <- 1 / 7
  pool <- function(col, rows) if (!any(rows)) 0 else sum(n[rows, col]) / max(sum(N[rows, col]), 1)
  rr <- if (any(resp)) resp else rep(TRUE, nrow(n))
  mu <- sort(pmin(pmax(c(pool("s1", rr), pool("s0", rr), pool("u1", N[, "u1"] > 0),
                         pool("u0", N[, "u0"] > 0)), 1e-6), 1 - 1e-6), decreasing = TRUE)
  list(z = z, mu = mu)
}

# EM to the posterior mode. n, N: donors x 4 matrices (columns s1, u1, s0, u0); a missing control
# has n = N = 0 and the donor's other control stands in for it in the initialisation and the
# one-sided rule.
mimosa2_fit <- function(n, N, tol = 1e-5, maxit = 100L) {
  n <- as.matrix(n); N <- as.matrix(N)
  stopifnot(identical(colnames(n), MIMOSA2_SAMPLES), identical(colnames(N), MIMOSA2_SAMPLES))
  init <- mimosa2_init(n, N)
  z <- init$z
  par <- mimosa2_pack(init$mu, rep(3.5, 4))
  # The one-sided rule: responder components get likelihood 0 when delta is negative.
  neg <- mimosa2_delta(n, N) < 0
  masked <- function(ll) { ll[neg, 1:4] <- -Inf; ll }
  neg_post <- function(par, z) {
    theta <- mimosa2_unpack(par)
    v <- -(sum(z * mimosa2_component_loglik(theta, n, N)) + mimosa2_log_prior(theta))
    if (is.finite(v)) v else 1e10
  }
  observed <- function(par, pi) {
    theta <- mimosa2_unpack(par)
    lp <- sweep(masked(mimosa2_component_loglik(theta, n, N)), 2, log(pi), "+")
    m <- apply(lp, 1, max)
    sum(m + log(rowSums(exp(lp - m)))) + mimosa2_log_prior(theta)
  }
  pi <- colMeans(z)
  obj <- -Inf
  converged <- FALSE
  it <- 0L
  while (it < maxit) {
    it <- it + 1L
    par <- stats::optim(par, neg_post, z = z, method = "BFGS", control = list(maxit = 100))$par
    lp <- sweep(masked(mimosa2_component_loglik(mimosa2_unpack(par), n, N)), 2, log(pi), "+")
    m <- apply(lp, 1, max)
    z <- exp(lp - m)
    z <- z / rowSums(z)
    pi <- colMeans(z)
    new <- observed(par, pi)
    if (is.finite(obj) && abs(new - obj) < tol * abs(obj)) { obj <- new; converged <- TRUE; break }
    obj <- new
  }
  theta <- mimosa2_unpack(par)
  p_response <- unname(rowSums(z[, 1:4, drop = FALSE]))
  names(p_response) <- rownames(n)
  structure(list(z = z, pi = stats::setNames(pi, MIMOSA2_LABELS), mu = theta$mu,
                 log_precision = theta$log_precision, p_response = p_response,
                 q = mimosa2_qvalue(p_response), converged = converged, iterations = it,
                 objective = obj),
            class = "cj_mimosa2_fit")
}

#' @export
print.cj_mimosa2_fit <- function(x, ...) {
  cat("cj_mimosa2_fit:", nrow(x$z), "donors,", if (x$converged) "converged" else "not converged",
      "after", x$iterations, "iteration(s)\n")
  invisible(x)
}

#' Baseline-adjusted positivity calls with a native MIMOSA2 model
#'
#' A native implementation of the MIMOSA2 design (ADR 0004,
#' `docs/adr/0004-mimosa2-native.md`), written from that specification and
#' sharing no code with the upstream package. It extends the mixture of
#' beta-binomials of Finak et al. (2014, MIMOSA). One fit covers one
#' stimulation, parent population and follow-up visit. Each donor enters with
#' four samples: stimulated and control at the baseline visit, and stimulated
#' and control at the follow-up visit. The model asks whether the
#' stimulated-minus-control proportion rose from baseline to follow-up.
#'
#' @details
#' **One subset per fit.** The model uses one count per sample, not the
#' Boolean patterns COMPASS uses. `positive = "any"` counts the cells positive
#' for at least one marker (all cells minus the all-negative pattern).
#' `positive = "<marker>"` counts the cells positive for that marker, whatever
#' the other markers do.
#'
#' **The model.** Four proportion levels, each a beta distribution with its own
#' mean and precision, ordered by construction. Eleven components say which of a
#' donor's four samples share a level; components 1 to 4 are responders. EM
#' finds the posterior mode under the priors of the ADR. `mimosa2_p_response` is
#' the posterior mass of the four responder components. Those components are
#' one-sided: a donor whose `delta` (`ps1 - pu1 - ps0 + pu0`, follow-up minus
#' baseline of the stimulated-minus-control proportion) is negative cannot be a
#' responder. `mimosa2_fdr` is the q-value of Newton et al. (2004). Its
#' multiple-testing family is the donors of one fit, so q-values are not
#' comparable across fits. `mimosa2_response` is `mimosa2_fdr < fdr`.
#'
#' **Deviations from upstream MIMOSA2** (the ADR gives the reasons): mixing
#' weights enter as `log(pi)`; every level has its own precision; level means
#' are ordered by a monotone parametrisation, not by sorting; the one-sided rule
#' is a hard zero; the optimiser is `stats::optim(method = "BFGS")`; the
#' initialisation is data-driven; the stopping rule is a relative change of the
#' log posterior below `tol`.
#'
#' **Missing samples.** With `missing_control = "drop"` (default) a donor needs
#' all four samples for a fit. Otherwise it is left out and listed in `excluded`
#' with the reason. With `"marginal"`, a donor with exactly one missing control
#' enters: that sample's count and total are set to 0, so its term drops out of
#' the likelihood, and the initialisation and the one-sided rule use the donor's
#' other control in its place. This assumes a background that is stable across
#' visits, and rests on whether the stimulated sample sits at the high or the
#' low level, an absolute rather than a within-donor contrast. **It is
#' experimental and not yet calibrated.** Donors missing a stimulated sample, or
#' both controls, stay excluded. `scores$control_source` is `"paired"`,
#' `"marginal: u1 missing"` or `"marginal: u0 missing"`.
#'
#' **Fit status.** The model has 8 level parameters and 11 components, so it
#' needs many donors per fit (dozens at least). `groups$status` is `"fitted"`,
#' `"not converged"` (EM hit `maxit`; the estimates are still returned),
#' `"failed: <message>"` (an error was caught), or `"no estimate"` when the
#' responder probabilities are not finite.
#'
#' @param x A `cj_data` (or a `cj_fit`, whose data is used).
#' @param positive `"any"` (default) or one marker name; see *Details*.
#' @param baseline_visit The baseline visit. Defaults to the first level of
#'   `samples$visit`. Every other visit is a follow-up and gets its own fits.
#' @param fdr Threshold for `mimosa2_response`, on `mimosa2_fdr`. 0.01 is the
#'   usual MIMOSA2 default.
#' @param tol,maxit EM stopping rule: stop when the log posterior changes by
#'   less than `tol` (relative), or after `maxit` iterations.
#' @param missing_control `"drop"` (default) or `"marginal"` (experimental); see
#'   *Details*.
#' @return A `cj_mimosa2` list with `fits` (named list of `cj_mimosa2_fit`:
#'   `z`, `pi`, `mu`, `log_precision`, `p_response`, `q`, `converged`,
#'   `iterations`, `objective`), `scores` (one row per donor and fit: the four
#'   positive counts `ns1`, `nu1`, `ns0`, `nu0`, the four totals `Ns1`, `Nu1`,
#'   `Ns0`, `Nu0`, `delta`, `mimosa2_p_response`, `mimosa2_fdr`,
#'   `mimosa2_response`, `control_source`), `groups` (one row per fit with its
#'   status) and `excluded` (donors left out, with the reason).
#' @references Finak G et al. (2014). Mixture models for single-cell assays with
#'   applications to vaccine studies. *Biostatistics* 15(1):87-101.
#'
#'   Newton MA et al. (2004). Detecting differential gene expression with a
#'   semiparametric hierarchical mixture method. *Biostatistics* 5(2):155-176.
#'
#'   Design and deviations: `docs/adr/0004-mimosa2-native.md` (ADR 0004).
#' @export
#' @examples
#' \donttest{
#' d <- cj_example_draws(n_donor = 30, markers = c("IFNG", "TNF"), visits = c("pre", "post"),
#'                       stimulations = "CMV", n_draw = 5)
#' counts <- t(vapply(seq_len(nrow(d$samples)), function(i)
#'   rmultinom(1, d$samples$N[i], colMeans(d$p[, i, ]))[, 1], numeric(4)))
#' storage.mode(counts) <- "integer"
#' dimnames(counts) <- list(d$samples$sample_id, dimnames(d$p)[[3]])
#' x <- cytoJointR:::new_cj_data(counts, d$samples, d$contrasts, d$bits)
#' mm <- cj_mimosa2(x)
#' mm$scores
#' }
cj_mimosa2 <- function(x, positive = "any", baseline_visit = NULL, fdr = 0.01, tol = 1e-5, maxit = 100L,
                       missing_control = c("drop", "marginal")) {
  missing_control <- match.arg(missing_control)
  if (inherits(x, "cj_fit")) x <- x$data
  if (!inherits(x, "cj_data")) stop("x must be a cj_data or cj_fit", call. = FALSE)
  if (!(length(positive) == 1L && positive %in% c("any", x$markers)))
    stop("positive must be \"any\" or one of the markers: ", paste(x$markers, collapse = ", "), call. = FALSE)
  visits <- levels(x$samples$visit)
  if (is.null(baseline_visit)) baseline_visit <- visits[1]
  if (!(length(baseline_visit) == 1L && baseline_visit %in% visits))
    stop("baseline_visit must be one of the visit levels: ", paste(visits, collapse = ", "), call. = FALSE)
  followups <- setdiff(visits, baseline_visit)
  if (!length(followups)) stop("MIMOSA2 needs a follow-up visit besides the baseline", call. = FALSE)
  marginal <- missing_control == "marginal"
  if (marginal)
    message("missing_control = \"marginal\" is experimental and not yet calibrated")

  # Positive cells of every sample: all but the all-negative pattern, or one marker's patterns.
  cols <- if (positive == "any") rowSums(x$bits) > 0 else x$bits[, positive] == 1L
  pos <- rowSums(x$counts[, cols, drop = FALSE])
  ct <- x$contrasts
  ct$visit <- as.character(ct$visit)
  base <- ct[ct$visit == baseline_visit, ]
  key <- function(d) paste(d$donor_id, d$stimulation, d$parent_population, sep = " | ")

  fits <- list(); groups <- list(); tidy <- list(); excluded <- list()
  for (v in followups) {
    fu <- ct[ct$visit == v, ]
    for (g in unique(paste(fu$stimulation, fu$parent_population, sep = " | "))) {
      cf <- fu[paste(fu$stimulation, fu$parent_population, sep = " | ") == g, ]
      cb <- base[match(key(cf), key(base)), ]
      miss_u1 <- is.na(cf$control_sample_id)
      miss_u0 <- is.na(cb$control_sample_id)
      reason <- if (marginal)
        ifelse(is.na(cb$contrast_id), "no stimulated sample at baseline",
        ifelse(miss_u1 & miss_u0, "no control sample at either visit", NA_character_))
      else
        ifelse(is.na(cb$contrast_id), "no stimulated sample at baseline",
        ifelse(is.na(cb$control_sample_id), "no control sample at baseline",
        ifelse(is.na(cf$control_sample_id), "no control sample at follow-up", NA_character_)))
      gid <- paste(g, v, sep = " | ")
      if (any(!is.na(reason)))
        excluded[[gid]] <- data.frame(donor_id = cf$donor_id, stimulation = cf$stimulation,
                                      parent_population = cf$parent_population, visit = v,
                                      reason = reason)[!is.na(reason), ]
      keep <- is.na(reason)
      cf <- cf[keep, ]; cb <- cb[keep, ]
      control_source <- ifelse(is.na(cf$control_sample_id), "marginal: u1 missing",
                        ifelse(is.na(cb$control_sample_id), "marginal: u0 missing", "paired"))
      ids <- list(s1 = cf$stim_sample_id, u1 = cf$control_sample_id,
                  s0 = cb$stim_sample_id, u0 = cb$control_sample_id)
      # A missing sample has count 0 and total 0, so its term drops out of the likelihood.
      Ntot <- vapply(ids, function(i) as.numeric(x$samples$N[match(i, x$samples$sample_id)]),
                     numeric(nrow(cf)))
      Ntot <- matrix(Ntot, nrow = nrow(cf), dimnames = list(cf$donor_id, names(ids)))
      n <- vapply(ids, function(i) unname(pos[i]), numeric(nrow(cf)))
      n <- matrix(n, nrow = nrow(cf), dimnames = list(cf$donor_id, names(ids)))
      Ntot[is.na(Ntot)] <- 0; n[is.na(n)] <- 0
      grp <- data.frame(group = gid, stimulation = cf$stimulation[1], parent_population = cf$parent_population[1],
                        visit = v, baseline_visit = baseline_visit, n_donors = nrow(cf))
      if (nrow(cf) < 2L) {
        groups[[gid]] <- cbind(grp, status = "failed: fewer than 2 donors with the samples needed")
        next
      }
      fit <- tryCatch(mimosa2_fit(n, Ntot, tol = tol, maxit = maxit),
                      error = function(e) conditionMessage(e))
      if (is.character(fit)) {
        groups[[gid]] <- cbind(grp, status = paste("failed:", fit))
        next
      }
      p_resp <- fit$p_response
      ok <- all(is.finite(p_resp))
      groups[[gid]] <- cbind(grp, status = if (!ok) "no estimate" else if (!fit$converged) "not converged" else "fitted")
      fits[[gid]] <- fit
      q <- if (ok) fit$q else rep(NA_real_, nrow(cf))
      tidy[[gid]] <- data.frame(
        donor_id = cf$donor_id, stimulation = cf$stimulation, parent_population = cf$parent_population,
        visit = v, baseline_visit = baseline_visit, positive = positive,
        ns1 = unname(n[, "s1"]), nu1 = unname(n[, "u1"]), ns0 = unname(n[, "s0"]), nu0 = unname(n[, "u0"]),
        Ns1 = Ntot[, "s1"], Nu1 = Ntot[, "u1"], Ns0 = Ntot[, "s0"], Nu0 = Ntot[, "u0"],
        delta = mimosa2_delta(n, Ntot),
        mimosa2_p_response = if (ok) unname(p_resp) else NA_real_,
        mimosa2_fdr = unname(q), mimosa2_response = unname(q) < fdr,
        control_source = control_source, row.names = NULL)
    }
  }
  bind <- function(l, empty) if (length(l)) do.call(rbind, unname(l)) else empty
  excluded <- bind(excluded, data.frame(donor_id = character(), stimulation = character(),
                                        parent_population = character(), visit = character(),
                                        reason = character()))
  if (nrow(excluded))
    message(nrow(excluded), " donor-fit(s) left out of MIMOSA2 (",
            if (marginal) "a stimulated sample and one control are needed" else "all four samples are needed",
            "); see $excluded")
  structure(list(fits = fits, scores = bind(tidy, NULL), groups = bind(groups, NULL), excluded = excluded),
            class = "cj_mimosa2")
}

#' @export
print.cj_mimosa2 <- function(x, ...) {
  cat("cj_mimosa2:", length(x$fits), "fit(s),", NROW(x$scores), "donor row(s),",
      sum(x$scores$mimosa2_response, na.rm = TRUE), "responder call(s);", nrow(x$excluded), "left out\n")
  invisible(x)
}
