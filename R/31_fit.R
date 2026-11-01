# Fit the M1 count model with cmdstanr and wrap the result in cj_fit.

m1_compile <- function() {
  dir <- tools::R_user_dir("cytoJointR", "cache")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  cmdstanr::cmdstan_model(system.file("stan", "m1.stan", package = "cytoJointR"), dir = dir)
}

# The list Stan receives as data.
m1_stan_data <- function(d, counts, pr) {
  list(S = nrow(d$X), K = nrow(d$B), J = ncol(d$B), M = ncol(d$Xs), P = ncol(d$Bpair),
       Pf = ncol(d$X), Ndon = length(d$donors), G = nrow(d$dv),
       use_donor = d$use_donor, use_dv = d$use_dv, centered = d$centered,
       y = unname(counts), X = unname(d$X), Xs = unname(d$Xs), W = unname(d$W),
       donor_of = d$donor_of, dv_of = d$dv_of,
       Bt = unname(t(d$B)), Bpt = unname(t(d$Bpair)),
       alpha_mean = pr$alpha_mean, alpha_sd = pr$alpha_sd, gamma_sd = pr$gamma_sd,
       theta0_sd = pr$theta0_sd, delta_sd = pr$delta_sd, tau_sd = pr$donor_sd,
       tau_u_sd = pr$tau_u_sd, lkj_eta = pr$lkj_eta)
}

#' Fit the M1 count model
#'
#' Bayesian donor-paired log-linear model of exact marker-pattern counts, fitted
#' by HMC in Stan. Marker `j` has natural parameter
#' `a_sj = alpha_j + fixed effects + b_i + r_im + u_iv`; pairwise interactions are
#' `theta0 + delta_m`. Donor effects `(b_i, r_i1..r_iM)` are multivariate normal
#' with LKJ(2) correlation; `u_iv` is a donor-by-visit basal deviation. All effects
#' have unit loadings across markers and are noncentered. A stimulated sample
#' enters the likelihood even when its visit has no control sample.
#'
#' @param x A `cj_data` object.
#' @param chains,iter_warmup,iter_sampling,seed Passed to CmdStan.
#' @param priors Named list that overrides [cj_default_priors()]. The defaults are
#'   illustrative smoke-test values, not validated.
#' @param donor_visit Include the donor-by-visit deviation `u_iv` (set `FALSE` for
#'   the ablation). Ignored, with a message, when the data hold one visit.
#' @param donor_effects Include donor effects `(b_i, r_im)`; `FALSE` gives a
#'   fixed-effects-only model.
#' @param centered Sample donor effects `(b, r)` and `u` on their natural scale
#'   instead of the noncentered default. Both give the same posterior. With many
#'   cells per sample (data-dominated donors) the noncentered form needs very small
#'   steps; in a J = 3, 24-donor, 10 000-cell benchmark the centered form with a dense
#'   metric (set automatically) took about a quarter of the time. `FALSE` keeps the
#'   SCOPE default.
#' @param ... Passed to the CmdStanModel `$sample()` method (for example
#'   `adapt_delta`). Defaults: `refresh = 0`, `parallel_chains = min(chains, 4)`,
#'   `metric = "dense_e"` if `centered` else `"diag_e"`, `sig_figs = 18` (so probabilities sum to one to machine precision).
#' @return A `cj_fit`. `$draws$p` holds probabilities for every sample
#'   (`draw_kind = "in_sample_donor_conditional"`), `$design` the design and
#'   parameter map, `$diagnostics` R-hat, ESS, divergences, treedepth and E-BFMI,
#'   `$engine_fit` the CmdStanMCMC object (with `log_lik` per sample for LOO).
#' @examples
#' \donttest{
#' sim <- cj_simulate_counts(n_donors = 6, seed = 1)
#' fit <- cj_fit(sim$data, chains = 2, iter_warmup = 300, iter_sampling = 300, seed = 1)
#' summary(fit)
#' }
#' @export
cj_fit <- function(x, chains = 4L, iter_warmup = 1000L, iter_sampling = 1000L, seed = 1L,
                   priors = list(), donor_visit = TRUE, donor_effects = TRUE, centered = FALSE, ...) {
  stopifnot(inherits(x, "cj_data"))
  d <- cj_design(x, donor_effects, donor_visit, centered)
  pr <- cj_resolve_priors(priors, d)
  args <- utils::modifyList(list(refresh = 0, parallel_chains = min(chains, 4L), sig_figs = 18,
                           metric = if (centered) "dense_e" else "diag_e"), list(...))
  engine <- do.call(m1_compile()$sample,
                    c(list(data = m1_stan_data(d, x$counts, pr), chains = chains,
                           iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                           seed = seed), args))

  pd <- posterior::as_draws_matrix(engine$draws("p"))
  idx <- do.call(rbind, lapply(regmatches(colnames(pd), gregexpr("[0-9]+", colnames(pd))), as.integer))
  n <- nrow(pd)
  p <- array(NA_real_, c(n, nrow(x$counts), ncol(x$counts)),
             dimnames = list(NULL, rownames(x$counts), colnames(x$counts)))
  p[cbind(rep(seq_len(n), ncol(pd)), rep(idx[, 1], each = n), rep(idx[, 2], each = n))] <- as.vector(pd)
  p <- pmax(p, 1e-300)
  p <- p / as.vector(apply(p, 1:2, sum))   # CSV round-trip leaves ~1e-16 drift
  draws <- new_cj_draws(p, x$samples, x$contrasts, x$bits,
                        list(engine_id = "cmdstan-hmc", model_id = "M1", seed = as.integer(seed),
                             n_draws = n, draw_kind = "in_sample_donor_conditional",
                             donor_effects = d$use_donor == 1L, donor_visit = d$use_dv == 1L,
                             centered = d$centered == 1L))
  new_cj_fit(x, draws, engine, d, pr, m1_diagnostics(engine, d))
}

m1_diagnostics <- function(engine, d) {
  ds <- engine$diagnostic_summary(quiet = TRUE)
  vars <- intersect(c("alpha", "gamma", "theta0", "delta", "tau", "rho", "tau_u"),
                    engine$metadata()$stan_variables)
  tab <- posterior::summarise_draws(engine$draws(vars), "mean", "sd",
                                    ~posterior::quantile2(.x, c(0.05, 0.95)),
                                    "rhat", "ess_bulk", "ess_tail") |> as.data.frame()
  tab$label <- d$param_map$label[match(tab$variable, d$param_map$stan)]
  list(n_divergent = sum(ds$num_divergent), n_max_treedepth = sum(ds$num_max_treedepth),
       ebfmi = ds$ebfmi, max_rhat = max(tab$rhat, na.rm = TRUE),
       min_ess_bulk = min(tab$ess_bulk, na.rm = TRUE),
       min_ess_tail = min(tab$ess_tail, na.rm = TRUE),
       time_s = engine$time()$total, summary = tab)
}

#' @export
print.cj_fit <- function(x, ...) {
  dg <- x$diagnostics
  cat("<cj_fit> M1 count model,", nrow(x$data$samples), "samples,",
      length(x$design$donors), "donors,", length(x$design$markers), "markers\n")
  cat("  stimulations:", paste(x$design$stimulations, collapse = ", "),
      "| visits:", paste(x$design$visits, collapse = ", "), "\n")
  cat("  donor effects:", x$design$use_donor == 1L, "| donor x visit:", x$design$use_dv == 1L, "\n")
  cat(sprintf("  draws: %d | time: %.1f s\n", x$draws$meta$n_draws, dg$time_s))
  cat(sprintf("  max R-hat %.3f | min ESS bulk %.0f, tail %.0f\n",
              dg$max_rhat, dg$min_ess_bulk, dg$min_ess_tail))
  cat(sprintf("  divergences %d | max-treedepth hits %d | E-BFMI min %.2f\n",
              dg$n_divergent, dg$n_max_treedepth, min(dg$ebfmi)))
  if (dg$n_divergent > 0 || dg$max_rhat > 1.01 || min(dg$ebfmi) < 0.3 || dg$min_ess_bulk < 400)
    cat("  WARNING: diagnostics outside usual thresholds (R-hat 1.01, ESS 400, E-BFMI 0.3)\n")
  invisible(x)
}

#' @export
summary.cj_fit <- function(object, ...) {
  structure(list(fit = object, table = object$diagnostics$summary[
    c("label", "variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk", "ess_tail")]),
    class = "summary.cj_fit")
}

#' @export
print.summary.cj_fit <- function(x, digits = 3, ...) {
  print(x$fit)
  cat("\nPosterior summaries (90% interval):\n")
  print(x$table[-2], digits = digits, row.names = FALSE)
  invisible(x)
}
