# Fit the M1 count model with cmdstanr and wrap the result in cj_fit.

# The executable is named by the hash of the Stan file, so a worktree with a different
# m1.stan never overwrites the executable another checkout is using.
# Parallel test workers share the cache, so nobody compiles into it directly:
# build in a per-process temp dir, copy beside the target, then rename into place.
# ponytail: rename is atomic only within one filesystem, and NFS semantics apply
# (a stale handle is possible); use a lock file if NFS races are ever seen.
m1_compile <- function() {
  stan <- system.file("stan", "m1.stan", package = "cytoJointR")
  dir <- tools::R_user_dir("cytoJointR", "cache")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  exe <- file.path(dir, paste0("m1_", substr(unname(tools::md5sum(stan)), 1, 12)))
  if (!file.exists(exe)) {
    build <- tempfile("m1_build_")
    dir.create(build)
    on.exit(unlink(build, recursive = TRUE))
    built <- cmdstanr::cmdstan_model(stan, dir = build)
    staged <- tempfile("m1_stage_", dir)
    stopifnot(file.copy(built$exe_file(), staged), file.rename(staged, exe))
  }
  cmdstanr::cmdstan_model(stan, exe_file = exe, compile = FALSE)
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
       overdispersed = d$overdispersed, tau_u_sd = pr$tau_u_sd, lkj_eta = pr$lkj_eta,
       inv_phi_sd = pr$inv_phi_sd)
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
#' `model = "m1dm"` (default) adds a sample-level residual: counts are
#' Dirichlet-multinomial, `y_s ~ DM(phi * mu_s)` with `mu_s = softmax(eta_s)` and
#' `inv_phi = 1 / phi ~ half-normal(0, inv_phi_sd)`. As `inv_phi -> 0` this is
#' `model = "m1"`, the multinomial model with no sample residual; keep it for
#' ablation.
#'
#' Why `"m1dm"` is the default: in the 2026-11-02 benchmark (545 replicates, 18
#' cells) the pre-registered rule found no calibrated model, so the choice is post
#' hoc. `"m1"` under-covers (CPS 0.24 to 0.79 in cells it failed), while `"m1dm"`
#' errs the other way: its CPS intervals are conservative. Both keep RMS coverage
#' inside the band in every cell. See `docs/adr/0002-default-model-m1dm.md`. CCS is
#' experimental under either model, because no model covered it under donor
#' heterogeneity.
#'
#' Two kinds of probability are available. `level = "sample"` (default) returns the
#' realised sample probabilities `p_s | y_s ~ Dirichlet(phi * mu_s + y_s)` in
#' `fit$draws`. `level = "structured"` returns `mu_s` there. Both are always kept:
#' `fit$draws_sample` and `fit$draws_structured`, and `fit$draws` points at the one
#' you chose. For `model = "m1"` they are identical.
#'
#' @param x A `cj_data` object.
#' @param model `"m1dm"` (Dirichlet-multinomial sample residual) or `"m1"` (exact
#'   multinomial; draws are identical to the version before `"m1dm"` existed).
#' @param level Which probabilities go in `fit$draws`: `"sample"` or `"structured"`.
#' @param chains,iter_warmup,iter_sampling,seed Passed to CmdStan.
#' @param priors Named list that overrides [cj_default_priors()]. The defaults are
#'   illustrative smoke-test values, not validated, except `inv_phi_sd`, which comes
#'   from a prior-predictive check (see `?cj_default_priors`).
#' @param donor_visit Include the donor-by-visit deviation `u_iv` (set `FALSE` for
#'   the ablation). Ignored, with a message, when the data hold one visit.
#' @param donor_effects Include donor effects `(b_i, r_im)`; `FALSE` gives a
#'   fixed-effects-only model.
#' @param centered Sample donor effects `(b, r)` and `u` on their natural scale
#'   instead of the noncentered default. Both give the same posterior. With many
#'   cells per sample (data-dominated donors) the noncentered form needs very small
#'   steps; in a J = 3, 24-donor, 10 000-cell benchmark the centered form with a dense
#'   metric (set automatically) took about a quarter of the time. But the centered
#'   form funnels when any donor variance is near zero: in the 2026-11-02 pilot,
#'   null and dependence-only truths (no donor variation in the induced effect) gave
#'   R-hat up to 1.4 and bulk ESS below 30, and fitting `donor_visit = TRUE` to data
#'   without donor-by-visit variation gave R-hat 1.23. `FALSE` (default) is the safe
#'   choice; use `TRUE` only as a speed option when donor variation is clearly present.
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
cj_fit <- function(x, model = c("m1dm", "m1"), chains = 4L, iter_warmup = 1000L,
                   iter_sampling = 1000L, seed = 1L, priors = list(), donor_visit = TRUE,
                   donor_effects = TRUE, centered = FALSE, level = c("sample", "structured"), ...) {
  stopifnot(inherits(x, "cj_data"))
  model <- match.arg(model)
  level <- match.arg(level)
  d <- cj_design(x, donor_effects, donor_visit, centered, overdispersed = model == "m1dm")
  pr <- cj_resolve_priors(priors, d)
  # cmdstanr names CSVs from the global RNG and the minute, so two fits that start from the
  # same seed collide, and the older fit's finalizer then deletes the newer fit's files.
  # A private output_dir is not auto-deleted, so the files live as long as the session.
  out <- tempfile("m1_out_")
  dir.create(out)
  args <- utils::modifyList(list(refresh = 0, parallel_chains = min(chains, 4L), sig_figs = 18,
                           metric = if (centered) "dense_e" else "diag_e", output_dir = out), list(...))
  engine <- do.call(m1_compile()$sample,
                    c(list(data = m1_stan_data(d, x$counts, pr), chains = chains,
                           iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                           seed = seed), args))

  # Stan variable -> [draw, sample, pattern] array
  as_array <- function(var) {
    pd <- posterior::as_draws_matrix(engine$draws(var))
    idx <- do.call(rbind, lapply(regmatches(colnames(pd), gregexpr("[0-9]+", colnames(pd))), as.integer))
    n <- nrow(pd)
    p <- array(NA_real_, c(n, nrow(x$counts), ncol(x$counts)),
               dimnames = list(NULL, rownames(x$counts), colnames(x$counts)))
    p[cbind(rep(seq_len(n), ncol(pd)), rep(idx[, 1], each = n), rep(idx[, 2], each = n))] <- as.vector(pd)
    p <- pmax(p, 1e-300)
    p / as.vector(apply(p, 1:2, sum))   # CSV round-trip leaves ~1e-16 drift
  }
  as_draws <- function(p, kind)
    new_cj_draws(p, x$samples, x$contrasts, x$bits,
                 list(engine_id = "cmdstan-hmc", model_id = if (d$overdispersed) "M1-DM" else "M1",
                      seed = as.integer(seed), n_draws = dim(p)[1], draw_kind = "in_sample_donor_conditional",
                      level = kind, donor_effects = d$use_donor == 1L, donor_visit = d$use_dv == 1L,
                      centered = d$centered == 1L))
  p_sample <- as_array("p")
  draws_sample <- as_draws(p_sample, "sample")
  # m1: mu == p exactly, so reuse the draws.
  draws_structured <- if (d$overdispersed) as_draws(as_array("mu"), "structured") else draws_sample
  fit <- new_cj_fit(x, if (level == "sample") draws_sample else draws_structured, engine, d, pr,
                    m1_diagnostics(engine, d))
  fit$draws_sample <- draws_sample
  fit$draws_structured <- draws_structured
  fit
}

m1_diagnostics <- function(engine, d) {
  ds <- engine$diagnostic_summary(quiet = TRUE)
  vars <- intersect(c("alpha", "gamma", "theta0", "delta", "tau", "rho", "tau_u", "inv_phi", "phi"),
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
  cat("<cj_fit>", if (isTRUE(x$design$overdispersed == 1L)) "M1-DM" else "M1", "count model,", nrow(x$data$samples), "samples,",
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
