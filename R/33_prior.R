# Priors of the M1 model and prior-predictive sampling.

#' Default priors of the M1 count model
#'
#' The values come from the bundle's `priors.example.json`. That file is marked
#' `ILLUSTRATIVE_SMOKE_INPUT_NOT_VALIDATED_DEFAULTS`: these are smoke-test inputs,
#' not validated defaults. Check them with [cj_prior_predict()] on your assay
#' before you trust a fit. `visit_sd`, `visit_stim_sd` and `tau_u_sd` are not in the
#' file; they are placeholders of the same kind.
#'
#' @return A named list: `baseline_any_positive` (sets the prior mean of `alpha`
#'   through the independence rule `alpha0 = logit(1 - (1 - q0)^(1/J))`, an
#'   initialisation aid only), `alpha_mean` (`NULL` or a length-J vector that
#'   overrides it), `alpha_sd`, `beta_sd` (stimulation main effects), `visit_sd`,
#'   `visit_stim_sd`, `theta0_sd`, `delta_sd`, `donor_sd` (half-normal scale of
#'   `tau`), `tau_u_sd` (half-normal scale of `tau_u`), `lkj_eta`.
#' @export
cj_default_priors <- function() {
  structure(list(baseline_any_positive = 0.005, alpha_mean = NULL, alpha_sd = 1,
                 beta_sd = 1.5, visit_sd = 1, visit_stim_sd = 1,
                 theta0_sd = 0.5, delta_sd = 0.3, donor_sd = 0.5, tau_u_sd = 0.5,
                 lkj_eta = 2),
            status = "ILLUSTRATIVE_SMOKE_INPUT_NOT_VALIDATED_DEFAULTS")
}

# Merge user priors over the defaults and expand them to what Stan receives.
cj_resolve_priors <- function(priors, d) {
  def <- cj_default_priors()
  bad <- setdiff(names(priors), names(def))
  if (length(bad)) stop("unknown priors: ", paste(bad, collapse = ", "), call. = FALSE)
  pr <- utils::modifyList(def, priors, keep.null = TRUE)
  J <- length(d$markers)
  if (is.null(pr$alpha_mean))
    pr$alpha_mean <- rep(stats::qlogis(1 - (1 - pr$baseline_any_positive)^(1 / J)), J)
  stopifnot(length(pr$alpha_mean) == J)
  pr$gamma_sd <- unname(unlist(pr[c("beta_sd", "visit_sd", "visit_stim_sd")])[d$col_info$prior])
  pr
}

# LKJ(eta) correlation matrix by the onion method (Lewandowski et al. 2009).
rlkj <- function(d, eta) {
  if (d == 1L) return(matrix(1))
  b <- eta + (d - 2) / 2
  r <- 2 * stats::rbeta(1, b, b) - 1
  S <- matrix(c(1, r, r, 1), 2)
  for (m in seq_len(d - 2L) + 1L) {
    b <- b - 0.5
    w <- stats::rnorm(m)
    w <- sqrt(stats::rbeta(1, m / 2, b)) * w / sqrt(sum(w^2))
    q <- t(chol(S)) %*% w
    S <- rbind(cbind(S, q), cbind(t(q), 1))
  }
  S
}

# One draw of the Stan parameters from the prior (noncentered scale).
m1_prior_draw <- function(d, pr) {
  J <- length(d$markers); M <- length(d$stimulations); P <- ncol(d$Bpair); Pf <- ncol(d$X)
  par <- list(
    alpha = stats::rnorm(J, pr$alpha_mean, pr$alpha_sd),
    gamma = sapply(seq_len(J), \(j) stats::rnorm(Pf, 0, pr$gamma_sd)) |> matrix(Pf, J),
    theta0 = stats::rnorm(P, 0, pr$theta0_sd),
    delta = matrix(stats::rnorm(M * P, 0, pr$delta_sd), M, P))
  if (d$use_donor) {
    Dd <- ncol(d$W)
    par$tau <- abs(stats::rnorm(Dd, 0, pr$donor_sd))
    par$L_Omega <- t(chol(rlkj(Dd, pr$lkj_eta)))
    par$z <- matrix(stats::rnorm(Dd * length(d$donors)), Dd)
  }
  if (d$use_dv) {
    par$tau_u <- abs(stats::rnorm(1, 0, pr$tau_u_sd))
    par$u_raw <- stats::rnorm(nrow(d$dv))
  }
  m1_zero_donor(par, d)
}

#' Prior-predictive pattern probabilities
#'
#' Samples the M1 priors and returns the implied pattern probabilities for every
#' sample of `x`, so you can inspect baseline prevalence, the whole pattern table
#' and high-degree patterns before fitting. Counts in `x` are not used.
#' The returned `cj_draws` has `meta$engine_id = "prior"`. Donor effects are drawn
#' for the donors in `x`, hence `draw_kind = "in_sample_donor_conditional"`.
#'
#' @param x A `cj_data` object that supplies the design.
#' @param n_draws Number of prior draws.
#' @param priors Named list that overrides [cj_default_priors()].
#' @inheritParams cj_design
#' @param seed Integer seed.
#' @return A `cj_draws` object.
#' @examples
#' sim <- cj_simulate_counts(n_donors = 4, seed = 1)
#' pp <- cj_prior_predict(sim$data, n_draws = 50, seed = 1)
#' # prior distribution of the control sample's any-positive probability
#' quantile(1 - pp$p[, 1, 1], c(0.05, 0.5, 0.95))
#' @export
cj_prior_predict <- function(x, n_draws = 200L, priors = list(), donor_effects = TRUE,
                             donor_visit = TRUE, seed = 1L) {
  d <- cj_design(x, donor_effects, donor_visit)
  pr <- cj_resolve_priors(priors, d)
  set.seed(seed)
  draws <- lapply(seq_len(n_draws), \(i) m1_softmax(m1_eta(d, m1_prior_draw(d, pr))))
  p <- aperm(simplify2array(draws), c(3L, 1L, 2L))
  dimnames(p) <- list(NULL, x$samples$sample_id, pattern_labels(x$bits))
  new_cj_draws(p, x$samples, x$contrasts, x$bits,
               list(engine_id = "prior", model_id = "M1", seed = as.integer(seed),
                    n_draws = as.integer(n_draws), draw_kind = "in_sample_donor_conditional",
                    priors = pr))
}
