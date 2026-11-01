# Shared helpers for the M1 model tests.

skip_if_no_stan <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if(is.null(tryCatch(cmdstanr::cmdstan_version(error_on_NA = FALSE),
                                     error = function(e) NULL)), "CmdStan not installed")
}

# True parameters under their Stan names, including rho from the truth's correlation.
truth_key <- function(sim) {
  tr <- sim$truth; d <- sim$design
  v <- cytoJointR:::m1_stan_vector(tr[c("alpha", "gamma", "theta0", "delta", "tau", "tau_u")], d)
  if (d$use_donor) {
    om <- tr$L_Omega %*% t(tr$L_Omega)
    ab <- utils::combn(nrow(om), 2L)
    v <- c(v, setNames(om[cbind(ab[1, ], ab[2, ])], sprintf("rho[%d]", seq_len(ncol(ab)))))
  }
  v
}

# Truth against the posterior mean and the 90% and 99% intervals.
recovery <- function(fit, sim) {
  tv <- truth_key(sim)
  dr <- fit$engine_fit$draws(names(tv), format = "draws_matrix")
  q <- apply(dr, 2L, stats::quantile, c(0.005, 0.05, 0.95, 0.995))
  data.frame(param = fit$design$param_map$label[match(names(tv), fit$design$param_map$stan)],
             truth = unname(tv), mean = colMeans(dr), q5 = q[2, ], q95 = q[3, ],
             in90 = tv >= q[2, ] & tv <= q[3, ], in99 = tv >= q[1, ] & tv <= q[4, ],
             row.names = NULL)
}

expect_diagnostics_clean <- function(fit) {
  dg <- fit$diagnostics
  testthat::expect_equal(dg$n_divergent, 0)
  testthat::expect_lt(dg$max_rhat, 1.02)
  testthat::expect_gt(dg$min_ess_bulk, 200)
  testthat::expect_gt(min(dg$ebfmi), 0.3)
}
