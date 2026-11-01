# Rung 1: J = 2, one stimulation, one visit, fixed effects only.

test_that("the R log-likelihood equals the Stan log-likelihood at a fixed parameter vector", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 6, N = 5000, donor_effects = FALSE, seed = 11)
  d <- sim$design
  pr <- cytoJointR:::cj_resolve_priors(list(), d)
  set.seed(1)
  par <- cytoJointR:::m1_prior_draw(d, pr)
  p <- cytoJointR:::m1_softmax(cytoJointR:::m1_eta(d, par))
  ll_r <- vapply(seq_len(nrow(p)), \(s) cj_log_multinomial(sim$data$counts[s, ], p[s, ]), 0)

  v <- cytoJointR:::m1_stan_vector(par, d)
  gq <- cytoJointR:::m1_compile()$generate_quantities(
    fitted_params = posterior::as_draws_matrix(matrix(v, 1, dimnames = list(NULL, names(v)))),
    data = cytoJointR:::m1_stan_data(d, sim$data$counts, pr), seed = 1, sig_figs = 18)
  expect_equal(as.numeric(gq$draws("log_lik", format = "matrix")), ll_r, tolerance = 1e-9)
})

test_that("rung 1 recovers fixed effects", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 10, N = 5000, donor_effects = FALSE, seed = 11)
  fit <- cj_fit(sim$data, chains = 4, seed = 1, donor_effects = FALSE, donor_visit = FALSE)
  rec <- recovery(fit, sim)
  print(rec, digits = 3)
  expect_s3_class(fit, "cj_fit")
  expect_equal(fit$draws$meta$draw_kind, "in_sample_donor_conditional")
  expect_identical(dim(fit$draws$p)[2:3], dim(sim$data$counts))
  expect_diagnostics_clean(fit)
  expect_true(all(rec$in99))
  expect_lt(fit$diagnostics$time_s, 60)
})
