# Rung 1: J = 2, one stimulation, one visit, fixed effects only.

test_that("the R log-likelihood equals the Stan log-likelihood at a fixed parameter vector", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 6, N = 5000, donor_effects = FALSE, seed = 11)
  ll <- loglik_r_vs_stan(sim)
  expect_equal(ll$stan, ll$r, tolerance = 1e-9)
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
