# Rung 2: donor effects (b_i, r_i) on top of rung 1.

test_that("the R and Stan log-likelihoods agree with donor effects", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 5, N = 5000, donor_visit = FALSE, seed = 21)
  ll <- loglik_r_vs_stan(sim)
  expect_equal(ll$stan, ll$r, tolerance = 1e-9)
})

test_that("rung 2 recovers fixed effects, donor scales and correlation", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 24, N = 5000, donor_visit = FALSE, seed = 21)
  fit <- cj_fit(sim$data, chains = 4, seed = 1,
                donor_visit = FALSE)
  rec <- recovery(fit, sim)
  print(rec, digits = 3)
  print(fit)
  expect_diagnostics_clean(fit)
  expect_true(all(rec$in99))
  expect_lt(fit$diagnostics$time_s, 60)

  # shrunken donor effects track the realised ones
  b <- fit$engine_fit$draws("z", format = "draws_matrix")
  expect_gt(stats::cor(colMeans(b)[seq(1, ncol(b), by = 2)], sim$truth$z[1, ]), 0.9)
})

test_that("the centered parameterisation gives the same likelihood and recovery", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 5, N = 5000, visits = c("v1", "v2"), donor_visit = TRUE,
                            seed = 22)
  ll <- loglik_r_vs_stan(sim, centered = TRUE)
  expect_equal(ll$stan, ll$r, tolerance = 1e-9)

  sim <- cj_simulate_counts(n_donors = 24, N = 5000, donor_visit = FALSE, seed = 21)
  fit <- cj_fit(sim$data, chains = 4, seed = 1, donor_visit = FALSE, centered = TRUE)
  rec <- recovery(fit, sim)
  print(rec[7:9, ], digits = 3)
  print(fit)
  expect_diagnostics_clean(fit)
  expect_true(all(rec$in99))
  expect_lt(fit$diagnostics$time_s, 60)
})
