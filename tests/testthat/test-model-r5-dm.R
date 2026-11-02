# Rung 5: Dirichlet-multinomial sample residual (model = "m1dm").

dm_sim <- function(phi, seed) {
  cj_simulate_counts(n_donors = 12, N = 5000, donor_visit = FALSE, truth = list(phi = phi), seed = seed)
}

test_that("the R Dirichlet-multinomial log-likelihood equals the Stan one", {
  skip_if_no_stan()
  sim <- dm_sim(500, 21)
  ll <- loglik_r_vs_stan(sim)
  expect_equal(ll$stan, ll$r, tolerance = 1e-9)
})

test_that("the simulator keeps mu and p, and phi = Inf is plain M1", {
  a <- cj_simulate_counts(n_donors = 4, N = 2000, seed = 5)
  expect_identical(a$mu, a$p)
  b <- cj_simulate_counts(n_donors = 4, N = 2000, truth = list(phi = Inf), seed = 5)
  expect_identical(a[names(a) != "design"], b[names(b) != "design"])
  c <- cj_simulate_counts(n_donors = 4, N = 2000, truth = list(phi = 50), seed = 5)
  expect_identical(c$mu, a$mu)
  expect_false(isTRUE(all.equal(c$p, c$mu)))
  expect_equal(unname(rowSums(c$p)), rep(1, nrow(c$p)))
})

test_that("m1dm recovers phi = 500", {
  skip_if_no_stan()
  sim <- dm_sim(500, 21)
  fit <- cj_fit(sim$data, chains = 4, seed = 1, donor_visit = FALSE)
  rec <- recovery(fit, sim)
  print(rec[grepl("phi", rec$param), ], digits = 3)
  expect_true(rec$in99[rec$param == "inv_phi"])
  expect_diagnostics_clean(fit, max_divergent = 4)
  expect_recovered(rec)
  expect_lt(fit$diagnostics$time_s, 60)
  # sample level is the default; structured is kept beside it
  expect_identical(fit$draws, fit$draws_sample)
  expect_identical(fit$draws$meta$level, "sample")
  expect_identical(fit$draws_structured$meta$level, "structured")
  expect_false(isTRUE(all.equal(fit$draws$p, fit$draws_structured$p)))
  expect_identical(cj_fit(sim$data, chains = 1, iter_warmup = 100, iter_sampling = 100, seed = 1,
                          donor_visit = FALSE, level = "structured")$draws$meta$level, "structured")
})

test_that("with phi = Inf data, inv_phi concentrates near 0", {
  skip_if_no_stan()
  sim <- dm_sim(Inf, 22)
  fit <- cj_fit(sim$data, chains = 4, seed = 1, donor_visit = FALSE)
  inv_phi <- as.vector(fit$engine_fit$draws("inv_phi", format = "draws_matrix"))
  expect_lt(stats::quantile(inv_phi, 0.95), 0.005)   # phi > 200 at N = 5000
  expect_diagnostics_clean(fit, max_divergent = 4)
})

test_that("model = \"m1\" is the multinomial model, without sample residual", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 10, N = 5000, donor_visit = FALSE, seed = 11)
  fit <- cj_fit(sim$data, model = "m1", chains = 4, seed = 1, donor_visit = FALSE)
  expect_false(any(grepl("^(inv_)?phi", dimnames(fit$engine_fit$draws())$variable)))
  expect_identical(fit$design$overdispersed, 0L)
  expect_identical(fit$draws$p, fit$draws_structured$p)
  expect_identical(fit$draws$meta$model_id, "M1")
  rec <- recovery(fit, sim)
  expect_recovered(rec)
  expect_diagnostics_clean(fit, max_divergent = 4)
})

test_that("cj_prior_predict draws inv_phi for m1dm only", {
  sim <- cj_simulate_counts(n_donors = 4, seed = 1)
  a <- cj_prior_predict(sim$data, n_draws = 20, seed = 2)
  b <- cj_prior_predict(sim$data, n_draws = 20, model = "m1", seed = 2)
  expect_length(a$meta$inv_phi, 20)
  expect_true(all(a$meta$inv_phi >= 0))
  expect_null(b$meta$inv_phi)
  expect_equal(unname(apply(a$p, 1:2, sum)), matrix(1, 20, nrow(sim$data$counts)))
})
