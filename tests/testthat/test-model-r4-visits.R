# Rung 4: categorical visits, visit-by-stimulation effects, donor-by-visit deviation u_iv, J = 3.

r4_sim <- function() {
  cj_simulate_counts(n_donors = 12, markers = c("IFNG", "TNF", "IL2"), stimulations = "CMV",
                     visits = c("v1", "v2"), N = 1500, donor_visit = TRUE, seed = 41)
}

test_that("the R and Stan log-likelihoods agree with visits and u_iv", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 4, markers = c("IFNG", "TNF", "IL2"),
                            stimulations = c("CMV", "SEB"), visits = c("v1", "v2"),
                            N = 3000, donor_visit = TRUE, seed = 42)
  ll <- loglik_r_vs_stan(sim)
  expect_equal(ll$stan, ll$r, tolerance = 1e-9)
})

test_that("rung 4 recovers visit effects, interactions, donor scales and tau_u", {
  skip_if_no_stan()
  sim <- r4_sim()
  fit <- cj_fit(sim$data, chains = 4, iter_warmup = 500, iter_sampling = 500, seed = 1)
  rec <- recovery(fit, sim)
  print(rec, digits = 3)
  print(fit)
  expect_diagnostics_clean(fit)
  expect_true(all(rec$in99))
  expect_true("tau_u" %in% rec$param)
  expect_equal(dim(fit$draws$p)[2:3], dim(sim$data$counts))
})

test_that("donor_visit = FALSE is the ablation without u_iv", {
  skip_if_no_stan()
  sim <- r4_sim()
  fit <- cj_fit(sim$data, chains = 4, iter_warmup = 300, iter_sampling = 300, seed = 1,
                donor_visit = FALSE)
  expect_false("tau_u" %in% fit$design$param_map$label)
  expect_false(fit$draws$meta$donor_visit)
  expect_true("log_lik" %in% fit$engine_fit$metadata()$stan_variables)
  print(fit)
})
