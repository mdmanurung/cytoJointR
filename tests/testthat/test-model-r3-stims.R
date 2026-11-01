# Rung 3: several stimulations sharing a control sample; missing controls.

test_that("rung 3 recovers two stimulations sharing one control per donor", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 20, stimulations = c("CMV", "SEB"), N = 5000,
                            donor_visit = FALSE, seed = 31)
  expect_equal(sum(sim$data$samples$condition == "control"), 20)   # one control, two contrasts
  expect_equal(length(unique(sim$data$contrasts$control_sample_id)), 20)
  fit <- cj_fit(sim$data, chains = 4, seed = 1, donor_visit = FALSE)
  rec <- recovery(fit, sim)
  print(rec, digits = 3)
  print(fit)
  expect_diagnostics_clean(fit)
  expect_true(all(rec$in99))
  expect_lt(fit$diagnostics$time_s, 60)
})

test_that("stimulated samples without a control still enter the fit", {
  skip_if_no_stan()
  sim <- cj_simulate_counts(n_donors = 20, stimulations = c("CMV", "SEB"), N = 5000,
                            donor_visit = FALSE, missing_control = 0.3, seed = 32)
  expect_equal(sum(is.na(sim$data$contrasts$control_sample_id)), 2 * 6)
  fit <- cj_fit(sim$data, chains = 4, seed = 1, donor_visit = FALSE)
  rec <- recovery(fit, sim)
  print(rec, digits = 3)
  expect_equal(dim(fit$draws$p)[2], nrow(sim$data$samples))       # p for every sample
  expect_equal(sum(fit$design$X[, "stim_CMV"]), 20)               # all stimulated rows kept
  expect_diagnostics_clean(fit)
  expect_true(all(rec$in99))
})
