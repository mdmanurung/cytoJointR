test_that("rlkj has the LKJ marginal: r ~ 2 Beta(b, b) - 1, b = eta - 1 + d / 2", {
  set.seed(1)
  r <- replicate(4000, cytoJointR:::rlkj(3L, 2)[1, 2])
  b <- 2 - 1 + 3 / 2
  expect_equal(var(r), 1 / (2 * b + 1), tolerance = 0.07)
  expect_true(all(replicate(50, min(eigen(cytoJointR:::rlkj(4L, 2))$values)) > 0))
})

test_that("cj_prior_predict returns valid cj_draws labelled as prior draws", {
  sim <- cj_simulate_counts(n_donors = 4, markers = c("A", "B", "C"),
                            visits = c("v1", "v2"), donor_visit = TRUE, seed = 1)
  pp <- cj_prior_predict(sim$data, n_draws = 30, seed = 3)
  expect_s3_class(pp, "cj_draws")
  expect_equal(pp$meta$engine_id, "prior")
  expect_equal(dim(pp$p), c(30, nrow(sim$data$samples), 8))
  pp2 <- cj_prior_predict(sim$data, n_draws = 30, seed = 3)
  expect_identical(pp$p, pp2$p)
})

test_that("priors are validated and defaults are marked illustrative", {
  expect_match(attr(cj_default_priors(), "status"), "NOT_VALIDATED")
  sim <- cj_simulate_counts(n_donors = 3, seed = 1)
  expect_error(cj_prior_predict(sim$data, priors = list(nope = 1)), "unknown priors")
})
