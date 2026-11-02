test_that("C0 returns valid cj_draws whose mean is (a0 + y) / (K * a0 + N)", {
  sim <- cj_simulate_counts(n_donors = 4, N = 1000, seed = 1)
  x <- sim$data
  dr <- cytoJointR:::cj_fit_dirichlet(x, a0 = 0.5, n_draws = 4000, seed = 1)
  expect_s3_class(dr, "cj_draws")
  expect_identical(dr$meta$engine_id, "conjugate")
  expect_identical(dr$meta$model_id, "C0")
  expect_identical(dr$meta$draw_kind, "in_sample_donor_conditional")
  expect_identical(dim(dr$p), c(4000L, dim(x$counts)))
  K <- ncol(x$counts)
  expect_equal(apply(dr$p, 2:3, mean), (x$counts + 0.5) / (K * 0.5 + 1000), tolerance = 0.02,
               ignore_attr = TRUE)
  expect_identical(cytoJointR:::cj_fit_dirichlet(x, n_draws = 10, seed = 3),
                   cytoJointR:::cj_fit_dirichlet(x, n_draws = 10, seed = 3))
  eps <- c(RMS = 0.005, poly_absolute = 0.005, poly_marginal = 0.005, CPS = 0.02,
           CCS = 0.005, breadth = 0.05, marker_magnitude = 0.005)
  sc <- cj_scores(dr, eps)   # the score layer runs on it unchanged
  expect_gt(nrow(sc), 0)
})
