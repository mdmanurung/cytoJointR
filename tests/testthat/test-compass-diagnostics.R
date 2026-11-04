# A COMPASSResult stand-in: 2 donors, 2 markers (subsets "11", "01", "10", null "00"),
# n iterations of 0/1 gamma drawn with the given per-subset probabilities.
fake_fit <- function(p, n = 400, seed = 1) {
  set.seed(seed)
  cats <- cbind(A = c(1, 0, 1, 0), B = c(1, 1, 0, 0), Counts = c(2, 1, 1, 0))
  g <- array(0L, c(2, 4, n), dimnames = list(c("d1", "d2"), NULL, NULL))
  for (i in 1:2) for (k in 1:4) g[i, k, ] <- rbinom(n, 1, p[i, k])
  mg <- apply(g, c(1, 2), mean)
  list(fit = list(gamma = g, mean_gamma = mg, categories = cats))
}

test_that("compass_rhat_ess handles constant and non-mixing runs", {
  expect_equal(unname(cytoJointR:::compass_rhat_ess(matrix(0, 50, 4))["rhat"]), 1)
  expect_true(is.na(cytoJointR:::compass_rhat_ess(matrix(0, 50, 4))["ess"]))
  # Each run stuck on its own value never mixes.
  stuck <- matrix(rep(c(0, 1, 0, 1), each = 50), 50, 4)
  expect_equal(unname(cytoJointR:::compass_rhat_ess(stuck)["rhat"]), Inf)
  set.seed(1)
  iid <- matrix(rnorm(2000), 500, 4)
  expect_lt(cytoJointR:::compass_rhat_ess(iid)[["rhat"]], 1.02)
  shifted <- iid + rep(c(0, 0, 1, 1), each = 500)
  expect_gt(cytoJointR:::compass_rhat_ess(shifted)[["rhat"]], 1.05)
})

test_that("compass_n_runs validates and raises diagnostics = TRUE to 4 runs", {
  expect_equal(cytoJointR:::compass_n_runs(1, FALSE), 1L)
  expect_equal(cytoJointR:::compass_n_runs(1, TRUE), 5L)
  expect_equal(cytoJointR:::compass_n_runs(7, TRUE), 7L)
  expect_error(cytoJointR:::compass_n_runs(0, FALSE), "n_runs")
  expect_warning(cytoJointR:::compass_n_runs(2, FALSE), "weak diagnostics")
})

test_that("unstable_call is raised by call disagreement or a high R-hat, not otherwise", {
  # d1: every subset responds with probability 0.9 in every run. d2: probability 0.5 in
  # runs 1-2 and 0.9 in runs 3-4 (the max mean gamma crosses 0.5 in some runs only).
  hi <- rbind(rep(0.9, 4), rep(0.9, 4))
  lo <- rbind(rep(0.9, 4), c(0.2, 0.2, 0.2, 0))
  fits <- list(fake_fit(lo, seed = 1), fake_fit(lo, seed = 2), fake_fit(hi, seed = 3), fake_fit(hi, seed = 4))
  d <- cytoJointR:::compass_run_diagnostics(fits)
  expect_equal(d$compass_runs_called, c(4, 2))
  expect_equal(d$compass_call_disagree, c(FALSE, TRUE))
  expect_equal(d$unstable_call, c(FALSE, TRUE))
  expect_gt(d$compass_gamma_rhat[2], 1.05)
  # R-hat flags only when asked to: rhat_max = Inf (default) leaves the call rule alone.
  expect_false(any(d$compass_rhat_high))
  strict <- cytoJointR:::compass_run_diagnostics(fits, rhat_max = 1.05)
  expect_equal(strict$compass_rhat_high[2], TRUE)
  expect_equal(strict$unstable_call, strict$compass_call_disagree | strict$compass_rhat_high)
  # d2: two runs call it, run 1 does not, so the reported call has support 2 of 4 (a tie).
  expect_equal(d$compass_call, c(TRUE, FALSE))
  expect_equal(d$compass_call_unsupported, c(FALSE, TRUE))
  # Agreeing runs: nobody flagged.
  same <- cytoJointR:::compass_run_diagnostics(lapply(1:4, function(s) fake_fit(hi, seed = s)))
  expect_false(any(same$unstable_call))
  expect_true(all(same$compass_gamma_rhat < 1.05))
})

test_that("draw-level FS and PFS average to COMPASS's FS and PFS", {
  skip_if_not_installed("COMPASS")
  f <- fake_fit(rbind(c(.7, .5, .3, 0), c(.2, .2, .9, 0)))
  dr <- cytoJointR:::compass_draws(f)
  expect_equal(colMeans(dr$FS), COMPASS::FunctionalityScore(structure(f, class = "COMPASSResult")))
  expect_equal(colMeans(dr$PFS), COMPASS::PolyfunctionalityScore(structure(f, class = "COMPASSResult")))
})

test_that("cj_compass(n_runs, diagnostics) adds convergence columns and keeps run 1", {
  skip_if_not_installed("COMPASS")
  d <- cj_example_draws(n_donor = 6, n_draw = 5, markers = c("IFNG", "TNF", "IL2"),
                        visits = "v1", stimulations = "CMV")
  x <- example_data(d)
  one <- cj_compass(x, iterations = 300, replications = 2, seed = 5)
  expect_null(one$diagnostics)
  expect_false("unstable_call" %in% names(one$scores))
  cc <- cj_compass(x, iterations = 300, replications = 2, seed = 5, n_runs = 4)
  # The first run is the single-run fit: same seed, same numbers.
  expect_equal(cc$scores$compass_FS, one$scores$compass_FS)
  expect_equal(cc$diagnostics$seeds, 5:8)
  need <- c("compass_gamma_rhat", "compass_gamma_ess", "compass_FS_rhat", "compass_PFS_rhat",
            "compass_runs_called", "compass_max_gamma", "unstable_call")
  expect_true(all(need %in% names(cc$scores)))
  expect_equal(nrow(cc$scores), nrow(one$scores))
  expect_equal(cc$diagnostics$summary$n_unstable, sum(cc$scores$unstable_call))
  expect_equal(cc$diagnostics$summary$n_donors, nrow(cc$scores))
  # Only the first fit is kept.
  expect_length(cc$fits, 1)
  # diagnostics = TRUE means 5 runs; the seeds are reproducible (runs 1-4 are those of n_runs = 4).
  dd <- cj_compass(x, iterations = 300, replications = 2, seed = 5, diagnostics = TRUE)
  expect_equal(dd$diagnostics$n_runs, 5L)
  expect_equal(dd$diagnostics$seeds, 5:9)
  expect_equal(dd$scores$compass_FS, cc$scores$compass_FS)
  # cj_run_compass passes the option on.
  sim <- cj_simulate_cells(donors = 4, cells_per_sample = 300, seed = 1)
  samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
  samples$visit <- factor("v1"); samples$parent_population <- "CD4"
  out <- cj_run_compass(sim$cells, c("IFNG", "TNF", "IL2"), samples, transform = "none",
                        control = "unstim", gate_args = list(method = "quantile"),
                        iterations = 200, replications = 2, n_runs = 3) |> suppressWarnings()
  expect_true("unstable_call" %in% names(out))
})
