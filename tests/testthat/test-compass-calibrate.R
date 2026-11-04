# Calibration of non-paired COMPASS calls (R/61_compass_calibrate.R). The map functions need no COMPASS.

test_that("the isotonic map is monotone and moves a liberal cut to the paired call", {
  set.seed(1)
  pop <- runif(200)
  paired <- as.numeric(pop > 0.7)          # paired call agrees only above 0.7
  m <- cytoJointR:::compass_cal_map(c(0.1, 0.6, 0.65, 0.75, 0.95), pop, paired, "isotonic")
  expect_false(is.unsorted(m))
  expect_true(all(m[1:3] < 0.5))
  expect_true(all(m[4:5] > 0.5))
})

test_that("the quantile map matches the paired quantile at the same rank", {
  pop <- seq(0, 1, length.out = 101)
  paired <- pop^3
  m <- cytoJointR:::compass_cal_map(c(0.25, 0.5, 1), pop, paired, "quantile")
  expect_equal(m, unname(quantile(paired, c(0.25, 0.5, 1), type = 1)) |> unname(), tolerance = 0.02)
  expect_false(is.unsorted(m))
})

test_that("the fallback map sends the fallback cut to 0.5 and keeps 0 and 1", {
  cut <- cytoJointR:::COMPASS_CAL_FALLBACK_CUT
  expect_equal(cytoJointR:::compass_cal_fallback(c(0, cut, 1)), c(0, 0.5, 1))
  expect_lt(cytoJointR:::compass_cal_fallback(0.8), 0.5)
})

test_that("compass_cal_auc ranks and is NA with one class", {
  expect_equal(cytoJointR:::compass_cal_auc(c(.1, .2, .8, .9), c(FALSE, FALSE, TRUE, TRUE)), 1)
  expect_equal(cytoJointR:::compass_cal_auc(c(.9, .8, .2, .1), c(FALSE, FALSE, TRUE, TRUE)), 0)
  expect_true(is.na(cytoJointR:::compass_cal_auc(c(.1, .2), c(TRUE, TRUE))))
})

test_that("compass_with_seed is reproducible and restores the caller's RNG", {
  set.seed(5); a <- runif(1); set.seed(5)
  x <- cytoJointR:::compass_with_seed(9, runif(3))
  expect_equal(x, cytoJointR:::compass_with_seed(9, runif(3)))
  expect_equal(runif(1), a)
})

cal_scores <- function(raw, src) {
  data.frame(compass_max_gamma = raw, control_source = src, stringsAsFactors = FALSE)
}

test_that("compass_calibrate_scores maps non-paired rows, keeps paired rows, and falls back", {
  set.seed(2)
  pop <- runif(40)
  ps <- data.frame(group = "g", cls = "population", gamma_pseudo = pop,
                   gamma_paired = as.numeric(pop > 0.7))
  sc <- cal_scores(c(0.9, 0.9, 0.6), c("paired", "population", "population"))
  r <- cytoJointR:::compass_calibrate_scores(sc, rep("g", 3), ps, "isotonic")
  expect_equal(r$scores$compass_gamma_cal[1], 0.9)            # paired row untouched
  expect_equal(r$scores$calibration, c("paired", "isotonic", "isotonic"))
  expect_equal(r$scores$compass_call, c(TRUE, TRUE, FALSE))
  expect_equal(r$info$n_pseudo, 40)
  # Calibration only tightens: a raw 0.4 stays uncalled even if the map would lift it.
  lift <- cytoJointR:::compass_calibrate_scores(cal_scores(0.4, "population"), "g",
            data.frame(group = "g", cls = "population", gamma_pseudo = c(rep(0.2, 20), rep(0.4, 20)),
                       gamma_paired = c(rep(0, 20), rep(1, 20))), "isotonic")
  expect_lte(lift$scores$compass_gamma_cal, 0.4)
  expect_false(lift$scores$compass_call)
  # Too few pseudo donors, or one side empty: the fallback cut applies.
  few <- cytoJointR:::compass_calibrate_scores(sc, rep("g", 3), ps[1:6, ], "isotonic")
  expect_equal(few$scores$calibration[2:3], c("fallback", "fallback"))
  one_side <- transform(ps, gamma_paired = 0.1)
  none <- cytoJointR:::compass_calibrate_scores(sc, rep("g", 3), one_side, "isotonic")
  expect_equal(none$scores$calibration[2], "fallback")
  expect_true(none$scores$compass_call[2])   # raw 0.9 is above the fallback cut
})

test_that("a pseudo-missing control is forced through the same resolver", {
  d <- cj_example_draws(n_donor = 4, n_draw = 5, markers = c("A", "B"), visits = c("v1", "v2"),
                        stimulations = "CMV")
  x <- cytoJointR:::new_cj_data(
    {m <- t(vapply(seq_len(nrow(d$samples)), function(i)
      rmultinom(1, d$samples$N[i], colMeans(d$p[, i, ]))[, 1], numeric(4)))
     storage.mode(m) <- "integer"; dimnames(m) <- list(d$samples$sample_id, dimnames(d$p)[[3]]); m},
    d$samples, d$contrasts, d$bits)
  id <- x$contrasts$contrast_id[x$contrasts$donor_id == "D1" & x$contrasts$visit == "v2"]
  r <- cytoJointR:::compass_resolve_controls(x, "donor", force_missing = id)$contrasts
  expect_equal(r$control_source[r$contrast_id == id], "donor_visit:v1")
  expect_equal(sum(r$control_source == "paired"), nrow(r) - 1L)
  p <- cytoJointR:::compass_resolve_controls(x, "population", force_missing = id)$contrasts
  expect_true(is.na(p$control_sample_id[p$contrast_id == id]))
  expect_equal(p$control_source[p$contrast_id == id], "population")
})

# 14 donors, one visit; the controls of D1 to D4 are missing.
cal_data <- function(n_donor = 14, n_miss = 4) {
  miss <- sprintf("D%d_v1_CD4_CMV", seq_len(n_miss))
  d <- cj_example_draws(n_donor = n_donor, n_draw = 5, markers = c("IFNG", "TNF", "IL2"),
                        visits = "v1", stimulations = "CMV", missing_control = miss)
  x <- example_data(d)
  keep <- !x$samples$sample_id %in% sub("CD4_CMV", "control", miss)
  cytoJointR:::new_cj_data(x$counts[keep, ], x$samples[keep, ], x$contrasts, x$bits)
}

test_that("missing_control_calibrate adds columns, leaves the default alone, and pseudo-fits paired donors", {
  skip_if_not_installed("COMPASS")
  x <- cal_data()
  base <- suppressMessages(cj_compass(x, missing_control = "population", iterations = 200,
                                      replications = 2, seed = 3))
  expect_true("compass_max_gamma" %in% names(base$scores))
  expect_false(any(c("compass_call", "compass_gamma_cal") %in% names(base$scores)))
  expect_null(base$calibration)
  cc <- suppressWarnings(suppressMessages(cj_compass(
    x, missing_control = "population", missing_control_calibrate = TRUE, max_pseudo = 5,
    iterations = 200, replications = 2, seed = 3)))
  # Raw gamma is unchanged by calibration; the paired rows keep it.
  expect_equal(cc$scores$compass_max_gamma, base$scores$compass_max_gamma)
  pr <- cc$scores$control_source == "paired"
  expect_equal(cc$scores$compass_gamma_cal[pr], cc$scores$compass_max_gamma[pr])
  expect_equal(cc$scores$compass_call, cc$scores$compass_gamma_cal > 0.5)
  # Only paired donors are pseudo-removed, at most max_pseudo, and their paired gamma is the real fit's.
  ps <- cc$calibration$pseudo
  expect_lte(nrow(ps), 5)
  expect_true(all(ps$donor_id %in% cc$scores$donor_id[pr]))
  real <- cc$scores$compass_max_gamma[match(ps$donor_id, cc$scores$donor_id)]
  expect_equal(ps$gamma_paired, real)
  # 5 paired donors are below the minimum, so the four missing rows use the fallback cut.
  expect_equal(unname(cc$scores$calibration[!pr]), rep("fallback", 4))
  expect_equal(cc$calibration$table$method, "fallback")
})

test_that("missing_control_calibrate with 'drop' or nothing missing is a no-op", {
  skip_if_not_installed("COMPASS")
  x <- cal_data(n_miss = 0)
  cc <- suppressMessages(cj_compass(x, missing_control = "population", missing_control_calibrate = TRUE,
                                    iterations = 100, replications = 2, seed = 3))
  expect_true(all(cc$scores$calibration == "paired"))
  expect_equal(cc$scores$compass_gamma_cal, cc$scores$compass_max_gamma)
})
