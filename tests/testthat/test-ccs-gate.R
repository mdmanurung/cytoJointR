# CCS support gate: a cj_fit supplies the counts; a stub with the same fields avoids Stan.
stub_fit <- function(n_cells = 1e5, polyfunctional_cells = NULL) {
  d <- cj_example_draws(n_draw = 30)
  d$samples$N <- n_cells
  dat <- example_data(d)
  if (!is.null(polyfunctional_cells)) {
    dat$counts[, rowSums(d$bits) >= 2] <- polyfunctional_cells
  }
  structure(list(draws = d, data = dat), class = "cj_fit")
}
ccs_rows <- function(x) x[x$score == "CCS", ]

test_that("a sparse fit is gated with a reason and keeps its rows", {
  fit <- stub_fit(polyfunctional_cells = 0L)
  sc <- suppressMessages(cj_scores(fit, EPS))
  ref <- suppressMessages(cj_scores(fit, EPS, ccs_gate = NULL))
  expect_equal(nrow(sc), nrow(ref))
  ccs <- ccs_rows(sc)
  expect_true(all(ccs$status == "not_estimable"))
  expect_true(all(is.na(ccs$estimate) & is.na(ccs$lower) & is.na(ccs$p_positive)))
  expect_match(ccs$reason, "control polyfunctional cells 0 < 20", fixed = TRUE)
  expect_match(ccs$reason, "stimulated polyfunctional cells 0 < 20", fixed = TRUE)
  expect_true(all(ccs$experimental))
  # Other scores are untouched.
  expect_equal(sc[sc$score != "CCS", setdiff(names(sc), "reason")],
               ref[ref$score != "CCS", setdiff(names(ref), "reason")])
  expect_true(all(is.na(sc$reason[sc$score != "CCS"])))
})

test_that("a dense fit passes the gate and CCS equals the ungated value", {
  fit <- stub_fit()
  sc <- suppressMessages(cj_scores(fit, EPS))
  ref <- suppressMessages(cj_scores(fit, EPS, ccs_gate = NULL))
  expect_false(any(ccs_rows(sc)$status == "not_estimable"))
  expect_false(anyNA(ccs_rows(sc)$estimate))
  expect_true(all(is.na(sc$reason)))
  expect_equal(sc$estimate, ref$estimate)
})

test_that("ccs_gate = NULL switches the gate off", {
  fit <- stub_fit(polyfunctional_cells = 0L)
  sc <- suppressMessages(cj_scores(fit, EPS, ccs_gate = NULL))
  expect_false(anyNA(ccs_rows(sc)$estimate))
  expect_true(all(is.na(sc$reason)))
  expect_false(anyNA(ccs_rows(cj_score_draws(fit, ccs_gate = NULL))$value))
})

test_that("the CV condition withholds CCS when the control mass is not pinned down", {
  fit <- stub_fit()
  gate <- list(min_poly_cells = 20, max_cv = 1e-4)
  ccs <- ccs_rows(suppressMessages(cj_scores(fit, EPS, ccs_gate = gate)))
  expect_true(all(ccs$status == "not_estimable"))
  expect_match(ccs$reason, "mass CV", fixed = TRUE)
  expect_no_match(ccs$reason, "polyfunctional cells", fixed = TRUE)
})

test_that("cj_score_draws and cj_patterns apply the same gate", {
  fit <- stub_fit(polyfunctional_cells = 0L)
  dr <- cj_score_draws(fit)
  expect_equal(nrow(dr), nrow(cj_score_draws(fit, ccs_gate = NULL)))
  expect_true(all(ccs_rows(dr)$status == "not_estimable" & is.na(ccs_rows(dr)$value)))
  expect_true(all(dr$status[dr$score != "CCS"] == "ok"))
  pt <- cj_patterns(fit, c("RMS", "CCS"), c(RMS = 0.005, CCS = 0.005))
  expect_true(all(ccs_rows(pt)$status == "not_estimable" & is.na(ccs_rows(pt)$estimate)))
  expect_true(all(pt$status[pt$score == "RMS"] == "ok"))
  expect_match(ccs_rows(pt)$reason, "control polyfunctional cells 0 < 20", fixed = TRUE)
})

test_that("draws without counts report support_unknown", {
  d <- cj_example_draws(n_draw = 20)
  ccs <- ccs_rows(suppressMessages(cj_scores(d, EPS)))
  expect_true(all(ccs$status == "support_unknown"))
  expect_true(all(is.na(ccs$estimate)))
})

test_that("a contrast without a control stays unavailable", {
  fit <- stub_fit(polyfunctional_cells = 0L)
  fit$draws <- cj_example_draws(n_draw = 30, missing_control = "D1_v1_CD4_CMV")
  fit$draws$samples$N <- 1e5
  ccs <- ccs_rows(suppressMessages(cj_scores(fit, EPS)))
  expect_equal(ccs$status[ccs$contrast_id == "D1_v1_CD4_CMV"], "unavailable")
})

test_that("a malformed ccs_gate is an error", {
  d <- cj_example_draws(n_draw = 10)
  expect_error(cj_scores(d, EPS, ccs_gate = list(min_poly_cells = 20)), "ccs_gate")
  expect_error(cj_score_draws(d, ccs_gate = 20), "ccs_gate")
})
