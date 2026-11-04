test_that("cj_mimosa2 passes each donor's four samples to the model by name", {
  b <- mimosa2_data()
  mm <- cj_mimosa2(b$x)
  expect_s3_class(mm, "cj_mimosa2")
  expect_s3_class(mm$fits[[1]], "cj_mimosa2_fit")
  expect_identical(mm$groups$status, "fitted")
  s <- mm$scores[mm$scores$donor_id == "D07", ]
  ifng <- b$x$counts[, "IFNG+ TNF-"]
  expect_equal(c(s$ns1, s$nu1, s$ns0, s$nu0),
               unname(ifng[c("D07_post_CMV", "D07_post_control", "D07_pre_CMV", "D07_pre_control")]))
  expect_equal(c(s$Ns1, s$Nu1, s$Ns0, s$Nu0), rep(20000, 4))
  expect_equal(s$delta, (s$ns1 - s$nu1 - s$ns0 + s$nu0) / 20000)
})

test_that("cj_mimosa2 separates simulated responders from non-responders", {
  b <- mimosa2_data()
  s <- cj_mimosa2(b$x)$scores
  truth <- s$donor_id %in% b$responders
  # A clear effect (0.002 to 0.006 at 20,000 cells): every responder outranks every non-responder.
  expect_gt(min(s$mimosa2_p_response[truth]), max(s$mimosa2_p_response[!truth]))
  expect_false(any(s$mimosa2_response[!truth]))
})

test_that("positive = <marker> counts that marker's cells; bad arguments fail", {
  b <- mimosa2_data()
  s <- cj_mimosa2(b$x, positive = "TNF")$scores
  expect_true(all(s$ns1 == 0 & s$nu1 == 0))
  expect_error(cj_mimosa2(b$x, positive = "IL2"), "positive must be")
  expect_error(cj_mimosa2(b$x, baseline_visit = "v9"), "baseline_visit must be")
})

test_that("a donor without all four samples is left out with the reason", {
  x <- mimosa2_data()$x
  x$contrasts$control_sample_id[x$contrasts$contrast_id == "D03_pre_CD4_CMV"] <- NA
  x$contrasts$control_sample_id[x$contrasts$contrast_id == "D05_post_CD4_CMV"] <- NA
  expect_message(mm <- cj_mimosa2(x), "2 donor-fit")
  expect_setequal(mm$excluded$donor_id, c("D03", "D05"))
  expect_identical(mm$excluded$reason[mm$excluded$donor_id == "D03"], "no control sample at baseline")
  expect_identical(mm$excluded$reason[mm$excluded$donor_id == "D05"], "no control sample at follow-up")
  expect_false(any(c("D03", "D05") %in% mm$scores$donor_id))
  expect_identical(mm$groups$n_donors, 38L)
})

# Hand-computed beta-binomial kernel: lgamma form of B(a + n, b + N - n) / B(a, b).
bb_hand <- function(n, N, mu, nu) {
  a <- mu * exp(nu); b <- (1 - mu) * exp(nu)
  lgamma(a + n) + lgamma(b + N - n) - lgamma(a + b + N) - lgamma(a) - lgamma(b) + lgamma(a + b)
}

test_that("a beta-binomial term with N = 0, n = 0 is exactly 0", {
  expect_identical(cytoJointR:::mimosa2_bb(0, 0, 0.003, exp(7)), 0)
  expect_identical(cytoJointR:::mimosa2_bb(c(0, 5), c(0, 100), 0.05, 20)[1], 0)
  expect_false(anyNA(cytoJointR:::mimosa2_bb(0, 0, 1e-8, exp(20))))
})

test_that("component log-likelihoods equal hand-computed sums of beta-binomial terms", {
  theta <- list(mu = c(0.01, 0.005, 0.002, 0.001), log_precision = c(5, 6, 7, 8))
  n <- matrix(c(30, 10, 12, 8), 1, dimnames = list("d", c("s1", "u1", "s0", "u0")))
  N <- matrix(c(5000, 4000, 6000, 5500), 1, dimnames = dimnames(n))
  ll <- cytoJointR:::mimosa2_component_loglik(theta, n, N)
  expect_identical(dim(ll), c(1L, 11L))
  mu <- theta$mu; nu <- theta$log_precision
  # R2: s1 at level 1, u1 at level 3, {s0, u0} at level 4.
  r2 <- bb_hand(30, 5000, mu[1], nu[1]) + bb_hand(10, 4000, mu[3], nu[3]) +
    bb_hand(12 + 8, 6000 + 5500, mu[4], nu[4])
  # NSR1: {u1, u0} at level 1, {s1, s0} at level 4.
  nsr1 <- bb_hand(10 + 8, 4000 + 5500, mu[1], nu[1]) + bb_hand(30 + 12, 5000 + 6000, mu[4], nu[4])
  expect_equal(unname(ll[1, "R2"]), r2, tolerance = 1e-10)
  expect_equal(unname(ll[1, "NSR1"]), nsr1, tolerance = 1e-10)
  # A missing control (n = N = 0) drops its term: R2 with u1 missing loses the u1 term.
  n2 <- n; N2 <- N; n2[, "u1"] <- 0; N2[, "u1"] <- 0
  ll2 <- cytoJointR:::mimosa2_component_loglik(theta, n2, N2)
  expect_equal(unname(ll2[1, "R2"]), r2 - bb_hand(10, 4000, mu[3], nu[3]), tolerance = 1e-10)
})

test_that("the level means stay strictly ordered after the fit", {
  b <- mimosa2_data()
  f <- mimosa2_fit_data(b$x)
  expect_s3_class(f, "cj_mimosa2_fit")
  expect_true(all(diff(f$mu) < 0))
  expect_true(all(f$log_precision > 0))
  expect_equal(sum(f$pi), 1)
  expect_equal(unname(rowSums(f$z)), rep(1, 40))
})

test_that("q-values are running means of 1 - P(response) in decreasing order of P", {
  p <- c(a = 0.5, b = 0.9, c = 0.1, d = 0.8)
  # Sorted: 0.9, 0.8, 0.5, 0.1 -> 1 - p = 0.1, 0.2, 0.5, 0.9 -> running means 0.1, 0.15, 0.8 / 3, 0.425.
  expect_equal(cytoJointR:::mimosa2_qvalue(p), c(a = 0.8 / 3, b = 0.1, c = 0.425, d = 0.15))
  expect_true(all(is.na(cytoJointR:::mimosa2_qvalue(c(0.5, NA)))))
})

test_that("missing_control = 'marginal' keeps a donor with one missing control and says so", {
  x <- mimosa2_data()$x
  x$contrasts$control_sample_id[x$contrasts$contrast_id == "D07_pre_CD4_CMV"] <- NA
  # D09 has no stimulated sample at baseline: no contrast row at the baseline visit.
  x$contrasts <- x$contrasts[x$contrasts$contrast_id != "D09_pre_CD4_CMV", ]
  expect_message(mm <- cj_mimosa2(x, missing_control = "marginal"), "experimental and not yet calibrated")
  s <- mm$scores
  expect_identical(s$control_source[s$donor_id == "D07"], "marginal: u0 missing")
  expect_identical(s$Nu0[s$donor_id == "D07"], 0)
  expect_identical(s$nu0[s$donor_id == "D07"], 0)
  expect_true(is.finite(s$mimosa2_p_response[s$donor_id == "D07"]))
  expect_true(all(s$control_source[s$donor_id != "D07"] == "paired"))
  expect_false("D09" %in% s$donor_id)
  expect_identical(mm$excluded$reason[mm$excluded$donor_id == "D09"], "no stimulated sample at baseline")
  expect_identical(mm$groups$n_donors, 39L)
  expect_identical(mm$groups$status, "fitted")
  # The default still drops the donor with a missing control.
  expect_false("D07" %in% suppressMessages(cj_mimosa2(x))$scores$donor_id)
})

test_that("a donor with both controls missing stays excluded under 'marginal'", {
  x <- mimosa2_data()$x
  x$contrasts$control_sample_id[x$contrasts$contrast_id %in% c("D03_pre_CD4_CMV", "D03_post_CD4_CMV")] <- NA
  mm <- suppressMessages(cj_mimosa2(x, missing_control = "marginal"))
  expect_identical(mm$excluded$reason[mm$excluded$donor_id == "D03"], "no control sample at either visit")
})
