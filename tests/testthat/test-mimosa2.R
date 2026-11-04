test_that("cj_mimosa2 passes each donor's four samples to MIMOSA2 by name", {
  skip_if_not_installed("MIMOSA2")
  b <- mimosa2_data()
  mm <- cj_mimosa2(b$x)
  expect_s3_class(mm, "cj_mimosa2")
  expect_s3_class(mm$fits[[1]], "MIMOSA2model")
  expect_identical(mm$groups$status, "fitted")
  s <- mm$scores[mm$scores$donor_id == "D07", ]
  ifng <- b$x$counts[, "IFNG+ TNF-"]
  expect_equal(c(s$ns1, s$nu1, s$ns0, s$nu0),
               unname(ifng[c("D07_post_CMV", "D07_post_control", "D07_pre_CMV", "D07_pre_control")]))
  expect_equal(c(s$Ns1, s$Nu1, s$Ns0, s$Nu0), rep(20000, 4))
  expect_equal(s$delta, (s$ns1 - s$nu1 - s$ns0 + s$nu0) / 20000)
})

test_that("cj_mimosa2 separates simulated responders from non-responders", {
  skip_if_not_installed("MIMOSA2")
  b <- mimosa2_data()
  s <- cj_mimosa2(b$x)$scores
  truth <- s$donor_id %in% b$responders
  # A clear effect (0.002 to 0.006 at 20,000 cells): every responder outranks every non-responder.
  expect_gt(min(s$mimosa2_p_response[truth]), max(s$mimosa2_p_response[!truth]))
  expect_false(any(s$mimosa2_response[!truth]))
})

test_that("positive = <marker> counts that marker's cells; bad arguments fail", {
  skip_if_not_installed("MIMOSA2")
  b <- mimosa2_data()
  s <- cj_mimosa2(b$x, positive = "TNF")$scores
  expect_true(all(s$ns1 == 0 & s$nu1 == 0))
  expect_error(cj_mimosa2(b$x, positive = "IL2"), "positive must be")
  expect_error(cj_mimosa2(b$x, baseline_visit = "v9"), "baseline_visit must be")
})

test_that("a donor without all four samples is left out with the reason", {
  skip_if_not_installed("MIMOSA2")
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
