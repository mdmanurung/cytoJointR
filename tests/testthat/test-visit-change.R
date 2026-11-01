test_that("visit change is the per-draw difference of net scores", {
  d <- cj_example_draws(n_draw = 30)
  vc <- cj_visit_change(d, EPS, "v1", "v2")
  sd <- cj_score_draws(d)
  one <- function(id) sd$value[sd$contrast_id == id & sd$score == "CCS"]
  delta <- one("D2_v2_SEB") - one("D2_v1_SEB")
  r <- vc$donor[vc$donor$donor_id == "D2" & vc$donor$stimulation == "SEB" &
                  vc$donor$score == "CCS", ]
  expect_equal(r$estimate, mean(delta))
  expect_equal(r$p_positive, mean(delta > EPS[["CCS"]]))
  # Equal donor weights: the population mean is the mean of donor means.
  pop <- vc$population[vc$population$stimulation == "SEB" & vc$population$score == "CCS", ]
  donors <- vc$donor[vc$donor$stimulation == "SEB" & vc$donor$score == "CCS", ]
  expect_equal(pop$estimate, mean(donors$estimate))
  expect_match(pop$target, "observed-donor average")
})

test_that("an unavailable contrast propagates to visit change", {
  d <- cj_example_draws(n_draw = 10, missing_control = "D1_v2_CMV")
  vc <- cj_visit_change(d, EPS, "v1", "v2")
  u <- vc$donor[vc$donor$donor_id == "D1" & vc$donor$stimulation == "CMV", ]
  expect_true(all(u$status == "unavailable" & is.na(u$estimate) & !is.na(u$reason)))
  ok <- vc$donor[vc$donor$donor_id == "D2" & vc$donor$stimulation == "CMV", ]
  expect_true(all(ok$status == "ok"))
  pop <- vc$population[vc$population$stimulation == "CMV", ]
  expect_true(all(pop$n_donors == 2 & pop$n_donors_unavailable == 1))
  # All donors unavailable -> population unavailable too.
  d2 <- cj_example_draws(n_draw = 10, n_donor = 1, missing_control = "D1_v1_CMV")
  expect_true(all(cj_visit_change(d2, EPS, "v1", "v2")$population$status[
    cj_visit_change(d2, EPS, "v1", "v2")$population$stimulation == "CMV"] == "unavailable"))
})

test_that("from and to must be distinct visit levels", {
  d <- cj_example_draws(n_draw = 5)
  expect_error(cj_visit_change(d, EPS, "v1", "v1"), "different")
  expect_error(cj_visit_change(d, EPS, "v1", "v9"), "visit levels")
})
