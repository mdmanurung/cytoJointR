# A stand-in cj_compass result: one scores row per contrast, with the COMPASS
# response probability set by hand so every agreement class is reachable.
stand_in_compass <- function(x, resp) {
  ct <- x$contrasts
  key <- paste(ct$donor_id, ct$visit)
  list(scores = data.frame(contrast_id = ct$contrast_id, donor_id = ct$donor_id, stimulation = ct$stimulation,
                           visit = ct$visit, parent_population = ct$parent_population,
                           compass_FS = 0.1, compass_PFS = 0.1,
                           compass_response_degree1 = unname(resp[key])))
}

test_that("cj_baseline_report joins both COMPASS visits to MIMOSA2 and classes every row", {
  skip_if_not_installed("MIMOSA2")
  b <- mimosa2_data()
  donors <- sprintf("D%02d", 1:40)
  # COMPASS at post: responders D01-D20 called, plus D21-D25 (called at pre too) and D26-D27 (post only).
  post <- setNames(ifelse(donors %in% c(sprintf("D%02d", 1:27)), 0.9, 0.1), paste(donors, "post"))
  pre <- setNames(ifelse(donors %in% sprintf("D%02d", 21:25), 0.9, 0.1), paste(donors, "pre"))
  comp <- stand_in_compass(b$x, c(post, pre))
  mm <- cj_mimosa2(b$x)
  rep <- cj_baseline_report(b$x, compass = comp, mimosa2 = mm)
  expect_identical(nrow(rep), 40L)
  expect_true(all(rep$visit == "post" & rep$baseline_visit == "pre"))
  cls <- setNames(as.character(rep$agreement), rep$donor_id)
  expect_true(all(cls[sprintf("D%02d", 1:20)] == "both"))
  expect_true(all(cls[sprintf("D%02d", 21:25)] == "pre_existing"))
  expect_true(all(cls[c("D26", "D27")] == "compass_only"))
  expect_true(all(cls[sprintf("D%02d", 28:40)] == "neither"))
  # Values come from the right visit, by name.
  expect_equal(rep$compass_response_baseline[rep$donor_id == "D21"], 0.9)
  expect_equal(rep$compass_response_followup[rep$donor_id == "D30"], 0.1)
  expect_equal(rep$mimosa2_p_response[rep$donor_id == "D07"],
               mm$scores$mimosa2_p_response[mm$scores$donor_id == "D07"])
  expect_identical(levels(rep$agreement), cytoJointR:::BASELINE_AGREEMENT_LEVELS)
})

test_that("a donor MIMOSA2 left out, or without a baseline COMPASS call, is not assessed", {
  skip_if_not_installed("MIMOSA2")
  b <- mimosa2_data()
  x <- b$x
  x$contrasts$control_sample_id[x$contrasts$contrast_id == "D03_pre_CD4_CMV"] <- NA
  resp <- setNames(rep(0.9, 80), paste(rep(sprintf("D%02d", 1:40), 2), rep(c("post", "pre"), each = 40)))
  resp["D30 pre"] <- NA # e.g. a baseline contrast COMPASS could not run
  comp <- stand_in_compass(x, resp)
  rep <- suppressMessages(cj_baseline_report(x, compass = comp, mimosa2 = cj_mimosa2(x)))
  expect_identical(as.character(rep$agreement[rep$donor_id == "D03"]), "not_assessed")
  expect_true(is.na(rep$mimosa2_response[rep$donor_id == "D03"]))
  # D30 is a non-responder called by COMPASS at post: without its baseline call the class is unknown.
  expect_identical(as.character(rep$agreement[rep$donor_id == "D30"]), "not_assessed")
})

test_that("cj_baseline_report runs real COMPASS when no result is given", {
  skip_if_not_installed("MIMOSA2")
  skip_if_not_installed("COMPASS")
  b <- mimosa2_data()
  rep <- suppressMessages(cj_baseline_report(b$x, compass_args = list(iterations = 400, replications = 2)))
  expect_identical(nrow(rep), 40L)
  expect_false(anyNA(rep$compass_response_followup))
  expect_false(anyNA(rep$compass_response_baseline))
})
