test_that("pattern contributions sum to the global score in every draw", {
  d <- cj_example_draws(n_draw = 15)
  eng <- cytoJointR:::score_engine(d, keep_patterns = TRUE)
  for (e in eng) for (sc in c("RMS", "CPS", "CCS"))
    expect_equal(unname(rowSums(e$contrib[, , sc])), unname(e$G[, sc]))
  # Posterior means reconcile too; interval bounds are never summed.
  pt <- cj_patterns(d, c("RMS", "CPS", "CCS"), EPS)
  sc <- cj_scores(d, EPS)
  tot <- aggregate(estimate ~ contrast_id + score, pt, sum)
  ref <- merge(tot, sc[is.na(sc$marker), c("contrast_id", "score", "estimate")],
               by = c("contrast_id", "score"))
  expect_equal(ref$estimate.x, ref$estimate.y)
})

test_that("a shared control sample gives each contrast the same draw", {
  d <- cj_example_draws(n_draw = 10, stimulations = c("CMV", "SEB"))
  # Make SEB identical to CMV: both contrasts must then agree draw by draw.
  d$p[, "D1_v1_SEB", ] <- d$p[, "D1_v1_CMV", ]
  sd <- cj_score_draws(d)
  pick <- function(id) sd[sd$contrast_id == id & sd$score == "RMS", "value"]
  expect_identical(pick("D1_v1_CMV"), pick("D1_v1_SEB"))
  # And the value is the stimulated minus control mass in the SAME draw k.
  pos <- rowSums(d$bits) > 0
  by_hand <- rowSums(d$p[, "D1_v1_CMV", pos]) - rowSums(d$p[, "D1_v1_control", pos])
  expect_equal(pick("D1_v1_CMV"), unname(by_hand))
})

test_that("an IPF failure is flagged and kept, not dropped", {
  d <- cj_example_draws(n_draw = 20)
  real <- cytoJointR::cj_counterfactual
  n <- 0
  testthat::local_mocked_bindings(cj_counterfactual = function(...) {
    n <<- n + 1
    if (n %% 5 == 0) stop("IPF failed to converge")
    real(...)
  })
  sd <- cj_score_draws(d)
  expect_equal(nrow(sd), 20 * 12 * 9)  # nothing dropped: 12 contrasts, 9 score columns
  failed <- sd[sd$status == "ipf_failed", ]
  expect_true(nrow(failed) > 0 && all(is.na(failed$value)))
  expect_true(all(sd$value[sd$status == "ok"] |> is.finite()))
  s <- cj_scores(d, EPS, max_failure = 0.5)
  expect_true(all(s$n_failed == 4))
  expect_true(all(s$status == "ok_with_failed_draws"))
  strict <- cj_scores(d, EPS)  # 20% failed > default 1% tolerance
  expect_true(all(strict$status == "failure_tolerance_exceeded" & is.na(strict$estimate)))
})

test_that("the three epsilon probabilities sum to one", {
  s <- cj_scores(cj_example_draws(n_draw = 100), EPS)
  expect_equal(s$p_positive + s$p_negative + s$p_negligible, rep(1, nrow(s)))
  expect_equal(unique(s$epsilon[s$score == "CPS"]), 0.02)
})

test_that("epsilon is required and named per score", {
  d <- cj_example_draws(n_draw = 10)
  expect_error(cj_scores(d), "epsilon")
  expect_error(cj_scores(d, 0.01), "named")
  expect_error(cj_scores(d, EPS[-1]), "lacks: RMS")
})

test_that("negative effects stay negative and are not clipped", {
  d <- cj_example_draws(n_draw = 100, effect = -0.5)
  sd <- cj_score_draws(d)
  expect_true(all(sd$value[sd$score == "RMS"] < 0))
  s <- cj_scores(d, EPS)
  rms <- s[s$score == "RMS", ]
  expect_true(all(rms$estimate < 0 & rms$p_negative == 1 & rms$p_positive == 0))
  pl <- cj_plot_forest(s, "RMS")
  expect_s3_class(pl, "ggplot")
  expect_true(all(ggplot2::layer_data(pl, 4)$x < 0))
})

test_that("a contrast without control is unavailable with a status, not bare NA", {
  d <- cj_example_draws(n_draw = 10, missing_control = "D1_v1_CMV")
  sd <- cj_score_draws(d)
  expect_true(all(sd$status[sd$contrast_id == "D1_v1_CMV"] == "unavailable"))
  s <- cj_scores(d, EPS)
  u <- s[s$contrast_id == "D1_v1_CMV", ]
  expect_true(all(u$status == "unavailable" & is.na(u$estimate)))
  expect_true(all(s$status[s$contrast_id != "D1_v1_CMV"] == "ok"))
})

test_that("cj_fit input adds raw counts to the pattern table", {
  d <- cj_example_draws(n_draw = 10)
  x <- example_data(d)
  fit <- cytoJointR:::new_cj_fit(x, d, NULL, NULL, NULL, NULL)
  pt <- cj_patterns(fit, "RMS", EPS["RMS"])
  r <- pt[pt$contrast_id == "D1_v1_CMV", ][1, ]
  expect_equal(r$n_S, unname(x$counts["D1_v1_CMV", r$pattern]))
  expect_equal(r$n_U, unname(x$counts["D1_v1_control", r$pattern]))
  expect_equal(r$N_S, 1000)
  expect_false("n_U" %in% names(cj_patterns(d, "RMS", EPS["RMS"])))
})

test_that("markers are labelled overlapping and plots build", {
  d <- cj_example_draws(n_draw = 20)
  expect_true(all(cj_markers(d, 0.005)$overlapping))
  pt <- cj_patterns(d, c("RMS", "CCS"), EPS)
  expect_s3_class(cj_plot_pattern_heatmap(pt, "p_positive", score = "RMS"), "ggplot")
  expect_s3_class(cj_plot_pattern_heatmap(pt, "estimate", score = "CCS"), "ggplot")
  expect_s3_class(cj_plot_contributions(pt, "RMS", "D1_v1_CMV"), "ggplot")
})
