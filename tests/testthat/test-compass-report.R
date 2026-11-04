# A stand-in for cj_run_compass(): the report only needs these columns, so the tests
# run without COMPASS. `resp` is the degree-1 response probability per contrast.
fake_compass <- function(x, resp) {
  ct <- x$contrasts[!is.na(x$contrasts$control_sample_id), ]
  data.frame(contrast_id = ct$contrast_id, donor_id = ct$donor_id, stimulation = ct$stimulation,
             visit = ct$visit, parent_population = ct$parent_population,
             compass_FS = resp / 7, compass_PFS = resp / 20,
             compass_response_degree1 = resp, compass_response_degree2 = resp / 2,
             control_source = "paired", control_fallback = FALSE, control_shift_check = NA_character_)
}

# Six donors: donors 1-3 respond strongly, 4-6 do not respond.
report_data <- function() {
  sim <- cj_simulate_counts(n_donors = 6, markers = c("IFNG", "TNF"), stimulations = "CMV",
                            N = 20000, seed = 3,
                            truth = list(alpha = c(IFNG = -6, TNF = -6), gamma = rbind(stim_CMV = c(0, 0)),
                                         theta0 = 0.5, delta = matrix(0, 1, 1), tau = c(0.3, 0.1),
                                         Omega = diag(2), tau_u = 0))
  sim$data
}

test_that("C0 report joins COMPASS and magnitudes by contrast and needs no Stan", {
  x <- report_data()
  n <- sum(!is.na(x$contrasts$control_sample_id))
  rep <- cj_compass_report(x, compass = fake_compass(x, rep(0.9, n)), n_draws = 500, seed = 1)
  expect_s3_class(rep, "cj_compass_report")
  expect_equal(rep$contrast_id, x$contrasts$contrast_id)
  expect_true(all(c("compass_FS", "compass_PFS", "compass_call", "rms", "rms_lower", "rms_upper",
                    "cps", "cps_lower", "cps_upper", "concordance", "discordant") %in% names(rep)))
  expect_true(all(rep$rms_lower <= rep$rms & rep$rms <= rep$rms_upper, na.rm = TRUE))
  expect_equal(unique(rep$magnitude_model), "c0")
  # Same draws, same numbers as the score layer.
  sc <- suppressMessages(cj_scores(cj_fit_dirichlet(x, n_draws = 500, seed = 1), REPORT_EPSILON,
                                   level = 0.9, ccs_gate = NULL))
  r <- sc[sc$score == "RMS" & is.na(sc$marker), ]
  expect_equal(rep$rms, r$estimate[match(rep$contrast_id, r$contrast_id)])
  # Per-marker table: one row per contrast and marker, no NICR without cells.
  mk <- attr(rep, "markers")
  expect_equal(nrow(mk), nrow(rep) * 2)
  expect_false("nicr" %in% names(mk))
})

test_that("concordance flags both kinds of disagreement and a missing control", {
  x <- report_data()
  n <- sum(!is.na(x$contrasts$control_sample_id))
  # Strong donors by construction: a huge stimulation effect for donors 1-3 only.
  x$counts[x$contrasts$stim_sample_id[1:3], ] <- x$counts[x$contrasts$stim_sample_id[1:3], ] +
    rep(c(0L, 400L, 400L, 400L), each = 3)
  resp <- c(0.9, 0.2, 0.9, 0.9, 0.2, 0.2)[seq_len(n)]
  rep <- cj_compass_report(x, compass = fake_compass(x, resp), n_draws = 500, seed = 1)
  mag <- rep$rms_lower > 1e-4 # the default floor is the RMS epsilon
  expect_equal(rep$magnitude_positive, mag)
  expect_true(all(mag[1:3]))
  expect_equal(as.character(rep$concordance[1:3]),
               c("response_with_magnitude", "magnitude_without_response", "response_with_magnitude"))
  # Donor 4 is called a responder, but its RMS interval includes 0.
  expect_false(mag[4])
  expect_equal(as.character(rep$concordance[4]), "response_without_magnitude")
  expect_equal(as.character(rep$concordance[5:6]), rep("no_response_no_magnitude", 2))
  expect_equal(rep$discordant, c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE))
  expect_equal(rep$compass_call, resp > 0.5)
  # response_threshold moves the call; magnitude_floor moves the magnitude test.
  r2 <- cj_compass_report(x, compass = fake_compass(x, resp), n_draws = 500, seed = 1,
                          response_threshold = 0.95)
  expect_true(all(!r2$compass_call))
  r3 <- cj_compass_report(x, compass = fake_compass(x, resp), n_draws = 500, seed = 1,
                          magnitude_floor = 1)
  expect_true(all(!r3$magnitude_positive))
})

test_that("a contrast without control is not assessed, never concordant by default", {
  x <- report_data()
  x$contrasts$control_sample_id[2] <- NA_character_
  n <- sum(!is.na(x$contrasts$control_sample_id))
  comp <- fake_compass(x, rep(0.9, n))
  # COMPASS with missing_control = "donor" returns a row for it; give it one.
  extra <- comp[1, ]
  extra$contrast_id <- x$contrasts$contrast_id[2]
  extra$donor_id <- x$contrasts$donor_id[2]
  rep <- cj_compass_report(x, compass = rbind(comp, extra), n_draws = 300, seed = 1)
  i <- match(x$contrasts$contrast_id[2], rep$contrast_id)
  expect_true(is.na(rep$rms[i]))
  expect_equal(as.character(rep$concordance[i]), "not_assessed")
  expect_true(is.na(rep$discordant[i]))
})

test_that("cell input adds NICR with bootstrap intervals", {
  sim <- cj_simulate_cells(donors = 5, cells_per_sample = 600, seed = 1)
  samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
  samples$visit <- factor("v1")
  samples$parent_population <- "CD4"
  x <- cj_from_cells(sim$cells, c("IFNG", "TNF", "IL2"), samples, transform = "none",
                     control = "unstim", method = "quantile")
  n <- sum(!is.na(x$contrasts$control_sample_id))
  rep <- cj_compass_report(sim$cells, c("IFNG", "TNF", "IL2"), samples, transform = "none",
                           control = "unstim", gate_args = list(method = "quantile"),
                           compass = fake_compass(x, rep(0.9, n)), n_boot = 20, n_draws = 300, seed = 1)
  mk <- attr(rep, "markers")
  expect_true(all(c("nicr", "nicr_lo", "nicr_hi", "poly_nicr", "poly_nicr_lo", "poly_nicr_hi") %in% names(mk)))
  expect_true(all(mk$nicr_lo <= mk$nicr_hi, na.rm = TRUE))
  # Same NICR point estimate as the fluorescence scores.
  fl <- cj_fluorescence_scores(x)
  expect_equal(mk$nicr, fl$nicr[match(paste(mk$contrast_id, mk$marker), paste(fl$contrast_id, fl$marker))])
})

test_that("m1dm option uses the structured fit", {
  skip_if_no_stan()
  x <- report_data()
  n <- sum(!is.na(x$contrasts$control_sample_id))
  rep <- cj_compass_report(x, magnitude_model = "m1dm", compass = fake_compass(x, rep(0.9, n)),
                           fit_args = list(chains = 2, iter_warmup = 200, iter_sampling = 200,
                                           refresh = 0, show_messages = FALSE))
  expect_equal(unique(rep$magnitude_model), "m1dm")
  expect_false(anyNA(rep$rms))
})

test_that("the forest plot maps magnitude to position and COMPASS response to colour", {
  x <- report_data()
  n <- sum(!is.na(x$contrasts$control_sample_id))
  rep <- cj_compass_report(x, compass = fake_compass(x, seq(0.1, 0.9, length.out = n)),
                           n_draws = 300, seed = 1)
  p <- cj_plot_compass_report(rep)
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), n)
  expect_equal(sort(p$data$estimate, decreasing = TRUE), p$data$estimate)
  expect_equal(p$data$response[match(rep$contrast_id, as.character(p$data$label))], rep$compass_response_degree1)
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("a real COMPASS run feeds the report", {
  skip_if_not_installed("COMPASS")
  sim <- cj_simulate_cells(donors = 6, cells_per_sample = 400, seed = 1)
  samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
  samples$visit <- factor("v1")
  samples$parent_population <- "CD4"
  rep <- cj_compass_report(sim$cells, c("IFNG", "TNF", "IL2"), samples, transform = "none",
                           control = "unstim", gate_args = list(method = "quantile"),
                           compass_args = list(iterations = 200, replications = 2), n_draws = 300)
  expect_false(all(is.na(rep$compass_FS)))
  expect_false(all(is.na(rep$rms)))
})
