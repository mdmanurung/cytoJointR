test_that("cj_compass runs real COMPASS and returns native results", {
  skip_if_not_installed("COMPASS")
  d <- cj_example_draws(n_donor = 6, n_draw = 5, markers = c("IFNG", "TNF", "IL2"),
                        visits = "v1", stimulations = c("CMV", "SEB"),
                        missing_control = "D6_v1_CD4_SEB")
  x <- example_data(d)
  cc <- cj_compass(x, iterations = 200, replications = 2)
  expect_s3_class(cc, "cj_compass")
  expect_true(all(vapply(cc$fits, inherits, TRUE, "COMPASSResult")))
  # COMPASS's own pattern-name format, null pattern last.
  expect_true("IFNG&!TNF&IL2" %in% colnames(cc$fits[[1]]$data$n_s))
  expect_equal(tail(colnames(cc$fits[[1]]$data$n_s), 1), "!IFNG&!TNF&!IL2")
  # The contrast without control is excluded and reported, not dropped silently.
  expect_equal(cc$excluded$contrast_id, "D6_v1_CD4_SEB")
  expect_false("D6" %in% cc$scores$donor_id[cc$scores$stimulation == "SEB"])
  # Tidy scores equal COMPASS's own numbers.
  f <- cc$fits[["CMV | v1 | CD4"]]
  tid <- cc$scores[cc$scores$stimulation == "CMV", ]
  expect_equal(tid$compass_FS, unname(COMPASS::FunctionalityScore(f)[tid$donor_id]))
  # Response() is a one-column matrix by donor: no NA, same numbers as COMPASS.
  expect_false(anyNA(cc$scores[c("compass_response_degree1", "compass_response_degree2")]))
  expect_equal(tid$compass_response_degree2,
               unname(COMPASS::Response(f, degree = 2)[tid$donor_id, 1]))
  expect_true(all(tid$compass_response_degree1 >= 0 & tid$compass_response_degree1 <= 1))
  # n_s and n_u are not swapped: stimulated above control gives FS well above
  # the control-vs-control fit.
  dd <- f$data
  null <- COMPASS::SimpleCOMPASS(dd$n_u, dd$n_u, dd$meta, "donor_id",
                                 iterations = 200, replications = 2)
  expect_gt(mean(tid$compass_FS), mean(COMPASS::FunctionalityScore(null)))
  # The adapter computes the posterior, so COMPASS's own accessors work.
  pd <- COMPASS::PosteriorDiff(f)
  expect_equal(nrow(pd), nrow(f$data$n_s))
  expect_equal(rownames(pd), rownames(f$data$n_s))
  expect_false(anyNA(pd))
  # The default filter's effect is reported.
  expect_true(all(c("n_removed", "removed") %in% names(cc$filters)))
  expect_equal(nrow(cc$filters), 2)
  # The plot helper reaches COMPASS's own heatmap.
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(cj_plot_compass(cc, row_annotation = "stimulation"))
})

# Three visits, one stimulation. D1 lacks the v2 control (v1 and v3 tie), D2 lacks v3,
# D3 has no control at any visit.
mc_data <- function() {
  miss <- c("D1_v2_CD4_CMV", "D2_v3_CD4_CMV", "D3_v1_CD4_CMV", "D3_v2_CD4_CMV", "D3_v3_CD4_CMV")
  d <- cj_example_draws(n_donor = 6, n_draw = 5, markers = c("IFNG", "TNF", "IL2"),
                        visits = c("v1", "v2", "v3"), stimulations = "CMV",
                        missing_control = miss)
  x <- example_data(d)
  # As in cells data: a visit without a contrast control has no control sample at all.
  gone <- c("D1_v2_control", "D2_v3_control", paste0("D3_", c("v1", "v2", "v3"), "_control"))
  keep <- !x$samples$sample_id %in% gone
  cytoJointR:::new_cj_data(x$counts[keep, ], x$samples[keep, ], x$contrasts, x$bits)
}

test_that("missing_control = 'drop' is the default and keeps the old behaviour", {
  skip_if_not_installed("COMPASS")
  x <- mc_data()
  cc <- cj_compass(x, iterations = 100, replications = 2)
  expect_equal(nrow(cc$excluded), 5)
  expect_true(all(cc$scores$control_source == "paired"))
  expect_false(any(cc$scores$donor_id == "D3"))
  expect_no_message(cj_compass(x, iterations = 100, replications = 2, missing_control = "drop"))
})

test_that("'population' runs a zero control row and labels it", {
  skip_if_not_installed("COMPASS")
  x <- mc_data()
  expect_message(cc <- cj_compass(x, missing_control = "population", iterations = 100,
                                  replications = 2), "5 of 18 contrast")
  expect_equal(nrow(cc$excluded), 0)
  expect_equal(sum(cc$scores$control_source == "population"), 5)
  expect_equal(sum(cc$scores$control_source == "paired"), 13)
  f <- cc$fits[["CMV | v2 | CD4"]]
  expect_equal(unname(f$data$counts_u["D1"]), 0)
  # Every COMPASS output path is finite for the zero-control donors.
  expect_false(anyNA(cc$scores[c("compass_FS", "compass_PFS", "compass_response_degree1",
                                 "compass_response_degree2")]))
  for (f in cc$fits) {
    expect_false(anyNA(COMPASS::PosteriorDiff(f)))
    expect_true(all(is.finite(COMPASS::FunctionalityScore(f))))
  }
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(cj_plot_compass(cc, row_annotation = "control_source"))
})

test_that("'donor' borrows the nearest other visit, ties to the earlier, and falls back", {
  skip_if_not_installed("COMPASS")
  x <- mc_data()
  cc <- suppressMessages(cj_compass(x, missing_control = "donor", iterations = 100, replications = 2))
  src <- function(id) cc$scores$control_source[cc$scores$contrast_id == id]
  expect_equal(src("D1_v2_CD4_CMV"), "donor_visit:v1")   # v1 and v3 tie
  expect_equal(src("D2_v3_CD4_CMV"), "donor_visit:v2")
  expect_equal(src("D3_v1_CD4_CMV"), "population")
  expect_equal(src("D4_v1_CD4_CMV"), "paired")
  expect_equal(cc$scores$control_fallback[cc$scores$donor_id == "D3"], rep(TRUE, 3))
  expect_false(any(cc$scores$control_fallback[cc$scores$donor_id != "D3"]))
  # The borrowed row holds the donor's own v1 control counts.
  f <- cc$fits[["CMV | v2 | CD4"]]
  n1 <- x$counts["D1_v1_control", ]
  expect_equal(unname(rowSums(f$data$n_u)["D1"]), sum(n1))
  names(n1) <- cytoJointR:::compass_pattern_names(x$bits)   # counts columns follow bits rows
  expect_equal(unname(f$data$n_u["D1", ]), unname(n1[colnames(f$data$n_u)]))
  expect_message(cj_compass(x, missing_control = "donor", iterations = 100, replications = 2),
                 "against the population control \\(3 of them a fallback")
})

test_that("the category filter ignores zero-control rows", {
  skip_if_not_installed("COMPASS")
  x <- mc_data()
  ps <- suppressMessages(cj_compass(x, missing_control = "population", iterations = 100, replications = 2))
  # Default filter reads stimulated counts only: same removals as on the stimulated matrix.
  s <- x$contrasts$stim_sample_id[x$contrasts$visit == "v2"]
  keep <- COMPASS_DEFAULT_FILTER(x$counts[s, ])
  fl <- ps$filters[ps$filters$visit == "v2", ]
  expect_equal(fl$n_removed, sum(!keep[colnames(x$counts)] & rowSums(x$bits) > 0))
  expect_equal(fl$n_zero_control, 2)   # D1 and D3 at v2
  # Control cells moved come only from rows that have a control.
  rem <- strsplit(fl$removed, "; ")[[1]]
  has <- !is.na(x$contrasts$control_sample_id[x$contrasts$visit == "v2"])
  cs <- x$contrasts$control_sample_id[x$contrasts$visit == "v2"][has]
  cols <- colnames(x$counts)[match(rem, cytoJointR:::compass_pattern_names(x$bits))]
  expect_equal(fl$cells_moved_control, sum(x$counts[cs, cols]))
})

test_that("a run with no control in any donor is excluded with a reason", {
  skip_if_not_installed("COMPASS")
  x <- mc_data()
  x$contrasts$control_sample_id[x$contrasts$visit == "v1"] <- NA
  cc <- suppressMessages(cj_compass(x, missing_control = "population", iterations = 100, replications = 2))
  expect_equal(sum(cc$excluded$visit == "v1"), 6)
  expect_match(cc$excluded$reason[1], "no donor")
  expect_false(any(cc$scores$visit == "v1"))
})

test_that("missing_control options say they are experimental, once per call, with FPR numbers", {
  skip_if_not_installed("COMPASS")
  x <- mc_data()
  m <- testthat::capture_messages(
    cj_compass(x, missing_control = "population", iterations = 100, replications = 2))
  expect_equal(sum(grepl("experimental", m)), 1)
  expect_true(any(grepl("0.19 to 0.39", m)))
  m <- testthat::capture_messages(suppressWarnings(
    cj_compass(x, missing_control = "donor", iterations = 100, replications = 2)))
  expect_equal(sum(grepl("experimental", m)), 1)
  expect_true(any(grepl("0.81 to 0.94", m)))
  expect_no_match(testthat::capture_messages(
    cj_compass(x, iterations = 100, replications = 2)), "experimental")
})

test_that("compass_shift_check separates a large background difference from none", {
  bits <- cj_bits(c("A", "B"))
  null <- which(rowSums(bits) == 0)
  mk <- function(frac_pos, n = 50000) {
    v <- rep(0L, nrow(bits)); v[-null] <- as.integer(round(n * frac_pos / 3)); v[null] <- n - sum(v); v
  }
  same <- rbind(mk(0.05), mk(0.052))
  expect_equal(cytoJointR:::compass_shift_check(same, bits), "ok")
  big <- rbind(mk(0.05), mk(0.10))
  expect_equal(cytoJointR:::compass_shift_check(big, bits), "differs")
  expect_equal(cytoJointR:::compass_shift_check(same[1, , drop = FALSE], bits), "untestable")
})

test_that("control_shift_check is set and warns when the donor's controls differ", {
  skip_if_not_installed("COMPASS")
  x <- mc_data()
  null <- which(rowSums(x$bits) == 0)
  # D1 lacks the v2 control. Same background at v1 and v3: ok. D2 lacks v3 and has v1, v2.
  expect_no_warning(cc <- suppressMessages(
    cj_compass(x, missing_control = "donor", iterations = 100, replications = 2)))
  chk <- function(cc, id) cc$scores$control_shift_check[cc$scores$contrast_id == id]
  expect_equal(chk(cc, "D1_v2_CD4_CMV"), "ok")
  expect_true(is.na(chk(cc, "D4_v1_CD4_CMV")))
  # Move 60% of D1's v3 control cells from the null pattern into a positive pattern.
  y <- x
  r <- "D1_v3_control"
  moved <- as.integer(0.6 * y$counts[r, null])
  y$counts[r, null] <- y$counts[r, null] - moved
  y$counts[r, setdiff(seq_len(ncol(y$counts)), null)[1]] <- y$counts[r, setdiff(seq_len(ncol(y$counts)), null)[1]] + moved
  expect_warning(cc <- suppressMessages(
    cj_compass(y, missing_control = "donor", iterations = 100, replications = 2)), "control_shift_check")
  expect_equal(chk(cc, "D1_v2_CD4_CMV"), "differs")
  # D2 lacks v3 and has v1, v2: untestable would need one control; D3 has none and is not checked.
  expect_true(is.na(chk(cc, "D3_v1_CD4_CMV")))
})

test_that("'donor' skips an empty control instead of failing", {
  skip_if_not_installed("COMPASS")
  x <- mc_data()
  # D2 lacks its v3 control; empty its v1 control so only v2 is usable.
  x$counts["D2_v1_control", ] <- 0L
  cc <- suppressWarnings(suppressMessages(
    cj_compass(x, missing_control = "donor", iterations = 100, replications = 2)))
  expect_equal(cc$scores$control_source[cc$scores$contrast_id == "D2_v3_CD4_CMV"], "donor_visit:v2")
  expect_equal(cc$scores$control_shift_check[cc$scores$contrast_id == "D2_v3_CD4_CMV"], "untestable")
})
