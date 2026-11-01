test_that("cj_compass runs real COMPASS and returns native results", {
  skip_if_not_installed("COMPASS")
  d <- cj_example_draws(n_donor = 6, n_draw = 5, markers = c("IFNG", "TNF", "IL2"),
                        visits = "v1", stimulations = c("CMV", "SEB"),
                        missing_control = "D6_v1_SEB")
  x <- example_data(d)
  cc <- cj_compass(x, iterations = 200, replications = 2)
  expect_s3_class(cc, "cj_compass")
  expect_true(all(vapply(cc$fits, inherits, TRUE, "COMPASSResult")))
  # COMPASS's own pattern-name format, null pattern last.
  expect_true("IFNG&!TNF&IL2" %in% colnames(cc$fits[[1]]$data$n_s))
  expect_equal(tail(colnames(cc$fits[[1]]$data$n_s), 1), "!IFNG&!TNF&!IL2")
  # The contrast without control is excluded and reported, not dropped silently.
  expect_equal(cc$excluded$contrast_id, "D6_v1_SEB")
  expect_false("D6" %in% cc$scores$donor_id[cc$scores$stimulation == "SEB"])
  # Tidy scores equal COMPASS's own numbers.
  f <- cc$fits[["CMV | v1 | CD4"]]
  tid <- cc$scores[cc$scores$stimulation == "CMV", ]
  expect_equal(tid$compass_FS, unname(COMPASS::FunctionalityScore(f)[tid$donor_id]))
  # The default filter's effect is reported.
  expect_true(all(c("n_removed", "removed") %in% names(cc$filters)))
  expect_equal(nrow(cc$filters), 2)
  # The plot helper reaches COMPASS's own heatmap.
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  expect_no_error(cj_plot_compass(cc, row_annotation = "stimulation"))
})
