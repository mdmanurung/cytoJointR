test_that("simulator returns complete truths and missing controls", {
  x <- cj_simulate_cells(donors=12,cells_per_sample=100,
                         markers=c("IFNG","TNF","IL2"),
                         missing_controls=.25,seed=4)
  expect_equal(nrow(x$truth),12L)
  expect_equal(sum(!x$truth$has_control),3L)
  out <- cj_compress(x$cells,c("IFNG","TNF","IL2"),
                     thresholds=c(IFNG=0,TNF=0,IL2=0))
  expect_equal(sum(out$summary$n),nrow(x$cells))
  expect_true(all(is.finite(as.matrix(x$truth[,4:8]))))
})
