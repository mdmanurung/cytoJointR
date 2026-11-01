test_that("Gaussian compressed likelihood equals event-level likelihood", {
  set.seed(31)
  for (positive in c(TRUE, FALSE)) {
    # Generate from a truncated normal by rejection for test purposes.
    raw <- rnorm(10000, mean=if (positive) 1 else -1, sd=1)
    x <- if (positive) raw[raw>0] else raw[raw<=0]
    a <- cj_log_trunc_gaussian_raw(x, 0.3, 1.3, 0, positive)
    b <- cj_log_trunc_gaussian_compressed(length(x),sum(x),sum(x*x),
                                          0.3,1.3,0,positive)
    expect_equal(a,b,tolerance=1e-8)
  }
})

test_that("compression preserves all events and their moments", {
  x <- data.frame(sample_id=rep(c("u1","s1"),each=4),
    donor_id="D1",condition=rep(c("unstim","stim"),each=4),
    IFNG=c(-2,1,-3,2,-1,-2,3,4), TNF=c(-2,-1,2,3,-1,2,3,-2))
  out <- cj_compress(x, c("IFNG","TNF"), thresholds=c(IFNG=0,TNF=0))
  expect_equal(sum(out$summary$n),8L)
  expect_equal(sum(out$summary$sx_IFNG),sum(x$IFNG))
  expect_equal(sum(out$summary$ss_TNF),sum(x$TNF^2))
  expect_equal(nrow(out$summary),8L)
})
