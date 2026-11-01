test_that("signed scores and exact contributions reconcile", {
  b <- cj_bits(c("IFNG", "TNF"))
  u <- c(9940,30,20,10)/10000
  s <- c(9700,110,70,120)/10000
  Iu <- Is <- matrix(0, 4, 2)
  Iu[b==1] <- 100; Is[b==1] <- 150
  w <- b / rep(colSums(b), each=nrow(b))
  x <- cj_score_reference(u,s,b,Iu,Is,w,c(IFNG=10,TNF=25))
  expect_equal(unname(x$global["RMS"]), 0.024)
  expect_equal(unname(x$global["CPS"]), 120/300-10/60)
  expect_equal(sum(x$patterns$RMS), unname(x$global["RMS"]))
  expect_equal(sum(x$patterns$CCS), unname(x$global["CCS"]))
  expect_equal(x$patterns$poly_absolute,
               x$patterns$poly_marginal+x$patterns$CCS)
  expect_equal(unname(x$global["poly_absolute"]),
               unname(x$global["poly_marginal"]+x$global["CCS"]))
  expect_equal(x$markers$NICR,
               x$markers$frequency_contribution+x$markers$intensity_contribution)
  expect_equal(x$markers$FIS, c(50,50))
  expect_equal(x$markers$FIS_standardized, c(5,2))
  expect_true(cj_score_reference(s,u,b)$global["RMS"] < 0)
  # No arbitrary small-denominator truncation in the new algebra.
  rare <- c(1-3e-12,1e-12,1e-12,1e-12)
  expect_true(is.finite(cj_score_reference(rare,rare,b)$global["CPS"]))
  expect_equal(unname(cj_scores_from_probs(u,s,b)["coordination"]),
               unname(x$global["CCS"]))
  expect_warning(cj_scores_from_probs(u,s,b,min_functional=1e-6),"deprecated")
})

test_that("fixed marker margins do not force zero magnitude", {
  b <- cj_bits(c("A","B"))
  u <- c(.64,.16,.16,.04); s <- c(.68,.12,.12,.08)
  expect_equal(drop(crossprod(b,u)),drop(crossprod(b,s)))
  x <- cj_score_reference(u,s,b)
  expect_equal(unname(x$global["RMS"]), -.04)
  expect_equal(unname(x$global["CCS"]), .04)
})

test_that("invalid partitions, intensity and draws are rejected", {
  b <- cj_bits(c("A","B")); u <- rep(0.25,4)
  expect_error(cj_score_reference(c(1,0,0,0),u,b), "positive")
  expect_error(cj_score_reference(u,u,b[c(1,1,3,4),]), "enumerate")
  expect_error(cj_score_reference(u,u,b,matrix(1,4,2),matrix(1,4,2)), "zero")
  expect_error(cj_summarize_signed(c(0,NA), 0.01), "finite")
  z <- cj_summarize_signed(c(-0.1,0,0.1), 0.02)
  expect_equal(unname(sum(z[c("p_positive","p_negative","p_negligible")])), 1)
})
