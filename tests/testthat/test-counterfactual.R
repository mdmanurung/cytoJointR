test_that("IPF preserves stimulated margins and baseline association", {
  bits <- cj_bits(c("IFNG","TNF"))
  pu <- c(.80,.06,.09,.05)
  ps <- c(.42,.24,.12,.22)
  q <- cj_counterfactual(pu,ps,bits)
  expect_equal(drop(crossprod(bits,q)),drop(crossprod(bits,ps)),
               tolerance=1e-8)
  # In the 2x2 setting the odds ratio is maintained.
  odds <- function(p) p[1]*p[4]/(p[2]*p[3])
  expect_equal(odds(q),odds(pu),tolerance=1e-7)
  expect_equal(sum(q),1)
})

test_that("coordination is zero if dependence stays constant", {
  bits <- cj_bits(c("IFNG","TNF"))
  pu <- c(.80,.06,.09,.05)
  ps <- cj_counterfactual(pu,c(.36,.24,.18,.22),bits)
  sc <- cj_scores_from_probs(pu,ps,bits)
  expect_lt(abs(unname(sc["coordination"])), 1e-8)
})
