# Old path: one cj_score_reference() call per draw (the scalar oracle).
old_engine <- function(d, i) {
  ct <- d$contrasts
  s <- match(ct$stim_sample_id[i], d$samples$sample_id)
  u <- match(ct$control_sample_id[i], d$samples$sample_id)
  nd <- dim(d$p)[1]
  cols <- c(GLOBAL_SCORES, paste0("marker_magnitude:", d$markers))
  G <- matrix(NA_real_, nd, length(cols), dimnames = list(NULL, cols))
  contrib <- array(NA_real_, c(nd, nrow(d$bits), 3L))
  status <- rep("ok", nd)
  for (k in seq_len(nd)) {
    p_u <- d$p[k, u, ]; p_s <- d$p[k, s, ]
    ref <- tryCatch(cj_score_reference(unname(p_u / sum(p_u)), unname(p_s / sum(p_s)), d$bits),
                    error = function(e) NULL)
    if (is.null(ref)) { status[k] <- "ipf_failed"; next }
    G[k, ] <- c(ref$global[GLOBAL_SCORES], ref$markers$magnitude)
    contrib[k, , ] <- as.matrix(ref$patterns[, CONTRIBUTION_SCORES])
  }
  list(G = G, contrib = contrib, status = status)
}

check_engine <- function(d) {
  eng <- cytoJointR:::compute_engine(d, keep_patterns = TRUE)
  n_checked <- 0L
  for (i in seq_along(eng)) {
    e <- eng[[i]]
    if (!e$available) { expect_true(is.na(d$contrasts$control_sample_id[i])); next }
    o <- old_engine(d, i)
    expect_identical(e$status, o$status)
    expect_equal(unname(e$G), unname(o$G), tolerance = 1e-10)
    expect_equal(unname(e$contrib), o$contrib, tolerance = 1e-10)
    n_checked <- n_checked + 1L
  }
  n_checked
}

test_that("vectorised engine equals the per-draw reference path", {
  expect_gt(check_engine(cj_example_draws(n_draw = 30)), 0)
  # Shared control per donor/visit across two stimulations, one unavailable contrast.
  d <- cj_example_draws(n_draw = 30, missing_control = "D2_v1_CD4_CMV", effect = -1)
  expect_gt(check_engine(d), 0)
  expect_true(any(!vapply(cytoJointR:::compute_engine(d, FALSE), `[[`, TRUE, "available")))
})

test_that("failed draws are flagged and kept, matching the reference", {
  d <- cj_example_draws(n_draw = 30)
  u <- match("D1_v1_control", d$samples$sample_id)
  s <- match("D1_v1_CMV", d$samples$sample_id)
  stopifnot(!is.na(u), !is.na(s))
  d$p[3, u, 2] <- 0                                   # zero baseline probability
  d$p[7, s, d$bits[, 1] == 0] <- 1e-300               # marker 1 margin is 1: not interior
  i <- match("D1_v1_CD4_CMV", d$contrasts$contrast_id)
  eng <- cytoJointR:::compute_engine(d, TRUE)[[i]]
  expect_equal(which(eng$status != "ok"), c(3L, 7L))
  expect_true(all(is.na(eng$G[c(3, 7), ])))
  check_engine(d)
  # The summary counts them rather than dropping them.
  sc <- cj_scores(d, EPS, max_failure = 0.5)
  expect_equal(sc$n_failed[sc$contrast_id == "D1_v1_CD4_CMV"][1], 2)
  expect_equal(sc$n_draws[sc$contrast_id == "D1_v1_CD4_CMV"][1], 30)
})

test_that("the engine is computed once per object", {
  d <- cj_example_draws(n_draw = 20)
  cytoJointR:::score_engine(d)
  first <- cytoJointR:::.engine_cache$val$eng
  cj_scores(d, EPS); cj_markers(d, 0.005)
  expect_identical(cytoJointR:::.engine_cache$val$eng, first)
  d2 <- cj_example_draws(n_draw = 20, seed = 2)
  expect_false(identical(cytoJointR:::score_engine(d2), first))
})
