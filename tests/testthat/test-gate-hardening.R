# cj_gate(): method = "auto", by = batch, QC flags, cofactor rule. The simulator
# gives a known truth, so every claim is checked against it.

gh_sim <- function(profile, markers = "M", batches = 1L, shift = 0, scale = 0, seed = 1, ...) {
  cj_simulate_gate_cells(profile, markers, samples_per_batch = 3, cells_per_sample = 2500,
                         batches = batches, batch_shift = shift, batch_scale = scale, seed = seed, ...)
}
gh_gate <- function(s, markers = "M", ...) {
  cj_gate(s$cells, markers, s$samples, transform = "none", control = "unstim", ...)
}
gh_kappa <- function(call, truth) {
  po <- mean(call == truth)
  pe <- mean(call) * mean(truth) + (1 - mean(call)) * (1 - mean(truth))
  (po - pe) / (1 - pe)
}

test_that("defaults are unchanged and the result gains qc, without batch fields", {
  s <- gh_sim("clean")
  g <- gh_gate(s)
  expect_identical(g$method, "tail")
  ctrl <- s$cells$M[s$cells$condition == "unstim"]
  pk <- gate_peak(ctrl)
  expect_equal(g$thresholds[["M"]], pk[["mode"]] + qnorm(0.999) * pk[["spread"]])
  expect_null(g$by)
  expect_null(g$batch_thresholds)
  expect_s3_class(g$qc, "data.frame")
  expect_false(any(g$qc$cut_in_broad_tail, g$qc$high_control_positive))
})

test_that("tail ratio is near 1 for a clean peak and large for a heavy upper tail", {
  set.seed(1)
  expect_lt(abs(gate_tail_ratio(rnorm(20000)) - 1), 0.15)
  s <- gh_sim("very_skewed")
  expect_gt(gate_tail_ratio(s$cells$M[s$cells$condition == "unstim"]), 3)
})

test_that("auto keeps the tail cut for a clean peak and takes the quantile for a heavy upper tail", {
  clean <- gh_gate(gh_sim("clean"), method = "auto")
  expect_identical(unname(clean$method_used), "tail")
  expect_identical(clean$auto$choice, "tail")
  # skewed negative: some seeds sit near the 1.5 limit, so require the choice on most seeds
  pick <- vapply(1:6, function(i) gh_gate(gh_sim("skewed", seed = i), method = "auto")$auto$choice, "")
  expect_gte(sum(pick == "quantile"), 3)
  s <- gh_sim("skewed", seed = 2, pos_mu = 14)
  g_auto <- gh_gate(s, method = "auto"); g_tail <- gh_gate(s)
  expect_identical(g_auto$auto$tail_ratio > 1.5, g_auto$auto$choice == "quantile")
  expect_equal(g_auto$auto$n, sum(s$cells$condition == "unstim"))
})

test_that("auto does not use the control quantile when a positive population sits in the control", {
  g <- gh_gate(gh_sim("wide_background"), method = "auto")
  expect_identical(g$auto$choice, "tail")
  expect_match(g$auto$reason, "positive population")
  g2 <- gh_gate(gh_sim("very_skewed"), method = "auto")
  expect_identical(g2$auto$choice, "tail")
})

test_that("auto keeps the tail cut with too few control events and says so", {
  s <- gh_sim("skewed")
  s$cells <- s$cells[s$cells$condition == "stim" | ave(seq_len(nrow(s$cells)), s$cells$sample_id, FUN = seq_along) <= 200, ]
  g <- gh_gate(s, method = "auto")
  expect_identical(g$auto$choice, "tail")
  expect_match(g$auto$reason, "too few control events")
})

test_that("auto beats tail against the truth on a skewed negative peak", {
  k <- vapply(1:4, function(i) {
    s <- gh_sim("skewed", seed = 10 + i)
    tr <- s$cells$truth_M
    c(tail = gh_kappa(s$cells$M > gh_gate(s)$thresholds[["M"]], tr),
      auto = gh_kappa(s$cells$M > gh_gate(s, method = "auto")$thresholds[["M"]], tr))
  }, c(tail = 0, auto = 0))
  expect_gte(mean(k["auto", ]), mean(k["tail", ]))
})

test_that("by = batch finds each batch's own cut and tracks batch offsets", {
  s <- gh_sim("clean", batches = 3, shift = 2, seed = 4)
  g <- gh_gate(s, by = "batch")
  expect_identical(g$by, "batch")
  expect_identical(rownames(g$batch_thresholds), c("B1", "B2", "B3"))
  # the batch with the larger offset has the larger cut
  off <- tapply(s$cells$M[s$cells$condition == "unstim"], s$cells$batch[s$cells$condition == "unstim"], median)
  expect_identical(order(g$batch_thresholds[, "M"]), order(unname(off)))
  expect_equal(g$thresholds[["M"]], gh_gate(s)$thresholds[["M"]])
  k <- vapply(1:4, function(i) {
    s <- gh_sim("clean", batches = 3, shift = 3, seed = 20 + i)
    g <- gh_gate(s, by = "batch")
    tr <- s$cells$truth_M
    c(pool = gh_kappa(s$cells$M > g$thresholds[["M"]], tr),
      batch = gh_kappa(s$cells$M > g$batch_thresholds[s$cells$batch, "M"], tr))
  }, c(pool = 0, batch = 0))
  expect_gt(mean(k["batch", ]), 0.97)
  expect_gt(mean(k["batch", ]), mean(k["pool", ]) + 0.05)
  expect_gt(g$batch_n_control[["B1"]], 0)
})

test_that("cj_from_cells applies the per-batch cuts", {
  s <- gh_sim("clean", batches = 3, shift = 2, seed = 4)
  samp <- s$samples
  samp$visit <- factor("v1"); samp$parent_population <- "CD4"
  g <- gh_gate(s, by = "batch", method = "quantile")
  x <- cj_from_cells(s$cells, samples = samp, gate = g)
  call <- s$cells$M > g$batch_thresholds[s$cells$batch, "M"]
  by_sample <- tapply(call, s$cells$sample_id, sum)
  expect_equal(unname(x$counts[names(by_sample), 2L]), unname(as.integer(by_sample)))
  # sum(X - t) is relative to each event's own batch cut
  i <- s$cells$sample_id == rownames(x$counts)[3]
  expect_equal(sum(x$intensity$sx[3, , "M"]), sum(s$cells$M[i] - g$batch_thresholds[s$cells$batch[i], "M"]))
})

test_that("batches with too few control events warn and use the pooled cut", {
  s <- gh_sim("clean", batches = 3, shift = 1, seed = 5)
  keep <- !(s$cells$batch == "B3" & s$cells$condition == "unstim" &
              ave(seq_len(nrow(s$cells)), s$cells$sample_id, FUN = seq_along) > 100)
  s$cells <- s$cells[keep, ]
  expect_warning(g <- gh_gate(s, by = "batch"), "B3")
  expect_equal(g$batch_thresholds["B3", "M"], g$thresholds[["M"]])
  expect_match(g$notes, "B3")
})

test_that("by input is checked", {
  s <- gh_sim("clean")
  expect_error(gh_gate(s, by = "nope"), "column of samples")
  expect_error(gh_gate(s, method = "manual", thresholds = c(M = 1), by = "batch"), "not used with")
  s2 <- s; s2$samples$batch[1] <- NA
  expect_error(gh_gate(s2, by = "batch"), "missing values")
})

test_that("QC flags: silent on a clean marker, firing on heavy-tailed and wide-background markers", {
  clean <- gh_gate(gh_sim("clean"))$qc
  expect_false(clean$high_control_positive)
  expect_false(clean$cut_in_broad_tail)
  wide <- gh_gate(gh_sim("wide_background"))$qc
  expect_true(wide$high_control_positive)
  expect_gt(wide$pct_control_pos_tailrule, 1)
  sk <- gh_gate(gh_sim("very_skewed"))$qc
  expect_true(sk$high_control_positive)
  expect_true(sk$cut_in_broad_tail)
  # thresholds are overridable
  lax <- gh_gate(gh_sim("wide_background"), qc_thresholds = c(max_control_pos = 0.5))$qc
  expect_false(lax$high_control_positive)
  expect_error(gh_gate(gh_sim("clean"), qc_thresholds = c(bad = 1)), "qc_thresholds")
})

test_that("batch_cut_spread fires with batch offsets and stays quiet without", {
  shifted <- gh_gate(gh_sim("clean", batches = 4, shift = 2, seed = 6), by = "batch")$qc
  expect_true(shifted$batch_cut_spread_flag)
  quiet <- gh_gate(gh_sim("clean", batches = 4, shift = 0, seed = 6), by = "batch")$qc
  expect_false(quiet$batch_cut_spread_flag)
  expect_true(is.na(gh_gate(gh_sim("clean"))$qc$batch_cut_spread))
})

test_that("print and plot show the flags", {
  g <- gh_gate(gh_sim("very_skewed", batches = 2, shift = 1))
  out <- capture.output(print(g))
  expect_true(any(grepl("high_control_positive", out)))
  gb <- gh_gate(gh_sim("clean", batches = 2, shift = 2), by = "batch")
  expect_true(any(grepl("batch_cut_spread", capture.output(print(gb)))))
  p <- cj_plot_gate(gb)
  expect_s3_class(p, "ggplot")
  expect_true(any(grepl("batch_cut_spread", levels(p$data$marker))))
  expect_gte(length(p$layers), 3L)
})

test_that("cofactor rule is scale-free and ignores a positive tail", {
  set.seed(2)
  n <- 20000
  x <- c(rnorm(n, 0, 300), 10^runif(n / 5, log10(3000), log10(3e5)))
  cells <- data.frame(sample_id = "c1", A = x, B = 2 * x)
  samples <- data.frame(sample_id = "c1", condition = "unstim")
  cf <- cj_cofactor(cells, c(A = "A"), samples, control = "unstim")
  expect_equal(unname(cf), 3 * 300, tolerance = 0.2)
  cells$A <- 10 * cells$A
  expect_equal(unname(cj_cofactor(cells, c(A = "A"), samples, control = "unstim")), 10 * unname(cf), tolerance = 0.05)
  g <- cj_gate(cells, c(A = "A", B = "B"), samples, transform = "asinh", cofactor = "auto", control = "unstim")
  expect_match(g$transform$label, "cofactor per marker")
  expect_error(cj_gate(cells, c(A = "A"), samples, transform = "none", cofactor = "auto", control = "unstim"), "asinh")
})

test_that("a cofactor named by marker is applied per marker", {
  cells <- data.frame(sample_id = "c1", A = c(0, 100, 1000, 5), B = c(0, 100, 1000, 5))
  m <- gate_matrix(cells, c(A = "A", B = "B"), gate_transform("asinh", c(A = 10, B = 100)))
  expect_equal(m[, "A"], asinh(cells$A / 10))
  expect_equal(m[, "B"], asinh(cells$B / 100))
  expect_error(gate_transform("asinh", c(10, 100)), "named by marker")
})

test_that("simulator returns a truth column per marker and respects the batch design", {
  s <- cj_simulate_gate_cells(c("clean", "t"), c("A", "B"), batches = 2, samples_per_batch = 2,
                              cells_per_sample = 500, seed = 3)
  expect_setequal(names(s$cells), c("sample_id", "donor_id", "condition", "batch", "A", "B", "truth_A", "truth_B"))
  expect_equal(nrow(s$samples), 2 * 2 * 2)
  expect_lt(mean(s$cells$truth_A[s$cells$condition == "unstim"]), 0.01)
  expect_gt(mean(s$cells$truth_A[s$cells$condition == "stim"]), 0.1)
})
