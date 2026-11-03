# Adoption API: gating, cj_from_cells, fluorescence scores, wrappers.
# Brute-force checks recompute every score from the raw cells of one sample.

MK <- c("IFNG", "TNF", "IL2")

sim_in <- function(donors = 8, cells = 600, seed = 3, populations = "CD4") {
  sim <- cj_simulate_cells(donors = donors, cells_per_sample = cells, seed = seed)
  samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
  samples$visit <- factor("v1")
  samples$parent_population <- populations
  list(cells = sim$cells, samples = samples)
}

build <- function(transform = "none", cofactor = NULL, method = "quantile", ...) {
  s <- sim_in(...)
  x <- cj_from_cells(s$cells, MK, s$samples, transform = transform, cofactor = cofactor,
                     control = "unstim", method = method)
  c(s, list(x = x))
}

recode_control_for_test <- function(samples) {
  samples$condition[samples$condition == "unstim"] <- "control"
  samples
}

# Transformed values of one sample's cells, as a matrix with marker columns.
cell_X <- function(s, id, tr = identity) tr(as.matrix(s$cells[s$cells$sample_id == id, MK]))

test_that("quantile gate is the pooled control quantile; stimulated cells play no role", {
  s <- sim_in()
  g <- cj_gate(s$cells, MK, s$samples, method = "quantile", transform = "none", control = "unstim")
  ctrl <- s$cells[s$cells$condition == "unstim", ]
  expect_equal(g$thresholds, vapply(MK, function(m) unname(quantile(ctrl[[m]], 0.999)), 1))
  expect_equal(g$n_control, nrow(ctrl))
  expect_s3_class(cj_plot_gate(g), "ggplot")
  s$cells[s$cells$condition == "stim", MK] <- s$cells[s$cells$condition == "stim", MK] + 5
  g2 <- cj_gate(s$cells, MK, s$samples, method = "quantile", transform = "none", control = "unstim")
  expect_identical(g2$thresholds, g$thresholds)
})

test_that("tail gate recovers the simulated gate at 0 (misclassification under 0.5%)", {
  s <- sim_in(donors = 12, cells = 1500, seed = 1)
  g <- cj_gate(s$cells, MK, s$samples, transform = "none", control = "unstim")
  expect_identical(g$method, "tail")
  X <- as.matrix(s$cells[MK])
  miss <- vapply(MK, function(m) mean((X[, m] > g$thresholds[[m]]) != (X[, m] > 0)), 1)
  expect_lt(max(miss), 0.005)
  expect_true(all(g$thresholds > -0.5 & g$thresholds < 1.3))
})

test_that("mindensity needs openCyto: absent gives a warning and the quantile", {
  s <- sim_in()
  local_mocked_bindings(has_pkg = function(pkg) FALSE)
  expect_warning(g <- cj_gate(s$cells, MK, s$samples, method = "mindensity", transform = "none",
                              control = "unstim"), "openCyto")
  expect_true(all(g$method_used == "quantile"))
  expect_match(g$notes, "quantile used")
})

test_that("mindensity runs through openCyto where it loads", {
  skip_if_not(identical(Sys.getenv("CJ_TEST_OPENCYTO"), "true"), "set CJ_TEST_OPENCYTO=true to load openCyto")
  skip_if_not_installed("openCyto")
  skip_if_not_installed("flowCore")
  s <- sim_in(donors = 12, cells = 1500, seed = 1)
  g <- cj_gate(s$cells, MK, s$samples, method = "mindensity", transform = "none", control = "unstim")
  expect_true(all(g$method_used == "mindensity"))
  expect_true(all(is.finite(g$thresholds)))
})

test_that("manual thresholds are kept, reordered by marker name, and a transform is required", {
  s <- sim_in()
  g <- cj_gate(s$cells, MK, s$samples, thresholds = c(IL2 = 0.5, IFNG = 0.1, TNF = 0.2),
               transform = "none", control = "unstim")
  expect_identical(g$method, "manual")
  expect_identical(g$thresholds, c(IFNG = 0.1, TNF = 0.2, IL2 = 0.5))
  expect_error(cj_gate(s$cells, MK, s$samples, thresholds = c(IFNG = 0), transform = "none"), "one value per marker")
  expect_error(cj_gate(s$cells, MK, s$samples, thresholds = c(IFNG = 0, TNF = 0, IL2 = 0)), "transform is required")
  expect_error(cj_gate(s$cells, MK, s$samples, transform = "asinh"), "[Tt]ransform = \"none\"|cofactor", ignore.case = TRUE)
  expect_error(cj_gate(s$cells, MK, s$samples, method = "quantile", thresholds = c(IFNG = 0, TNF = 0, IL2 = 0),
                       transform = "none"), "only with method")
  # markers map names to columns
  cl <- setNames(s$cells, sub("^IFNG$", "BV786-A", names(s$cells)))
  g2 <- cj_gate(cl, c(IFNG = "BV786-A", TNF = "TNF", IL2 = "IL2"), s$samples, thresholds = c(IFNG = 0, TNF = 0, IL2 = 0),
                transform = "none", control = "unstim")
  expect_identical(names(g2$thresholds), MK)
})

test_that("the transform is applied before thresholding, and the gate fixes it", {
  s <- sim_in()
  g <- cj_gate(s$cells, MK, s$samples, method = "quantile", transform = "asinh", cofactor = 5, control = "unstim")
  ctrl <- s$cells[s$cells$condition == "unstim", MK]
  expect_equal(unname(g$thresholds[["TNF"]]), unname(quantile(asinh(ctrl$TNF / 5), 0.999)))
  expect_error(cj_from_cells(s$cells, MK, s$samples, gate = g, transform = "none"), "fixes the transform")
  x <- cj_from_cells(s$cells, MK, s$samples, gate = g)
  expect_equal(x$events$X[, "IL2"], asinh(s$cells$IL2 / 5))
  expect_identical(x$intensity$transform, "asinh(x / 5)")
  # non-finite transformed values stop the run
  bad <- s$cells; bad$IL2[1] <- NA
  expect_error(cj_gate(bad, MK, s$samples, method = "quantile", transform = "none", control = "unstim"), "not finite")
})

test_that("cj_from_cells reproduces cj_from_events counts and brute-force intensity sums", {
  b <- build(transform = "asinh", cofactor = 5)
  x <- b$x
  th <- x$intensity$thresholds
  ev_cells <- b$cells[c('sample_id', MK)]; ev_cells[MK] <- asinh(ev_cells[MK] / 5)
  ref <- cj_from_events(ev_cells, MK, recode_control_for_test(b$samples), thresholds = th)
  expect_identical(x$counts, ref$counts)
  # sx, ss for every sample, pattern, marker: loop over cells
  for (id in x$samples$sample_id) {
    X <- cell_X(b, id, function(m) asinh(m / 5))
    code <- drop((X > rep(th, each = nrow(X))) %*% c(1, 2, 4)) + 1
    for (k in 1:8) for (j in seq_along(MK)) {
      v <- X[code == k, j] - th[[j]]
      expect_equal(x$intensity$sx[id, k, j], sum(v), tolerance = 1e-9)
      expect_equal(x$intensity$ss[id, k, j], sum(v^2), tolerance = 1e-9)
    }
  }
  expect_error(cj_from_cells(b$cells[-1], MK, b$samples, transform = "none", control = "unstim"), "sample_id")
  expect_identical(x$samples$condition[grepl("unstim", x$samples$sample_id)][1], "control")
})

# Independent per-cell scores of one contrast: list of per-marker values.
brute <- function(b, tr, row, th) {
  XS <- cell_X(b, row$stim_sample_id, tr); XU <- cell_X(b, row$control_sample_id, tr)
  J <- length(MK)
  one <- function(X, j, min_deg, only = NULL) {
    pos <- sweep(X, 2L, th, `>`)
    sel <- pos[, j] & rowSums(pos) >= min_deg
    if (!is.null(only)) sel <- sel & rowSums(pos) == only
    sum((X[, j] - th[[j]])[sel]) / nrow(X)
  }
  lapply(seq_len(J), function(j) {
    deg_terms <- vapply(1:J, function(d) {
      pos <- function(X) sweep(X, 2L, th, `>`)
      f <- function(X) { P <- pos(X); sel <- P[, j] & rowSums(P) == d; sum((X[, j] - th[[j]])[sel]) / nrow(X) }
      f(XS) - f(XU)
    }, 1)
    list(nicr = one(XS, j, 1) - one(XU, j, 1), poly = one(XS, j, 2) - one(XU, j, 2),
         b_S = one(XS, j, 1), b_U = one(XU, j, 1), pb_S = one(XS, j, 2), pb_U = one(XU, j, 2),
         mono = deg_terms[1], by_degree = deg_terms)
  })
}

test_that("NICR, Poly-NICR, burdens and IW-FS/PFS equal a per-cell computation", {
  for (tr_name in c("none", "asinh")) {
    b <- if (tr_name == "none") build() else build(transform = "asinh", cofactor = 5)
    tr <- if (tr_name == "none") identity else function(m) asinh(m / 5)
    th <- b$x$intensity$thresholds
    fs <- suppressMessages(cj_fluorescence_scores(b$x))
    ct <- b$x$contrasts[!is.na(b$x$contrasts$control_sample_id), ]
    expect_equal(nrow(fs), nrow(ct) * 3L)
    J <- 3L
    for (r in seq_len(nrow(ct))) {
      br <- brute(b, tr, ct[r, ], th)
      for (j in seq_len(J)) {
        row <- fs[fs$contrast_id == ct$contrast_id[r] & fs$marker == MK[j], ]
        expect_equal(row$nicr, br[[j]]$nicr, tolerance = 1e-10)
        expect_equal(row$poly_nicr, br[[j]]$poly, tolerance = 1e-10)
        expect_equal(row$b_stim, br[[j]]$b_S, tolerance = 1e-10)
        expect_equal(row$b_ctrl, br[[j]]$b_U, tolerance = 1e-10)
        expect_equal(row$poly_b_stim, br[[j]]$pb_S, tolerance = 1e-10)
        expect_equal(row$iw_fs, br[[j]]$nicr / (2^J - 1), tolerance = 1e-10)
        # PFS weights d / choose(J, d) / J * 2 / (J + 1), applied to the degree-d net burden
        w <- (1:J) / choose(J, 1:J) / J * 2 / (J + 1)
        expect_equal(row$iw_pfs, sum(w * br[[j]]$by_degree), tolerance = 1e-10)
        # NICR is the sum of its degree parts
        expect_equal(row$nicr, sum(br[[j]]$by_degree), tolerance = 1e-10)
      }
    }
  }
})

test_that("burden of the polyfunctional cells never exceeds the total; Poly-NICR <= NICR when degree-1 NICR >= 0", {
  b <- build()
  fs <- suppressMessages(cj_fluorescence_scores(b$x))
  expect_true(all(fs$poly_b_stim <= fs$b_stim + 1e-12))
  expect_true(all(fs$poly_b_ctrl <= fs$b_ctrl + 1e-12))
  th <- b$x$intensity$thresholds
  ct <- b$x$contrasts[!is.na(b$x$contrasts$control_sample_id), ]
  mono <- unlist(lapply(seq_len(nrow(ct)), function(r) vapply(brute(b, identity, ct[r, ], th), `[[`, 1, "mono")))
  expect_true(any(mono >= 0))
  keep <- mono >= 0
  expect_true(all(fs$poly_nicr[keep] <= fs$nicr[keep] + 1e-12))
  expect_equal(fs$nicr - fs$poly_nicr, mono, tolerance = 1e-10)
})

test_that("frequency + intensity + unmatched equals NICR exactly, and matches a per-cell split", {
  b <- build(donors = 8, cells = 2000)
  fs <- suppressMessages(cj_fluorescence_scores(b$x))
  expect_equal(fs$nicr_freq + fs$nicr_int + fs$nicr_unmatched, fs$nicr, tolerance = 1e-12)
  th <- b$x$intensity$thresholds
  ct <- b$x$contrasts[!is.na(b$x$contrasts$control_sample_id), ]
  r <- 1L; j <- 2L
  XS <- cell_X(b, ct$stim_sample_id[r]); XU <- cell_X(b, ct$control_sample_id[r])
  code <- function(X) drop((X > rep(th, each = nrow(X))) %*% c(1, 2, 4)) + 1
  cS <- code(XS); cU <- code(XU)
  frq <- int <- unm <- 0
  for (k in 1:8) if ((k - 1) %/% 2^(j - 1) %% 2 == 1) {
    nS <- sum(cS == k); nU <- sum(cU == k)
    pS <- nS / nrow(XS); pU <- nU / nrow(XU)
    IS <- if (nS) mean(XS[cS == k, j] - th[[j]]) else 0
    IU <- if (nU) mean(XU[cU == k, j] - th[[j]]) else 0
    if (nS && nU) { frq <- frq + (pS - pU) * (IS + IU) / 2; int <- int + (IS - IU) * (pS + pU) / 2
    } else unm <- unm + pS * IS - pU * IU
  }
  row <- fs[fs$contrast_id == ct$contrast_id[r] & fs$marker == MK[j], ]
  expect_equal(c(row$nicr_freq, row$nicr_int, row$nicr_unmatched), c(frq, int, unm), tolerance = 1e-10)
  # a pattern seen in one condition only lands in `unmatched`
  expect_true(any(abs(fs$nicr_unmatched) > 0))
})

test_that("FIS: weights sum to one; single-positive reference equals the per-cell shift; sparse support is unavailable", {
  b <- build(donors = 8, cells = 2000)
  bits <- b$x$bits
  w <- cytoJointR:::fis_weights(c(5, 3, 2, 0, 4, 0, 0, 1), bits, "pooled")
  expect_equal(colSums(w), c(IFNG = 1, TNF = 1, IL2 = 1), ignore_attr = TRUE)
  expect_true(all(w[bits[, "IFNG"] == 0L, 1] == 0))
  fs1 <- suppressMessages(cj_fluorescence_scores(b$x, reference = "single"))
  th <- b$x$intensity$thresholds
  ct <- b$x$contrasts[!is.na(b$x$contrasts$control_sample_id), ]
  r <- 1L; j <- 1L
  XS <- cell_X(b, ct$stim_sample_id[r]); XU <- cell_X(b, ct$control_sample_id[r])
  only <- function(X) { P <- sweep(X, 2L, th, `>`); X[P[, j] & rowSums(P) == 1L, j] - th[[j]] }
  row <- fs1[fs1$contrast_id == ct$contrast_id[r] & fs1$marker == MK[j], ]
  expect_equal(row$fis, mean(only(XS)) - mean(only(XU)), tolerance = 1e-10)
  expect_true(any(is.na(fs1$fis)) || any(!is.na(fs1$fis)))
  fsp <- suppressMessages(cj_fluorescence_scores(b$x, reference = "pooled", fis_min_n = 1e6))
  expect_true(all(is.na(fsp$fis)))
  expect_true(all(grepl("fewer than", fsp$fis_reason)))
  expect_error(cj_fluorescence_scores(b$x, reference = "bogus"), "reference must be")
})

test_that("bootstrap: reproducible, keeps the RNG stream, intervals bracket and need events", {
  b <- build(donors = 6, cells = 400)
  set.seed(5); a <- runif(1)
  set.seed(5)
  f1 <- suppressMessages(cj_fluorescence_scores(b$x, n_boot = 40, seed = 11))
  expect_equal(runif(1), a)
  f2 <- suppressMessages(cj_fluorescence_scores(b$x, n_boot = 40, seed = 11))
  expect_identical(f1, f2)
  expect_true(all(f1$nicr_lo <= f1$nicr_hi))
  expect_true(all(f1$nicr_lo <= f1$nicr_hi & f1$poly_nicr_lo <= f1$poly_nicr_hi, na.rm = TRUE))
  # the bootstrap is centred near the estimate and is not degenerate
  expect_gt(mean(f1$nicr_hi - f1$nicr_lo), 0)
  expect_gt(mean(f1$nicr >= f1$nicr_lo & f1$nicr <= f1$nicr_hi), 0.7)
  # resampled counts keep the sample size
  s <- cytoJointR:::boot_stats(1:50, rep(1:8, length.out = 50), matrix(rnorm(100), 50, 2), 8L)
  expect_equal(sum(s$n), 50)
  nokeep <- cj_from_cells(b$cells, MK, b$samples, transform = "none", control = "unstim",
                          method = "quantile", keep_events = FALSE)
  expect_error(cj_fluorescence_scores(nokeep, n_boot = 5), "needs events")
  expect_error(cj_fluorescence_scores(cj_example()), "cj_from_cells")
})

test_that("contrasts without a control are omitted with a message and listed", {
  b <- build(donors = 12, cells = 300, seed = 1)
  expect_message(fs <- cj_fluorescence_scores(b$x), "without a control")
  ex <- attr(fs, "excluded")
  expect_gt(nrow(ex), 0)
  expect_false(any(ex$contrast_id %in% fs$contrast_id))
})

test_that("cj_subset keeps the intensity summaries of the kept samples", {
  s <- sim_in(donors = 6, cells = 300)
  s2 <- s; s2$cells$sample_id <- paste0(s2$cells$sample_id, "_cd8")
  s2$samples$sample_id <- paste0(s2$samples$sample_id, "_cd8"); s2$samples$parent_population <- "CD8"
  cells <- rbind(s$cells, s2$cells); samples <- rbind(s$samples, s2$samples)
  x <- cj_from_cells(cells, MK, samples, transform = "none", control = "unstim", method = "quantile")
  sub <- cj_subset(x, "CD8")
  expect_identical(dimnames(sub$intensity$sx)[[1]], sub$samples$sample_id)
  expect_equal(sub$intensity$sx, x$intensity$sx[sub$samples$sample_id, , ])
  expect_identical(max(sub$events$sample), nrow(sub$samples))
  expect_equal(nrow(sub$events$X), sum(sub$counts))
  all_fs <- suppressMessages(cj_fluorescence_scores(x))
  expect_equal(suppressMessages(cj_fluorescence_scores(sub))$nicr,
               all_fs$nicr[all_fs$parent_population == "CD8"])
})

test_that("wrappers fail clearly when the optional package is absent", {
  b <- build()
  local_mocked_bindings(has_pkg = function(pkg) FALSE)
  expect_error(cj_run_compass(b$x), "needs the COMPASS")
  expect_error(cj_run_polyicsflow(b$x), "needs the polyICSFlow")
  expect_message(qs <- cj_quickstart(b$cells, MK, b$samples, gate = b$x$gate, run = c("compass", "polyicsflow")),
                 "not installed")
  expect_named(qs, c("gate", "data", "description", "fluorescence"))
})

test_that("cj_quickstart returns gate, data, description and fluorescence scores", {
  b <- build()
  qs <- cj_quickstart(b$cells, MK, b$samples, transform = "none", control = "unstim", method = "quantile",
                      run = character(), n_boot = 5)
  expect_s3_class(qs$gate, "cj_gate")
  expect_s3_class(qs$data, "cj_data")
  expect_true("nicr_lo" %in% names(qs$fluorescence))
  expect_error(cj_quickstart(b$cells, MK, b$samples, transform = "none", control = "unstim",
                             method = "quantile", run = character(), fit = TRUE), "epsilon")
})

test_that("polyICSFlow event table agrees with the pattern counts", {
  b <- build()
  ev <- cytoJointR:::polyics_events(b$x)
  expect_equal(nrow(ev$positivity), sum(b$x$counts))
  expect_equal(as.vector(table(factor(ev$sample, levels = seq_len(nrow(b$x$counts))))), unname(b$x$samples$N))
  expect_identical(cytoJointR:::polyics_labels(b$x$bits)[2], "IFNG+TNF-IL2-")
})

test_that("cj_run_compass runs end to end from cells", {
  skip_if_not_installed("COMPASS")
  b <- build(donors = 8, cells = 1500, seed = 3)
  tab <- suppressMessages(cj_run_compass(b$cells, MK, b$samples, transform = "none", control = "unstim",
                                         gate_args = list(method = "quantile"), iterations = 200, replications = 2,
                                         seed = 1))
  expect_true(all(c("compass_FS", "compass_PFS", "compass_response_degree2") %in% names(tab)))
  expect_s3_class(attr(tab, "compass"), "cj_compass")
  expect_equal(sort(tab$donor_id), sort(setdiff(b$x$contrasts$donor_id[!is.na(b$x$contrasts$control_sample_id)], character())))
  expect_true(all(tab$compass_FS >= 0 & tab$compass_FS <= 1))
  # the same table from the cj_data
  tab2 <- suppressMessages(cj_run_compass(b$x, iterations = 200, replications = 2, seed = 1))
  expect_equal(nrow(tab2), nrow(tab))
})

test_that("cj_run_polyicsflow runs end to end where polyICSFlow loads", {
  skip_if_not_installed("polyICSFlow")
  b <- build(donors = 6, cells = 1000)
  out <- suppressMessages(cj_run_polyicsflow(b$x))
  expect_s3_class(out, "cj_polyicsflow")
  expect_true(all(c("polyics_net_any", "polyics_net_poly") %in% names(out$net)))
  expect_true(all(out$net$polyics_net_any >= out$net$polyics_net_poly - 1e-9, na.rm = TRUE))
})
