# Calibration of non-paired COMPASS calls (R/61_compass_calibrate.R). The map functions need no COMPASS.

test_that("the isotonic map is monotone and moves a liberal cut to the paired call", {
  set.seed(1)
  pop <- runif(200)
  paired <- as.numeric(pop > 0.7)          # paired call agrees only above 0.7
  m <- cytoJointR:::compass_cal_map(c(0.1, 0.6, 0.65, 0.75, 0.95), pop, paired, "isotonic")
  expect_false(is.unsorted(m))
  expect_true(all(m[1:3] < 0.5))
  expect_true(all(m[4:5] > 0.5))
})

test_that("the quantile map matches the paired quantile at the same rank", {
  pop <- seq(0, 1, length.out = 101)
  paired <- pop^3
  m <- cytoJointR:::compass_cal_map(c(0.25, 0.5, 1), pop, paired, "quantile")
  expect_equal(m, unname(quantile(paired, c(0.25, 0.5, 1), type = 1)) |> unname(), tolerance = 0.02)
  expect_false(is.unsorted(m))
})

test_that("the fallback map sends the fallback cut to 0.5 and keeps 0 and 1", {
  cut <- cytoJointR:::COMPASS_CAL_FALLBACK_CUT
  expect_equal(cytoJointR:::compass_cal_fallback(c(0, cut, 1)), c(0, 0.5, 1))
  expect_lt(cytoJointR:::compass_cal_fallback(0.8), 0.5)
})

test_that("compass_cal_auc ranks and is NA with one class", {
  expect_equal(cytoJointR:::compass_cal_auc(c(.1, .2, .8, .9), c(FALSE, FALSE, TRUE, TRUE)), 1)
  expect_equal(cytoJointR:::compass_cal_auc(c(.9, .8, .2, .1), c(FALSE, FALSE, TRUE, TRUE)), 0)
  expect_true(is.na(cytoJointR:::compass_cal_auc(c(.1, .2), c(TRUE, TRUE))))
})

test_that("compass_with_seed is reproducible and restores the caller's RNG", {
  set.seed(5); a <- runif(1); set.seed(5)
  x <- cytoJointR:::compass_with_seed(9, runif(3))
  expect_equal(x, cytoJointR:::compass_with_seed(9, runif(3)))
  expect_equal(runif(1), a)
})

cal_scores <- function(raw, src) {
  data.frame(compass_max_gamma = raw, control_source = src, stringsAsFactors = FALSE)
}

test_that("compass_calibrate_scores maps non-paired rows, keeps paired rows, and falls back", {
  set.seed(2)
  pop <- runif(40)
  ps <- data.frame(group = "g", cls = "population", gamma_pseudo = pop,
                   gamma_paired = as.numeric(pop > 0.7))
  sc <- cal_scores(c(0.9, 0.9, 0.6), c("paired", "population", "population"))
  r <- cytoJointR:::compass_calibrate_scores(sc, rep("g", 3), ps, "isotonic")
  expect_equal(r$scores$compass_gamma_cal[1], 0.9)            # paired row untouched
  expect_equal(r$scores$calibration, c("paired", "isotonic", "isotonic"))
  expect_equal(r$scores$compass_call, c(TRUE, TRUE, FALSE))
  expect_equal(r$info$n_pseudo, 40)
  # Calibration only tightens: a raw 0.4 stays uncalled even if the map would lift it.
  lift <- cytoJointR:::compass_calibrate_scores(cal_scores(0.4, "population"), "g",
            data.frame(group = "g", cls = "population", gamma_pseudo = c(rep(0.2, 20), rep(0.4, 20)),
                       gamma_paired = c(rep(0, 20), rep(1, 20))), "isotonic")
  expect_lte(lift$scores$compass_gamma_cal, 0.4)
  expect_false(lift$scores$compass_call)
  # Too few pseudo donors, or one side empty: the fallback cut applies.
  few <- cytoJointR:::compass_calibrate_scores(sc, rep("g", 3), ps[1:6, ], "isotonic")
  expect_equal(few$scores$calibration[2:3], c("fallback", "fallback"))
  one_side <- transform(ps, gamma_paired = 0.1)
  none <- cytoJointR:::compass_calibrate_scores(sc, rep("g", 3), one_side, "isotonic")
  expect_equal(none$scores$calibration[2], "fallback")
  expect_true(none$scores$compass_call[2])   # raw 0.9 is above the fallback cut
})

test_that("a pseudo-missing control is forced through the same resolver", {
  d <- cj_example_draws(n_donor = 4, n_draw = 5, markers = c("A", "B"), visits = c("v1", "v2"),
                        stimulations = "CMV")
  x <- cytoJointR:::new_cj_data(
    {m <- t(vapply(seq_len(nrow(d$samples)), function(i)
      rmultinom(1, d$samples$N[i], colMeans(d$p[, i, ]))[, 1], numeric(4)))
     storage.mode(m) <- "integer"; dimnames(m) <- list(d$samples$sample_id, dimnames(d$p)[[3]]); m},
    d$samples, d$contrasts, d$bits)
  id <- x$contrasts$contrast_id[x$contrasts$donor_id == "D1" & x$contrasts$visit == "v2"]
  r <- cytoJointR:::compass_resolve_controls(x, "donor", force_missing = id)$contrasts
  expect_equal(r$control_source[r$contrast_id == id], "donor_visit:v1")
  expect_equal(sum(r$control_source == "paired"), nrow(r) - 1L)
  p <- cytoJointR:::compass_resolve_controls(x, "population", force_missing = id)$contrasts
  expect_true(is.na(p$control_sample_id[p$contrast_id == id]))
  expect_equal(p$control_source[p$contrast_id == id], "population")
})

# 14 donors, one visit; the controls of D1 to D4 are missing.
cal_data <- function(n_donor = 14, n_miss = 4) {
  miss <- sprintf("D%d_v1_CD4_CMV", seq_len(n_miss))
  d <- cj_example_draws(n_donor = n_donor, n_draw = 5, markers = c("IFNG", "TNF", "IL2"),
                        visits = "v1", stimulations = "CMV", missing_control = miss)
  x <- example_data(d)
  keep <- !x$samples$sample_id %in% sub("CD4_CMV", "control", miss)
  cytoJointR:::new_cj_data(x$counts[keep, ], x$samples[keep, ], x$contrasts, x$bits)
}

test_that("missing_control_calibrate adds columns, leaves the default alone, and pseudo-fits paired donors", {
  skip_if_not_installed("COMPASS")
  x <- cal_data()
  base <- suppressMessages(cj_compass(x, missing_control = "population", iterations = 200,
                                      replications = 2, seed = 3))
  expect_equal(tail(names(base$scores), 1), "compass_max_gamma")   # appended last: positional users see no change
  expect_false(any(c("compass_call", "compass_gamma_cal") %in% names(base$scores)))
  expect_null(base$calibration)
  cc <- suppressWarnings(suppressMessages(cj_compass(
    x, missing_control = "population", missing_control_calibrate = TRUE, max_pseudo = 5,
    iterations = 200, replications = 2, seed = 3)))
  # Raw gamma is unchanged by calibration; the paired rows keep it.
  expect_equal(cc$scores$compass_max_gamma, base$scores$compass_max_gamma)
  pr <- cc$scores$control_source == "paired"
  expect_equal(cc$scores$compass_gamma_cal[pr], cc$scores$compass_max_gamma[pr])
  expect_equal(cc$scores$compass_call, cc$scores$compass_gamma_cal > 0.5)
  # Only paired donors are pseudo-removed, at most max_pseudo, and their paired gamma is the real fit's.
  ps <- cc$calibration$pseudo
  expect_equal(nrow(ps), 5)          # max_pseudo = 5 of 10 paired donors: one fold of 5
  expect_true(all(ps$fold == 1))
  expect_equal(cc$calibration$n_failed, 0)
  expect_true(all(ps$donor_id %in% cc$scores$donor_id[pr]))
  real <- cc$scores$compass_max_gamma[match(ps$donor_id, cc$scores$donor_id)]
  expect_equal(ps$gamma_paired, real)
  # 5 paired donors are below the minimum, so the four missing rows use the fallback cut.
  expect_equal(unname(cc$scores$calibration[!pr]), rep("fallback", 4))
  expect_equal(cc$calibration$table$method, "fallback")
})

test_that("missing_control_calibrate with 'drop' or nothing missing is a no-op", {
  skip_if_not_installed("COMPASS")
  x <- cal_data(n_miss = 0)
  cc <- suppressMessages(cj_compass(x, missing_control = "population", missing_control_calibrate = TRUE,
                                    iterations = 100, replications = 2, seed = 3))
  expect_true(all(cc$scores$calibration == "paired"))
  expect_equal(cc$scores$compass_gamma_cal, cc$scores$compass_max_gamma)
})

test_that("a tied pseudo-missing gamma gives the paired mean for both maps", {
  pop <- rep(0.3, 20); paired <- c(rep(0, 10), rep(1, 10))
  for (m in c("isotonic", "quantile"))
    expect_equal(cytoJointR:::compass_cal_map(c(0.2, 0.9), pop, paired, m), c(0.5, 0.5))
})

test_that("non-finite pseudo gammas are dropped before counting and mapping", {
  set.seed(3)
  pop <- runif(40)
  ps <- data.frame(group = "g", cls = "population", gamma_pseudo = pop, gamma_paired = as.numeric(pop > 0.7))
  ps$gamma_pseudo[1:3] <- NA
  ps$gamma_paired[4:5] <- NaN
  sc <- cal_scores(c(0.9, 0.6), c("population", "population"))
  r <- expect_no_error(cytoJointR:::compass_calibrate_scores(sc, rep("g", 2), ps, "isotonic"))
  expect_equal(r$info$n_pseudo, 35)
  expect_false(anyNA(r$scores$compass_gamma_cal))
  # Only NA pseudo rows left: fallback, not an error.
  ps$gamma_pseudo <- NA_real_
  r <- cytoJointR:::compass_calibrate_scores(sc, rep("g", 2), ps, "isotonic")
  expect_equal(unique(r$scores$calibration), "fallback")
})

test_that("two groups get their own maps and the AUC is flagged only where the map was used", {
  set.seed(4)
  pop <- runif(40)
  good <- data.frame(group = "g1", cls = "population", gamma_pseudo = pop, gamma_paired = as.numeric(pop > 0.7))
  bad <- data.frame(group = "g2", cls = "population", gamma_pseudo = pop, gamma_paired = as.numeric(pop < 0.7))
  few <- data.frame(group = "g3", cls = "population", gamma_pseudo = pop[1:6], gamma_paired = as.numeric(pop[1:6] < 0.7))
  sc <- cal_scores(c(0.9, 0.9, 0.9, 0.9), c("population", "population", "population", "paired"))
  gid <- c("g1", "g2", "g3", "g1")
  r <- cytoJointR:::compass_calibrate_scores(sc, gid, rbind(good, bad, few), "isotonic")
  expect_equal(r$scores$calibration, c("isotonic", "isotonic", "fallback", "paired"))
  expect_true(r$scores$compass_call[1])
  i <- r$info[match(c("g1", "g2", "g3"), r$info$group), ]
  expect_equal(i$used_map, c(TRUE, TRUE, FALSE))
  expect_lt(i$pseudo_auc[2], 0.8)
})

test_that("cj_compass warns on a poor pseudo AUC, reports failed folds and ignores calibrate under 'drop'", {
  skip_if_not_installed("COMPASS")
  x <- cal_data()
  set.seed(5)
  pop <- runif(20)
  fake <- function(...) {
    ps <- data.frame(group = "CMV | v1 | CD4", donor_id = sprintf("D%d", 5:24), cls = "population", fold = 1L,
                     gamma_pseudo = pop, gamma_paired = as.numeric(pop < 0.5))
    attr(ps, "n_failed") <- 2L
    ps
  }
  testthat::local_mocked_bindings(compass_pseudo_fits = fake, .package = "cytoJointR")
  expect_warning(
    expect_message(
      cc <- cj_compass(x, missing_control = "population", missing_control_calibrate = TRUE,
                       iterations = 100, replications = 2, seed = 3),
      "2 pseudo-missing refit"),
    "does not rank donors")
  expect_equal(cc$calibration$n_failed, 2)
  expect_true(all(cc$calibration$table$used_map))
  expect_warning(cd <- cj_compass(x, missing_control = "drop", missing_control_calibrate = TRUE,
                                  iterations = 100, replications = 2, seed = 3), "ignored")
  expect_null(cd$calibration)
  expect_false("compass_call" %in% names(cd$scores))
})

test_that("cj_compass leaves .Random.seed as it was, also when it was absent", {
  skip_if_not_installed("COMPASS")
  x <- cal_data()
  run <- function() suppressWarnings(suppressMessages(cj_compass(
    x, missing_control = "population", missing_control_calibrate = TRUE, max_pseudo = 5,
    iterations = 100, replications = 2, seed = 3)))
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(11)
  before <- get(".Random.seed", envir = globalenv())
  run()
  expect_identical(get(".Random.seed", envir = globalenv()), before)
  rm(".Random.seed", envir = globalenv())
  run()
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

# 30 donors, 3 markers, one stimulation; half respond with 3% of stimulated cells moved to A+B.
# The controls of D1 to D6 are missing.
sim_cal_data <- function(n_donor = 30, n_miss = 6, seed = 1) {
  set.seed(seed)
  markers <- c("A", "B", "C")
  bits <- cj_bits(markers)
  pat <- apply(bits, 1, \(b) paste(markers[b == 1L], collapse = "+"))
  p0 <- apply(bits, 1, \(b) prod(ifelse(b == 1L, 0.03, 0.97)))
  q <- as.numeric(pat == "A+B")
  ids <- sprintf("D%02d", seq_len(n_donor))
  resp <- ids %in% ids[seq(2, n_donor, by = 2)]
  grid <- expand.grid(donor_id = ids, condition = c("control", "CMV"), stringsAsFactors = FALSE)
  grid$sample_id <- paste(grid$donor_id, "v1", grid$condition, sep = "_")
  grid$N <- 30000L
  grid$counts <- lapply(seq_len(nrow(grid)), \(i) {
    w <- if (grid$condition[i] == "CMV" && grid$donor_id[i] %in% ids[resp]) 0.03 else 0
    rmultinom(1, grid$N[i], (1 - w) * p0 + w * q)[, 1]
  })
  gone <- paste0(ids[seq_len(n_miss)], "_v1_control")
  grid <- grid[!grid$sample_id %in% gone, ]
  counts <- do.call(rbind, grid$counts)
  dimnames(counts) <- list(grid$sample_id, cytoJointR:::pattern_labels(bits))
  storage.mode(counts) <- "integer"
  meta <- data.frame(sample_id = grid$sample_id, donor_id = grid$donor_id, visit = factor("v1"),
                     parent_population = "CD4", condition = grid$condition, N = grid$N)
  list(x = cytoJointR:::new_cj_data(counts, meta, cytoJointR:::derive_contrasts(meta), bits), resp = setNames(resp, ids))
}

test_that("the isotonic path runs end to end with 24 paired donors and keeps rows aligned", {
  skip_if_not_installed("COMPASS")
  d <- sim_cal_data()
  cc <- suppressWarnings(suppressMessages(cj_compass(
    d$x, missing_control = "population", missing_control_calibrate = TRUE, max_pseudo = 20,
    iterations = 400, replications = 2, seed = 3)))
  expect_equal(nrow(cc$scores), 30)
  expect_equal(sum(cc$scores$control_source != "paired"), 6)
  expect_equal(cc$calibration$table$method, "isotonic")
  expect_equal(nrow(cc$calibration$pseudo), 20)
  ps <- cc$calibration$pseudo
  # Pseudo donors are paired donors; their paired gamma is the real fit's value for that donor.
  f <- cc$fits[[1]]
  mg <- f$fit$mean_gamma
  real <- apply(mg[, -ncol(mg), drop = FALSE], 1, max)
  expect_equal(ps$gamma_paired, unname(real[ps$donor_id]))
  expect_equal(cc$scores$compass_max_gamma, unname(real[cc$scores$donor_id]))
  # Responders are enough that the pseudo set holds both calls, so the map is used.
  expect_true(any(ps$gamma_paired > 0.5) && any(ps$gamma_paired <= 0.5))
  np <- cc$scores$control_source != "paired"
  expect_true(all(cc$scores$compass_gamma_cal[np] <= cc$scores$compass_max_gamma[np]))
  expect_equal(cc$scores$compass_call, cc$scores$compass_gamma_cal > 0.5)
})

test_that("a user seed does not hijack the calibration arguments", {
  skip_if_not_installed("COMPASS")
  x <- cal_data()
  cc <- suppressWarnings(suppressMessages(cj_compass(
    x, missing_control = "population", missing_control_calibrate = TRUE, max_pseudo = 5, pseudo_reps = 2,
    iterations = 100, replications = 2, seed = 3)))
  expect_equal(sort(unique(cc$calibration$pseudo$rep)), 1:2)
  expect_equal(nrow(cc$calibration$pseudo), 10)
})
