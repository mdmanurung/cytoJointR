eb_sim <- function(seed = 4, n = 24) {
  cj_simulate_counts(n_donors = n, markers = c("IFNG", "TNF", "IL2"), stimulations = "CMV", N = 10000, seed = seed,
                     truth = list(alpha = c(IFNG = -8.3, TNF = -8, IL2 = -7.7),
                                  gamma = rbind(stim_CMV = c(0.2, 0.4, 0.6)), theta0 = rep(0.6, 3),
                                  delta = matrix(0, 1, 3), tau = c(0.5, 0.2), Omega = diag(2), tau_u = 0))
}

test_that("C0-EB returns valid draws and is fast without Stan", {
  x <- eb_sim()$data
  t <- system.time(d <- cytoJointR:::cj_fit_dirichlet_eb(x, n_draws = 1000, seed = 1))[["elapsed"]]
  expect_s3_class(d, "cj_draws")
  expect_lt(t, 10)
  expect_equal(unname(apply(d$p, 1:2, sum)), matrix(1, 1000, nrow(x$samples)), tolerance = 1e-8)
  expect_true(all(d$p > 0))
  expect_equal(d$meta$model_id, "C0-EB")
  expect_true(all(d$meta$hyper$converged))
  # reproducible
  d2 <- cytoJointR:::cj_fit_dirichlet_eb(x, n_draws = 1000, seed = 1)
  expect_identical(d$p, d2$p)
})

test_that("C0-EB posterior means follow the counts when the signal is strong", {
  x <- eb_sim()$data
  d <- cytoJointR:::cj_fit_dirichlet_eb(x, n_draws = 2000, seed = 1)
  pos <- 1 - colMeans(d$p[, , 1])
  obs <- 1 - x$counts[, 1] / rowSums(x$counts)
  expect_gt(cor(pos, obs), 0.95)
})

test_that("C0-EB ranks weak responses better than C0 on a seeded rare-response simulation", {
  rho <- vapply(1:6, function(s) {
    sim <- eb_sim(seed = s)
    truth <- vapply(seq_len(nrow(sim$data$contrasts)), function(i) {
      ct <- sim$data$contrasts[i, ]
      sum(sim$p[ct$stim_sample_id, -1]) - sum(sim$p[ct$control_sample_id, -1])
    }, numeric(1))
    est <- function(d) {
      sc <- suppressMessages(cj_scores(d, EPS, ccs_gate = NULL))
      sc <- sc[sc$score == "RMS" & is.na(sc$marker), ]
      sc$estimate[match(sim$data$contrasts$contrast_id, sc$contrast_id)]
    }
    c(c0 = cor(est(cytoJointR:::cj_fit_dirichlet(sim$data, n_draw = 1000, seed = 1)), truth, method = "spearman"),
      eb = cor(est(cytoJointR:::cj_fit_dirichlet_eb(sim$data, n_draws = 1000, seed = 1)), truth, method = "spearman"))
  }, c(c0 = 0, eb = 0))
  expect_gt(mean(rho["eb", ]), mean(rho["c0", ]) + 0.05)
})

test_that("C0-EB handles several stimulations, visits and a missing control", {
  d0 <- cj_example_draws(n_donor = 6, n_draw = 5, markers = c("IFNG", "TNF", "IL2"), visits = c("v1", "v2"),
                         stimulations = c("CMV", "SEB"), missing_control = "D6_v1_CD4_SEB")
  x <- example_data(d0)
  d <- cytoJointR:::cj_fit_dirichlet_eb(x, n_draws = 200, seed = 1)
  expect_equal(dim(d$p)[2], nrow(x$samples))
  sc <- suppressMessages(cj_scores(d, EPS, ccs_gate = NULL))
  expect_true(any(sc$status == "unavailable")) # the contrast without a control stays unavailable
  expect_false(anyNA(sc$estimate[sc$status == "ok"]))
})

test_that("cj_compass_report accepts magnitude_model = 'c0eb'", {
  x <- eb_sim()$data
  ct <- x$contrasts[!is.na(x$contrasts$control_sample_id), ]
  comp <- data.frame(contrast_id = ct$contrast_id, donor_id = ct$donor_id, stimulation = ct$stimulation,
                     visit = ct$visit, parent_population = ct$parent_population, compass_FS = 0.1, compass_PFS = 0.1,
                     compass_response_degree1 = 0.2, compass_response_degree2 = 0.1, control_source = "paired",
                     control_fallback = FALSE, control_shift_check = NA_character_)
  rep <- cj_compass_report(x, magnitude_model = "c0eb", compass = comp, n_draws = 500, seed = 1)
  expect_equal(unique(rep$magnitude_model), "c0eb")
  expect_true(all(rep$rms_lower <= rep$rms & rep$rms <= rep$rms_upper))
})

two_stim <- function(seed = 5, n = 16) {
  cj_simulate_counts(n_donors = n, markers = c("IFNG", "TNF"), stimulations = c("CMV", "SEB"), N = 20000, seed = seed,
                     truth = list(alpha = c(IFNG = -6, TNF = -6),
                                  gamma = rbind(stim_CMV = c(1.6, 1.6), stim_SEB = c(0, 0)),
                                  theta0 = 0.5, delta = matrix(0, 2, 1), tau = c(0.4, 0.1, 0.1),
                                  Omega = diag(3), tau_u = 0))
}

test_that("each stimulation gets its own hyperparameters (strong CMV, null SEB)", {
  x <- two_stim()$data
  d <- cytoJointR:::cj_fit_dirichlet_eb(x, n_draws = 1000, seed = 1)
  h <- d$meta$hyper
  mean_cmv <- h$mean[h$stimulation %in% "CMV"]; mean_seb <- h$mean[h$stimulation %in% "SEB"]
  expect_gt(mean_cmv - mean_seb, 1) # fold change about 5 against 1
  expect_lt(abs(mean_seb), 0.5)
  sc <- suppressMessages(cj_scores(d, EPS, ccs_gate = NULL))
  sc <- sc[sc$score == "RMS" & is.na(sc$marker) & sc$status == "ok", ]
  seb <- sc[grepl("SEB", sc$contrast_id), ]
  cmv <- sc[grepl("CMV", sc$contrast_id), ]
  expect_gte(mean(seb$lower <= 0 & seb$upper >= 0), 0.75) # the null stimulation's interval covers 0
  expect_true(all(cmv$lower > 0))
})

test_that("extra control replicates share one baseline draw and no sample is left at 1/K", {
  x <- eb_sim()$data
  # a second control sample for donor 1: same counts as the first, new id
  ctrl1 <- x$contrasts$control_sample_id[1]
  s1 <- x$samples[x$samples$sample_id == ctrl1, ]
  s1$sample_id <- "extra_control"
  cnt <- rbind(x$counts, extra_control = x$counts[ctrl1, ])
  samples <- rbind(x$samples, s1)
  x2 <- cytoJointR:::new_cj_data(cnt, samples, x$contrasts, x$bits)
  d <- cytoJointR:::cj_fit_dirichlet_eb(x2, n_draws = 500, seed = 1)
  i <- match(c(x$contrasts$control_sample_id[1], "extra_control"), x2$samples$sample_id)
  q <- 1 - d$p[, i, 1]
  expect_equal(q[, 1], q[, 2]) # one baseline per block
  K <- ncol(x$counts)
  expect_false(isTRUE(all.equal(unname(d$p[1, i[2], ]), rep(1 / K, K))))
  expect_lt(max(q[, 2]), 0.01) # a rare baseline, not the uniform pattern distribution
})

test_that("a partial epsilon is filled from the defaults, and ignored arguments warn", {
  x <- eb_sim()$data
  ct <- x$contrasts[!is.na(x$contrasts$control_sample_id), ]
  comp <- data.frame(contrast_id = ct$contrast_id, donor_id = ct$donor_id, stimulation = ct$stimulation,
                     visit = ct$visit, parent_population = ct$parent_population, compass_FS = 0.1, compass_PFS = 0.1,
                     compass_response_degree1 = 0.2, compass_response_degree2 = 0.1, control_source = "paired",
                     control_fallback = FALSE, control_shift_check = NA_character_)
  rep <- cj_compass_report(x, compass = comp, epsilon = c(RMS = 5e-4), n_draws = 300, seed = 1)
  expect_equal(attr(rep, "settings")$epsilon[["CPS"]], REPORT_EPSILON[["CPS"]])
  expect_equal(attr(rep, "settings")$epsilon[["RMS"]], 5e-4)
  expect_warning(cj_compass_report(x, magnitude_model = "c0", compass = comp, fit_args = list(sd_min = 0.1), n_draws = 300),
                 "ignored")
})

test_that("a parent population with few donors warns", {
  x <- eb_sim(n = 4)$data
  expect_warning(cytoJointR:::cj_fit_dirichlet_eb(x, n_draws = 200, seed = 1), "poorly determined")
})
