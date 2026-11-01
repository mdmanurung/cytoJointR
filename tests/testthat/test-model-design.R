test_that("design has named columns, dropped empty cells and a parameter map", {
  sim <- cj_simulate_counts(n_donors = 4, markers = c("IFNG", "TNF", "IL2"),
                            stimulations = c("CMV", "SEB"), visits = c("v1", "v2"), seed = 1)
  d <- cj_design(sim$data)
  expect_equal(colnames(d$X), c("stim_CMV", "stim_SEB", "visit_v2",
                                "visit_v2:stim_CMV", "visit_v2:stim_SEB"))
  expect_equal(colnames(d$W), c("basal", "induced_CMV", "induced_SEB"))
  expect_equal(colnames(d$Bpair), c("IFNG:TNF", "IFNG:IL2", "TNF:IL2"))
  expect_equal(nrow(d$dv), 8)                                  # 4 donors x 2 visits
  expect_true(all(d$param_map$stan %in% d$param_map$stan[!duplicated(d$param_map$stan)]))
  expect_equal(d$param_map$label[d$param_map$stan == "gamma[4,2]"], "gamma[visit_v2:stim_CMV, TNF]")

  # a stimulation absent at the second visit loses its visit-by-stimulation column
  keep <- !(sim$data$samples$condition == "SEB" & sim$data$samples$visit == "v2")
  x <- sim$data
  x$samples <- x$samples[keep, ]; x$counts <- x$counts[keep, ]
  x$contrasts <- x$contrasts[x$contrasts$stim_sample_id %in% x$samples$sample_id, ]
  d2 <- cj_design(x)
  expect_equal(d2$dropped, "visit_v2:stim_SEB")
  expect_false("visit_v2:stim_SEB" %in% colnames(d2$X))
})

test_that("donor_visit is switched off for a single visit", {
  sim <- cj_simulate_counts(n_donors = 3, seed = 1)
  expect_message(d <- cj_design(sim$data, donor_visit = TRUE), "single visit")
  expect_equal(d$use_dv, 0L)
})

test_that("the R linear predictor matches a direct evaluation of the model", {
  sim <- cj_simulate_counts(n_donors = 3, markers = c("A", "B", "C"), donor_visit = TRUE,
                            visits = c("v1", "v2"), seed = 2)
  d <- sim$design; tr <- sim$truth
  eta <- cytoJointR:::m1_eta(d, tr)
  s <- 5L; k <- 8L                                            # all three markers positive
  de <- tr$tau * (tr$L_Omega %*% tr$z[, d$donor_of[s]])
  a <- tr$alpha + as.vector(d$X[s, ] %*% tr$gamma) +
    sum(d$W[s, ] * de) + tr$tau_u * tr$u_raw[d$dv_of[s]]
  th <- tr$theta0 + as.vector(d$Xs[s, ] %*% tr$delta)
  expect_equal(unname(eta[s, k]), sum(a) + sum(th), tolerance = 1e-12)
  expect_equal(unname(eta[, 1]), rep(0, nrow(eta)))           # eta_s0 = 0
})
