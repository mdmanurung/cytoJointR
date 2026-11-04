# Hand-built data: 40 donors, visits pre and post, one stimulation. Positive cells sit in
# the IFNG-only pattern. The first `n_resp` donors respond at post; the rest do not.
mimosa2_data <- function(n_donor = 40, n_resp = 20, N = 20000L, seed = 3) {
  set.seed(seed)
  bits <- cj_bits(c("IFNG", "TNF"))
  lab <- pattern_labels(bits)
  ifng_only <- which(bits[, "IFNG"] == 1L & bits[, "TNF"] == 0L)
  neg <- which(rowSums(bits) == 0L)
  donors <- sprintf("D%02d", seq_len(n_donor))
  samples <- expand.grid(condition = c("control", "CMV"), visit = c("pre", "post"), donor_id = donors,
                         stringsAsFactors = FALSE)
  samples$visit <- factor(samples$visit, levels = c("pre", "post"))
  samples$parent_population <- "CD4"
  samples$sample_id <- paste(samples$donor_id, samples$visit, samples$condition, sep = "_")
  samples$N <- N
  resp <- samples$donor_id %in% donors[seq_len(n_resp)]
  p <- ifelse(samples$condition == "control", 0.001,
              ifelse(samples$visit == "post" & resp, 0.006, 0.002))
  counts <- matrix(0L, nrow(samples), length(lab), dimnames = list(samples$sample_id, lab))
  counts[, ifng_only] <- as.integer(rbinom(nrow(samples), N, p))
  counts[, neg] <- N - counts[, ifng_only]
  stim <- samples[samples$condition == "CMV", ]
  contrasts <- data.frame(contrast_id = paste(stim$donor_id, stim$visit, "CD4_CMV", sep = "_"),
                          donor_id = stim$donor_id, visit = stim$visit, parent_population = "CD4",
                          stimulation = "CMV", stim_sample_id = stim$sample_id,
                          control_sample_id = sub("_CMV$", "_control", stim$sample_id))
  samples <- samples[c("sample_id", "donor_id", "visit", "parent_population", "condition", "N")]
  list(x = cytoJointR:::new_cj_data(counts, samples, contrasts, bits), responders = donors[seq_len(n_resp)])
}

test_that("cj_mimosa2 passes each donor's four samples to MIMOSA2 by name", {
  skip_if_not_installed("MIMOSA2")
  b <- mimosa2_data()
  mm <- cj_mimosa2(b$x)
  expect_s3_class(mm, "cj_mimosa2")
  expect_s3_class(mm$fits[[1]], "MIMOSA2model")
  expect_identical(mm$groups$status, "fitted")
  s <- mm$scores[mm$scores$donor_id == "D07", ]
  ifng <- b$x$counts[, "IFNG+ TNF-"]
  expect_equal(c(s$ns1, s$nu1, s$ns0, s$nu0),
               unname(ifng[c("D07_post_CMV", "D07_post_control", "D07_pre_CMV", "D07_pre_control")]))
  expect_equal(c(s$Ns1, s$Nu1, s$Ns0, s$Nu0), rep(20000, 4))
  expect_equal(s$delta, (s$ns1 - s$nu1 - s$ns0 + s$nu0) / 20000)
})

test_that("cj_mimosa2 separates simulated responders from non-responders", {
  skip_if_not_installed("MIMOSA2")
  b <- mimosa2_data()
  s <- cj_mimosa2(b$x)$scores
  truth <- s$donor_id %in% b$responders
  # A clear effect (0.002 to 0.006 at 20,000 cells): every responder outranks every non-responder.
  expect_gt(min(s$mimosa2_p_response[truth]), max(s$mimosa2_p_response[!truth]))
  expect_false(any(s$mimosa2_response[!truth]))
})

test_that("positive = <marker> counts that marker's cells; bad arguments fail", {
  skip_if_not_installed("MIMOSA2")
  b <- mimosa2_data()
  s <- cj_mimosa2(b$x, positive = "TNF")$scores
  expect_true(all(s$ns1 == 0 & s$nu1 == 0))
  expect_error(cj_mimosa2(b$x, positive = "IL2"), "positive must be")
  expect_error(cj_mimosa2(b$x, baseline_visit = "v9"), "baseline_visit must be")
})

test_that("a donor without all four samples is left out with the reason", {
  skip_if_not_installed("MIMOSA2")
  x <- mimosa2_data()$x
  x$contrasts$control_sample_id[x$contrasts$contrast_id == "D03_pre_CD4_CMV"] <- NA
  x$contrasts$control_sample_id[x$contrasts$contrast_id == "D05_post_CD4_CMV"] <- NA
  expect_message(mm <- cj_mimosa2(x), "2 donor-fit")
  expect_setequal(mm$excluded$donor_id, c("D03", "D05"))
  expect_identical(mm$excluded$reason[mm$excluded$donor_id == "D03"], "no control sample at baseline")
  expect_identical(mm$excluded$reason[mm$excluded$donor_id == "D05"], "no control sample at follow-up")
  expect_false(any(c("D03", "D05") %in% mm$scores$donor_id))
  expect_identical(mm$groups$n_donors, 38L)
})
