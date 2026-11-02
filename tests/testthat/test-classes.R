toy_parts <- function() {
  bits <- cj_bits(c("IFNG", "TNF"))
  visit <- factor(c("v1", "v1", "v1"), levels = "v1")
  samples <- data.frame(
    sample_id = c("D1_ctrl", "D1_cmv", "D1_seb"), donor_id = "D1", visit = visit,
    parent_population = "CD4", condition = c("control", "CMV", "SEB"))
  counts <- matrix(c(900L, 50L, 30L, 20L,
                     800L, 90L, 60L, 50L,
                     500L, 200L, 150L, 150L), nrow = 3, byrow = TRUE,
                   dimnames = list(samples$sample_id, cytoJointR:::pattern_labels(bits)))
  samples$N <- as.numeric(rowSums(counts))
  # One control sample shared by two stimulations.
  contrasts <- data.frame(
    contrast_id = c("D1_v1_CD4_CMV", "D1_v1_CD4_SEB"), donor_id = "D1", visit = visit[1:2],
    parent_population = "CD4", stimulation = c("CMV", "SEB"),
    stim_sample_id = c("D1_cmv", "D1_seb"), control_sample_id = "D1_ctrl")
  list(counts = counts, samples = samples, contrasts = contrasts, bits = bits)
}

test_that("cj_data accepts a shared control and rejects reordered counts", {
  x <- toy_parts()
  expect_s3_class(do.call(cytoJointR:::new_cj_data, x), "cj_data")
  x$counts <- x$counts[3:1, ]
  expect_error(do.call(cytoJointR:::new_cj_data, x), "rownames")
})

test_that("a contrast without a control sample is valid (net scores unavailable)", {
  x <- toy_parts()
  x$contrasts$control_sample_id[2] <- NA
  expect_s3_class(do.call(cytoJointR:::new_cj_data, x), "cj_data")
})

test_that("cj_draws rejects probabilities that do not sum to one", {
  x <- toy_parts()
  p <- array(0.25, c(2, 3, 4), dimnames = list(NULL, x$samples$sample_id,
                                               colnames(x$counts)))
  meta <- list(engine_id = "test", model_id = "M1", seed = 1L, n_draws = 2L,
               draw_kind = "in_sample_donor_conditional")
  expect_s3_class(cytoJointR:::new_cj_draws(p, x$samples, x$contrasts, x$bits, meta),
                  "cj_draws")
  p[1, 1, 1] <- 0.5
  expect_error(cytoJointR:::new_cj_draws(p, x$samples, x$contrasts, x$bits, meta),
               "sum to one")
})
