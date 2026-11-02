test_that("net percent is signed, display column is the only clipped one", {
  d <- cj_describe(cj_example())$patterns
  k <- d[d$contrast_id == "D1_v1_CD4_CMV" & d$pattern == "IFNG- TNF-", ]
  expect_equal(k$net_percent, 800 / 1000 * 100 - 900 / 1000 * 100)
  expect_true(k$net_percent < 0)
  expect_equal(k$net_percent_display, 0)
  expect_true(k$display_only)
  expect_true(all(d$net_percent_display[!d$display_only] == d$net_percent[!d$display_only]))
})

test_that("exactly and at-least aggregates are consistent", {
  d <- cj_describe(cj_example())
  e <- d$exactly[d$exactly$contrast_id == "D2_v1_CD4_CMV", ]
  a <- d$at_least[d$at_least$contrast_id == "D2_v1_CD4_CMV", ]
  expect_equal(sum(e$net_percent), 0)  # both sides sum to 100 percent
  expect_equal(a$net_percent[a$degree == 1], -e$net_percent[e$degree == 0])
  expect_equal(d$patterns$functionality[1:4], c("none", "mono", "mono", "poly"))
})

test_that("a missing control gives NA net values and plots are ggplot objects", {
  x <- cj_example()
  x$contrasts$control_sample_id[1] <- NA
  d <- cj_describe(x)$patterns
  expect_true(all(is.na(d$net_percent[d$contrast_id == x$contrasts$contrast_id[1]])))
  expect_message(g <- cj_plot_degree_bars(x), "no control")
  expect_s3_class(g, "ggplot")
  expect_s3_class(cj_plot_pattern_heatmap(x, x$contrasts$contrast_id[2]), "ggplot")
  cells <- data.frame(sample_id = "s", A = c(0, 2, 2, 0, 3), B = c(0, 0, 2, 3, 3))
  expect_s3_class(cj_plot_pattern_scatter(cells, c("A", "B"), c(A = 1, B = 1)), "ggplot")
})

test_that("summary reports contrasts lacking a control sample", {
  x <- cj_example()
  x$contrasts$control_sample_id[1] <- NA
  expect_output(print(x), "without a control sample: 1")
})

test_that("cj_describe works on data with no contrasts", {
  x <- cj_example()
  keep <- x$samples$condition == "control"
  samples <- x$samples[keep, ]
  y <- cj_from_counts(
    data.frame(sample_id = rep(samples$sample_id, each = 4), IFNG = c(0, 1, 0, 1), TNF = c(0, 0, 1, 1),
               count = as.vector(t(x$counts[samples$sample_id, ]))),
    samples, c("IFNG", "TNF"))
  expect_equal(nrow(y$contrasts), 0)
  d <- cj_describe(y)
  expect_equal(nrow(d$patterns), 0)
  expect_equal(nrow(d$exactly), 0)
  expect_equal(nrow(d$at_least), 0)
})
