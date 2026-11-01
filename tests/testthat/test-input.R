mk <- function() {
  s <- expand.grid(sample_id = c("c", "s1", "s2"), A = 0:1, B = 0:1, stringsAsFactors = FALSE)
  s$count <- c(90, 80, 70, 5, 10, 15, 3, 6, 9, 2, 4, 6)
  samples <- data.frame(sample_id = c("c", "s1", "s2"), donor_id = "D1", visit = factor("v1"),
                        parent_population = "CD4", condition = c("control", "X", "Y"))
  list(counts = s, samples = samples)
}
build <- function(x, ...) cj_from_counts(x$counts, x$samples, c("A", "B"))

test_that("one control serves two stimulations and is stored once", {
  x <- build(mk())
  expect_s3_class(x, "cj_data")
  expect_equal(nrow(x$samples), 3L)
  expect_equal(x$contrasts$control_sample_id, c("c", "c"))
  expect_equal(unname(x$counts["c", ]), c(90L, 5L, 3L, 2L))
})

test_that("a missing control gives NA; two controls error", {
  m <- mk(); m$samples$donor_id[1] <- "D2"
  m$counts$sample_id[m$counts$sample_id == "c"] <- "c"
  expect_true(all(is.na(build(m)$contrasts$control_sample_id)))
  m <- mk(); m$samples$condition[2] <- "control"
  expect_error(build(m), "more than one control")
})

test_that("bad counts and metadata are rejected", {
  m <- mk(); m$counts$count[1] <- 0.5; expect_error(build(m), "nonnegative integers")
  m <- mk(); m$counts$count[1] <- -1; expect_error(build(m), "nonnegative integers")
  m <- mk(); m$counts$count[1] <- NA; expect_error(build(m), "nonnegative integers")
  m <- mk(); m$counts <- rbind(m$counts, m$counts[1, ]); expect_error(build(m), "duplicate pattern")
  m <- mk(); m$samples <- rbind(m$samples, transform(m$samples[1, ], donor_id = "D9"))
  expect_error(build(m), "varies within a sample_id")
  m <- mk(); m$samples$condition[2] <- "Control"; expect_error(build(m), "named like a control")
  m <- mk(); m$counts$pct <- 1; expect_error(build(m), "percentage")
  m <- mk(); m$samples$visit <- "v1"; expect_error(build(m), "factor")
})

test_that("N reconstructs the all-negative count only when the residual is valid", {
  m <- mk(); m$counts <- m$counts[!(m$counts$A == 0 & m$counts$B == 0), ]
  m$samples$N <- c(100, 100, 100)
  x <- build(m)
  expect_equal(unname(x$counts[, 1]), c(90L, 80L, 70L))
  expect_equal(x$provenance$all_negative_reconstructed, c("c", "s1", "s2"))
  m$samples$N <- c(100, 100, 5)
  expect_error(build(m), "negative all-negative")
  m$samples$N <- 0
  expect_error(build(m), "N > 0 but no|do not sum|negative")
})

test_that("all samples at N = 0 is an error", {
  m <- mk(); m$counts$count <- 0
  expect_error(build(m), "N = 0")
})

test_that("events, polyICSFlow and the pattern tabulator agree with counts", {
  cells <- data.frame(sample_id = rep(c("c", "s"), each = 4), A = c(0, 0, 2, 2, 2, 2, 0, 2),
                      B = c(0, 2, 0, 2, 2, 2, 0, 0))
  samples <- data.frame(sample_id = c("c", "s"), donor_id = "D1", visit = factor("v1"),
                        parent_population = "CD4", condition = c("control", "X"))
  x <- cj_from_events(cells, c("A", "B"), samples, thresholds = c(A = 1, B = 1))
  expect_equal(unname(x$counts["c", ]), c(1L, 1L, 1L, 1L))
  expect_equal(unname(x$counts["s", ]), c(1L, 1L, 0L, 2L))
  p <- data.frame(MarkerComb = factor(c("A-B-", "A+B+")), n = c(3L, 1L))
  y <- cj_from_polyicsflow(p, samples[1, ], c("A", "B"))
  expect_equal(unname(y$counts[1, ]), c(3L, 0L, 0L, 1L))
  expect_error(cj_from_polyicsflow(transform(p, MarkerComb = "B+A+"), samples[1, ], c("A", "B")),
               "MarkerComb")
  # Parent events 2..5; marker A positive at events 2,3; B at 3,4.
  expect_equal(cytoJointR:::tabulate_patterns(c(FALSE, TRUE, TRUE, TRUE, TRUE),
                 list(c(FALSE, TRUE, TRUE, FALSE, FALSE), c(FALSE, FALSE, TRUE, TRUE, FALSE))),
               c(1L, 1L, 1L, 1L))
})

test_that("cj_from_gatingset runs on a real GatingSet", {
  # A subprocess probe: flowWorkspace can segfault on load and would kill this session.
  ok <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
          c("-e", shQuote("quit(status = !requireNamespace('flowWorkspace', quietly = TRUE))")),
          stdout = FALSE, stderr = FALSE, timeout = 120)) == 0
  skip_if_not(ok, "flowWorkspace does not load in this environment")
  skip_if_not_installed("flowWorkspace")
  skip_if_not_installed("flowCore")
  set.seed(1)
  fr <- function() flowCore::flowFrame(cbind(A = rnorm(200), B = rnorm(200)))
  fs <- flowCore::flowSet(list(c = fr(), s = fr()))
  gs <- flowWorkspace::GatingSet(fs)
  g <- function(ch) flowCore::rectangleGate(setNames(list(c(0, Inf)), ch))
  flowWorkspace::gs_pop_add(gs, g("A"), name = "Apos")
  flowWorkspace::gs_pop_add(gs, g("B"), name = "Bpos")
  flowWorkspace::recompute(gs)
  samples <- data.frame(sample_id = c("c", "s"), donor_id = "D1", visit = factor("v1"),
                        parent_population = "root", condition = c("control", "X"))
  x <- cj_from_gatingset(gs, "root", c(A = "Apos", B = "Bpos"), samples)
  expect_equal(unname(rowSums(x$counts)), c(200, 200))
})
