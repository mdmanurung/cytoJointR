fwd_files <- function() system.file("extdata", c("CytoTrol_CytoTrol_1.fcs", "CytoTrol_CytoTrol_2.fcs"),
                                    package = "flowWorkspaceData")

test_that("cj_read_fcs reads compensated events by description or channel name", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspaceData")
  f <- fwd_files()
  cells <- cj_read_fcs(f, c(CD38 = "CD38 APC", DR = "V545-A"))
  expect_named(cells, c("sample_id", "CD38", "DR"))
  expect_setequal(unique(cells$sample_id), c("CytoTrol_CytoTrol_1", "CytoTrol_CytoTrol_2"))
  # Same values as flowCore by hand: read, compensate with the file's SPILL, take the channel.
  fr <- flowCore::read.FCS(f[1], transformation = FALSE, truncate_max_range = FALSE)
  ref <- flowCore::exprs(flowCore::compensate(fr, flowCore::spillover(fr)$SPILL))
  one <- cells[cells$sample_id == "CytoTrol_CytoTrol_1", ]
  expect_equal(one$CD38, unname(ref[, "R660-A"]))
  expect_equal(one$DR, unname(ref[, "V545-A"]))
  expect_identical(attr(cells, "compensation")$keyword, c("SPILL", "SPILL"))
  raw <- cj_read_fcs(f[1], c(CD38 = "CD38 APC"), compensate = FALSE)
  expect_equal(raw$CD38, unname(flowCore::exprs(fr)[, "R660-A"]))
  expect_true(is.na(attr(raw, "compensation")$keyword))
})

test_that("cj_read_fcs fails on an unknown or ambiguous marker and on bad ids", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspaceData")
  f <- fwd_files()
  expect_error(cj_read_fcs(f, "IFNG"), "matches 0 channels")
  expect_error(cj_read_fcs(f, "CD38"), "matches 0 channels") # descriptions match exactly, not by prefix
  expect_error(cj_read_fcs(f, "CD38 APC", sample_id = c("a", "a")), "one unique id per file")
})

test_that("cj_read_fcs output feeds cj_from_cells", {
  skip_if_not_installed("flowCore")
  skip_if_not_installed("flowWorkspaceData")
  cells <- cj_read_fcs(fwd_files(), c(CD38 = "CD38 APC", DR = "HLA-DR V500"))
  samples <- data.frame(sample_id = c("CytoTrol_CytoTrol_1", "CytoTrol_CytoTrol_2"), donor_id = "D1",
                        visit = factor("v1"), parent_population = "all", condition = c("control", "stim"))
  x <- cj_from_cells(cells, c("CD38", "DR"), samples, transform = "asinh", cofactor = 150, method = "quantile")
  expect_s3_class(x, "cj_data")
  expect_equal(unname(rowSums(x$counts)), as.numeric(table(cells$sample_id)[x$samples$sample_id]))
})

test_that("cj_cells_from_gatingset returns a population's events by marker description", {
  # A subprocess probe: flowWorkspace can segfault on load and would kill this session.
  ok <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
          c("-e", shQuote("quit(status = !requireNamespace('flowWorkspace', quietly = TRUE))")),
          stdout = FALSE, stderr = FALSE, timeout = 120)) == 0
  skip_if_not(ok, "flowWorkspace does not load in this environment")
  skip_if_not_installed("flowWorkspaceData")
  gs <- flowWorkspace::load_gs(system.file("extdata", "gs_manual", package = "flowWorkspaceData"))
  cells <- cj_cells_from_gatingset(gs, "CD4", c(CD38 = "CD38 APC", DR = "HLA-DR V500"))
  expect_named(cells, c("sample_id", "CD38", "DR"))
  expect_identical(nrow(cells), as.integer(flowWorkspace::gh_pop_get_count(gs[[1]], "CD4")))
  ref <- flowCore::exprs(flowWorkspace::gs_pop_get_data(gs, "CD4")[[1]])
  expect_equal(cells$CD38, unname(ref[, "<R660-A>"]))
})
