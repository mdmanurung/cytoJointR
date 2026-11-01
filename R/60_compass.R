# Adapter to the real COMPASS package (docs/adr/0001). Native COMPASSResult
# objects are returned untouched; COMPASS FS/PFS/response probabilities are
# COMPASS quantities, not cytoJointR sign probabilities.

# COMPASS default category filter (COMPASS::COMPASS), applied to stimulated counts.
COMPASS_DEFAULT_FILTER <- function(x) colSums(x > 5) > 2

# "IFNG&!TNF&IL2": COMPASS's own pattern-name format, markers in bits order.
compass_pattern_names <- function(bits) {
  apply(bits, 1L, function(b) paste0(ifelse(b == 1L, "", "!"), colnames(bits), collapse = "&"))
}

compass_require <- function() {
  if (!requireNamespace("COMPASS", quietly = TRUE))
    stop("The COMPASS package is required; install it to use cj_compass()", call. = FALSE)
}

#' Fit COMPASS to a cytoJointR data set
#'
#' Runs one `COMPASS::SimpleCOMPASS()` per stimulation, visit and parent
#' population, using each donor's stimulated and control sample as the
#' `n_s` and `n_u` rows. The native `COMPASSResult` objects are returned in
#' `fits`, so every COMPASS function (`COMPASS::scores()`,
#' `COMPASS::PosteriorDiff()`, `plot()`, ...) works on them.
#'
#' Contrasts without a control sample cannot enter COMPASS; they are listed in
#' `excluded`. `SimpleCOMPASS()` has no category filter, so the adapter applies
#' COMPASS's own default from `COMPASS::COMPASS()`,
#' `function(x) colSums(x > 5) > 2` on the stimulated counts, and moves removed
#' patterns' cells to the all-negative pattern as COMPASS does. What was removed
#' is reported in `filters`. Pass `category_filter = NULL` to keep every pattern.
#'
#' COMPASS FS, PFS and response probabilities have COMPASS's definitions and
#' range. They are not cytoJointR sign probabilities and should not be
#' compared with `p_positive` as if they were.
#'
#' @param x A `cj_data` (or a `cj_fit`, whose `data` is used).
#' @param category_filter COMPASS category filter; the default is COMPASS's own.
#' @param verbose Show COMPASS's console output.
#' @param ... Passed to `COMPASS::SimpleCOMPASS()`, e.g. `iterations`,
#'   `replications`, `seed`. COMPASS's defaults apply when omitted.
#' @return A `cj_compass` list with `fits` (named list of `COMPASSResult`),
#'   `scores` (tidy table with `compass_FS`, `compass_PFS`,
#'   `compass_response_degree1`, `compass_response_degree2`), `groups` (one
#'   row per COMPASS run with its status), `filters` (what the category filter
#'   removed per run) and `excluded` (contrasts left out, with the reason).
#' @export
#' @examplesIf requireNamespace("COMPASS", quietly = TRUE)
#' \donttest{
#' d <- cj_example_draws(n_donor = 4, markers = c("IFNG", "TNF"), visits = "v1",
#'                       stimulations = "CMV", n_draw = 5)
#' counts <- t(vapply(seq_len(nrow(d$samples)), function(i)
#'   rmultinom(1, d$samples$N[i], colMeans(d$p[, i, ]))[, 1], numeric(4)))
#' storage.mode(counts) <- "integer"
#' dimnames(counts) <- list(d$samples$sample_id, dimnames(d$p)[[3]])
#' x <- cytoJointR:::new_cj_data(counts, d$samples, d$contrasts, d$bits)
#' cc <- cj_compass(x, iterations = 200, replications = 2)
#' cc$scores
#' }
cj_compass <- function(x, category_filter = COMPASS_DEFAULT_FILTER, verbose = FALSE, ...) {
  compass_require()
  if (inherits(x, "cj_fit")) x <- x$data
  if (!inherits(x, "cj_data")) stop("x must be a cj_data or cj_fit", call. = FALSE)
  if (grepl("[&!+*?()\\[\\]{}|^$\\\\]", paste(x$markers, collapse = ""), perl = TRUE))
    stop("COMPASS parses marker names with '&', '!' and a regular expression; ",
         "rename markers to use only letters, digits, '_', '.' or '-'", call. = FALSE)
  ct <- x$contrasts
  none <- is.na(ct$control_sample_id)
  excluded <- ct[none, c("contrast_id", "donor_id", "visit", "parent_population", "stimulation")]
  if (nrow(excluded)) excluded$reason <- "no control sample"
  ct <- ct[!none, ]

  # COMPASS wants the all-negative pattern last; the rest by degree.
  degree <- rowSums(x$bits)
  ord <- c(which(degree > 0)[order(degree[degree > 0])], which(degree == 0))
  nm <- compass_pattern_names(x$bits)[ord]
  null_col <- length(nm)
  gid <- paste(ct$stimulation, ct$visit, ct$parent_population, sep = " | ")
  fits <- list(); groups <- list(); filters <- list(); tidy <- list()
  for (g in unique(gid)) {
    cg <- ct[gid == g, ]
    cg <- cg[order(cg$donor_id), ]
    n_s <- x$counts[cg$stim_sample_id, ord, drop = FALSE]
    n_u <- x$counts[cg$control_sample_id, ord, drop = FALSE]
    dimnames(n_s) <- dimnames(n_u) <- list(cg$donor_id, nm)
    # Same filter, same null-category bookkeeping as COMPASS::COMPASS().
    del <- integer()
    if (!is.null(category_filter)) {
      del <- setdiff(which(!match.fun(category_filter)(n_s)), null_col)
      if (length(del)) {
        n_s[, null_col] <- n_s[, null_col] + rowSums(n_s[, del, drop = FALSE])
        n_u[, null_col] <- n_u[, null_col] + rowSums(n_u[, del, drop = FALSE])
        n_s <- n_s[, -del, drop = FALSE]; n_u <- n_u[, -del, drop = FALSE]
      }
    }
    key <- cg[1, c("stimulation", "visit", "parent_population")]
    filters[[g]] <- cbind(group = g, key, n_patterns = length(nm), n_removed = length(del),
                          removed = paste(nm[del], collapse = "; "),
                          cells_moved_stimulated = sum(x$counts[cg$stim_sample_id, ord[del]]),
                          cells_moved_control = sum(x$counts[cg$control_sample_id, ord[del]]),
                          row.names = NULL)
    meta <- data.frame(donor_id = cg$donor_id, stimulation = cg$stimulation,
                       visit = cg$visit, parent_population = cg$parent_population)
    run <- function() COMPASS::SimpleCOMPASS(n_s, n_u, meta, "donor_id", ...)
    fit <- tryCatch({
      res <- NULL
      if (verbose) res <- run() else utils::capture.output(res <- suppressMessages(run()))
      res
    }, error = function(e) conditionMessage(e))
    ok <- !is.character(fit)
    groups[[g]] <- cbind(group = g, key, n_donors = nrow(cg),
                         status = if (ok) "fitted" else paste("failed:", fit), row.names = NULL)
    if (!ok) next
    fits[[g]] <- fit
    resp <- function(deg) if (length(x$markers) >= deg)
      suppressMessages(COMPASS::Response(fit, degree = deg)) else NA_real_
    tidy[[g]] <- data.frame(
      contrast_id = cg$contrast_id, donor_id = cg$donor_id, stimulation = cg$stimulation,
      visit = cg$visit, parent_population = cg$parent_population,
      compass_FS = unname(COMPASS::FunctionalityScore(fit)[cg$donor_id]),
      compass_PFS = unname(COMPASS::PolyfunctionalityScore(fit)[cg$donor_id]),
      compass_response_degree1 = unname(resp(1)[cg$donor_id]),
      compass_response_degree2 = unname(resp(2)[cg$donor_id]))
  }
  bind <- function(l) if (length(l)) do.call(rbind, unname(l))
  structure(list(fits = fits, scores = bind(tidy), groups = bind(groups),
                 filters = bind(filters), excluded = excluded,
                 note = paste("COMPASS FS, PFS and response probabilities are COMPASS",
                              "quantities, not cytoJointR sign probabilities.")),
            class = "cj_compass")
}

#' Draw COMPASS's own heatmap for one run
#'
#' Calls `plot()` on a native `COMPASSResult`, which draws COMPASS's mean-gamma
#' heatmap.
#'
#' @param x Output of [cj_compass()].
#' @param group Name of a run in `names(x$fits)`; default the first.
#' @param ... Passed to COMPASS's `plot()` method.
#' @return What COMPASS's `plot()` returns, invisibly.
#' @export
cj_plot_compass <- function(x, group = names(x$fits)[1], ...) {
  compass_require()
  if (!group %in% names(x$fits))
    stop("group must be one of: ", paste(names(x$fits), collapse = ", "), call. = FALSE)
  invisible(plot(x$fits[[group]], ...))
}
