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

# Decide the control sample of every stimulated sample. Returns the contrasts to run
# (with control_source, control_fallback and a control_sample_id that is NA for
# "population") and the contrasts left out with the reason.
compass_resolve_controls <- function(x, missing_control) {
  ct <- x$contrasts
  none <- is.na(ct$control_sample_id)
  ct$control_source <- ifelse(none, NA_character_, "paired")
  ct$control_fallback <- FALSE
  if (missing_control != "drop") {
    ctrl <- x$samples[x$samples$condition == CONTROL_CONDITION, ]
    lev <- levels(x$samples$visit)
    for (i in which(none)) {
      cand <- ctrl[ctrl$donor_id == ct$donor_id[i] & ctrl$parent_population == ct$parent_population[i] &
                    ctrl$visit != ct$visit[i], ]
      if (missing_control == "donor" && nrow(cand)) {
        pos <- match(as.character(cand$visit), lev)
        # Nearest visit by level order; a tie goes to the earlier visit.
        pick <- cand[order(abs(pos - match(as.character(ct$visit[i]), lev)), pos)[1], ]
        ct$control_sample_id[i] <- pick$sample_id
        ct$control_source[i] <- paste0("donor_visit:", pick$visit)
      } else {
        ct$control_source[i] <- "population"
        ct$control_fallback[i] <- missing_control == "donor"
      }
    }
  }
  drop <- is.na(ct$control_source)
  # A run in which every donor lacks a control has no population control to borrow.
  gid <- paste(ct$stimulation, ct$visit, ct$parent_population, sep = " | ")
  nocontrol <- !drop & unname(tapply(is.na(ct$control_sample_id), gid, all)[gid])
  keep_cols <- c("contrast_id", "donor_id", "visit", "parent_population", "stimulation")
  excluded <- ct[drop | nocontrol, keep_cols]
  if (nrow(excluded))
    excluded$reason <- ifelse(nocontrol[drop | nocontrol], "no donor in this run has a control sample",
                              "no control sample")
  list(contrasts = ct[!drop & !nocontrol, ], excluded = excluded)
}

#' Fit COMPASS to a cytoJointR data set
#'
#' Runs one `COMPASS::SimpleCOMPASS()` per stimulation, visit and parent
#' population, using each donor's stimulated and control sample as the
#' `n_s` and `n_u` rows. The native `COMPASSResult` objects are returned in
#' `fits`, so every COMPASS function (`COMPASS::scores()`,
#' `COMPASS::PosteriorDiff()`, `plot()`, ...) works on them.
#'
#' By default (`missing_control = "drop"`) contrasts without a control sample
#' cannot enter COMPASS; they are listed in `excluded`. The other two settings
#' keep them (see below). `SimpleCOMPASS()` has no category filter, so the adapter applies
#' COMPASS's own default from `COMPASS::COMPASS()`,
#' `function(x) colSums(x > 5) > 2` on the stimulated counts, and moves removed
#' patterns' cells to the all-negative pattern as COMPASS does. What was removed
#' is reported in `filters`. Pass `category_filter = NULL` to keep every pattern.
#'
#' COMPASS FS, PFS and response probabilities have COMPASS's definitions and
#' range. `compass_response_degree*` is `COMPASS::Response()`: the mean, over
#' the cell subsets with at least that many markers, of the posterior
#' probability that the subset responds. They are not cytoJointR sign probabilities and should not be
#' compared with `p_positive` as if they were.
#'
#' @section Missing control samples:
#' COMPASS integrates out each donor's control proportions, so an all-zero
#' `n_u` row adds nothing to the likelihood and that donor's stimulated counts
#' are compared with the population control distribution that the other donors
#' define. `missing_control` chooses what a stimulated sample without a control
#' of its own visit gets:
#' * `"drop"`: nothing; the contrast is excluded (default, unchanged).
#' * `"donor"`: the same donor's control from the nearest other visit with the
#'   same parent population. Nearest is by order of the visit levels; a tie goes
#'   to the earlier visit. A donor with no control at all falls back to
#'   `"population"`, flagged by `control_fallback = TRUE`.
#' * `"population"`: an all-zero control row.
#'
#' Every `scores` row has `control_source`: `"paired"`, `"donor_visit:<visit>"`
#' or `"population"`. A `"population"` row is scored against the other donors'
#' controls, so it is only as good as the assumption that this donor's
#' unstimulated pattern distribution resembles theirs (`docs/adr/0003`). A
#' borrowed control carries any real visit effect on the background. A run in
#' which no donor has a control is skipped and listed in `excluded`.
#' `COMPASS::PosteriorDiff()` on a `"population"` row is relative to that
#' population control. The category filter reads only stimulated counts, so
#' zero rows do not change which patterns it removes.
#'
#' @param x A `cj_data` (or a `cj_fit`, whose `data` is used).
#' @param category_filter COMPASS category filter; the default is COMPASS's own.
#' @param verbose Show COMPASS's console output.
#' @param missing_control What to do with a stimulated sample that has no control
#'   sample at its visit: `"drop"`, `"donor"` or `"population"`. See *Missing
#'   control samples*.
#' @param ... Passed to `COMPASS::SimpleCOMPASS()`, e.g. `iterations`,
#'   `replications`, `seed`. COMPASS's defaults apply when omitted.
#' @return A `cj_compass` list with `fits` (named list of `COMPASSResult`),
#'   `scores` (tidy table with `compass_FS`, `compass_PFS`,
#'   `compass_response_degree1`, `compass_response_degree2`, `control_source`,
#'   `control_fallback`), `groups` (one
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
cj_compass <- function(x, category_filter = COMPASS_DEFAULT_FILTER, verbose = FALSE,
                       missing_control = c("drop", "donor", "population"), ...) {
  compass_require()
  missing_control <- match.arg(missing_control)
  if (inherits(x, "cj_fit")) x <- x$data
  if (!inherits(x, "cj_data")) stop("x must be a cj_data or cj_fit", call. = FALSE)
  if (grepl("[&!+*?()\\[\\]{}|^$\\\\]", paste(x$markers, collapse = ""), perl = TRUE))
    stop("COMPASS parses marker names with '&', '!' and a regular expression; ",
         "rename markers to use only letters, digits, '_', '.' or '-'", call. = FALSE)
  ct <- compass_resolve_controls(x, missing_control)
  excluded <- ct$excluded
  ct <- ct$contrasts
  n_alt <- sum(ct$control_source != "paired")
  if (n_alt)
    message(n_alt, " of ", nrow(ct), " contrast(s) ran without a paired control: ",
            sum(startsWith(ct$control_source, "donor_visit")), " with the donor's control from another visit, ",
            sum(ct$control_source == "population"), " against the population control (",
            sum(ct$control_fallback), " of them a fallback from \"donor\")")

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
    # A sample without a control gets an all-zero row: COMPASS then uses the population control.
    n_u <- x$counts[match(cg$control_sample_id, rownames(x$counts)), ord, drop = FALSE]
    n_u[is.na(cg$control_sample_id), ] <- 0L
    n_u_raw <- n_u
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
                          cells_moved_control = sum(n_u_raw[, del]),
                          n_zero_control = sum(is.na(cg$control_sample_id)),
                          row.names = NULL)
    meta <- data.frame(donor_id = cg$donor_id, stimulation = cg$stimulation,
                       visit = cg$visit, parent_population = cg$parent_population,
                       control_source = cg$control_source)
    # SimpleCOMPASS() skips the posterior that PosteriorDiff() and friends read;
    # COMPASS() computes it with the same unexported helper.
    run <- function() {
      f <- COMPASS::SimpleCOMPASS(n_s, n_u, meta, "donor_id", ...)
      f$fit$posterior <- utils::getFromNamespace("compute_posterior", "COMPASS")(f)
      f
    }
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
    # Response() returns a one-column matrix (rows = donors), so index its column
    # by donor; indexing the matrix itself by name gives NA.
    resp <- function(deg) if (length(x$markers) >= deg)
      suppressMessages(COMPASS::Response(fit, degree = deg))[cg$donor_id, 1] else
        rep(NA_real_, nrow(cg))
    tidy[[g]] <- data.frame(
      contrast_id = cg$contrast_id, donor_id = cg$donor_id, stimulation = cg$stimulation,
      visit = cg$visit, parent_population = cg$parent_population,
      compass_FS = unname(COMPASS::FunctionalityScore(fit)[cg$donor_id]),
      compass_PFS = unname(COMPASS::PolyfunctionalityScore(fit)[cg$donor_id]),
      compass_response_degree1 = unname(resp(1)),
      compass_response_degree2 = unname(resp(2)),
      control_source = cg$control_source, control_fallback = cg$control_fallback)
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
