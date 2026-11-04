# Joint report: COMPASS at the baseline and follow-up visits beside the
# baseline-adjusted MIMOSA2 call, one row per follow-up contrast.

BASELINE_AGREEMENT_LEVELS <- c("both", "pre_existing", "compass_only", "mimosa2_only", "neither",
                               "not_assessed")

#' COMPASS beside the baseline-adjusted MIMOSA2 call
#'
#' One row per donor, stimulation, parent population and follow-up visit. Each
#' row has COMPASS at the baseline visit and at the follow-up visit (FS, PFS,
#' response probability and call) and the MIMOSA2 call of [cj_mimosa2()], which
#' asks whether the stimulated-minus-control signal rose from baseline.
#'
#' @details
#' The COMPASS call is `compass_response_degree1 > response_threshold`, as in
#' [cj_compass_report()]. `agreement` compares it at the follow-up visit with
#' `mimosa2_response`:
#' * `"both"`: both call a response.
#' * `"pre_existing"`: COMPASS calls a response at the follow-up visit and at
#'   baseline, and MIMOSA2 does not. The response was there before; it did not
#'   rise.
#' * `"compass_only"`: COMPASS calls a response at the follow-up visit only, and
#'   MIMOSA2 does not.
#' * `"mimosa2_only"`: MIMOSA2 calls a rise, and COMPASS calls no response at the
#'   follow-up visit.
#' * `"neither"`: neither calls a response.
#' * `"not_assessed"`: a call is missing, for example a donor that MIMOSA2 left
#'   out because a sample is missing (see `mimosa2$excluded`), or a COMPASS
#'   follow-up call without the baseline call needed to tell `"pre_existing"`
#'   from `"compass_only"`.
#'
#' COMPASS sees every Boolean subset; MIMOSA2 sees one count (`positive` in
#' [cj_mimosa2()]). A `"compass_only"` row can be a response in a subset that
#' the MIMOSA2 count dilutes.
#'
#' @param x A `cj_data` (or a `cj_fit`, whose data is used).
#' @param compass A [cj_compass()] result for `x`. `NULL` runs
#'   `cj_compass(x, ...)` with `compass_args`.
#' @param mimosa2 A [cj_mimosa2()] result for `x`. `NULL` runs
#'   `cj_mimosa2(x, ...)` with `mimosa2_args`.
#' @param response_threshold COMPASS response probability above which a donor
#'   is called, as in [cj_compass_report()].
#' @param compass_args,mimosa2_args Lists of arguments for the two runs when
#'   `compass` or `mimosa2` is `NULL`.
#' @return A data frame: the keys, `baseline_visit`, `compass_FS_baseline`,
#'   `compass_PFS_baseline`, `compass_response_baseline`,
#'   `compass_call_baseline`, the same four with `_followup`, MIMOSA2's
#'   `delta`, `mimosa2_p_response`, `mimosa2_fdr`, `mimosa2_response`, and
#'   `agreement` (a factor with the levels above).
#' @export
cj_baseline_report <- function(x, compass = NULL, mimosa2 = NULL, response_threshold = 0.5,
                               compass_args = list(), mimosa2_args = list()) {
  if (inherits(x, "cj_fit")) x <- x$data
  if (!inherits(x, "cj_data")) stop("x must be a cj_data or cj_fit", call. = FALSE)
  if (is.null(mimosa2)) mimosa2 <- do.call(cj_mimosa2, c(list(x), mimosa2_args))
  if (is.null(compass)) compass <- do.call(cj_compass, c(list(x), compass_args))
  baseline <- if (!is.null(mimosa2_args$baseline_visit)) mimosa2_args$baseline_visit else
    if (NROW(mimosa2$scores)) mimosa2$scores$baseline_visit[1] else levels(x$samples$visit)[1]

  keys <- c("donor_id", "stimulation", "parent_population")
  comp <- compass$scores
  comp$visit <- as.character(comp$visit)
  if (anyDuplicated(comp[c(keys, "visit")]))
    stop("compass$scores has more than one row per donor, stimulation, parent population and visit", call. = FALSE)
  pick <- function(visit_col, suffix) {
    out <- comp[c(keys, "visit", "compass_FS", "compass_PFS", "compass_response_degree1")]
    names(out) <- c(keys, visit_col, paste0(c("compass_FS", "compass_PFS", "compass_response"), suffix))
    out[[paste0("compass_call", suffix)]] <- out[[paste0("compass_response", suffix)]] > response_threshold
    out
  }

  # Rows: the follow-up contrasts of x; COMPASS at follow-up and baseline; MIMOSA2.
  rows <- x$contrasts[as.character(x$contrasts$visit) != baseline, c(keys, "visit")]
  rows$visit <- as.character(rows$visit)
  rows$baseline_visit <- baseline
  out <- merge(rows, pick("visit", "_followup"), by = c(keys, "visit"), all.x = TRUE, sort = FALSE)
  out <- merge(out, pick("baseline_visit", "_baseline"), by = c(keys, "baseline_visit"), all.x = TRUE, sort = FALSE)
  m2 <- mimosa2$scores
  if (NROW(m2)) {
    m2 <- m2[c(keys, "visit", "delta", "mimosa2_p_response", "mimosa2_fdr", "mimosa2_response")]
    if (anyDuplicated(m2[c(keys, "visit")]))
      stop("mimosa2$scores has more than one row per donor, stimulation, parent population and visit", call. = FALSE)
    out <- merge(out, m2, by = c(keys, "visit"), all.x = TRUE, sort = FALSE)
  } else {
    out[c("delta", "mimosa2_p_response", "mimosa2_fdr")] <- NA_real_
    out$mimosa2_response <- NA
  }
  stopifnot("one report row per follow-up contrast" = nrow(out) == nrow(rows))

  cf <- out$compass_call_followup; cb <- out$compass_call_baseline; mr <- out$mimosa2_response
  out$agreement <- factor(
    ifelse(is.na(cf) | is.na(mr), "not_assessed",
    ifelse(cf & mr, "both",
    ifelse(mr, "mimosa2_only",
    ifelse(!cf, "neither",
    ifelse(is.na(cb), "not_assessed",
    ifelse(cb, "pre_existing", "compass_only")))))),
    levels = BASELINE_AGREEMENT_LEVELS)
  cols <- c(keys, "visit", "baseline_visit",
            paste0(rep(c("compass_FS", "compass_PFS", "compass_response", "compass_call"), 2),
                   rep(c("_baseline", "_followup"), each = 4)),
            "delta", "mimosa2_p_response", "mimosa2_fdr", "mimosa2_response", "agreement")
  out <- out[order(out$stimulation, out$parent_population, out$visit, out$donor_id), cols]
  rownames(out) <- NULL
  out
}
