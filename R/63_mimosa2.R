# Adapter to MIMOSA2 (github.com/RGLab/MIMOSA2): baseline-adjusted positivity calls
# for one cell subset. Native MIMOSA2model objects are returned untouched.

mimosa2_require <- function() {
  if (!requireNamespace("MIMOSA2", quietly = TRUE))
    stop("The MIMOSA2 package is required; install it with remotes::install_github(\"RGLab/MIMOSA2\")",
         call. = FALSE)
}

#' Baseline-adjusted positivity calls with MIMOSA2
#'
#' Runs one `MIMOSA2::MIMOSA2()` per stimulation, parent population and
#' follow-up visit. Each donor enters with four samples: stimulated and control
#' at the baseline visit, and stimulated and control at the follow-up visit.
#' MIMOSA2 then asks whether the stimulated-minus-control proportion rose from
#' baseline to follow-up. The native `MIMOSA2model` objects are returned in
#' `fits`, so `MIMOSA2::getZ()`, `MIMOSA2::getFDR()` and `MIMOSA2::Boxplot()`
#' work on them.
#'
#' @details
#' **One subset per fit.** MIMOSA2 models one count per sample, not the
#' Boolean patterns COMPASS uses. `positive = "any"` counts the cells positive
#' for at least one marker (all cells minus the all-negative pattern).
#' `positive = "<marker>"` counts the cells positive for that marker, whatever
#' the other markers do.
#'
#' **All four samples are needed.** A donor that lacks any of them for a fit is
#' left out and listed in `excluded` with the reason. MIMOSA2 does not relax the
#' paired-control requirement; `cj_compass(missing_control = )` is the option
#' for stimulated samples without a control.
#'
#' **What a call means.** `mimosa2_p_response` is the posterior mass of
#' MIMOSA2's four responder components. Those components are one-sided: a donor
#' whose `delta` (`ps1 - pu1 - ps0 + pu0`, follow-up minus baseline of the
#' stimulated-minus-control proportion) is negative cannot be a responder.
#' `mimosa2_fdr` is MIMOSA2's q-value by the method of Newton et al. (2004). Its
#' multiple-testing family is the donors of one fit, so q-values are not
#' comparable across fits. `mimosa2_response` is `mimosa2_fdr < fdr`.
#'
#' **Fit status.** MIMOSA2 is an EM fit of an 11-component mixture with 8
#' parameters, so it needs many donors per fit (its own example uses 100). It
#' stops silently when its optimiser fails. `groups$status` is `"fitted"`,
#' `"failed: <message>"`, or `"no estimate"` when the responder probabilities
#' are not finite.
#'
#' @param x A `cj_data` (or a `cj_fit`, whose data is used).
#' @param positive `"any"` (default) or one marker name; see *Details*.
#' @param baseline_visit The baseline visit. Defaults to the first level of
#'   `samples$visit`. Every other visit is a follow-up and gets its own fits.
#' @param fdr Threshold for `mimosa2_response`, on `mimosa2_fdr`. 0.01 is
#'   MIMOSA2's `getResponse()` default.
#' @param tol,maxit Passed to `MIMOSA2::MIMOSA2()`. `tol = 1e-5` is the value of
#'   MIMOSA2's README; its own default is `1e-8`.
#' @return A `cj_mimosa2` list with `fits` (named list of `MIMOSA2model`),
#'   `scores` (one row per donor and fit: the four positive counts `ns1`, `nu1`,
#'   `ns0`, `nu0`, the four totals `Ns1`, `Nu1`, `Ns0`, `Nu0`, `delta`,
#'   `mimosa2_p_response`, `mimosa2_fdr`, `mimosa2_response`), `groups` (one
#'   row per fit with its status) and `excluded` (donors left out, with the
#'   reason).
#' @references Finak G et al. (2014). Mixture models for single-cell assays with
#'   applications to vaccine studies. *Biostatistics* 15(1):87-101.
#' @export
#' @examplesIf requireNamespace("MIMOSA2", quietly = TRUE)
#' \donttest{
#' d <- cj_example_draws(n_donor = 30, markers = c("IFNG", "TNF"), visits = c("pre", "post"),
#'                       stimulations = "CMV", n_draw = 5)
#' counts <- t(vapply(seq_len(nrow(d$samples)), function(i)
#'   rmultinom(1, d$samples$N[i], colMeans(d$p[, i, ]))[, 1], numeric(4)))
#' storage.mode(counts) <- "integer"
#' dimnames(counts) <- list(d$samples$sample_id, dimnames(d$p)[[3]])
#' x <- cytoJointR:::new_cj_data(counts, d$samples, d$contrasts, d$bits)
#' mm <- cj_mimosa2(x)
#' mm$scores
#' }
cj_mimosa2 <- function(x, positive = "any", baseline_visit = NULL, fdr = 0.01, tol = 1e-5, maxit = 100L) {
  mimosa2_require()
  if (inherits(x, "cj_fit")) x <- x$data
  if (!inherits(x, "cj_data")) stop("x must be a cj_data or cj_fit", call. = FALSE)
  if (!(length(positive) == 1L && positive %in% c("any", x$markers)))
    stop("positive must be \"any\" or one of the markers: ", paste(x$markers, collapse = ", "), call. = FALSE)
  visits <- levels(x$samples$visit)
  if (is.null(baseline_visit)) baseline_visit <- visits[1]
  if (!(length(baseline_visit) == 1L && baseline_visit %in% visits))
    stop("baseline_visit must be one of the visit levels: ", paste(visits, collapse = ", "), call. = FALSE)
  followups <- setdiff(visits, baseline_visit)
  if (!length(followups)) stop("MIMOSA2 needs a follow-up visit besides the baseline", call. = FALSE)

  # Positive cells of every sample: all but the all-negative pattern, or one marker's patterns.
  cols <- if (positive == "any") rowSums(x$bits) > 0 else x$bits[, positive] == 1L
  pos <- rowSums(x$counts[, cols, drop = FALSE])
  ct <- x$contrasts
  ct$visit <- as.character(ct$visit)
  base <- ct[ct$visit == baseline_visit, ]
  key <- function(d) paste(d$donor_id, d$stimulation, d$parent_population, sep = " | ")

  fits <- list(); groups <- list(); tidy <- list(); excluded <- list()
  for (v in followups) {
    fu <- ct[ct$visit == v, ]
    for (g in unique(paste(fu$stimulation, fu$parent_population, sep = " | "))) {
      cf <- fu[paste(fu$stimulation, fu$parent_population, sep = " | ") == g, ]
      cb <- base[match(key(cf), key(base)), ]
      reason <- ifelse(is.na(cb$contrast_id), "no stimulated sample at baseline",
                ifelse(is.na(cb$control_sample_id), "no control sample at baseline",
                ifelse(is.na(cf$control_sample_id), "no control sample at follow-up", NA_character_)))
      gid <- paste(g, v, sep = " | ")
      if (any(!is.na(reason)))
        excluded[[gid]] <- data.frame(donor_id = cf$donor_id, stimulation = cf$stimulation,
                                      parent_population = cf$parent_population, visit = v,
                                      reason = reason)[!is.na(reason), ]
      keep <- is.na(reason)
      cf <- cf[keep, ]; cb <- cb[keep, ]
      ids <- list(ns1 = cf$stim_sample_id, nu1 = cf$control_sample_id,
                  ns0 = cb$stim_sample_id, nu0 = cb$control_sample_id)
      Ntot <- vapply(ids, function(i) as.numeric(x$samples$N[match(i, x$samples$sample_id)]),
                     numeric(nrow(cf)))
      Ntot <- matrix(Ntot, nrow = nrow(cf), dimnames = list(cf$donor_id, names(ids)))
      n <- lapply(ids, function(i) unname(pos[i]))
      grp <- data.frame(group = gid, stimulation = cf$stimulation[1], parent_population = cf$parent_population[1],
                        visit = v, baseline_visit = baseline_visit, n_donors = nrow(cf))
      if (nrow(cf) < 2L) {
        groups[[gid]] <- cbind(grp, status = "failed: fewer than 2 donors with all four samples")
        next
      }
      fit <- tryCatch(MIMOSA2::MIMOSA2(Ntot = Ntot, ns1 = n$ns1, nu1 = n$nu1, ns0 = n$ns0, nu0 = n$nu0,
                                       tol = tol, maxit = maxit),
                      error = function(e) conditionMessage(e))
      if (is.character(fit)) {
        groups[[gid]] <- cbind(grp, status = paste("failed:", fit))
        next
      }
      p_resp <- MIMOSA2::getZ(fit)[, "P.R"]
      ok <- all(is.finite(p_resp))
      groups[[gid]] <- cbind(grp, status = if (ok) "fitted" else "no estimate")
      fits[[gid]] <- fit
      q <- if (ok) MIMOSA2::getFDR(fit) else NA_real_
      tidy[[gid]] <- data.frame(
        donor_id = cf$donor_id, stimulation = cf$stimulation, parent_population = cf$parent_population,
        visit = v, baseline_visit = baseline_visit, positive = positive,
        ns1 = n$ns1, nu1 = n$nu1, ns0 = n$ns0, nu0 = n$nu0,
        Ns1 = Ntot[, "ns1"], Nu1 = Ntot[, "nu1"], Ns0 = Ntot[, "ns0"], Nu0 = Ntot[, "nu0"],
        delta = n$ns1 / Ntot[, "ns1"] - n$nu1 / Ntot[, "nu1"] - n$ns0 / Ntot[, "ns0"] + n$nu0 / Ntot[, "nu0"],
        mimosa2_p_response = if (ok) unname(p_resp) else NA_real_,
        mimosa2_fdr = q, mimosa2_response = q < fdr, row.names = NULL)
    }
  }
  bind <- function(l, empty) if (length(l)) do.call(rbind, unname(l)) else empty
  excluded <- bind(excluded, data.frame(donor_id = character(), stimulation = character(),
                                        parent_population = character(), visit = character(),
                                        reason = character()))
  if (nrow(excluded))
    message(nrow(excluded), " donor-fit(s) left out of MIMOSA2 (all four samples are needed); see $excluded")
  structure(list(fits = fits, scores = bind(tidy, NULL), groups = bind(groups, NULL), excluded = excluded),
            class = "cj_mimosa2")
}

#' @export
print.cj_mimosa2 <- function(x, ...) {
  cat("cj_mimosa2:", length(x$fits), "fit(s),", NROW(x$scores), "donor row(s),",
      sum(x$scores$mimosa2_response, na.rm = TRUE), "responder call(s);", nrow(x$excluded), "left out\n")
  invisible(x)
}
