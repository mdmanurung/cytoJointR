# Descriptive fluorescence-aware scores (NICR, Poly-NICR, FIS, intensity-weighted
# FS/PFS, frequency/intensity decomposition) as plug-in estimates from events.
# No model: p_ck = n_ck / N_c and I_ckj = mean(X_j - t_j) among the cells of
# pattern k in condition c, so p_ck * I_ckj = sum(X_j - t_j) / N_c exactly.

FLUOR_SCORES <- c("nicr", "poly_nicr", "fis", "iw_fs", "iw_pfs", "nicr_freq", "nicr_int")

# COMPASS PolyfunctionalityScore weight of a degree-d subset out of J markers.
compass_pfs_weight <- function(d, J) ifelse(d == 0, 0, d / choose(J, d) / J * 2 / (J + 1))

# Reference weights w[k, j] for FIS: fixed, common to both conditions, sum to 1
# over the marker-positive patterns of marker j. "pooled" uses the pooled counts
# of the stimulation group (all its stimulated and control samples); "single"
# puts all weight on the pattern where only marker j is positive. A matrix
# (patterns x markers, rows named like colnames(x$counts)) is used as given.
fis_weights <- function(n_pool, bits, reference) {
  J <- ncol(bits)
  if (is.matrix(reference)) return(reference)
  w <- matrix(0, nrow(bits), J, dimnames = list(NULL, colnames(bits)))
  for (j in seq_len(J)) {
    pos <- bits[, j] == 1L
    w[, j] <- if (reference == "single") as.numeric(pos & rowSums(bits) == 1L) else
      ifelse(pos, n_pool, 0)
    if (sum(w[, j]) > 0) w[, j] <- w[, j] / sum(w[, j])
  }
  w
}

# All scores for one contrast. n: pattern counts (K) and sx: sums of X - t
# (K x J) for the stimulated (S) and control (U) sample. Returns one row per marker.
fluor_contrast <- function(nS, sxS, nU, sxU, bits, w, fis_min_n) {
  K <- nrow(bits); J <- ncol(bits)
  d <- rowSums(bits)
  NS <- sum(nS); NU <- sum(nU)
  out <- data.frame(marker = colnames(bits), n_stim = NS, n_ctrl = NU, b_stim = NA_real_,
                    b_ctrl = NA_real_, nicr = NA_real_, poly_b_stim = NA_real_,
                    poly_b_ctrl = NA_real_, poly_nicr = NA_real_, nicr_freq = NA_real_,
                    nicr_int = NA_real_, nicr_unmatched = NA_real_, fis = NA_real_,
                    fis_reason = NA_character_, iw_fs = NA_real_, iw_pfs = NA_real_,
                    stringsAsFactors = FALSE)
  if (NS == 0 || NU == 0) {
    out$fis_reason <- "no events in a sample"
    return(out)
  }
  pS <- nS / NS; pU <- nU / NU
  wp <- compass_pfs_weight(d, J)
  for (j in seq_len(J)) {
    pos <- bits[, j] == 1L
    poly <- pos & d >= 2L
    bS <- sxS[, j] / NS  # p_Sk * I_Skj
    bU <- sxU[, j] / NU
    dB <- bS - bU
    out$b_stim[j] <- sum(bS[pos]); out$b_ctrl[j] <- sum(bU[pos])
    out$nicr[j] <- sum(dB[pos])
    out$poly_b_stim[j] <- sum(bS[poly]); out$poly_b_ctrl[j] <- sum(bU[poly])
    out$poly_nicr[j] <- sum(dB[poly])
    # Exact split: dB = dp (I_S + I_U) / 2 + dI (p_S + p_U) / 2 where both
    # conditions have cells; the rest (a pattern seen in one condition) is "unmatched".
    both <- pos & nS > 0 & nU > 0
    IS <- ifelse(nS > 0, sxS[, j] / pmax(nS, 1), 0)
    IU <- ifelse(nU > 0, sxU[, j] / pmax(nU, 1), 0)
    out$nicr_freq[j] <- sum(((pS - pU) * (IS + IU) / 2)[both])
    out$nicr_int[j] <- sum(((IS - IU) * (pS + pU) / 2)[both])
    out$nicr_unmatched[j] <- sum(dB[pos & !both])
    out$iw_fs[j] <- sum(dB[pos]) / (2^J - 1)
    out$iw_pfs[j] <- sum((dB * wp)[pos])
    use <- w[, j] > 0
    if (!any(use)) {
      out$fis_reason[j] <- "no reference pattern"
    } else if (any(nS[use] < fis_min_n | nU[use] < fis_min_n)) {
      lack <- which(use & (nS < fis_min_n | nU < fis_min_n))
      out$fis_reason[j] <- paste0("fewer than ", fis_min_n, " events in: ",
                                  paste(utils::head(rownames(sxS)[lack], 3), collapse = "; "))
    } else out$fis[j] <- sum(w[use, j] * (IS - IU)[use])
  }
  out
}

# Per-sample pattern counts and X - t sums from resampled events.
boot_stats <- function(rows, code, Xt, K) {
  i <- rows[sample.int(length(rows), replace = TRUE)]
  n <- tabulate(code[i], nbins = K)
  sx <- matrix(0, K, ncol(Xt))
  rs <- rowsum(Xt[i, , drop = FALSE], code[i], reorder = TRUE)
  sx[as.integer(rownames(rs)), ] <- rs
  list(n = n, sx = sx)
}

#' Fluorescence-aware net scores per donor, stimulation, visit and marker
#'
#' Descriptive, model-free scores that use the intensity of the positive cells as
#' well as their frequency. Everything is a plug-in estimate from the events of
#' one sample; there is no truncated-Gaussian model and no Stan fit. Needs a
#' `cj_data` from [cj_from_cells()].
#'
#' @details
#' Notation. `X_j` is the transformed fluorescence of marker `j` and `t_j` its
#' gate threshold. For a condition `c` (stimulated `S` or control `U`) and a
#' pattern `k`: `p_ck = n_ck / N_c` is the share of parent events in the pattern
#' and `I_ckj = mean(X_j - t_j)` over the pattern's events. For a marker-positive
#' pattern `I_ckj > 0`. Then `p_ck * I_ckj = sum(X_j - t_j) / N_c`. All values
#' are in transformed fluorescence units per parent cell.
#'
#' * `b_stim`, `b_ctrl`: burden `B_cj = sum over patterns positive for j of
#'   p_ck I_ckj`, which equals the mean of `(X_j - t_j) * 1[X_j > t_j]` over all
#'   parent events.
#' * `nicr` = `B_Sj - B_Uj`, the net integrated cytokine response of marker `j`.
#'   Signed, never clipped. Markers are never summed: raw channels are on
#'   different scales.
#' * `poly_nicr`: the same sum over patterns with degree at least 2. It is not
#'   bounded by `nicr` in general: the degree-1 part of `nicr` can be negative.
#'   `poly_b_stim <= b_stim` and `poly_b_ctrl <= b_ctrl` always hold.
#' * `nicr_freq`, `nicr_int`, `nicr_unmatched`: exact split of `nicr`, per
#'   pattern, `dB_k = dp_k (I_S + I_U)/2 + dI_k (p_S + p_U)/2`, summed over
#'   marker-positive patterns that have events in both conditions. A pattern seen
#'   in one condition only has no conditional mean in the other; its whole change
#'   goes to `nicr_unmatched`. `nicr_freq + nicr_int + nicr_unmatched = nicr`
#'   exactly. It is an algebraic allocation, not a causal one.
#' * `fis`: fixed-reference intensity shift `sum_k w_kj (I_Skj - I_Ukj)`. The
#'   weights `w_kj` are nonnegative, sum to 1 over the patterns positive for `j`,
#'   and are the same in both conditions and every bootstrap replicate. FIS is
#'   `NA` with the cause in `fis_reason` unless every pattern with positive
#'   weight has at least `fis_min_n` events in both samples. Control samples
#'   rarely hold polyfunctional events, so with `reference = "pooled"` FIS is
#'   often unavailable; `reference = "single"` (weight on the pattern where only
#'   `j` is positive) is supported far more often.
#' * `iw_fs`, `iw_pfs` (experimental): intensity-weighted analogues of COMPASS
#'   FS and PFS, per marker. COMPASS weights the posterior probability that a
#'   cell subset responds, `gamma_k`, by `1/(2^J - 1)` for FS and by
#'   `(d_k / choose(J, d_k)) / J * 2 / (J + 1)` for PFS (`d_k` = number of
#'   markers in subset `k`). Here `gamma_k` is replaced by the net integrated
#'   intensity of marker `j` in subset `k`, `dB_kj = p_Sk I_Skj - p_Uk I_Ukj`,
#'   summed over subsets positive for `j`:
#'   `iw_fs_j = sum dB_kj / (2^J - 1)` and
#'   `iw_pfs_j = sum dB_kj * d_k / choose(J, d_k) / J * 2 / (J + 1)`.
#'   `iw_fs_j` is therefore `nicr_j / (2^J - 1)` (a rescaling). `iw_pfs_j` gives
#'   higher-degree subsets more weight. Both are signed, in fluorescence units
#'   per parent cell, not comparable in size with COMPASS FS/PFS, and have not
#'   been validated beyond their algebra.
#'
#' Contrasts without a control sample are omitted and listed in
#' `attr(, "excluded")`, with a message.
#'
#' @section Bootstrap:
#' With `n_boot > 0` each replicate resamples the events within every sample
#' with replacement (a control shared by several stimulations is resampled once)
#' and recomputes all scores. The reported interval is the percentile interval
#' of the replicates (`<score>_lo`, `<score>_hi`) for `nicr`, `poly_nicr`,
#' `fis`, `iw_fs`, `iw_pfs`, `nicr_freq` and `nicr_int`. It covers
#' within-sample cell sampling only, not donor-to-donor variation: donors are
#' the replicates for any group comparison. It needs `keep_events = TRUE` in
#' [cj_from_cells()]. Replicates where a score is `NA` are dropped for that
#' score.
#'
#' @param x A `cj_data` from [cj_from_cells()].
#' @param n_boot Number of bootstrap replicates; 0 for none.
#' @param level Interval level, for example 0.95.
#' @param seed Seed for the bootstrap; the caller's RNG stream is left unchanged.
#' @param reference FIS weights: `"pooled"`, `"single"`, or a matrix with one
#'   row per pattern (rownames `colnames(x$counts)`) and one column per marker.
#' @param fis_min_n Events needed in each sample for every reference pattern.
#' @return A data frame with one row per contrast and marker: `contrast_id`,
#'   `donor_id`, `visit`, `parent_population`, `stimulation`, `marker`, `n_stim`,
#'   `n_ctrl`, the columns above, and the interval columns. `attr(, "excluded")`
#'   lists omitted contrasts.
#' @export
#' @examples
#' sim <- cj_simulate_cells(donors = 6, cells_per_sample = 300, seed = 1)
#' samples <- unique(sim$cells[c("sample_id", "donor_id", "condition")])
#' samples$visit <- factor("v1")
#' samples$parent_population <- "CD4"
#' x <- cj_from_cells(sim$cells, c("IFNG", "TNF", "IL2"), samples, transform = "none",
#'                    control = "unstim", method = "quantile")
#' fs <- cj_fluorescence_scores(x, n_boot = 20, seed = 1)
#' head(fs[c("donor_id", "marker", "nicr", "poly_nicr", "nicr_lo", "nicr_hi")])
cj_fluorescence_scores <- function(x, n_boot = 0L, level = 0.95, seed = NULL,
                                   reference = "pooled", fis_min_n = 1L) {
  if (!inherits(x, "cj_data") || is.null(x$intensity))
    stop("x must come from cj_from_cells(); counts alone carry no fluorescence", call. = FALSE)
  if (!is.matrix(reference) && !reference %in% c("pooled", "single"))
    stop("reference must be \"pooled\", \"single\" or a pattern-by-marker matrix", call. = FALSE)
  if (n_boot > 0 && is.null(x$events))
    stop("the bootstrap needs events; build x with cj_from_cells(keep_events = TRUE)", call. = FALSE)
  ct <- x$contrasts
  none <- is.na(ct$control_sample_id)
  excluded <- ct[none, c("contrast_id", "donor_id", "visit", "parent_population", "stimulation")]
  if (nrow(excluded)) {
    message(nrow(excluded), " contrast(s) without a control sample omitted: ",
            paste(utils::head(excluded$contrast_id, 5), collapse = ", "),
            if (nrow(excluded) > 5) ", ...")
    excluded$reason <- "no control sample"
  }
  ct <- ct[!none, , drop = FALSE]
  if (!nrow(ct)) stop("no contrast has a control sample", call. = FALSE)
  bits <- x$bits; K <- nrow(bits); J <- ncol(bits)
  iS <- match(ct$stim_sample_id, x$samples$sample_id)
  iU <- match(ct$control_sample_id, x$samples$sample_id)

  # Reference weights per stimulation x parent population, fixed for all replicates.
  grp <- paste(ct$stimulation, ct$parent_population, sep = "\r")
  W <- lapply(split(seq_len(nrow(ct)), grp), function(r)
    fis_weights(colSums(x$counts[unique(c(iS[r], iU[r])), , drop = FALSE]), bits, reference))
  gw <- function(r) W[[grp[r]]]

  slice <- function(a, i) matrix(a[i, , ], K, J, dimnames = dimnames(a)[2:3])
  one <- function(nmat, sxarr) {
    do.call(rbind, lapply(seq_len(nrow(ct)), function(r)
      fluor_contrast(nmat[iS[r], ], slice(sxarr, iS[r]), nmat[iU[r], ], slice(sxarr, iU[r]),
                     bits, gw(r), fis_min_n)))
  }
  est <- one(x$counts, x$intensity$sx)
  key <- ct[rep(seq_len(nrow(ct)), each = J),
            c("contrast_id", "donor_id", "visit", "parent_population", "stimulation")]
  est <- cbind(key, est, row.names = NULL)

  if (n_boot > 0) {
    ev <- x$events
    code <- drop((x$events$X > rep(x$intensity$thresholds, each = nrow(ev$X))) %*% 2^(seq_len(J) - 1L)) + 1L
    Xt <- sweep(ev$X, 2L, x$intensity$thresholds)
    rows <- split(seq_along(ev$sample), factor(ev$sample, levels = seq_len(nrow(x$samples))))
    used <- sort(unique(c(iS, iU)))
    run <- function() {
      draws <- array(NA_real_, c(n_boot, nrow(est), length(FLUOR_SCORES)),
                     dimnames = list(NULL, NULL, FLUOR_SCORES))
      for (b in seq_len(n_boot)) {
        nmat <- x$counts; nmat[] <- 0L
        sxarr <- x$intensity$sx; sxarr[] <- 0
        for (i in used) if (length(rows[[i]])) {
          st <- boot_stats(rows[[i]], code, Xt, K)
          nmat[i, ] <- st$n; sxarr[i, , ] <- st$sx
        }
        r <- one(nmat, sxarr)
        for (s in FLUOR_SCORES) draws[b, , s] <- r[[s]]
      }
      draws
    }
    draws <- if (is.null(seed)) run() else withr_seed(seed, run())
    a <- (1 - level) / 2
    for (s in FLUOR_SCORES) {
      q <- apply(draws[, , s, drop = FALSE], 2L, function(v)
        if (all(is.na(v))) c(NA_real_, NA_real_) else stats::quantile(v, c(a, 1 - a), na.rm = TRUE, names = FALSE))
      est[[paste0(s, "_lo")]] <- q[1, ]
      est[[paste0(s, "_hi")]] <- q[2, ]
    }
  }
  attr(est, "excluded") <- excluded
  est
}
