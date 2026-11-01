# Deterministic reference algebra only: supplied probabilities and conditional
# excesses are NOT a fitted model and carry no uncertainty by themselves.
# First marker is the low bit, matching cj_bits().
#' Documented in man/ (hand-written Rd).
#' @noRd
#' @export
cj_score_reference <- function(p_u, p_s, bits, excess_u = NULL,
                               excess_s = NULL, reference_weights = NULL,
                               reference_scale = NULL) {
  if (!is.matrix(bits) || !is.numeric(bits) || anyNA(bits) ||
      !all(bits %in% 0:1) || ncol(bits) < 1L ||
      nrow(bits) != 2^ncol(bits) || anyDuplicated(as.data.frame(bits)))
    stop("bits must enumerate every exact binary pattern once")
  K <- nrow(bits); J <- ncol(bits)
  for (p in list(p_u, p_s)) {
    if (!is.numeric(p) || length(p) != K || any(!is.finite(p)) ||
        any(p <= 0) || abs(sum(p) - 1) > 1e-10)
      stop("Supply normalized, strictly positive fitted pattern probabilities; no pseudocounts are added")
  }
  if (is.null(colnames(bits))) colnames(bits) <- paste0("marker", seq_len(J))
  if (anyDuplicated(colnames(bits))) stop("Marker names must be unique")
  p_cf <- cj_counterfactual(p_u, p_s, bits)
  degree <- rowSums(bits); positive <- degree >= 1L; poly <- degree >= 2L
  q_u <- sum(p_u[positive]); q_s <- sum(p_s[positive]); dp <- p_s - p_u
  # No denominator truncation: small positive q remains visible in diagnostics.
  cps <- numeric(K)
  cps[poly] <- p_s[poly]/q_s - p_u[poly]/q_u
  rms <- ifelse(positive, dp, 0)
  ccs <- ifelse(poly, p_s - p_cf, 0)
  poly_absolute <- ifelse(poly, dp, 0)
  poly_marginal <- ifelse(poly, p_cf - p_u, 0)
  breadth <- degree * (p_s/q_s - p_u/q_u)
  pattern_label <- apply(bits, 1L, function(b)
    paste0(colnames(bits), ifelse(b == 1L, "+", "-"), collapse = " "))
  # pattern_id follows the actual bits, even for a permuted table.
  patterns <- data.frame(
    pattern_id = drop(bits %*% 2^(seq_len(J)-1L)), pattern = pattern_label,
    degree = degree, p_u = p_u, p_s = p_s, p_cf = p_cf,
    net_pattern = dp, RMS = rms, CPS = cps, CCS = ccs,
    poly_absolute = poly_absolute, poly_marginal = poly_marginal,
    coordination_residual = p_s - p_cf, breadth = breadth)
  global <- c(RMS = sum(rms), poly_absolute = sum(poly_absolute),
              poly_marginal = sum(poly_marginal),
              CPS = sum(cps), CCS = sum(ccs), breadth = sum(breadth))
  markers <- data.frame(marker = colnames(bits),
                        magnitude = drop(crossprod(bits, dp)))
  intensity <- NULL
  if (xor(is.null(excess_u), is.null(excess_s)))
    stop("Supply both condition-specific excess matrices")
  if (is.null(excess_u) && !is.null(reference_weights))
    stop("FIS reference weights require both excess matrices")
  if (!is.null(reference_scale) && is.null(reference_weights))
    stop("Standardized FIS requires fixed reference weights")
  if (!is.null(reference_scale)) {
    if (!is.numeric(reference_scale) || length(reference_scale) != J ||
        !identical(names(reference_scale), colnames(bits)) ||
        any(!is.finite(reference_scale)) || any(reference_scale <= 0))
      stop("Reference scales must be finite positive values named in marker order")
  }
  if (!is.null(excess_u)) {
    for (x in list(excess_u, excess_s)) {
      if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(K,J)) ||
          any(!is.finite(x)) || any(x < 0) || any(x[bits == 0] != 0))
        stop("Excess must be finite nonnegative K by J matrices, zero at marker-negative patterns")
    }
    burden_u <- sweep(excess_u, 1L, p_u, `*`)
    burden_s <- sweep(excess_s, 1L, p_s, `*`)
    net <- burden_s - burden_u
    frequency <- sweep((excess_s + excess_u)/2, 1L, dp, `*`)
    strength <- sweep(excess_s - excess_u, 1L, (p_s + p_u)/2, `*`)
    markers$burden_u <- colSums(burden_u)
    markers$burden_s <- colSums(burden_s)
    markers$NICR <- colSums(net)
    markers$Poly_NICR <- colSums(net[poly, , drop = FALSE])
    markers$frequency_contribution <- colSums(frequency)
    markers$intensity_contribution <- colSums(strength)
    markers$FIS <- NA_real_
    markers$FIS_standardized <- NA_real_
    fis_contribution <- matrix(NA_real_, K, J)
    fis_standardized_contribution <- matrix(NA_real_, K, J)
    if (!is.null(reference_weights)) {
      w <- reference_weights
      if (!is.matrix(w) || !is.numeric(w) || !identical(dim(w), c(K,J)) ||
          any(!is.finite(w)) || any(w < 0) || any(w[bits == 0] != 0) ||
          any(abs(colSums(w) - 1) > 1e-10))
        stop("Fixed FIS weights must be nonnegative K by J, zero for negative states, columns sum to one")
      fis_contribution <- w * (excess_s - excess_u)
      markers$FIS <- colSums(fis_contribution)
      if (!is.null(reference_scale)) {
        fis_standardized_contribution <- sweep(fis_contribution, 2L,
                                               reference_scale, `/`)
        markers$FIS_standardized <- colSums(fis_standardized_contribution)
      }
    }
    intensity <- data.frame(
      pattern_id = rep(patterns$pattern_id, times = J),
      pattern = rep(pattern_label, times = J),
      marker = rep(colnames(bits), each = K),
      degree = rep(degree, times = J),
      marker_positive = as.vector(bits) == 1L,
      excess_u = as.vector(excess_u), excess_s = as.vector(excess_s),
      net_burden = as.vector(net),
      frequency_contribution = as.vector(frequency),
      intensity_contribution = as.vector(strength),
      FIS_contribution = as.vector(fis_contribution),
      FIS_standardized_contribution = as.vector(fis_standardized_contribution))
    intensity <- intensity[intensity$marker_positive, , drop = FALSE]
  }
  list(global = global, patterns = patterns, markers = markers,
       intensity = intensity,
       diagnostics = list(q_u = q_u, q_s = q_s,
         ipf_margin_error = max(abs(crossprod(bits, p_cf - p_s))),
         status = "DETERMINISTIC_REFERENCE_NOT_FITTED_INFERENCE"))
}

# No automatic omission of invalid draws. They indicate a model/draw/score
# failure requiring inspection; filtering them can distort the target law.
#' Documented in man/ (hand-written Rd).
#' @noRd
#' @export
cj_summarize_signed <- function(draws, epsilon, level = 0.95) {
  if (!is.numeric(draws) || !length(draws) || any(!is.finite(draws)))
    stop("Draws must be a nonempty finite numeric vector; resolve missing/failed draws explicitly")
  if (length(epsilon) != 1L || !is.finite(epsilon) || epsilon < 0)
    stop("Supply one prespecified nonnegative epsilon on the score scale")
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must be between zero and one")
  bounds <- stats::quantile(draws, c((1-level)/2, (1+level)/2), names = FALSE)
  c(estimate = mean(draws), lower = bounds[1], upper = bounds[2],
    epsilon = epsilon,
    p_positive = mean(draws > epsilon),
    p_negative = mean(draws < -epsilon),
    p_negligible = mean(abs(draws) <= epsilon),
    n_draws = length(draws), n_invalid = 0)
}
