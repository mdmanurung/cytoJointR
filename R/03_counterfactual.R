# Model-based IPF counterfactual: matching stimulated margins while preserving
# all baseline log-linear interactions. Baseline probabilities must be positive.
#' @rdname cytoJointR-reference
#' @export
cj_counterfactual <- function(p_u, p_s, bits, tol=1e-10, max_iter=10000L) {
  stopifnot(is.matrix(bits), nrow(bits)==length(p_u),
            length(p_u)==length(p_s), all(p_u>0), all(p_s>0),
            all(is.finite(p_u)), all(is.finite(p_s)))
  p_u <- p_u/sum(p_u); p_s <- p_s/sum(p_s)
  target <- drop(crossprod(bits, p_s))
  if (any(target<=0 | target>=1)) stop("Target margins must be interior")
  q <- p_u
  for (iter in seq_len(max_iter)) {
    for (j in seq_len(ncol(bits))) {
      plus <- bits[,j]==1L
      a <- sum(q[plus])
      if (a<=0 || a>=1) stop("Baseline support inadequate")
      q[plus] <- q[plus] * target[j]/a
      q[!plus] <- q[!plus] * (1-target[j])/(1-a)
      q <- q/sum(q)
    }
    if (max(abs(drop(crossprod(bits,q)) - target)) < tol) return(q)
  }
  stop("IPF failed to converge")
}

# Compatibility name for earlier reference callers. All score algebra now lives
# in cj_score_reference(); no separate small-denominator policy is applied here.
#' @rdname cytoJointR-reference
#' @export
cj_scores_from_probs <- function(p_u, p_s, bits, min_functional = NULL) {
  if (!is.null(min_functional))
    warning("min_functional is deprecated and ignored; inspect q_u/q_s and support diagnostics")
  g <- cj_score_reference(p_u, p_s, bits)$global
  c(magnitude = g[["RMS"]], poly_absolute = g[["poly_absolute"]],
    poly_conditional = g[["CPS"]], coordination = g[["CCS"]],
    breadth = g[["breadth"]])
}
