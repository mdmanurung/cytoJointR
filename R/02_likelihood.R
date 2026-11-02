# Verification reference: Gaussian conditionally truncated at a fixed gate.
# One marker and one gate-defined population (positive OR negative).
# The conditional density normalization is essential; omitting it breaks the joint model.
#' @rdname cytoJointR-reference
#' @export
cj_log_trunc_gaussian_raw <- function(x, mu, sigma, threshold, positive) {
  stopifnot(length(sigma)==1L, sigma>0, all(is.finite(x)))
  if (positive && any(x <= threshold)) stop("Values violate positive gate")
  if (!positive && any(x > threshold)) stop("Values violate negative gate")
  lp <- if (positive) pnorm(threshold, mean=mu, sd=sigma,
                            lower.tail=FALSE, log.p=TRUE)
        else pnorm(threshold, mean=mu, sd=sigma,
                   lower.tail=TRUE, log.p=TRUE)
  sum(dnorm(x, mean=mu, sd=sigma, log=TRUE)) - length(x)*lp
}

#' @rdname cytoJointR-reference
#' @export
cj_log_trunc_gaussian_compressed <- function(n, sx, ss, mu, sigma,
                                              threshold, positive) {
  stopifnot(length(n)==1L, n>=0L, sigma>0)
  if (n == 0L) return(0)
  lp <- if (positive) pnorm(threshold, mean=mu, sd=sigma,
                            lower.tail=FALSE, log.p=TRUE)
        else pnorm(threshold, mean=mu, sd=sigma,
                   lower.tail=TRUE, log.p=TRUE)
  -0.5*n*log(2*pi) - n*log(sigma) -
    (ss - 2*mu*sx + n*mu^2)/(2*sigma^2) - n*lp
}

#' @rdname cytoJointR-reference
#' @export
cj_log_multinomial <- function(n, p) {
  stopifnot(length(n)==length(p), all(n>=0),
            all(p>0), abs(sum(p)-1)<1e-10)
  lgamma(sum(n)+1) - sum(lgamma(n+1)) + sum(n*log(p))
}
