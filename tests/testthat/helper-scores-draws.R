# Test helpers for the score, visit-change and COMPASS tests.
EPS <- c(RMS = 0.005, poly_absolute = 0.005, poly_marginal = 0.005, CPS = 0.02,
         CCS = 0.005, breadth = 0.05, marker_magnitude = 0.005)

# A cj_data whose counts are multinomial draws from the mean probabilities of d.
example_data <- function(d, seed = 2) {
  set.seed(seed)
  K <- dim(d$p)[3]
  counts <- t(vapply(seq_len(nrow(d$samples)), function(i)
    rmultinom(1, d$samples$N[i], colMeans(d$p[, i, ]))[, 1], numeric(K)))
  storage.mode(counts) <- "integer"
  dimnames(counts) <- list(d$samples$sample_id, dimnames(d$p)[[3]])
  cytoJointR:::new_cj_data(counts, d$samples, d$contrasts, d$bits)
}
