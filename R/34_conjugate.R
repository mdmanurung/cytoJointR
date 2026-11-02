# C0 baseline: independent conjugate Dirichlet per sample, no Stan, no pooling.

# p_s ~ Dirichlet(a0 + y_s) for every sample. Returns a cj_draws, so the score
# layer runs on it unchanged. Draws are conditional on the observed sample, as in
# cj_fit(); the mean is (a0 + y) / (K * a0 + N). Sets the global RNG seed.
cj_fit_dirichlet <- function(x, a0 = 0.5, n_draws = 4000L, seed = 1L) {
  stopifnot(inherits(x, "cj_data"), a0 > 0)
  set.seed(seed)
  dims <- c(n_draws, dim(x$counts))
  g <- array(stats::rgamma(prod(dims), shape = rep(as.vector(x$counts) + a0, each = n_draws)),
             dims, dimnames = list(NULL, rownames(x$counts), colnames(x$counts)))
  p <- g / as.vector(apply(g, 1:2, sum))
  p <- pmax(p, 1e-300)   # tiny shape can underflow to 0; keep full support
  p <- p / as.vector(apply(p, 1:2, sum))
  new_cj_draws(p, x$samples, x$contrasts, x$bits,
               list(engine_id = "conjugate", model_id = "C0", seed = as.integer(seed),
                    n_draws = as.integer(n_draws), draw_kind = "in_sample_donor_conditional",
                    a0 = a0))
}
