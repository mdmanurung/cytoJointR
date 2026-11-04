# C0-EB: empirical-Bayes pooled conjugate-style draws, no Stan.
#
# C0 treats every sample alone, so a donor's RMS (the net any-positive fraction) is one noisy
# stimulated-minus-control difference. In rare responses the difference carries almost no information
# about the donor's baseline, yet the true RMS mostly follows the baseline. C0-EB pools donors on the
# any-positive scale and keeps C0's Dirichlet within the positive patterns.
#
# Any-positive model, per donor x visit x parent population ("block"):
#   log q_b      ~ N(mu, sd_q^2)        control any-positive fraction (the baseline)
#   log rho_bm   ~ N(lambda_m, sd_r^2)  fold change of stimulation m over the control
#   u_bc ~ Binomial(N_c, q_b) for every control c of the block;  s_bm ~ Binomial(N_stim, q_b * rho_bm)
# With null_mass = TRUE, rho is a mixture: a point mass at 1 (no response) with weight pi_m fitted per
# stimulation, and the log-normal otherwise.
# mu, sd_q, lambda_m and sd_r are fitted per parent population (and stimulation) by marginal maximum
# likelihood on a grid, so they come from the data. The stimulated counts also inform q_b, which
# is why the baseline is estimated better than from the control alone. Within the positive patterns each
# sample gets Dirichlet(a0 + counts), as C0. With pool_patterns = TRUE the prior is a0 + kappa * pooled share,
# kappa from a Dirichlet-multinomial marginal likelihood; dev simulations showed it shrinks CPS too far
# (coverage 0.2-0.5), so it is off by default.

EB_Q_POINTS <- 300L # grid points for log q
EB_R_POINTS <- 140L # grid points for log rho
EB_LOG_RHO <- c(-5, 6) # range of log rho
EB_SD_MIN <- 0.2 # lower bound of the fitted standard deviations (log scale): the smallest dev value with rare_weak coverage >= 0.90 (0.92, 0.91); see runs/2026-11-08_compass-magnitudes
EB_SD_MAX <- 3
EB_KAPPA_MAX <- 1e4

# log-likelihood of s positives out of n at success probability exp(lp), -Inf at or above 1.
eb_binom_ll <- function(s, n, lp) {
  ll <- s * lp + (n - s) * log1p(-exp(pmin(lp, -1e-12)))
  ll[!is.finite(ll) | lp >= 0] <- -Inf
  ll
}

# Data for the marginal likelihood of one block (a donor x visit x parent population).
# `ctrls` are all control samples of the block: they share q, so their likelihoods multiply.
# `u` (length Q) is that product; `m` holds one Q x (R + null) matrix per stimulated sample. With a null
# component, the last column is rho = 1 exactly. All are scaled by their maximum, a constant that drops
# out of the maximisation and of the draws.
eb_block_data <- function(pos, N, ctrls, stims, lq, lr, null_mass) {
  out <- list(ctrls = ctrls, stims = stims, u = rep(1, length(lq)))
  scale_exp <- function(ll) exp(ll - max(ll[is.finite(ll)], -Inf))
  for (i in ctrls) out$u <- out$u * scale_exp(eb_binom_ll(pos[i], N[i], lq))
  out$m <- lapply(stims, function(i) {
    lp <- outer(lq, lr, `+`)
    if (null_mass) lp <- cbind(lp, lq)
    scale_exp(eb_binom_ll(pos[i], N[i], lp))
  })
  out
}

# Hyperparameters from theta = (mu, log sd_q, then per stimulation: lambda, log sd_r [, logit pi]).
eb_unpack <- function(theta, n_stim, null_mass, sd_min) {
  stride <- 2L + null_mass
  at <- \(m, k) theta[2L + (m - 1L) * stride + k]
  list(mu = theta[1], sd_q = min(max(exp(theta[2]), sd_min), EB_SD_MAX),
       lambda = vapply(seq_len(n_stim), at, numeric(1), k = 1L),
       sd_r = pmin(pmax(exp(vapply(seq_len(n_stim), at, numeric(1), k = 2L)), sd_min), EB_SD_MAX),
       pi = if (null_mass) stats::plogis(vapply(seq_len(n_stim), at, numeric(1), k = 3L)) else rep(0, n_stim))
}

# Prior weights of rho on the grid (and, with a null component, the point mass at rho = 1 as last entry).
eb_rho_weights <- function(h, m, lr, hr, null_mass) {
  w <- stats::dnorm(lr, h$lambda[m], h$sd_r[m]) * hr
  if (null_mass) c((1 - h$pi[m]) * w, h$pi[m]) else w
}

# Negative marginal log-likelihood of the hyperparameters for one parent population.
eb_negll <- function(theta, blocks, stim_of, n_stim, lq, lr, hq, hr, sd_min, null_mass) {
  h <- eb_unpack(theta, n_stim, null_mass, sd_min)
  fq <- stats::dnorm(lq, h$mu, h$sd_q) * hq
  g <- lapply(seq_len(n_stim), eb_rho_weights, h = h, lr = lr, hr = hr, null_mass = null_mass)
  ll <- 0
  for (bi in seq_along(blocks)) {
    b <- blocks[[bi]]
    v <- b$u
    for (k in seq_along(b$m)) v <- v * drop(b$m[[k]] %*% g[[stim_of[[bi]][k]]])
    ll <- ll + log(sum(fq * v) + 1e-300)
  }
  -ll
}

# Dirichlet-multinomial marginal log-likelihood of the within-positive counts at concentration kappa.
eb_kappa_ll <- function(log_kappa, Y, w) {
  a <- exp(log_kappa) * w
  A <- sum(a)
  n <- rowSums(Y)
  sum(lgamma(A) - lgamma(A + n)) + sum(rowSums(lgamma(sweep(Y, 2L, a, `+`)) - rep(lgamma(a), each = nrow(Y))))
}

#' @noRd
cj_fit_dirichlet_eb <- function(x, a0 = 0.5, n_draws = 4000L, seed = 1L, pool_patterns = FALSE,
                                sd_min = EB_SD_MIN, null_mass = FALSE, min_blocks = 8L) {
  stopifnot(inherits(x, "cj_data"), a0 > 0)
  set.seed(seed)
  sm <- x$samples; cnt <- x$counts; bits <- x$bits
  K <- ncol(cnt); S <- nrow(cnt)
  neg <- which(rowSums(bits) == 0L)
  N <- rowSums(cnt)
  pos_n <- N - cnt[, neg]
  is_ctrl <- sm$condition == CONTROL_CONDITION
  p <- array(NA_real_, c(n_draws, S, K), dimnames = list(NULL, sm$sample_id, colnames(cnt)))
  hyper <- list()

  # Any-positive fractions: draws of q for every control and q * rho for each stimulated sample ----
  q_draw <- matrix(NA_real_, n_draws, S)
  for (pp in unique(sm$parent_population)) {
    idx <- which(sm$parent_population == pp)
    block_key <- paste(sm$donor_id[idx], sm$visit[idx], sep = "\r")
    n_blocks <- length(unique(block_key))
    if (n_blocks < min_blocks)
      warning("parent population ", pp, " has ", n_blocks, " donor-visit blocks; the across-donor ",
              "hyperparameters are poorly determined below about ", min_blocks, " (use magnitude_model = \"c0\")",
              call. = FALSE)
    blocks <- list(); stim_of <- list()
    stim_names <- sort(unique(sm$condition[idx][!is_ctrl[idx]]))
    q_hat <- (pos_n[idx] + 0.5) / (N[idx] + 1)
    lq_rng <- c(log(0.5 / max(N[idx])) - 3, min(log(max(q_hat) * 1.5) + 3, log(0.99)))
    lq <- seq(lq_rng[1], lq_rng[2], length.out = EB_Q_POINTS)
    lr <- seq(EB_LOG_RHO[1], EB_LOG_RHO[2], length.out = EB_R_POINTS)
    hq <- diff(lq[1:2]); hr <- diff(lr[1:2])
    for (k in unique(block_key)) {
      ii <- idx[block_key == k]
      ctrls <- ii[is_ctrl[ii]]; st <- ii[!is_ctrl[ii]]
      blocks[[length(blocks) + 1L]] <- eb_block_data(pos_n, N, ctrls, st, lq, lr, null_mass)
      stim_of[[length(stim_of) + 1L]] <- match(sm$condition[st], stim_names)
    }
    n_stim <- length(stim_names)
    # Start: population baseline and mean fold change from the pooled fractions.
    mu0 <- log(mean(q_hat[is_ctrl[idx]], na.rm = TRUE))
    if (!is.finite(mu0)) mu0 <- log(mean(q_hat))
    th0 <- c(mu0, log(0.5), unlist(lapply(stim_names, function(s) {
      r <- mean(q_hat[sm$condition[idx] == s]) / exp(mu0)
      c(log(max(r, 0.2)), log(0.3), if (null_mass) stats::qlogis(0.3))
    })))
    stride_lo <- c(-4, log(sd_min), if (null_mass) -7)
    stride_hi <- c(5, log(EB_SD_MAX), if (null_mass) 3)
    nll <- \(th) eb_negll(th, blocks, stim_of, n_stim, lq, lr, hq, hr, sd_min, null_mass)
    fit <- stats::optim(th0, nll, method = "L-BFGS-B",
                        lower = c(-Inf, log(sd_min), rep(stride_lo, n_stim)),
                        upper = c(Inf, log(EB_SD_MAX), rep(stride_hi, n_stim)))
    h <- eb_unpack(fit$par, n_stim, null_mass, sd_min)
    fq <- stats::dnorm(lq, h$mu, h$sd_q) * hq
    hyper[[pp]] <- data.frame(parent_population = pp, stimulation = c(NA, stim_names),
                              mean = c(h$mu, h$lambda), sd = c(h$sd_q, h$sd_r), null_mass = c(NA, h$pi),
                              converged = fit$convergence == 0L)
    gm <- lapply(seq_len(n_stim), eb_rho_weights, h = h, lr = lr, hr = hr, null_mass = null_mass)
    # Posterior draws per block: q from its marginal, then rho_m given q. Every control replicate of the
    # block gets the same q draw: they measure one baseline.
    for (bi in seq_along(blocks)) {
      b <- blocks[[bi]]
      v <- b$u
      for (k in seq_along(b$m)) v <- v * drop(b$m[[k]] %*% gm[[stim_of[[bi]][k]]])
      w <- fq * v
      qi <- sample.int(length(lq), n_draws, replace = TRUE, prob = w / sum(w))
      lq_d <- lq[qi] + stats::runif(n_draws, -hq / 2, hq / 2)
      for (i in b$ctrls) q_draw[, i] <- exp(lq_d)
      for (k in seq_along(b$stims)) {
        pr <- b$m[[k]][qi, , drop = FALSE] * rep(gm[[stim_of[[bi]][k]]], each = n_draws)
        cum <- t(apply(pr, 1L, cumsum))
        u <- stats::runif(n_draws) * cum[, ncol(cum)]
        ri <- pmin(rowSums(cum < u) + 1L, ncol(pr))
        lr_d <- ifelse(ri > length(lr), 0, lr[pmin(ri, length(lr))] + stats::runif(n_draws, -hr / 2, hr / 2))
        q_draw[, b$stims[k]] <- exp(lq_d + lr_d)
      }
    }
  }
  q_draw <- pmin(q_draw, 1 - 1e-9)
  # Every sample is in exactly one block, so every sample has draws; a block without a control only has
  # no control column, and its contrast is unavailable downstream.
  if (anyNA(q_draw)) stop("internal error: a sample received no any-positive draws", call. = FALSE)

  # Pattern vector within the positives: Dirichlet(a0 [+ kappa * pooled share] + counts) ----
  pos_idx <- setdiff(seq_len(K), neg)
  for (i in seq_len(S)) {
    y <- cnt[i, pos_idx]
    alpha <- a0 + y
    if (pool_patterns && length(pos_idx) > 1L) {
      same <- which(sm$parent_population == sm$parent_population[i] & sm$condition == sm$condition[i])
      Y <- cnt[same, pos_idx, drop = FALSE]
      Y <- Y[rowSums(Y) > 0, , drop = FALSE]
      if (nrow(Y) >= 2L) {
        w <- (colSums(Y) + a0) / sum(colSums(Y) + a0)
        lk <- stats::optimise(eb_kappa_ll, c(log(0.1), log(EB_KAPPA_MAX)), Y = Y, w = w, maximum = TRUE)$maximum
        alpha <- a0 + exp(lk) * w + y
      }
    }
    g <- matrix(stats::rgamma(n_draws * length(pos_idx), shape = rep(alpha, each = n_draws)), n_draws)
    w_d <- g / rowSums(g)
    p[, i, pos_idx] <- q_draw[, i] * w_d
    p[, i, neg] <- 1 - q_draw[, i]
  }
  p <- pmax(p, 1e-300)
  p <- p / as.vector(apply(p, 1:2, sum))
  new_cj_draws(p, x$samples, x$contrasts, x$bits,
               list(engine_id = "conjugate-eb", model_id = "C0-EB", seed = as.integer(seed),
                    n_draws = as.integer(n_draws), draw_kind = "in_sample_donor_conditional",
                    a0 = a0, pool_patterns = pool_patterns, null_mass = null_mass,
                    hyper = do.call(rbind, hyper)))
}
