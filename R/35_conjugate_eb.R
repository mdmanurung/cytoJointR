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
#   u_b ~ Binomial(N_control, q_b);  s_bm ~ Binomial(N_stim, q_b * rho_bm)
# mu, sd_q, lambda_m and sd_r are fitted per parent population (and stimulation) by marginal maximum
# likelihood on a grid, so they come from the data. The stimulated counts also inform q_b, which
# is why the baseline is estimated better than from the control alone. Within the positive patterns each
# sample gets Dirichlet(a0 + kappa * pooled share + counts); kappa comes from a Dirichlet-multinomial
# marginal likelihood (kappa = 0 recovers C0's within-positive part).

EB_Q_POINTS <- 300L # grid points for log q
EB_R_POINTS <- 140L # grid points for log rho
EB_LOG_RHO <- c(-5, 6) # range of log rho
EB_SD_MIN <- 0.05 # lower bound of the fitted standard deviations (log scale)
EB_SD_MAX <- 3
EB_KAPPA_MAX <- 1e4

# log-likelihood of s positives out of n at success probability exp(lp), -Inf at or above 1.
eb_binom_ll <- function(s, n, lp) {
  ll <- s * lp + (n - s) * log1p(-exp(pmin(lp, -1e-12)))
  ll[!is.finite(ll) | lp >= 0] <- -Inf
  ll
}

# Per parent population: data for the marginal likelihood of one block.
# Returns the grids, the control vector `u` (Q) and one Q x R matrix `m` per stimulated sample,
# all scaled by their maximum (a constant that drops out of the maximisation and of the draws).
eb_block_data <- function(pos, N, ctrl, stims, lq, lr) {
  out <- list(ctrl = ctrl, stims = stims)
  scale_exp <- function(ll) exp(ll - max(ll[is.finite(ll)], -Inf))
  out$u <- if (!is.na(ctrl)) scale_exp(eb_binom_ll(pos[ctrl], N[ctrl], lq)) else rep(1, length(lq))
  out$m <- lapply(stims, function(i) {
    lp <- outer(lq, lr, `+`)
    scale_exp(eb_binom_ll(pos[i], N[i], lp))
  })
  out
}

# Negative marginal log-likelihood of the hyperparameters for one parent population.
# theta = (mu, log sd_q, lambda_1, log sd_r1, lambda_2, ...) over the stimulations in `stim_names`.
eb_negll <- function(theta, blocks, stim_of, n_stim, lq, lr, hq, hr, sd_min) {
  sd_q <- min(max(exp(theta[2]), sd_min), EB_SD_MAX)
  fq <- stats::dnorm(lq, theta[1], sd_q) * hq
  g <- vapply(seq_len(n_stim), function(m) {
    sd_r <- min(max(exp(theta[2 + 2 * m]), sd_min), EB_SD_MAX)
    stats::dnorm(lr, theta[1 + 2 * m], sd_r) * hr
  }, numeric(length(lr)))
  g <- matrix(g, ncol = n_stim)
  ll <- 0
  for (b in blocks) {
    v <- b$u
    for (k in seq_along(b$m)) v <- v * drop(b$m[[k]] %*% g[, stim_of[[k]]])
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
cj_fit_dirichlet_eb <- function(x, a0 = 0.5, n_draws = 4000L, seed = 1L, pool_patterns = TRUE,
                                sd_min = EB_SD_MIN) {
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

  # Any-positive fractions: draws of q for the control and q * rho for each stimulated sample ----
  q_draw <- matrix(NA_real_, n_draws, S)
  for (pp in unique(sm$parent_population)) {
    idx <- which(sm$parent_population == pp)
    block_key <- paste(sm$donor_id[idx], sm$visit[idx], sep = "\r")
    blocks <- list(); stim_of <- list(); members <- list()
    stim_names <- sort(unique(sm$condition[idx][!is_ctrl[idx]]))
    q_hat <- (pos_n[idx] + 0.5) / (N[idx] + 1)
    lq_rng <- c(log(0.5 / max(N[idx])) - 3, min(log(max(q_hat) * 1.5) + 3, log(0.99)))
    lq <- seq(lq_rng[1], lq_rng[2], length.out = EB_Q_POINTS)
    lr <- seq(EB_LOG_RHO[1], EB_LOG_RHO[2], length.out = EB_R_POINTS)
    hq <- diff(lq[1:2]); hr <- diff(lr[1:2])
    for (k in unique(block_key)) {
      ii <- idx[block_key == k]
      ctrl <- ii[is_ctrl[ii]]; st <- ii[!is_ctrl[ii]]
      ctrl <- if (length(ctrl)) ctrl[1] else NA_integer_
      b <- eb_block_data(pos_n, N, ctrl, st, lq, lr)
      blocks[[length(blocks) + 1L]] <- b
      stim_of[[length(stim_of) + 1L]] <- match(sm$condition[st], stim_names)
      members[[length(members) + 1L]] <- ii
    }
    n_stim <- length(stim_names)
    # Start: population baseline and mean fold change from the pooled fractions.
    mu0 <- log(mean(q_hat[is_ctrl[idx]], na.rm = TRUE))
    if (!is.finite(mu0)) mu0 <- log(mean(q_hat))
    th0 <- c(mu0, log(0.5), unlist(lapply(stim_names, function(s) {
      r <- mean(q_hat[sm$condition[idx] == s]) / exp(mu0)
      c(log(max(r, 0.2)), log(0.3))
    })))
    nll <- \(th) eb_negll(th, blocks, stim_of, n_stim, lq, lr, hq, hr, sd_min)
    fit <- stats::optim(th0, nll, method = "L-BFGS-B",
                        lower = c(-Inf, log(sd_min), rep(c(-4, log(sd_min)), n_stim)),
                        upper = c(Inf, log(EB_SD_MAX), rep(c(5, log(EB_SD_MAX)), n_stim)))
    th <- fit$par
    sd_q <- min(max(exp(th[2]), sd_min), EB_SD_MAX)
    fq <- stats::dnorm(lq, th[1], sd_q) * hq
    hyper[[pp]] <- data.frame(parent_population = pp, stimulation = c(NA, stim_names),
                              mean = c(th[1], th[1 + 2 * seq_len(n_stim)]),
                              sd = c(sd_q, pmin(pmax(exp(th[2 + 2 * seq_len(n_stim)]), sd_min), EB_SD_MAX)),
                              converged = fit$convergence == 0L)
    gm <- lapply(seq_len(n_stim), function(m) {
      sd_r <- min(max(exp(th[2 + 2 * m]), sd_min), EB_SD_MAX)
      stats::dnorm(lr, th[1 + 2 * m], sd_r) * hr
    })
    # Posterior draws per block: q from its marginal, then rho_m given q.
    for (bi in seq_along(blocks)) {
      b <- blocks[[bi]]
      v <- b$u
      mq <- lapply(seq_along(b$m), function(k) drop(b$m[[k]] %*% gm[[stim_of[[bi]][k]]]))
      for (k in seq_along(mq)) v <- v * mq[[k]]
      w <- fq * v
      qi <- sample.int(length(lq), n_draws, replace = TRUE, prob = w / sum(w))
      lq_d <- lq[qi] + stats::runif(n_draws, -hq / 2, hq / 2)
      if (!is.na(b$ctrl)) q_draw[, b$ctrl] <- exp(lq_d)
      for (k in seq_along(b$stims)) {
        pr <- b$m[[k]][qi, , drop = FALSE] * rep(gm[[stim_of[[bi]][k]]], each = n_draws)
        cum <- t(apply(pr, 1L, cumsum))
        u <- stats::runif(n_draws) * cum[, ncol(cum)]
        ri <- pmin(rowSums(cum < u) + 1L, length(lr))
        lr_d <- lr[ri] + stats::runif(n_draws, -hr / 2, hr / 2)
        q_draw[, b$stims[k]] <- exp(lq_d + lr_d)
      }
    }
  }
  q_draw <- pmin(q_draw, 1 - 1e-9)
  # a block with no control: its control columns stay NA and the contrast is unavailable downstream

  # Pattern vector within the positives: Dirichlet(a0 + kappa * pooled share + counts) ----
  pos_idx <- setdiff(seq_len(K), neg)
  for (i in seq_len(S)) {
    if (anyNA(q_draw[, i])) next
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
  # Samples left without draws (a block with no control and no stimulated sample cannot occur);
  # a stimulated sample whose block has no control still has draws, its control does not exist.
  p[is.na(p)] <- 1 / K
  p <- pmax(p, 1e-300)
  p <- p / as.vector(apply(p, 1:2, sum))
  new_cj_draws(p, x$samples, x$contrasts, x$bits,
               list(engine_id = "conjugate-eb", model_id = "C0-EB", seed = as.integer(seed),
                    n_draws = as.integer(n_draws), draw_kind = "in_sample_donor_conditional",
                    a0 = a0, pool_patterns = pool_patterns, hyper = do.call(rbind, hyper)))
}
