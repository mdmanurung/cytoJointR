# Generator A: counts simulated exactly from M1 with known truth.

# Truth used when the caller gives none. Values are fixed functions of the
# index, so a simulation differs between runs only through seed-driven draws.
m1_default_truth <- function(d) {
  J <- length(d$markers); M <- length(d$stimulations); P <- ncol(d$Bpair)
  ci <- d$col_info
  gamma <- sapply(seq_len(J), \(j) {
    m <- match(ci$stimulation, d$stimulations)
    ifelse(ci$type == "stimulation", 2.5 + 0.2 * (j - 1) - 0.4 * (m - 1),
           ifelse(ci$type == "visit", 0.3, -0.6))
  }) |> matrix(nrow(ci), J, dimnames = list(ci$column, d$markers))
  Dd <- M + 1L
  Omega <- matrix(0.3, Dd, Dd); diag(Omega) <- 1
  list(alpha = setNames(-4.5 + 0.3 * (seq_len(J) - 1), d$markers), gamma = gamma,
       theta0 = setNames(0.6 + 0.1 * (seq_len(P) - 1), colnames(d$Bpair)),
       delta = matrix(rep(c(0.4, -0.4)[(seq_len(M) - 1) %% 2 + 1], P), M, P,
                      dimnames = list(d$stimulations, colnames(d$Bpair))),
       tau = setNames(c(0.5, rep(0.3, M)), colnames(d$W)), Omega = Omega, tau_u = 0.3)
}

#' Simulate counts from the M1 model (generator A)
#'
#' Draws data from the model the fit assumes, with known parameters. Each donor
#' has, at every visit, one control sample (unless dropped) and one sample per
#' stimulation. Donor effects are drawn from the truth's `tau` and `Omega`.
#'
#' @param n_donors Number of donors.
#' @param markers Marker names (2 to 5).
#' @param stimulations Stimulation names.
#' @param visits Visit labels (first is the reference visit).
#' @param N Cells per sample.
#' @param truth Named list that overrides the default truth: `alpha` (length J),
#'   `gamma` (fixed-effect columns by J, see `cj_design()$col_info`), `theta0`,
#'   `delta` (stimulations by pairs), `tau` (length 1 + M), `Omega`, `tau_u`.
#' @param donor_effects,donor_visit Simulate donor effects and donor-by-visit
#'   deviations (a switched-off effect is absent from the truth).
#' @param missing_control Fraction of donor-visits whose control sample is dropped.
#' @param seed Integer seed.
#' @return A list: `data` (`cj_data`), `truth` (Stan-scale parameters `alpha`,
#'   `gamma`, `theta0`, `delta`, `tau`, `L_Omega`, `z`, `tau_u`, `u_raw`, plus the
#'   realised donor effects `donor_effect` (`b`, `r`) and `u`), `p` (true
#'   per-sample probabilities, samples by patterns), `design`.
#' @examples
#' sim <- cj_simulate_counts(n_donors = 4, markers = c("IFNG", "TNF"), seed = 1)
#' sim$data$samples
#' round(sim$p, 3)
#' @export
cj_simulate_counts <- function(n_donors = 12L, markers = c("IFNG", "TNF"),
                               stimulations = "CMV", visits = "v1", N = 10000L,
                               truth = list(), donor_effects = TRUE, donor_visit = TRUE,
                               missing_control = 0, seed = 1L) {
  set.seed(seed)
  bits <- cj_bits(markers)
  donors <- sprintf("D%02d", seq_len(n_donors))
  grid <- expand.grid(condition = c(CONTROL_CONDITION, stimulations), visit = visits,
                      donor_id = donors, stringsAsFactors = FALSE)
  grid <- grid[order(grid$donor_id, grid$visit, match(grid$condition, c(CONTROL_CONDITION, stimulations))), ]
  dv <- unique(grid[c("donor_id", "visit")])
  drop <- dv[sample.int(nrow(dv), floor(nrow(dv) * missing_control)), ]
  grid <- grid[!(grid$condition == CONTROL_CONDITION &
                 paste(grid$donor_id, grid$visit) %in% paste(drop$donor_id, drop$visit)), ]
  samples <- data.frame(
    sample_id = paste(grid$donor_id, grid$visit, grid$condition, sep = "_"),
    donor_id = grid$donor_id, visit = factor(grid$visit, levels = visits),
    parent_population = "CD4", condition = grid$condition, N = as.numeric(N))
  d <- m1_design(samples, bits, donor_effects, donor_visit)

  tr <- utils::modifyList(m1_default_truth(d), truth)
  stopifnot(identical(dim(tr$gamma), c(ncol(d$X), length(markers))))
  par <- list(alpha = unname(tr$alpha), gamma = unname(tr$gamma), theta0 = unname(tr$theta0),
              delta = unname(tr$delta))
  if (d$use_donor) {
    par$tau <- unname(tr$tau)
    par$L_Omega <- t(chol(tr$Omega))
    par$z <- matrix(stats::rnorm(ncol(d$W) * length(d$donors)), ncol(d$W))
  }
  if (d$use_dv) {
    par$tau_u <- tr$tau_u
    par$u_raw <- stats::rnorm(nrow(d$dv))
  }
  par <- m1_zero_donor(par, d)
  p <- m1_softmax(m1_eta(d, par))
  counts <- t(vapply(seq_len(nrow(p)), \(s) stats::rmultinom(1, N, p[s, ])[, 1], integer(nrow(bits))))
  dimnames(counts) <- list(samples$sample_id, pattern_labels(bits))

  contrasts <- derive_contrasts(samples)
  data <- new_cj_data(counts, samples, contrasts, bits,
                      provenance = list(generator = "A (exact M1)", seed = seed))

  truth_out <- par
  if (d$use_donor) {
    de <- diag(par$tau, length(par$tau)) %*% par$L_Omega %*% par$z
    dimnames(de) <- list(colnames(d$W), d$donors)
    truth_out$donor_effect <- de
  }
  if (d$use_dv) truth_out$u <- setNames(par$tau_u * par$u_raw, paste(d$dv$donor_id, d$dv$visit, sep = "@"))
  dimnames(truth_out$gamma) <- dimnames(tr$gamma)
  list(data = data, truth = truth_out, p = `dimnames<-`(p, dimnames(counts)), design = d)
}
