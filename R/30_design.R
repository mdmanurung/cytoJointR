# M1 design: every fixed-effect matrix is built here, with named columns, and
# passed to Stan as data. The R mirror of the Stan linear predictor (m1_eta)
# lives here too; simulation, prior prediction and the likelihood check use it.

#' Design matrices and parameter map of the M1 count model
#'
#' Builds, in R, every matrix the Stan model receives as data. Fixed-effect
#' columns are named, so the map from Stan parameters to meaning is inspectable
#' and is stored in `fit$design`.
#'
#' Fixed effects per marker `j` (reference visit = first visit present):
#' stimulation (`stim_<m>`), visit (`visit_<v>`) and visit by stimulation
#' (`visit_<v>:stim_<m>`). Visit-by-stimulation columns with no stimulated sample
#' are dropped and listed in `$dropped`. Donor loadings `W` are `(1, stim
#' indicators)`: basal `b_i` and induced `r_im`, shared across markers and visits.
#'
#' @param x A `cj_data` object.
#' @param donor_effects Include donor effects `(b_i, r_i1..r_iM)`.
#' @param donor_visit Include the donor-by-visit basal deviation `u_iv`. Switched
#'   off with a message when the data hold a single visit.
#' @param centered Sample donor effects on their natural scale instead of the
#'   noncentered (standardised) scale. See [cj_fit()].
#' @param overdispersed Add the Dirichlet-multinomial sample residual (`inv_phi`,
#'   `phi` in `param_map`). `cj_fit()` sets it from its `model` argument.
#' @return A list of class `cj_design`: matrices `X`, `Xs`, `W`, `B`, `Bpair`,
#'   index vectors, flags, `col_info`, `dropped` and `param_map` (Stan name to label).
#' @export
cj_design <- function(x, donor_effects = TRUE, donor_visit = TRUE, centered = FALSE,
                      overdispersed = FALSE) {
  stopifnot(inherits(x, "cj_data"))
  m1_design(x$samples, x$bits, donor_effects, donor_visit, centered, overdispersed)
}

m1_design <- function(samples, bits, donor_effects = TRUE, donor_visit = TRUE, centered = FALSE,
                      overdispersed = FALSE) {
  pops <- unique(as.character(samples$parent_population))
  if (length(pops) > 1L)
    stop("samples hold ", length(pops), " parent populations (", paste(pops, collapse = ", "),
         "); v0.1 fits one. Use cj_subset(x, parent_population) first.", call. = FALSE)
  stims <- sort(setdiff(unique(as.character(samples$condition)), CONTROL_CONDITION))
  if (!length(stims)) stop("no stimulated samples: nothing to model", call. = FALSE)
  visits <- levels(droplevels(samples$visit))
  if (donor_visit && length(visits) < 2L) {
    message("single visit: donor_visit switched off")
    donor_visit <- FALSE
  }
  cond <- as.character(samples$condition)
  vis <- as.character(samples$visit)
  Xs <- vapply(stims, \(m) as.numeric(cond == m), numeric(nrow(samples)))
  Xs <- matrix(Xs, nrow(samples), dimnames = list(samples$sample_id, paste0("stim_", stims)))

  cols <- list(); info <- list(); dropped <- character()
  add <- function(name, v, type, stim, visit, prior) {
    cols[[name]] <<- v
    info[[name]] <<- data.frame(column = name, type = type, stimulation = stim,
                                visit = visit, prior = prior)
  }
  for (m in stims) add(paste0("stim_", m), as.numeric(cond == m), "stimulation", m, NA, "beta_sd")
  for (v in visits[-1])
    add(paste0("visit_", v), as.numeric(vis == v), "visit", NA, v, "visit_sd")
  for (v in visits[-1]) for (m in stims) {
    nm <- paste0("visit_", v, ":stim_", m)
    z <- as.numeric(vis == v & cond == m)
    if (sum(z) == 0) dropped <- c(dropped, nm)
    else add(nm, z, "visit_by_stimulation", m, v, "visit_stim_sd")
  }
  X <- matrix(unlist(cols), nrow(samples), dimnames = list(samples$sample_id, names(cols)))
  W <- cbind(basal = 1, `colnames<-`(Xs, paste0("induced_", stims)))

  pairs <- utils::combn(ncol(bits), 2L)
  Bpair <- apply(pairs, 2L, \(p) bits[, p[1]] * bits[, p[2]])
  colnames(Bpair) <- apply(pairs, 2L, \(p) paste(colnames(bits)[p], collapse = ":"))

  donors <- unique(as.character(samples$donor_id))
  dv_key <- paste(samples$donor_id, samples$visit, sep = "@")
  dv <- unique(data.frame(donor_id = as.character(samples$donor_id), visit = vis))
  d <- list(X = X, Xs = Xs, W = W, B = bits, Bpair = Bpair,
            markers = colnames(bits), stimulations = stims, visits = visits,
            col_info = do.call(rbind, unname(info)), dropped = dropped,
            donors = donors, donor_of = match(as.character(samples$donor_id), donors),
            dv = dv, dv_of = match(dv_key, paste(dv$donor_id, dv$visit, sep = "@")),
            use_donor = as.integer(donor_effects), use_dv = as.integer(donor_visit), centered = as.integer(centered),
            overdispersed = as.integer(overdispersed), sample_id = samples$sample_id)
  d$param_map <- m1_param_map(d)
  structure(d, class = "cj_design")
}

# Stan parameter name -> readable label, in Stan order (column-major matrices).
m1_param_map <- function(d) {
  J <- length(d$markers); M <- length(d$stimulations); Pf <- ncol(d$X); P <- ncol(d$Bpair)
  grid <- function(n1, n2) expand.grid(i = seq_len(n1), j = seq_len(n2))
  g <- grid(Pf, J); h <- grid(M, P)
  map <- rbind(
    data.frame(stan = sprintf("alpha[%d]", seq_len(J)), label = paste0("alpha[", d$markers, "]")),
    data.frame(stan = sprintf("gamma[%d,%d]", g$i, g$j),
               label = sprintf("gamma[%s, %s]", colnames(d$X)[g$i], d$markers[g$j])),
    data.frame(stan = sprintf("theta0[%d]", seq_len(P)), label = paste0("theta0[", colnames(d$Bpair), "]")),
    data.frame(stan = sprintf("delta[%d,%d]", h$i, h$j),
               label = sprintf("delta[%s, %s]", d$stimulations[h$i], colnames(d$Bpair)[h$j])))
  if (d$use_donor) {
    Dd <- ncol(d$W); ab <- utils::combn(Dd, 2L)
    map <- rbind(map,
      data.frame(stan = sprintf("tau[%d]", seq_len(Dd)), label = paste0("tau[", colnames(d$W), "]")),
      data.frame(stan = sprintf("rho[%d]", seq_len(ncol(ab))),
                 label = sprintf("rho[%s, %s]", colnames(d$W)[ab[1, ]], colnames(d$W)[ab[2, ]])))
  }
  if (d$use_dv) map <- rbind(map, data.frame(stan = "tau_u[1]", label = "tau_u"))
  if (d$overdispersed)
    map <- rbind(map, data.frame(stan = c("inv_phi[1]", "phi[1]"), label = c("inv_phi", "phi")))
  map
}

# Zero-size parameter blocks stay zero-size, as in the Stan model.
m1_zero_donor <- function(par, d) {
  if (!d$use_donor) par[c("tau", "L_Omega", "z")] <- list(numeric(0), matrix(0, 0, 0), matrix(0, 0, 0))
  if (!d$use_dv) par[c("tau_u", "u_raw")] <- list(numeric(0), numeric(0))
  par
}

# R mirror of eta_matrix() in m1.stan. `par` holds the Stan parameters
# (noncentered): alpha, gamma, theta0, delta, tau, L_Omega, z, tau_u, u_raw.
m1_eta <- function(d, par) {
  S <- nrow(d$X)
  shift <- numeric(S)
  if (d$use_donor) {
    DE <- if (d$centered) par$z else diag(par$tau, length(par$tau)) %*% par$L_Omega %*% par$z
    shift <- shift + rowSums(d$W * t(DE[, d$donor_of, drop = FALSE]))
  }
  if (d$use_dv) shift <- shift + (if (d$centered) par$u_raw else par$tau_u * par$u_raw)[d$dv_of]
  A <- sweep(d$X %*% par$gamma, 2L, par$alpha, "+") + shift
  TH <- sweep(d$Xs %*% par$delta, 2L, par$theta0, "+")
  eta <- A %*% t(d$B) + TH %*% t(d$Bpair)
  dimnames(eta) <- list(d$sample_id, pattern_labels(d$B))
  eta
}

# Row-wise softmax; probabilities are floored so cj_draws full support holds.
m1_softmax <- function(eta) {
  p <- exp(eta - apply(eta, 1L, max))
  p <- pmax(p / rowSums(p), 1e-300)
  p / rowSums(p)
}

# Stan parameters flattened in declaration order, column-major, with Stan names.
m1_stan_vector <- function(par, d) {
  flat <- function(nm, v) {
    if (!length(v)) return(NULL)
    v <- as.matrix(v)
    idx <- if (ncol(v) == 1L && nm %in% c("alpha", "theta0", "tau", "tau_u", "u_raw", "inv_phi"))
      sprintf("%s[%d]", nm, seq_along(v))
    else {
      g <- expand.grid(i = seq_len(nrow(v)), j = seq_len(ncol(v)))
      sprintf("%s[%d,%d]", nm, g$i, g$j)
    }
    setNames(as.vector(v), idx)
  }
  c(flat("alpha", par$alpha), flat("gamma", par$gamma), flat("theta0", par$theta0),
    flat("delta", par$delta), flat("tau", par$tau), flat("L_Omega", par$L_Omega),
    flat("z", par$z), flat("tau_u", par$tau_u), flat("u_raw", par$u_raw),
    flat("inv_phi", par$inv_phi))
}
