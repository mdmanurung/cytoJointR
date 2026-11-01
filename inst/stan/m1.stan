// cytoJointR M1: donor-paired log-linear model of exact marker-pattern counts.
// eta_sk = sum_j B_kj a_sj + sum_{j<l} B_kj B_kl theta_sjl, p_s = softmax(eta_s).
// a_sj  = alpha_j + X_s gamma_j + (donor effect, loading W_s) + (donor x visit effect).
// theta_s = theta0 + Xs_s delta.   All design matrices are built in R (cj_design()).
functions {
  matrix eta_matrix(vector alpha, matrix gamma, vector theta0, matrix delta,
                    vector tau, matrix L_Omega, matrix z, vector tau_u, vector u_raw,
                    matrix X, matrix Xs, matrix W, array[] int donor_of,
                    array[] int dv_of, matrix Bt, matrix Bpt) {
    int S = rows(X);
    int J = rows(Bt);
    vector[S] shift = rep_vector(0, S);
    if (num_elements(z) > 0) {
      matrix[cols(W), cols(z)] DE = diag_pre_multiply(tau, L_Omega) * z;
      shift += columns_dot_product(W', DE[:, donor_of])';
    }
    if (num_elements(u_raw) > 0) {
      vector[num_elements(u_raw)] u = tau_u[1] * u_raw;
      shift += u[dv_of];
    }
    matrix[S, J] A = rep_matrix(alpha', S) + X * gamma + shift * rep_row_vector(1, J);
    matrix[S, cols(delta)] TH = rep_matrix(theta0', S) + Xs * delta;
    return A * Bt + TH * Bpt;
  }
}
data {
  int<lower=1> S;                        // samples
  int<lower=2> K;                        // patterns, 2^J
  int<lower=2> J;                        // markers
  int<lower=1> M;                        // stimulations
  int<lower=1> P;                        // marker pairs
  int<lower=1> Pf;                       // fixed-effect columns (no intercept)
  int<lower=1> Ndon;
  int<lower=1> G;                        // donor x visit groups
  int<lower=0, upper=1> use_donor;
  int<lower=0, upper=1> use_dv;
  array[S, K] int<lower=0> y;
  matrix[S, Pf] X;
  matrix[S, M] Xs;                       // stimulation indicators
  matrix[S, 1 + M] W;                    // donor loadings: (1, Xs)
  array[S] int<lower=1, upper=Ndon> donor_of;
  array[S] int<lower=1, upper=G> dv_of;
  matrix[J, K] Bt;
  matrix[P, K] Bpt;
  // priors
  vector[J] alpha_mean;
  real<lower=0> alpha_sd;
  vector<lower=0>[Pf] gamma_sd;
  real<lower=0> theta0_sd;
  real<lower=0> delta_sd;
  real<lower=0> tau_sd;
  real<lower=0> tau_u_sd;
  real<lower=0> lkj_eta;
}
transformed data {
  int Dd = use_donor ? 1 + M : 0;
  int Nz = use_donor ? Ndon : 0;
  int Gu = use_dv ? G : 0;
}
parameters {
  vector[J] alpha;
  matrix[Pf, J] gamma;
  vector[P] theta0;
  matrix[M, P] delta;
  vector<lower=0>[Dd] tau;
  cholesky_factor_corr[Dd] L_Omega;
  matrix[Dd, Nz] z;
  vector<lower=0>[use_dv] tau_u;
  vector[Gu] u_raw;
}
model {
  alpha ~ normal(alpha_mean, alpha_sd);
  for (j in 1:J) gamma[, j] ~ normal(0, gamma_sd);
  theta0 ~ normal(0, theta0_sd);
  to_vector(delta) ~ normal(0, delta_sd);
  if (use_donor) {
    tau ~ normal(0, tau_sd);
    L_Omega ~ lkj_corr_cholesky(lkj_eta);
    to_vector(z) ~ std_normal();
  }
  if (use_dv) {
    tau_u ~ normal(0, tau_u_sd);
    u_raw ~ std_normal();
  }
  matrix[S, K] eta = eta_matrix(alpha, gamma, theta0, delta, tau, L_Omega, z, tau_u,
                                u_raw, X, Xs, W, donor_of, dv_of, Bt, Bpt);
  for (s in 1:S) y[s] ~ multinomial_logit(eta[s]');
}
generated quantities {
  array[S] vector[K] p;
  vector[S] log_lik;
  vector[Dd * (Dd - 1) %/% 2] rho;
  {
    matrix[S, K] eta = eta_matrix(alpha, gamma, theta0, delta, tau, L_Omega, z, tau_u,
                                  u_raw, X, Xs, W, donor_of, dv_of, Bt, Bpt);
    for (s in 1:S) {
      p[s] = softmax(eta[s]');
      log_lik[s] = multinomial_logit_lpmf(y[s] | eta[s]');
    }
    if (Dd > 1) {
      matrix[Dd, Dd] Omega = multiply_lower_tri_self_transpose(L_Omega);
      int c = 0;
      for (a in 1:(Dd - 1)) for (b in (a + 1):Dd) { c += 1; rho[c] = Omega[a, b]; }
    }
  }
}
