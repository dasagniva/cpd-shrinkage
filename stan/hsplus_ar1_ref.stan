// =============================================================================
// stan/hsplus_ar1_ref.stan
// Reference Model: Significance-Guided Global-Local Shrinkage under AR(1) Noise
// Implements Section 3 (Equations 1-6) with independent outlier scales (rho=0, delta=0).
// Non-centred parameterisations for increments omega and outliers psi.
// Prewhitened likelihood representation.
//
// The R wrapper fits this model to robust_standardize()d data, matching the Gibbs
// sampler (R/31_gibbs_sampler.R): mu_1 ~ N(0, 100^2), omega_t ~ N(0, sigma^2 tau^2
// kappa_t^2 eta_t^2 lambda_t^2), psi_t ~ N(0, tau_psi^2 eta_psi_t^2 lambda_psi_t^2).
// use_slab = 1 switches both to the regularised (Piironen-Vehtari) slab, the
// fallback in Section 5 P3 if divergences persist; the Gibbs sampler has no slab,
// so Stan-Gibbs agreement checks must use use_slab = 0.
// =============================================================================

data {
  int<lower=2> n;                      // Length of time series (n <= 300 for validation)
  vector[n] y;                         // Observed (standardised) time series
  vector<lower=0>[n] kappa;            // NSP-informed local scale multiplier (Eq. 5)
  real<lower=0> s_tau;                 // Prior scale for tau_omega (default 1 / sqrt(n))
  real<lower=0> a_phi;                 // Prior beta shape 1 for (phi + 1)/2 (default 10)
  real<lower=0> b_phi;                 // Prior beta shape 2 for (phi + 1)/2 (default 2)
  real<lower=0> mu1_sd;                // Prior sd of mu_1 (default 100, standardised units)
  int<lower=0, upper=1> use_slab;      // 1 = regularised horseshoe slab
  real<lower=0> slab_scale;            // Slab scale (used only if use_slab = 1)
  real<lower=0> slab_df;               // Slab degrees of freedom (used only if use_slab = 1)
  real<lower=0> sigma_min;             // Lower bound for sigma
  real<lower=sigma_min> sigma_max;     // Upper bound for sigma
}

transformed data {
  int n_omega = n - 1;
}

parameters {
  real mu1;
  real<lower=-0.95, upper=0.95> phi;   // Truncated away from the random-walk boundary (Section 3.4)
  real<lower=sigma_min, upper=sigma_max> sigma;

  real<lower=0> tau_omega;
  vector<lower=0>[n_omega] lambda_omega;
  vector<lower=0>[n_omega] eta_omega;
  vector[n_omega] z_omega;

  real<lower=0> tau_psi;
  vector<lower=0>[n] lambda_psi;
  vector<lower=0>[n] eta_psi;
  vector[n] z_psi;

  real<lower=0> c2_tilde;              // Slab auxiliary (inert when use_slab = 0)
}

transformed parameters {
  vector[n] mu;
  vector[n] psi;
  vector[n_omega] omega;
  vector<lower=0>[n] nu;               // Outlier scale nu_t (Eq. 6)

  {
    real c2 = square(slab_scale) * c2_tilde;
    for (t in 1:n_omega) {
      real lam = kappa[t + 1] * eta_omega[t] * lambda_omega[t];
      real lam_eff = use_slab == 1 ? sqrt(c2 * square(lam) / (c2 + square(tau_omega * lam))) : lam;
      omega[t] = sigma * tau_omega * lam_eff * z_omega[t];
    }
    for (t in 1:n) {
      real lam_p = eta_psi[t] * lambda_psi[t];
      real lam_p_eff = use_slab == 1 ? sqrt(c2 * square(lam_p) / (c2 + square(tau_psi * lam_p))) : lam_p;
      nu[t] = tau_psi * lam_p_eff;
      psi[t] = nu[t] * z_psi[t];
    }
  }

  mu[1] = mu1;
  for (t in 2:n) {
    mu[t] = mu[t - 1] + omega[t - 1];
  }
}

model {
  // (phi + 1)/2 ~ Beta(a_phi, b_phi), truncated by the parameter bounds
  target += beta_lpdf((phi + 1.0) / 2.0 | a_phi, b_phi);
  // Jeffreys prior pi(sigma^2) propto 1/sigma^2  =>  pi(sigma) propto 1/sigma
  target += -log(sigma);

  mu1 ~ normal(0, mu1_sd);

  tau_omega ~ cauchy(0, s_tau);
  tau_psi ~ cauchy(0, 1.0 / sqrt(n));

  eta_omega ~ cauchy(0, 1);
  lambda_omega ~ cauchy(0, 1);
  z_omega ~ std_normal();

  eta_psi ~ cauchy(0, 1);
  lambda_psi ~ cauchy(0, 1);
  z_psi ~ std_normal();

  c2_tilde ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);

  // Prewhitened likelihood: P_phi (y - mu - psi) ~ N(0, sigma^2 I_n) (Section 3.1)
  {
    vector[n] e = y - mu - psi;
    vector[n] u;
    u[1] = sqrt(1.0 - square(phi)) * e[1];
    u[2:n] = e[2:n] - phi * e[1:(n - 1)];
    u ~ normal(0, sigma);
    // log |det P_phi| = 0.5 * log(1 - phi^2)
    target += 0.5 * log1m(square(phi));
  }
}

generated quantities {
  vector<lower=0, upper=1>[n] o_t;     // Outlier score nu_t^2 / (nu_t^2 + sigma^2) (Eq. 8)
  real sigma_marg = sigma / sqrt(1.0 - square(phi));
  for (t in 1:n) {
    o_t[t] = square(nu[t]) / (square(nu[t]) + square(sigma));
  }
}
