// =============================================================================
// stan/hsplus_ar1_ref.stan
// Reference Model: Significance-Guided Global-Local Shrinkage under AR(1) Noise
// Implements Section 3 (Equations 1-6) with independent outlier scales (rho=0, delta=0).
// Non-centred parameterisations for increments omega and outliers psi.
// Prewhitened likelihood representation.
// =============================================================================

data {
  int<lower=1> n;                      // Length of time series (n <= 300 for validation)
  vector[n] y;                         // Observed time series
  vector<lower=0>[n] kappa;            // NSP-informed local scale multiplier (Eq. 5)
  real<lower=0> s_tau;                 // Prior scale for tau_omega (default 1 / sqrt(n))
  real<lower=0> a_phi;                 // Prior beta shape 1 for (phi + 1)/2 (default 10)
  real<lower=0> b_phi;                 // Prior beta shape 2 for (phi + 1)/2 (default 2)
  real<lower=0> slab_scale;            // Regularised horseshoe slab scale (default 5.0)
  real<lower=0> slab_df;               // Slab degrees of freedom (default 4.0)
  real<lower=0> sigma_min;             // Lower bound for sigma (default 1e-4)
  real<lower=sigma_min> sigma_max;     // Upper bound for sigma (default 1e4)
}

transformed data {
  int n_omega = n - 1;                 // Number of increments
}

parameters {
  // Base mean level
  real mu1;

  // AR(1) persistence truncated away from random-walk boundary |phi| <= 0.95 (Section 3.4)
  real<lower=-0.95, upper=0.95> phi;

  // Innovation scale (Jeffreys prior pi(sigma^2) proportional to 1/sigma^2 => pi(sigma) proportional to 1/sigma)
  real<lower=sigma_min, upper=sigma_max> sigma;

  // Global scale for trend increments
  real<lower=0> tau_omega;

  // Horseshoe+ auxiliary scales for trend increments omega_t (t = 2, ..., n)
  vector<lower=0>[n_omega] lambda_omega;  // Local half-Cauchy scale
  vector<lower=0>[n_omega] eta_omega;     // Secondary half-Cauchy scale (Horseshoe+)
  vector[n_omega] z_omega;                // Standard normal primitives for non-centred draw

  // Global scale for outlier process psi
  real<lower=0> tau_psi;

  // Horseshoe scales for outlier process psi_t (t = 1, ..., n; independent scales)
  vector<lower=0>[n] lambda_psi;          // Local outlier half-Cauchy scale
  vector<lower=0>[n] eta_psi;             // Secondary half-Cauchy scale
  vector[n] z_psi;                        // Standard normal primitives for non-centred draw

  // Regularised horseshoe slab auxiliary
  real<lower=0> c2_tilde;
}

transformed parameters {
  vector[n] mu;
  vector[n] psi;
  vector[n_omega] omega;
  real<lower=0> c2;

  // Regularised horseshoe slab (Piironen & Vehtari 2017) to prevent divergences
  c2 = slab_scale * slab_scale * c2_tilde;

  // Construct increments omega_t with Horseshoe+ local shrinkage (Eq. 3-4)
  {
    vector[n_omega] eff_lambda;
    for (t in 1:n_omega) {
      // Horseshoe+ scale modulated by kappa_t (Eq. 4-5)
      real lam = fmin(kappa[t + 1] * eta_omega[t] * lambda_omega[t], 1e4);
      real tau_lam = fmin(tau_omega * lam, 1e4);
      // Regularised shrinkage scale
      eff_lambda[t] = sqrt(fmax(0.0, (c2 * square(lam)) / (c2 + square(tau_lam) + 1e-12)));
    }
    // omega_t | lambda, tau, sigma ~ N(0, sigma^2 * tau_omega^2 * eff_lambda^2)
    omega = sigma * tau_omega * (eff_lambda .* z_omega);
  }

  // Piecewise constant mean trajectory mu_t = mu_{t-1} + omega_{t-1}
  mu[1] = mu1;
  for (t in 2:n) {
    mu[t] = mu[t - 1] + omega[t - 1];
  }

  // Construct sparse outliers psi_t with independent Horseshoe+ scales (Eq. 6 with rho=0, delta=0)
  {
    for (t in 1:n) {
      real lam_p = fmin(eta_psi[t] * lambda_psi[t], 1e4);
      real tau_lam_p = fmin(tau_psi * lam_p, 1e4);
      real eff_lam_p = sqrt(fmax(0.0, (c2 * square(lam_p)) / (c2 + square(tau_lam_p) + 1e-12)));
      psi[t] = tau_psi * eff_lam_p * z_psi[t];
    }
  }
}

model {
  // Priors
  // 1. phi: (phi + 1)/2 ~ Beta(a_phi, b_phi) truncated to |phi| <= 0.95 (Section 3.4)
  target += beta_lpdf((phi + 1.0) / 2.0 | a_phi, b_phi);

  // 2. sigma: Jeffreys prior pi(sigma^2) proportional to 1/sigma^2 => target += -log(sigma)
  target += -log(sigma);

  // 3. Base mean level prior
  mu1 ~ normal(0, 10.0 * sigma);

  // 4. Global scales: tau_omega ~ C+(0, s_tau), tau_psi ~ C+(0, 1/sqrt(n))
  tau_omega ~ cauchy(0, s_tau);
  tau_psi ~ cauchy(0, 1.0 / sqrt(n));

  // 5. Horseshoe+ hierarchy for omega
  eta_omega ~ cauchy(0, 1);
  lambda_omega ~ cauchy(0, 1);
  z_omega ~ std_normal();

  // 6. Horseshoe+ hierarchy for psi
  eta_psi ~ cauchy(0, 1);
  lambda_psi ~ cauchy(0, 1);
  z_psi ~ std_normal();

  // 7. Regularised slab parameter
  c2_tilde ~ inv_gamma(0.5 * slab_df, 0.5 * slab_df);

  // Likelihood: prewhitened form P_phi * (y - mu - psi) ~ N(0, sigma^2 * I_n) (Section 3.1)
  {
    vector[n] e = y - mu - psi;
    vector[n] u;

    // First observation from stationary distribution: Var(eps_1) = sigma^2 / (1 - phi^2)
    u[1] = sqrt(1.0 - square(phi)) * e[1];
    // AR(1) innovations u_t = e_t - phi * e_{t-1} ~ N(0, sigma^2)
    u[2:n] = e[2:n] - phi * e[1:(n - 1)];

    // Standard normal likelihood for innovations
    u ~ normal(0, sigma);

    // Jacobian adjustment for prewhitening transformation: log |det(P_phi)| = 0.5 * log(1 - phi^2)
    target += 0.5 * log1m(square(phi));
  }
}

generated quantities {
  // Posterior scores from Section 3.5
  vector<lower=0>[n] o_t;              // Outlier score nu_t^2 / (nu_t^2 + sigma^2) (Eq. 8)
  vector[n_omega] abs_omega = abs(omega);
  vector[n] abs_psi = abs(psi);
  real sigma_marg = sigma / sqrt(1.0 - square(phi)); // Marginal noise SD

  for (t in 1:n) {
    real nu2 = square(psi[t] / (z_psi[t] != 0 ? z_psi[t] : 1e-6));
    o_t[t] = nu2 / (nu2 + square(sigma));
  }
}
