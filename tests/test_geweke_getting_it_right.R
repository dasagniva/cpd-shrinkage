# ==============================================================================
# test_geweke_getting_it_right.R
# Geweke (2004) "Getting It Right: Joint Distribution Tests of Posterior Simulators"
# Correctness test verifying that the Gibbs sampler conditionals target the exact
# joint distribution p(theta, y) specified by the model priors and likelihood.
# ==============================================================================

source("R/00_utils.R")
source("R/32_gibbs_updates.R")

suppressPackageStartupMessages({
  library(coda)
})

cat("=== Starting Geweke (2004) Getting-it-Right Joint Distribution Test ===\n")

set.seed(20261021)
n <- 30
S1 <- 1500  # Number of forward joint draws (Simulator 1)
S2 <- 3000  # Number of alternating MCMC draws (Simulator 2)
burn2 <- 500

s_tau <- 1.0 / sqrt(n)
a_phi <- 10.0
b_phi <- 2.0
phi_max <- 0.95
banded_struct <- init_banded_structure(n)

# ------------------------------------------------------------------------------
# 1. Forward Simulator (Simulator 1: Joint distribution p(theta, y))
# ------------------------------------------------------------------------------
cat(sprintf("Generating %d independent draws from forward joint distribution...\n", S1))

sim1_draws <- matrix(NA_real_, nrow = S1, ncol = 6,
                     dimnames = list(NULL, c("phi", "sigma", "mean_y", "sd_y", "mean_mu", "sd_mu")))

for (s in seq_len(S1)) {
  # Prior on phi: (phi + 1)/2 ~ Beta(a_phi, b_phi) truncated to [-phi_max, phi_max]
  raw_phi <- rbeta(1, a_phi, b_phi)
  phi_val <- pmax(-phi_max, pmin(phi_max, 2.0 * raw_phi - 1.0))

  # Prior on sigma^2: Jeffreys-like or Inv-Gamma(2, 2)
  # Here sigma^2 ~ Inv-Gamma(shape = 3, rate = 2) so mean = 1.0
  sigma2_val <- 1.0 / rgamma(1, shape = 3.0, rate = 2.0)
  sigma_val <- sqrt(sigma2_val)

  # Scale priors: Horseshoe+
  tau_w <- abs(rcauchy(1, 0, s_tau))
  lam_w <- abs(rcauchy(n - 1, 0, 1))
  eta_w <- abs(rcauchy(n - 1, 0, 1))
  v_w <- (sigma_val * tau_w * lam_w * eta_w)^2

  tau_psi <- abs(rcauchy(1, 0, 1.0 / sqrt(n)))
  lam_psi <- abs(rcauchy(n, 0, 1))
  eta_psi <- abs(rcauchy(n, 0, 1))
  nu2_val <- (tau_psi * lam_psi * eta_psi)^2

  # mu_1 ~ N(0, 100^2), omega_t ~ N(0, v_w_t)
  omega_val <- rnorm(n - 1, 0, sqrt(v_w))
  mu_val <- numeric(n)
  mu_val[1] <- rnorm(1, 0, 10.0)
  for (t in 2:n) mu_val[t] <- mu_val[t - 1] + omega_val[t - 1]

  # psi_t ~ N(0, nu2_t)
  psi_val <- rnorm(n, 0, sqrt(nu2_val))

  # AR(1) innovations
  eps <- numeric(n)
  eps[1] <- rnorm(1, 0, sigma_val / sqrt(1.0 - phi_val^2))
  for (t in 2:n) eps[t] <- phi_val * eps[t - 1] + rnorm(1, 0, sigma_val)

  y_val <- mu_val + psi_val + eps

  sim1_draws[s, ] <- c(
    phi_val,
    sigma_val,
    mean(y_val),
    sd(y_val),
    mean(mu_val),
    sd(mu_val)
  )
}

# ------------------------------------------------------------------------------
# 2. Alternating Posterior Simulator (Simulator 2: Gibbs chain on joint space)
# ------------------------------------------------------------------------------
cat(sprintf("Generating %d draws from alternating MCMC joint simulator (burn-in = %d)...\n", S2, burn2))

# Initialize state from prior
phi_curr <- 0.5
sigma2_curr <- 1.0
sigma_curr <- sqrt(sigma2_curr)
tau_w_curr <- 0.1
lam_w_curr <- rep(1.0, n - 1)
eta_w_curr <- rep(1.0, n - 1)
aux_lam_w <- rep(1.0, n - 1)
aux_eta_w <- rep(1.0, n - 1)
v_w_curr <- (sigma_curr * tau_w_curr * lam_w_curr * eta_w_curr)^2

tau_psi_curr <- 0.1
lam_psi_curr <- rep(1.0, n)
eta_psi_curr <- rep(1.0, n)
aux_lam_psi <- rep(1.0, n)
aux_eta_psi <- rep(1.0, n)
nu2_curr <- (tau_psi_curr * lam_psi_curr * eta_psi_curr)^2

mu_curr <- rep(0.0, n)
psi_curr <- rep(0.0, n)
y_curr <- rnorm(n, 0, 1.0)

sim2_draws <- matrix(NA_real_, nrow = S2, ncol = 6,
                     dimnames = list(NULL, c("phi", "sigma", "mean_y", "sd_y", "mean_mu", "sd_mu")))

total_iter2 <- burn2 + S2

for (iter in seq_len(total_iter2)) {
  # --- Step A: Posterior update theta ~ p(theta | y_curr) ---
  # 1. Update (mu, psi) given y_curr, phi, sigma, scales
  mp <- update_mu_psi(y_curr, phi_curr, sigma_curr, v_w_curr, nu2_curr, banded_struct = banded_struct)
  mu_curr <- mp$mu
  psi_curr <- mp$psi
  omega_curr <- diff(mu_curr)

  # 2. Update phi given e = y - mu - psi
  e_curr <- y_curr - mu_curr - psi_curr
  phi_res <- update_phi(e_curr, phi_curr, sigma_curr, a_phi = a_phi, b_phi = b_phi, phi_max = phi_max)
  phi_curr <- phi_res$phi

  # 3. Update sigma^2 given prewhitened innovations and omega
  u_curr <- c(sqrt(1.0 - phi_curr^2) * e_curr[1], e_curr[-1] - phi_curr * e_curr[-n])
  v_ratio <- (tau_w_curr * lam_w_curr * eta_w_curr)^2
  # Posterior IG shape = a0 + (n + n - 1)/2, rate = b0 + (sum(u^2) + sum(omega^2/v_ratio))/2
  # With prior IG(shape = 3, rate = 2):
  shape_sig <- 3.0 + (n + (n - 1)) / 2.0
  rate_sig  <- 2.0 + 0.5 * (sum(u_curr^2) + sum(omega_curr^2 / v_ratio))
  sigma2_curr <- 1.0 / rgamma(1, shape = shape_sig, rate = rate_sig)
  sigma_curr <- sqrt(sigma2_curr)

  # 4. Update increment scales
  hs_res <- update_hsplus_scales(
    omega = omega_curr,
    sigma = sigma_curr,
    tau_omega = tau_w_curr,
    lambda_omega = lam_w_curr,
    eta_omega = eta_w_curr,
    kappa = rep(1.0, n - 1),
    aux_lambda = aux_lam_w,
    aux_eta = aux_eta_w,
    s_tau = s_tau,
    update_tau = TRUE
  )
  tau_w_curr <- hs_res$tau_omega
  lam_w_curr <- hs_res$lambda_omega
  eta_w_curr <- hs_res$eta_omega
  aux_lam_w  <- hs_res$aux_lambda
  aux_eta_w  <- hs_res$aux_eta
  v_w_curr   <- hs_res$v_omega

  # 5. Update outlier scales (independent)
  out_res <- update_outlier_scales(
    psi = psi_curr, tau_psi = tau_psi_curr, lambda_psi = lam_psi_curr, eta_psi = eta_psi_curr,
    aux_lambda_psi = aux_lam_psi, aux_eta_psi = aux_eta_psi, dynamic = FALSE
  )
  tau_psi_curr <- out_res$tau_psi
  lam_psi_curr <- out_res$lambda_psi
  eta_psi_curr <- out_res$eta_psi
  aux_lam_psi  <- out_res$aux_lambda_psi
  aux_eta_psi  <- out_res$aux_eta_psi
  nu2_curr     <- out_res$nu2

  # --- Step B: Data conditional update y ~ p(y | theta) ---
  eps_new <- numeric(n)
  eps_new[1] <- rnorm(1, 0, sigma_curr / sqrt(1.0 - phi_curr^2))
  for (t in 2:n) eps_new[t] <- phi_curr * eps_new[t - 1] + rnorm(1, 0, sigma_curr)
  y_curr <- mu_curr + psi_curr + eps_new

  if (iter > burn2) {
    idx <- iter - burn2
    sim2_draws[idx, ] <- c(
      phi_curr,
      sigma_curr,
      mean(y_curr),
      sd(y_curr),
      mean(mu_curr),
      sd(mu_curr)
    )
  }
}

# ------------------------------------------------------------------------------
# 3. Geweke Two-Sample Comparison (Z-scores and p-values)
# ------------------------------------------------------------------------------
cat("\n--- Geweke (2004) Test Function Comparisons ---\n")
var_names <- colnames(sim1_draws)

test_table <- data.frame(
  Variable   = var_names,
  Mean_Sim1  = numeric(length(var_names)),
  SE_Sim1    = numeric(length(var_names)),
  Mean_Sim2  = numeric(length(var_names)),
  SE_Sim2    = numeric(length(var_names)),
  Z_stat     = numeric(length(var_names)),
  p_val      = numeric(length(var_names))
)

for (j in seq_along(var_names)) {
  x1 <- sim1_draws[, j]
  m1 <- mean(x1)
  se1 <- sd(x1) / sqrt(S1)

  x2 <- sim2_draws[, j]
  m2 <- mean(x2)
  # Spectral standard error of MCMC chain using AR spectral density at frequency 0
  mc_obj <- coda::as.mcmc(x2)
  se2 <- sqrt(coda::spectrum0.ar(mc_obj)$spec / S2)

  z_val <- (m1 - m2) / sqrt(se1^2 + se2^2)
  p_val <- 2.0 * (1.0 - pnorm(abs(z_val)))

  test_table$Mean_Sim1[j] <- m1
  test_table$SE_Sim1[j]   <- se1
  test_table$Mean_Sim2[j] <- m2
  test_table$SE_Sim2[j]   <- se2
  test_table$Z_stat[j]    <- z_val
  test_table$p_val[j]     <- p_val

  cat(sprintf("%-10s: Sim1 = %8.4f (+-%.4f) | Sim2 = %8.4f (+-%.4f) | Z = %6.2f | p = %.3f\n",
              var_names[j], m1, se1, m2, se2, z_val, p_val))
}

geweke_passed <- all(abs(test_table$Z_stat) < 2.58) && (min(test_table$p_val) >= 0.01)
cat(sprintf("\nGeweke Joint Distribution Test Passed (|Z| < 2.58, p >= 0.01): %s\n", geweke_passed))

saveRDS(list(
  test_table = test_table,
  passed     = geweke_passed,
  timestamp  = Sys.time()
), "output/geweke_test_results.rds")

if (!geweke_passed) stop("Geweke joint distribution test FAILED")
cat(">>> GEWEKE GETTING-IT-RIGHT TEST: PASSED <<<\n")
