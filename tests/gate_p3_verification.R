# ==============================================================================
# gate_p3_verification.R
# Acceptance Gate Verification for Phase P3:
# 1. Parameter recovery on 20 simulated series (n <= 300)
# 2. Simulation-Based Calibration (SBC) rank histogram checks on the reference model
# ==============================================================================

source("R/00_utils.R")
source("R/10_dgp.R")
source("R/20_frontend_nsp.R")
source("R/21_frontend_capa.R")
source("R/22_prior_map.R")

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(SBC)
})

log_msg("=== Starting Phase P3 Acceptance Gate Verification ===")

# Compile Stan reference model
stan_file <- "stan/hsplus_ar1_ref.stan"
log_msg(sprintf("Compiling Stan reference model: %s ...", stan_file))
mod <- cmdstan_model(stan_file)
log_msg("Compilation complete.")

# ------------------------------------------------------------------------------
# Part 1: Parameter Recovery on 20 Simulated Series
# ------------------------------------------------------------------------------
n_recovery <- 20
n_len <- 100
phi_grid <- seq(0.2, 0.8, length.out = n_recovery)
sigma_grid <- seq(0.6, 1.4, length.out = n_recovery)

log_msg(sprintf("Running Part 1: Parameter recovery on %d simulated series (n = %d)...", n_recovery, n_len))

recovery_records <- data.frame(
  iter          = integer(n_recovery),
  true_phi      = numeric(n_recovery),
  est_phi       = numeric(n_recovery),
  phi_q05       = numeric(n_recovery),
  phi_q95       = numeric(n_recovery),
  phi_covered   = logical(n_recovery),
  true_sigma    = numeric(n_recovery),
  est_sigma     = numeric(n_recovery),
  sigma_q05     = numeric(n_recovery),
  sigma_q95     = numeric(n_recovery),
  sigma_covered = logical(n_recovery),
  max_rhat      = numeric(n_recovery),
  min_ess       = numeric(n_recovery)
)

for (i in seq_len(n_recovery)) {
  true_phi <- phi_grid[i]
  true_sigma <- sigma_grid[i]

  sim <- simulate_series(
    n = n_len,
    phi = true_phi,
    sigma = true_sigma,
    noise = "gaussian",
    k0 = 2,
    jump = 2.5,
    outlier_regime = "isolated",
    p_outlier = 0.02,
    collective = 0,
    seed = 20263000L + i
  )

  # Compute front-end prior maps
  nsp_res <- run_frontend_nsp(sim$y, alpha = 0.05)
  capa_res <- run_frontend_capa(sim$y, prewhiten = TRUE)
  pmap <- build_prior_maps(n = n_len, nsp_intervals = nsp_res$intervals, capa_results = capa_res, kappa_in = 2.0)

  stan_data <- list(
    n          = n_len,
    y          = as.vector(sim$y),
    kappa      = as.vector(pmap$kappa),
    s_tau      = 1.0 / sqrt(n_len),
    a_phi      = 10.0,
    b_phi      = 2.0,
    mu1_sd     = 100.0,
    use_slab   = 0L,
    slab_scale = 5.0,
    slab_df    = 4.0,
    sigma_min  = 1e-4,
    sigma_max  = 1e4
  )

  fit <- mod$sample(
    data            = stan_data,
    seed            = 20263000L + i,
    chains          = 2,
    parallel_chains = 2,
    iter_warmup     = 500,
    iter_sampling   = 500,
    adapt_delta     = 0.99,
    max_treedepth   = 12,
    init            = 0.5,
    show_messages   = FALSE,
    refresh         = 0
  )

  summ <- fit$summary(c("phi", "sigma", "tau_omega", "tau_psi"))

  phi_row <- summ[summ$variable == "phi", ]
  sigma_row <- summ[summ$variable == "sigma", ]

  p_cov <- (true_phi >= phi_row$q5 && true_phi <= phi_row$q95)
  s_cov <- (true_sigma >= sigma_row$q5 && true_sigma <= sigma_row$q95)

  recovery_records$iter[i]          <- i
  recovery_records$true_phi[i]      <- true_phi
  recovery_records$est_phi[i]       <- phi_row$mean
  recovery_records$phi_q05[i]       <- phi_row$q5
  recovery_records$phi_q95[i]       <- phi_row$q95
  recovery_records$phi_covered[i]   <- p_cov

  recovery_records$true_sigma[i]    <- true_sigma
  recovery_records$est_sigma[i]     <- sigma_row$mean
  recovery_records$sigma_q05[i]     <- sigma_row$q5
  recovery_records$sigma_q95[i]     <- sigma_row$q95
  recovery_records$sigma_covered[i] <- s_cov

  recovery_records$max_rhat[i]      <- max(summ$rhat, na.rm = TRUE)
  recovery_records$min_ess[i]       <- min(summ$ess_bulk, na.rm = TRUE)

  cat(sprintf("[%02d/%d] phi: true=%.2f est=%.2f [%.2f, %.2f] cov=%s | sigma: true=%.2f est=%.2f [%.2f, %.2f] cov=%s | Rhat=%.2f ESS=%.0f\n",
              i, n_recovery,
              true_phi, phi_row$mean, phi_row$q5, phi_row$q95, ifelse(p_cov, "OK", "NO"),
              true_sigma, sigma_row$mean, sigma_row$q5, sigma_row$q95, ifelse(s_cov, "OK", "NO"),
              max(summ$rhat, na.rm = TRUE), min(summ$ess_bulk, na.rm = TRUE)))
}

phi_coverage <- mean(recovery_records$phi_covered)
sigma_coverage <- mean(recovery_records$sigma_covered)
mean_phi_mae <- mean(abs(recovery_records$est_phi - recovery_records$true_phi))
mean_sigma_mae <- mean(abs(recovery_records$est_sigma - recovery_records$true_sigma))

cat("\n--- Part 1 Recovery Summary ---\n")
cat(sprintf("phi 90%% Credible Interval Coverage  : %.1f%% (MAE = %.3f)\n", 100 * phi_coverage, mean_phi_mae))
cat(sprintf("sigma 90%% Credible Interval Coverage: %.1f%% (MAE = %.3f)\n", 100 * sigma_coverage, mean_sigma_mae))
cat(sprintf("Median R-hat across all runs         : %.3f\n", median(recovery_records$max_rhat)))
cat(sprintf("85th percentile R-hat across runs    : %.3f\n", quantile(recovery_records$max_rhat, 0.85)))
cat(sprintf("Max R-hat across all runs            : %.3f\n", max(recovery_records$max_rhat)))

# ------------------------------------------------------------------------------
# Part 2: Simulation-Based Calibration (SBC) Rank Histogram Checks
# ------------------------------------------------------------------------------
log_msg("\nRunning Part 2: Simulation-Based Calibration (SBC)...")
n_sbc <- 30
n_sbc_len <- 50
sbc_ranks <- list(phi = integer(n_sbc), sigma = integer(n_sbc), tau_omega = integer(n_sbc))

for (k in seq_len(n_sbc)) {
  set.seed(20264000L + k)

  # Draw parameters from prior
  # 1. phi: (phi + 1)/2 ~ Beta(10, 2)
  phi_raw <- rbeta(1, 10.0, 2.0)
  prior_phi <- pmax(-0.95, pmin(0.95, 2.0 * phi_raw - 1.0))

  # 2. sigma: log-uniform draw in [0.5, 2.0]
  prior_sigma <- exp(runif(1, log(0.5), log(2.0)))

  # 3. Global scales
  s_tau <- 1.0 / sqrt(n_sbc_len)
  prior_tau_omega <- abs(rcauchy(1, 0, s_tau))
  prior_tau_psi <- abs(rcauchy(1, 0, 1.0 / sqrt(n_sbc_len)))

  # 4. Horseshoe+ hierarchy for increments and outliers (no slab: use_slab = 0)
  eta_omega <- abs(rcauchy(n_sbc_len - 1, 0, 1))
  lambda_omega <- abs(rcauchy(n_sbc_len - 1, 0, 1))
  omega <- prior_sigma * prior_tau_omega * eta_omega * lambda_omega * rnorm(n_sbc_len - 1)

  mu <- numeric(n_sbc_len)
  mu[1] <- rnorm(1, 0, 10.0)
  for (t in 2:n_sbc_len) mu[t] <- mu[t - 1] + omega[t - 1]

  # psi_t ~ N(0, tau_psi^2 eta^2 lambda^2): not scaled by sigma, as in the Stan model
  eta_psi <- abs(rcauchy(n_sbc_len, 0, 1))
  lambda_psi <- abs(rcauchy(n_sbc_len, 0, 1))
  psi <- prior_tau_psi * eta_psi * lambda_psi * rnorm(n_sbc_len)

  # 5. Prewhitened AR(1) innovations
  eps <- numeric(n_sbc_len)
  eps[1] <- rnorm(1, 0, prior_sigma / sqrt(1.0 - prior_phi^2))
  for (t in 2:n_sbc_len) eps[t] <- prior_phi * eps[t - 1] + rnorm(1, 0, prior_sigma)

  y_sbc <- mu + psi + eps

  stan_data_sbc <- list(
    n          = n_sbc_len,
    y          = as.vector(y_sbc),
    kappa      = rep(1.0, n_sbc_len),
    s_tau      = s_tau,
    a_phi      = 10.0,
    b_phi      = 2.0,
    mu1_sd     = 10.0,
    use_slab   = 0L,
    slab_scale = 5.0,
    slab_df    = 4.0,
    sigma_min  = 0.5,
    sigma_max  = 2.0
  )

  fit_sbc <- mod$sample(
    data            = stan_data_sbc,
    seed            = 20264000L + k,
    chains          = 1,
    iter_warmup     = 400,
    iter_sampling   = 400,
    adapt_delta     = 0.99,
    max_treedepth   = 12,
    init            = 0.5,
    show_messages   = FALSE,
    refresh         = 0
  )

  draws <- fit_sbc$draws(c("phi", "sigma", "tau_omega"))
  keep <- seq(1, dim(draws)[1], by = 4) # thin to reduce rank autocorrelation
  phi_draws <- as.vector(draws[keep, , "phi"])
  sigma_draws <- as.vector(draws[keep, , "sigma"])
  tau_draws <- as.vector(draws[keep, , "tau_omega"])

  sbc_ranks$phi[k]       <- sum(phi_draws < prior_phi)
  sbc_ranks$sigma[k]     <- sum(sigma_draws < prior_sigma)
  sbc_ranks$tau_omega[k] <- sum(tau_draws < prior_tau_omega)

  cat(sprintf("[SBC %02d/%d] phi: true=%.2f est=%.2f rk=%d | sigma: true=%.2f est=%.2f rk=%d | tau: true=%.2f est=%.2f rk=%d\n",
              k, n_sbc,
              prior_phi, mean(phi_draws), sbc_ranks$phi[k],
              prior_sigma, mean(sigma_draws), sbc_ranks$sigma[k],
              prior_tau_omega, mean(tau_draws), sbc_ranks$tau_omega[k]))
}

# Uniformity tests for ranks (normalized to [0, 1])
S_draws <- 100
p_ks_phi <- ks.test((sbc_ranks$phi + 0.5) / (S_draws + 1), "punif")$p.value
p_ks_sigma <- ks.test((sbc_ranks$sigma + 0.5) / (S_draws + 1), "punif")$p.value
p_ks_tau <- ks.test((sbc_ranks$tau_omega + 0.5) / (S_draws + 1), "punif")$p.value

cat("\n--- Part 2 SBC Uniformity Test Results ---\n")
cat(sprintf("phi       rank uniformity KS test p-value: %.3f\n", p_ks_phi))
cat(sprintf("sigma     rank uniformity KS test p-value: %.3f\n", p_ks_sigma))
cat(sprintf("tau_omega rank uniformity KS test p-value: %.3f\n", p_ks_tau))

# ------------------------------------------------------------------------------
# Gate Conclusion
# ------------------------------------------------------------------------------
cat("\n=== Phase P3 Gate Conclusion ===\n")
recovery_pass <- (phi_coverage >= 0.75 && sigma_coverage >= 0.75 && median(recovery_records$max_rhat) < 1.05 && quantile(recovery_records$max_rhat, 0.85) < 1.05)
sbc_pass <- (p_ks_phi >= 0.01 && p_ks_sigma >= 0.01 && p_ks_tau >= 0.01)

if (recovery_pass && sbc_pass) {
  cat("[PASS] Parameter recovery on 20 series and SBC rank tests passed.\n")
} else {
  cat(sprintf("[FAIL] recovery_pass = %s, sbc_pass = %s\n", recovery_pass, sbc_pass))
}

saveRDS(list(
  recovery_records = recovery_records,
  sbc_ranks = sbc_ranks,
  p_ks = list(phi = p_ks_phi, sigma = p_ks_sigma, tau_omega = p_ks_tau),
  recovery_pass = recovery_pass,
  sbc_pass = sbc_pass
), file = "output/gate_p3_results.rds")
log_msg("Saved results to output/gate_p3_results.rds")
if (!(recovery_pass && sbc_pass)) stop("Phase P3 acceptance gate FAILED")
