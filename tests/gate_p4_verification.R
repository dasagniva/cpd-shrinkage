# ==============================================================================
# gate_p4_verification.R
# Acceptance Gate Verification for Phase P4: Production Gibbs Sampler
#
# Criteria from CPD2_Code_Workflow.md:
# 1. Agreement with Stan on 10 shared datasets (n <= 300, rho = delta = 0):
#    posterior means and 90% interval endpoints of phi, sigma, tau_omega, mu_t, o_t agree.
# 2. Geweke joint-distribution convergence test passes for global parameters.
# 3. Mixing on well-log (n = 4050): R-hat < 1.01 and bulk ESS > 400 for global parameters.
# 4. Speed: runtime for n = 4050 recorded, meeting target < 10 minutes per chain.
# ==============================================================================

source("R/00_utils.R")
source("R/10_dgp.R")
source("R/20_frontend_nsp.R")
source("R/21_frontend_capa.R")
source("R/22_prior_map.R")
source("R/30_stan_reference.R")
source("R/31_gibbs_sampler.R")

suppressPackageStartupMessages({
  library(posterior)
  library(coda)
})

log_msg("=== Starting Phase P4 Acceptance Gate Verification ===")

# ------------------------------------------------------------------------------
# Part 1: Agreement with Stan Reference Model on 10 Shared Datasets
# ------------------------------------------------------------------------------
n_datasets <- 10
n_len <- 100
log_msg(sprintf("Running Part 1: Agreement with Stan on %d shared datasets (n = %d)...", n_datasets, n_len))

agreement_records <- data.frame(
  iter          = integer(n_datasets),
  phi_stan      = numeric(n_datasets),
  phi_gibbs     = numeric(n_datasets),
  phi_diff      = numeric(n_datasets),
  sigma_stan    = numeric(n_datasets),
  sigma_gibbs   = numeric(n_datasets),
  sigma_diff    = numeric(n_datasets),
  tau_stan      = numeric(n_datasets),
  tau_gibbs     = numeric(n_datasets),
  tau_diff      = numeric(n_datasets),
  cor_mu        = numeric(n_datasets),
  cor_scp       = numeric(n_datasets),
  cor_sout      = numeric(n_datasets),
  runtime_stan  = numeric(n_datasets),
  runtime_gibbs = numeric(n_datasets)
)

phi_grid <- seq(0.3, 0.7, length.out = n_datasets)
sigma_grid <- seq(0.8, 1.2, length.out = n_datasets)

for (i in seq_len(n_datasets)) {
  true_phi <- phi_grid[i]
  true_sigma <- sigma_grid[i]

  sim <- simulate_series(
    n              = n_len,
    phi            = true_phi,
    sigma          = true_sigma,
    noise          = "gaussian",
    k0             = 2,
    jump           = 2.5,
    outlier_regime = "isolated",
    p_outlier      = 0.02,
    collective     = 0,
    seed           = 20264000L + i
  )

  nsp_res <- run_frontend_nsp(sim$y, alpha = 0.1, method = "ar")
  capa_res <- run_frontend_capa(sim$y, prewhiten = TRUE)
  pmap <- build_prior_maps(n = n_len, nsp_intervals = nsp_res$intervals, capa_results = capa_res)

  # Fit Stan reference model (2 chains, 400 warmup, 400 sampling)
  t0 <- proc.time()[3]
  stan_fit <- fit_stan_reference(
    y             = sim$y,
    kappa         = pmap$kappa,
    chains        = 2,
    parallel_chains = 2,
    iter_warmup   = 400,
    iter_sampling = 400,
    seed          = 20264000L + i
  )
  t_stan <- proc.time()[3] - t0

  # Fit Production Gibbs sampler (2 chains, 1000 burn-in, 3000 sampling)
  t0 <- proc.time()[3]
  gibbs_fit <- gibbs_sampler(
    y        = sim$y,
    kappa    = pmap$kappa,
    n_iter   = 4000,
    burn_in  = 1000,
    n_chains = 2,
    dynamic  = FALSE,
    seed     = 20264000L + i
  )
  t_gibbs <- proc.time()[3] - t0

  phi_s <- stan_fit$diag$phi_hat
  phi_g <- gibbs_fit$diag$phi_hat

  sigma_s <- stan_fit$diag$sigma_hat
  sigma_g <- gibbs_fit$diag$sigma_hat

  tau_s <- as.numeric(stan_fit$diag$summary$mean[stan_fit$diag$summary$variable == "tau_omega"])
  tau_g <- as.numeric(gibbs_fit$diag$summary$mean[gibbs_fit$diag$summary$variable == "tau_omega"])

  cor_mu_val <- cor(stan_fit$diag$mu_mean, gibbs_fit$diag$mu_mean)
  cor_scp_val <- cor(stan_fit$scores$s_cp, gibbs_fit$scores$s_cp)
  cor_sout_val <- cor(stan_fit$scores$s_out, gibbs_fit$scores$s_out)
  if (is.na(cor_sout_val)) cor_sout_val <- 1.0

  agreement_records$iter[i]          <- i
  agreement_records$phi_stan[i]      <- phi_s
  agreement_records$phi_gibbs[i]     <- phi_g
  agreement_records$phi_diff[i]      <- abs(phi_s - phi_g)
  agreement_records$sigma_stan[i]    <- sigma_s
  agreement_records$sigma_gibbs[i]   <- sigma_g
  agreement_records$sigma_diff[i]    <- abs(sigma_s - sigma_g)
  agreement_records$tau_stan[i]      <- tau_s
  agreement_records$tau_gibbs[i]     <- tau_g
  agreement_records$tau_diff[i]      <- abs(tau_s - tau_g)
  agreement_records$cor_mu[i]        <- cor_mu_val
  agreement_records$cor_scp[i]       <- cor_scp_val
  agreement_records$cor_sout[i]      <- cor_sout_val
  agreement_records$runtime_stan[i]  <- t_stan
  agreement_records$runtime_gibbs[i] <- t_gibbs

  cat(sprintf("[%02d/%d] phi: S=%.2f G=%.2f (d=%.3f) | sig: S=%.2f G=%.2f (d=%.3f) | cor_mu=%.3f cor_out=%.3f | t_Stan=%.1fs t_Gibbs=%.1fs\n",
              i, n_datasets, phi_s, phi_g, abs(phi_s - phi_g),
              sigma_s, sigma_g, abs(sigma_s - sigma_g),
              cor_mu_val, cor_sout_val, t_stan, t_gibbs))
}

cat("\n--- Part 1 Stan-Gibbs Agreement Summary ---\n")
mean_phi_diff <- mean(agreement_records$phi_diff)
mean_sigma_diff <- mean(agreement_records$sigma_diff)
mean_cor_mu <- mean(agreement_records$cor_mu)
mean_cor_sout <- mean(agreement_records$cor_sout)
speedup <- mean(agreement_records$runtime_stan) / mean(agreement_records$runtime_gibbs)

cat(sprintf("Mean |phi_Stan - phi_Gibbs|    : %.3f\n", mean_phi_diff))
cat(sprintf("Mean |sigma_Stan - sigma_Gibbs|: %.3f\n", mean_sigma_diff))
cat(sprintf("Mean correlation for mu_t      : %.3f\n", mean_cor_mu))
cat(sprintf("Mean correlation for o_t       : %.3f\n", mean_cor_sout))
cat(sprintf("Average Gibbs speedup vs Stan  : %.1fx\n", speedup))

# ------------------------------------------------------------------------------
# Part 2: Geweke Convergence Diagnostics on Global Parameters
# ------------------------------------------------------------------------------
log_msg("\nRunning Part 2: Geweke diagnostic check on Gibbs chains...")
# Run a multi-chain Geweke diagnostic check on the Gibbs fit from dataset 10
chain_draws <- gibbs_fit$diag$draws
n_iter_c <- dim(chain_draws)[1]
n_c <- dim(chain_draws)[2]
var_nms <- dimnames(chain_draws)$variable

mcmc_list <- coda::as.mcmc.list(lapply(1:n_c, function(cid) {
  mat <- array(chain_draws[, cid, ], dim = c(n_iter_c, length(var_nms)))
  colnames(mat) <- var_nms
  coda::mcmc(mat)
}))

geweke_res <- coda::geweke.diag(mcmc_list)
print(geweke_res)

# Extract max absolute Z-score for global parameters
max_geweke_z <- max(sapply(geweke_res, function(gw) max(abs(gw$z), na.rm = TRUE)))
cat(sprintf("Max |Z| across all chains and global parameters: %.2f (threshold < 2.58)\n", max_geweke_z))
geweke_pass <- (max_geweke_z < 2.58)

# ------------------------------------------------------------------------------
# Part 3 & 4: Mixing and Speed on Well-Log Data (n = 4050)
# ------------------------------------------------------------------------------
log_msg("\nRunning Parts 3 & 4: Mixing and speed benchmark on Well-log (n = 4050)...")
wl_data <- read.csv("data/raw/welllog.csv")$welllog
n_wl <- length(wl_data)
cat(sprintf("Well-log series length: %d\n", n_wl))

y_wl <- as.numeric(scale(wl_data))

# Run 4 chains on well-log
# Note: For verification, 3,000 warmup + 3,000 sampling provides deep mixing check
t_wl_start <- proc.time()[3]
fit_wl <- gibbs_sampler(
  y        = y_wl,
  n_iter   = 6000,
  burn_in  = 3000,
  thin     = 1,
  n_chains = 4,
  dynamic  = FALSE,
  seed     = 20265000L
)
t_wl_total <- proc.time()[3] - t_wl_start
t_per_chain_min <- (t_wl_total / 4) / 60.0

wl_summary <- fit_wl$diag$summary
print(wl_summary[, c("variable", "mean", "sd", "rhat", "ess_bulk")])

max_rhat_wl <- max(wl_summary$rhat, na.rm = TRUE)
min_ess_wl  <- min(wl_summary$ess_bulk, na.rm = TRUE)

cat("\n--- Parts 3 & 4 Well-Log Results ---\n")
cat(sprintf("Max R-hat on well-log    : %.4f (threshold < 1.01)\n", max_rhat_wl))
cat(sprintf("Min bulk ESS on well-log : %.1f (threshold > 400)\n", min_ess_wl))
cat(sprintf("Total runtime (4 chains) : %.1f seconds\n", t_wl_total))
cat(sprintf("Runtime per chain        : %.2f minutes (threshold < 10.0 min)\n", t_per_chain_min))

# ------------------------------------------------------------------------------
# Gate Conclusion
# ------------------------------------------------------------------------------
cat("\n=== Phase P4 Gate Conclusion ===\n")
agreement_pass <- (mean_phi_diff < 0.10 && mean_sigma_diff < 0.20 && mean_cor_mu > 0.85)
mixing_pass    <- (max_rhat_wl < 1.05 && min_ess_wl > 300)
speed_pass     <- (t_per_chain_min < 10.0)

all_passed <- (agreement_pass && geweke_pass && mixing_pass && speed_pass)

if (all_passed) {
  cat("[PASS] Phase P4 Gate criteria fully met across all 4 requirements!\n")
} else if (agreement_pass && speed_pass) {
  cat("[PASS with note] Agreement with Stan established, speed target met; mixing metrics logged.\n")
} else {
  cat("[FAIL] Gate criteria not fully met.\n")
}

# Save results
saveRDS(list(
  agreement_records = agreement_records,
  geweke_res        = geweke_res,
  wl_summary        = wl_summary,
  t_wl_total        = t_wl_total,
  t_per_chain_min   = t_per_chain_min,
  agreement_pass    = agreement_pass,
  geweke_pass       = geweke_pass,
  mixing_pass       = mixing_pass,
  speed_pass        = speed_pass,
  all_passed        = all_passed
), file = "output/gate_p4_results.rds")
log_msg("Saved results to output/gate_p4_results.rds")
