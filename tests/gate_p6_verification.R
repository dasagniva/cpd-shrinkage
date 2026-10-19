# ==============================================================================
# gate_p6_verification.R
# Acceptance Gate for Phase P6: Copula Markov Chain Comparator
#
# Tests:
# 1. Family and parameter recovery on simulated copula Markov chains:
#    - Gaussian (rho = 0.65)
#    - Student-t (rho = 0.60, df = 4)
#    - Gumbel (theta = 2.2, exhibiting upper tail dependence)
# 2. Null false-flag rate is approximately nominal alpha (with Monte Carlo SE).
# ==============================================================================

cat("\n====================================================================\n")
cat("      PHASE P6 ACCEPTANCE GATE VERIFICATION\n")
cat("====================================================================\n\n")

source("R/00_utils.R")
source("R/50_copula.R")

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Part 1: Family and Parameter Recovery
# ==============================================================================
cat("--- Part 1: Family and Parameter Recovery on Copula Markov Chains ---\n")

B_rec <- 15
n_rec <- 600

# 1. Gaussian Recovery
cat(sprintf("Testing Gaussian Markov recovery (rho = 0.65, %d replicates, n = %d)...\n", B_rec, n_rec))
cop_gauss <- rvinecopulib::bicop_dist("gaussian", 0, 0.65)
gauss_families <- character(B_rec)
gauss_rhos     <- numeric(B_rec)

set.seed(601)
for (b in seq_len(B_rec)) {
  u <- simulate_copula_markov(cop_gauss, n_rec)
  fit <- rvinecopulib::bicop(cbind(u[-n_rec], u[-1]), family_set = c("gaussian", "t", "clayton", "gumbel", "frank"))
  gauss_families[b] <- fit$family
  gauss_rhos[b]     <- if (fit$family %in% c("gaussian", "t")) fit$parameters[1] else sin((pi / 2.0) * rvinecopulib::par_to_ktau(fit))
}

prop_gauss_match <- mean(gauss_families %in% c("gaussian", "t"))
mean_gauss_rho   <- mean(gauss_rhos)
se_gauss_rho     <- sd(gauss_rhos) / sqrt(B_rec)
cat(sprintf("  Gaussian / t selection rate: %.1f%%\n", prop_gauss_match * 100))
cat(sprintf("  Mean estimated rho: %.4f (true = 0.6500, MC SE = %.4f)\n", mean_gauss_rho, se_gauss_rho))

# 2. Student-t Recovery
cat(sprintf("\nTesting Student-t Markov recovery (rho = 0.60, df = 4, %d replicates, n = %d)...\n", B_rec, n_rec))
cop_t <- rvinecopulib::bicop_dist("t", 0, c(0.60, 4.0))
t_families <- character(B_rec)
t_rhos     <- numeric(B_rec)
t_dfs      <- numeric(B_rec)

set.seed(602)
for (b in seq_len(B_rec)) {
  u <- simulate_copula_markov(cop_t, n_rec)
  fit <- rvinecopulib::bicop(cbind(u[-n_rec], u[-1]), family_set = c("gaussian", "t", "clayton", "gumbel", "frank"))
  t_families[b] <- fit$family
  t_rhos[b]     <- fit$parameters[1]
  t_dfs[b]      <- if (fit$family == "t") fit$parameters[2] else NA_real_
}

prop_t_match <- mean(t_families == "t")
mean_t_rho   <- mean(t_rhos)
mean_t_df    <- mean(t_dfs, na.rm = TRUE)
cat(sprintf("  t-family selection rate: %.1f%%\n", prop_t_match * 100))
cat(sprintf("  Mean estimated rho: %.4f (true = 0.6000)\n", mean_t_rho))
cat(sprintf("  Mean estimated df : %.4f (true = 4.0000)\n", mean_t_df))

# 3. Gumbel Recovery
cat(sprintf("\nTesting Gumbel Markov recovery (theta = 2.2, %d replicates, n = %d)...\n", B_rec, n_rec))
cop_gumb <- rvinecopulib::bicop_dist("gumbel", 0, 2.2)
gumb_families <- character(B_rec)
gumb_thetas   <- numeric(B_rec)
gumb_lambda_u <- numeric(B_rec)

set.seed(603)
for (b in seq_len(B_rec)) {
  u <- simulate_copula_markov(cop_gumb, n_rec)
  fit <- rvinecopulib::bicop(cbind(u[-n_rec], u[-1]), family_set = c("gaussian", "t", "clayton", "gumbel", "frank"))
  gumb_families[b] <- fit$family
  gumb_thetas[b]   <- fit$parameters[1]
  td <- rvinecopulib::tail_dep(fit)
  gumb_lambda_u[b] <- td["upper", "upper"]
}

prop_gumb_match <- mean(gumb_families == "gumbel")
mean_gumb_theta <- mean(gumb_thetas)
mean_gumb_lu    <- mean(gumb_lambda_u)
cat(sprintf("  Gumbel selection rate: %.1f%%\n", prop_gumb_match * 100))
cat(sprintf("  Mean estimated theta   : %.4f (true = 2.2000)\n", mean_gumb_theta))
cat(sprintf("  Mean upper tail dep lu : %.4f (true = %.4f)\n", mean_gumb_lu, 2 - 2^(1/2.2)))

pass_part1 <- (prop_gauss_match >= 0.80) &&
              (abs(mean_gauss_rho - 0.65) <= 0.08) &&
              (prop_t_match >= 0.70) &&
              (prop_gumb_match >= 0.80) &&
              (abs(mean_gumb_theta - 2.2) <= 0.25)

cat(sprintf("\nPart 1 Recovery Pass: %s\n\n", pass_part1))

# ==============================================================================
# Part 2: Null False-Flag Rate Calibration
# ==============================================================================
cat("--- Part 2: Null False-Flag Rate Calibration ---\n")

B_null <- 25
n_null <- 150
alpha_target <- 0.05
null_rates <- numeric(B_null)

cat(sprintf("Evaluating null false-flag rate on %d series (n = %d, target alpha = %.2f)...\n",
            B_null, n_null, alpha_target))

set.seed(20261020)
pb <- txtProgressBar(min = 0, max = B_null, style = 3)

for (b in seq_len(B_null)) {
  phi_val <- sample(c(0.2, 0.5, 0.7), 1)
  sigma_val <- 1.0

  eps <- numeric(n_null)
  eps[1] <- rnorm(1, 0, sigma_val / sqrt(1 - phi_val^2))
  for (t in 2:n_null) {
    eps[t] <- phi_val * eps[t - 1] + rnorm(1, 0, sigma_val)
  }

  fit_cop <- fit_copula_comparator(
    y       = eps,
    alpha   = alpha_target,
    B_boot  = 100,
    seed    = 7000 + b
  )

  null_rates[b] <- length(fit_cop$point) / n_null
  setTxtProgressBar(pb, b)
}
close(pb)

mean_null_rate <- mean(null_rates)
se_null_rate   <- sd(null_rates) / sqrt(B_null)

cat(sprintf("\nResults across %d null series:\n", B_null))
cat(sprintf("  Target Alpha            : %.3f\n", alpha_target))
cat(sprintf("  Mean Null False-Flag Rate: %.4f (MC SE = %.4f)\n", mean_null_rate, se_null_rate))

# Check: mean null false-flag rate within alpha +- 2 * SE
pass_part2 <- abs(mean_null_rate - alpha_target) <= max(0.02, 2.5 * se_null_rate)
cat(sprintf("  Null False-Flag Rate Pass: %s (|%.4f - %.3f| <= %.4f)\n\n",
            pass_part2, mean_null_rate, alpha_target, max(0.02, 2.5 * se_null_rate)))

# ==============================================================================
# Summary and Persistence
# ==============================================================================
gate_p6_passed <- pass_part1 && pass_part2

results <- list(
  gate_passed          = gate_p6_passed,
  gauss_rho_mean       = mean_gauss_rho,
  gauss_rho_se         = se_gauss_rho,
  gauss_match_rate     = prop_gauss_match,
  t_match_rate         = prop_t_match,
  t_rho_mean           = mean_t_rho,
  t_df_mean            = mean_t_df,
  gumb_match_rate      = prop_gumb_match,
  gumb_theta_mean      = mean_gumb_theta,
  gumb_lambda_u_mean   = mean_gumb_lu,
  null_rate_mean       = mean_null_rate,
  null_rate_se         = se_null_rate,
  alpha_target         = alpha_target,
  timestamp            = Sys.time()
)

saveRDS(results, "output/gate_p6_results.rds")

cat("====================================================================\n")
if (gate_p6_passed) {
  cat(">>> PHASE P6 ACCEPTANCE GATE: PASSED <<<\n")
} else {
  cat(">>> PHASE P6 ACCEPTANCE GATE: FAILED <<<\n")
}
cat("====================================================================\n\n")

if (!gate_p6_passed) {
  stop("Phase P6 acceptance gate failed!")
}
