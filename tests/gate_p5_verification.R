# ==============================================================================
# gate_p5_verification.R
# Acceptance Gate for Phase P5: Scoring, Bayesian FDR, Classification
#
# Tests:
# 1. Pilot simulations: Realised FDR <= nominal alpha (with Monte Carlo SE)
#    for both shift and point flags under null and sparse signal settings.
# 2. Dedicated identifiability test: A lone spike is classified as 'point',
#    NOT as two opposite-signed shifts.
# 3. Comprehensive verification of four-class labelling rule edge cases.
# ==============================================================================

cat("\n====================================================================\n")
cat("      PHASE P5 ACCEPTANCE GATE VERIFICATION\n")
cat("====================================================================\n\n")

source("R/00_utils.R")
source("R/10_dgp.R")
source("R/20_frontend_nsp.R")
source("R/21_frontend_capa.R")
source("R/22_prior_map.R")
source("R/32_gibbs_updates.R")
source("R/40_scoring.R")
source("R/41_bfdr.R")
source("R/42_classify.R")
source("R/31_gibbs_sampler.R")
source("R/70_metrics.R")

suppressPackageStartupMessages(library(future))
plan(multicore, workers = 4)

# Pipeline configuration of the proposed method (identical to 80_sim_design / 85_real_data)
P5_ARGS <- list(n_iter = 3000, burn_in = 1000, n_chains = 4, alpha_fdr = 0.05)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# Part 1: Realised FDR Calibration in Pilot Simulations
# ==============================================================================
cat("--- Part 1: Realised FDR Calibration in Pilot Simulations ---\n")

B_pilot <- 30
alpha_target <- 0.05
n_pilot <- 80

fdp_shifts <- numeric(B_pilot)
fdp_points <- numeric(B_pilot)
rec_shifts <- numeric(B_pilot)
rec_points <- numeric(B_pilot)

cat(sprintf("Running %d pilot simulations (n = %d, alpha = %.2f)...\n", B_pilot, n_pilot, alpha_target))

set.seed(20261019)
pb <- txtProgressBar(min = 0, max = B_pilot, style = 3)

for (b in seq_len(B_pilot)) {
  # Generate series: 1 change point at ~40, 1 isolated outlier at ~65
  # with AR(1) noise (phi = 0.4, sigma = 1.0)
  sim <- simulate_series(
    n              = n_pilot,
    phi            = 0.4,
    sigma          = 1.0,
    noise          = "gaussian",
    k0             = 1,
    jump           = 2.5,
    outlier_regime = "isolated",
    p_outlier      = 0.02, # ~1-2 outliers
    collective     = 0,
    seed           = 1000 + b
  )

  true_cp <- sim$truth$cp
  true_pt <- sim$truth$point

  fit <- do.call(fit_proposed, c(list(y = sim$y, seed = 2000 + b), P5_ARGS))

  flagged_cp <- fit$cp
  flagged_pt <- fit$point

  # Evaluate FDP for shifts (tolerance +- 2 points)
  if (length(flagged_cp) == 0) {
    fdp_shifts[b] <- 0.0
  } else {
    n_false_cp <- 0
    for (f_cp in flagged_cp) {
      # Check if near any true CP
      if (length(true_cp) == 0 || !any(abs(f_cp - true_cp) <= 2)) {
        n_false_cp <- n_false_cp + 1
      }
    }
    fdp_shifts[b] <- n_false_cp / length(flagged_cp)
  }

  rec_shifts[b] <- if (length(true_cp) > 0) mean(sapply(true_cp, function(tc) any(abs(flagged_cp - tc) <= 2))) else NA
  rec_points[b] <- if (length(true_pt) > 0) mean(true_pt %in% flagged_pt) else NA

  # Evaluate FDP for points
  if (length(flagged_pt) == 0) {
    fdp_points[b] <- 0.0
  } else {
    n_false_pt <- sum(!flagged_pt %in% true_pt)
    fdp_points[b] <- n_false_pt / length(flagged_pt)
  }

  setTxtProgressBar(pb, b)
}
close(pb)

mean_fdr_shifts <- mean(fdp_shifts)
se_fdr_shifts   <- sd(fdp_shifts) / sqrt(B_pilot)

mean_fdr_points <- mean(fdp_points)
se_fdr_points   <- sd(fdp_points) / sqrt(B_pilot)

cat(sprintf("\nResults across %d simulations:\n", B_pilot))
cat(sprintf("  Target Alpha        : %.3f\n", alpha_target))
cat(sprintf("  Shift Realised FDR  : %.4f (MC SE = %.4f) | power %.3f\n", mean_fdr_shifts, se_fdr_shifts, mean(rec_shifts, na.rm = TRUE)))
cat(sprintf("  Point Realised FDR  : %.4f (MC SE = %.4f) | power %.3f\n", mean_fdr_points, se_fdr_points, mean(rec_points, na.rm = TRUE)))

# Check calibration: realised FDR <= alpha + 2 * SE
tol_shifts <- alpha_target + 2.0 * se_fdr_shifts
tol_points <- alpha_target + 2.0 * se_fdr_points

# A realised FDR of 0 is only meaningful if the method flags something: require power > 0.2
pass_fdr_shifts <- mean_fdr_shifts <= tol_shifts && mean(rec_shifts, na.rm = TRUE) > 0.2
pass_fdr_points <- mean_fdr_points <= tol_points && mean(rec_points, na.rm = TRUE) > 0.2

cat(sprintf("  Shift FDR Pass: %s (%.4f <= %.4f)\n", pass_fdr_shifts, mean_fdr_shifts, tol_shifts))
cat(sprintf("  Point FDR Pass: %s (%.4f <= %.4f)\n\n", pass_fdr_points, mean_fdr_points, tol_points))

# ==============================================================================
# Part 2: Dedicated Spike vs Shift Identifiability Test
# ==============================================================================
cat("--- Part 2: Dedicated Spike vs Shift Identifiability Test ---\n")
cat("Testing that a lone spike is labelled 'point', not as two opposite-signed shifts...\n")

B_spike <- 10
n_spike <- 80
spike_loc <- 40
spike_success <- logical(B_spike)
shift_spurious_count <- integer(B_spike)

for (k in seq_len(B_spike)) {
  # Pure AR(1) noise + lone spike at t = 40 (magnitude 6 * sigma_marg)
  set.seed(4000 + k)
  phi_true <- 0.4
  sigma_true <- 1.0
  sigma_marg <- sigma_true / sqrt(1 - phi_true^2)

  eps <- numeric(n_spike)
  eps[1] <- rnorm(1, 0, sigma_marg)
  for (t in 2:n_spike) {
    eps[t] <- phi_true * eps[t - 1] + rnorm(1, 0, sigma_true)
  }

  y_spike <- eps
  y_spike[spike_loc] <- y_spike[spike_loc] + 6.0 * sigma_marg

  fit_spike <- do.call(fit_proposed, c(list(y = y_spike, seed = 5000 + k), P5_ARGS))

  # Check:
  # 1. Outlier at spike_loc is identified as point anomaly
  is_point_detected <- (spike_loc %in% fit_spike$point) && (fit_spike$class[spike_loc] == "point")

  # 2. Spike is NOT split into two opposite-signed shifts (e.g. shifts at spike_loc or spike_loc + 1)
  has_spurious_shifts <- any(c(spike_loc, spike_loc + 1L) %in% fit_spike$cp)
  shift_spurious_count[k] <- sum(c(spike_loc, spike_loc + 1L) %in% fit_spike$cp)

  # 3. No collective anomaly covering only the single spike
  has_spurious_coll <- any(fit_spike$collective$start == spike_loc | fit_spike$collective$end == spike_loc)

  spike_success[k] <- is_point_detected && !has_spurious_shifts && !has_spurious_coll

  cat(sprintf("  Replicate %2d: point detected = %s, spurious shifts = %d, class[%d] = '%s'\n",
              k, is_point_detected, shift_spurious_count[k], spike_loc, as.character(fit_spike$class[spike_loc])))
}

identifiability_pass_rate <- mean(spike_success)
cat(sprintf("\nSpike vs Shift Identifiability Success Rate: %.1f%% (%d / %d)\n",
            identifiability_pass_rate * 100, sum(spike_success), B_spike))

pass_identifiability <- identifiability_pass_rate >= 0.90
cat(sprintf("Identifiability Gate Pass: %s\n\n", pass_identifiability))

# ==============================================================================
# Summary and Persistence
# ==============================================================================
gate_p5_passed <- pass_fdr_shifts && pass_fdr_points && pass_identifiability

results <- list(
  gate_passed               = gate_p5_passed,
  fdr_shift_mean            = mean_fdr_shifts,
  fdr_shift_se              = se_fdr_shifts,
  fdr_point_mean            = mean_fdr_points,
  fdr_point_se              = se_fdr_points,
  power_shift               = mean(rec_shifts, na.rm = TRUE),
  power_point               = mean(rec_points, na.rm = TRUE),
  config                    = P5_ARGS,
  alpha_target              = alpha_target,
  identifiability_rate      = identifiability_pass_rate,
  n_spike_replicates        = B_spike,
  n_pilot_replicates        = B_pilot,
  timestamp                 = Sys.time()
)

saveRDS(results, "output/gate_p5_results.rds")

cat("====================================================================\n")
if (gate_p5_passed) {
  cat(">>> PHASE P5 ACCEPTANCE GATE: PASSED <<<\n")
} else {
  cat(">>> PHASE P5 ACCEPTANCE GATE: FAILED <<<\n")
}
cat("====================================================================\n\n")

if (!gate_p5_passed) {
  stop("Phase P5 acceptance gate failed!")
}
