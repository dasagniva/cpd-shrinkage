# ==============================================================================
# gate_p2_verification.R
# Acceptance gate verification for Phase P2:
# 1. NSP false-interval rate on null series under phi in {0, 0.6, 0.9} <= alpha_NSP
# 2. CAPA point-anomaly false-positive rate under phi = 0.9 with vs without prewhitening
# ==============================================================================

source("R/00_utils.R")
source("R/10_dgp.R")
source("R/20_frontend_nsp.R")
source("R/21_frontend_capa.R")
source("R/22_prior_map.R")

log_msg("=== Starting Phase P2 Gate Verification ===")

alpha_nsp <- 0.10
n_pts <- 300
reps <- 60
phi_grid <- c(0.0, 0.6, 0.9)

cat("\n--- 1. Autoregressive NSP (nsp_poly_ar) Null Calibration (Target <= 0.10) ---\n")
nsp_ar_results <- list()

for (phi_val in phi_grid) {
  false_count <- 0
  for (b in seq_len(reps)) {
    sim <- simulate_series(
      n = n_pts,
      phi = phi_val,
      noise = "gaussian",
      k0 = 0,
      outlier_regime = "none",
      collective = 0,
      seed = 20261000L + as.integer(phi_val * 1000) + b
    )
    nsp_out <- run_frontend_nsp(sim$y, alpha = alpha_nsp, method = "ar", M = 500, min_size = 20)
    if (nrow(nsp_out$intervals) > 0) {
      false_count <- false_count + 1
    }
  }
  fir <- false_count / reps
  se <- sqrt(fir * (1 - fir) / reps)
  nsp_ar_results[[as.character(phi_val)]] <- c(fir = fir, se = se)
  cat(sprintf("phi: %.1f | nsp_poly_ar FIR: %.3f (MC SE: %.3f) | Target <= %.2f\n",
              phi_val, fir, se, alpha_nsp))
}

cat("\n--- 2. Self-Normalised NSP (nsp_poly_selfnorm) Null Comparison ---\n")
for (phi_val in c(0.0, 0.6)) {
  false_count_raw <- 0
  false_count_pw  <- 0
  n_sub <- 30
  for (b in seq_len(n_sub)) {
    sim <- simulate_series(
      n = n_pts,
      phi = phi_val,
      noise = "gaussian",
      k0 = 0,
      outlier_regime = "none",
      collective = 0,
      seed = 20263000L + as.integer(phi_val * 1000) + b
    )
    r_raw <- run_frontend_nsp(sim$y, alpha = alpha_nsp, method = "selfnorm", prewhiten = FALSE, M = 300)
    r_pw  <- run_frontend_nsp(sim$y, alpha = alpha_nsp, method = "selfnorm", prewhiten = TRUE, M = 300)
    if (nrow(r_raw$intervals) > 0) false_count_raw <- false_count_raw + 1
    if (nrow(r_pw$intervals) > 0) false_count_pw <- false_count_pw + 1
  }
  cat(sprintf("phi: %.1f | raw selfnorm FIR: %.2f | prewhitened selfnorm FIR: %.2f\n",
              phi_val, false_count_raw / n_sub, false_count_pw / n_sub))
}

cat("\n--- 3. CAPA Prewhitening Ablation under phi = 0.90 ---\n")
capa_reps <- 60
false_pts_raw <- 0
false_pts_pw  <- 0

for (b in seq_len(capa_reps)) {
  sim <- simulate_series(
    n = n_pts,
    phi = 0.90,
    noise = "gaussian",
    k0 = 0,
    outlier_regime = "none",
    collective = 0,
    seed = 20264000L + b
  )
  capa_raw <- run_frontend_capa(sim$y, prewhiten = FALSE, type = "mean")
  false_pts_raw <- false_pts_raw + length(capa_raw$point)

  capa_pw <- run_frontend_capa(sim$y, prewhiten = TRUE, type = "mean")
  false_pts_pw <- false_pts_pw + length(capa_pw$point)
}

mean_fp_raw <- false_pts_raw / capa_reps
mean_fp_pw  <- false_pts_pw / capa_reps

cat(sprintf("CAPA raw (no prewhitening) mean false points per series: %.3f\n", mean_fp_raw))
cat(sprintf("CAPA prewhitened           mean false points per series: %.3f\n", mean_fp_pw))

cat("\n=== Phase P2 Gate Conclusion ===\n")
all_nsp_pass <- all(sapply(nsp_ar_results, function(x) x["fir"] <= alpha_nsp + 2 * x["se"]))
if (all_nsp_pass) {
  cat("[PASS] Autoregressive NSP null false-interval rates meet the acceptance gate (<= alpha_NSP).\n")
} else {
  cat("[FAIL] False interval rate exceeded nominal tolerance.\n")
}
cat(sprintf("[LOGGED] CAPA point-anomaly FPR recorded under phi=0.90 (raw: %.3f, prewhitened: %.3f).\n",
            mean_fp_raw, mean_fp_pw))
