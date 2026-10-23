# ==============================================================================
# run_welllog_long.R
# Long Well-Log MCMC Run (50,000 iterations, 25,000 burn-in, 4 chains)
# Purpose: Diagnose whether convergence failure on well-log is slow mixing vs
#          multimodality (chains stuck on different feature explanations).
# ==============================================================================

source("R/00_utils.R")
source("R/10_dgp.R")
source("R/20_frontend_nsp.R")
source("R/21_frontend_capa.R")
source("R/22_prior_map.R")
source("R/31_gibbs_sampler.R")
source("R/32_gibbs_updates.R")
source("R/40_scoring.R")
source("R/41_bfdr.R")
source("R/42_classify.R")
source("R/85_real_data.R")

suppressPackageStartupMessages({
  library(posterior)
  library(coda)
  library(future)
})

workers_avail <- min(4, future::availableCores())
plan(multicore, workers = workers_avail)
log_msg(sprintf("=== Starting Long Well-Log Run: 4 chains x 50,000 iter (25,000 burn-in) using %d workers ===", workers_avail))

wl <- load_real_series("welllog")
cat(sprintf("Loaded well-log series: n = %d, range = [%.1f, %.1f]\n", wl$n, min(wl$y), max(wl$y)))

t0 <- proc.time()[3]
N_ITER <- 50000
BURN_IN <- 25000
fit_wl <- fit_proposed(
  wl$y,
  n_iter   = N_ITER,
  burn_in  = BURN_IN,
  n_chains = 4,
  seed     = 20265000L
)
t_wl_total <- proc.time()[3] - t0
t_per_chain_min <- (t_wl_total / 60) / ceiling(4 / workers_avail)

cat(sprintf("\n=== Execution Complete in %.2f minutes (%.2f min / chain batch) ===\n\n",
            t_wl_total / 60, t_per_chain_min))

# ------------------------------------------------------------------------------
# Diagnostics & Convergence Summary
# ------------------------------------------------------------------------------
wl_summary <- fit_wl$diag$summary
cat("--- Overall Posterior Summary (Draws Array) ---\n")
print(as.data.frame(wl_summary[, c("variable", "mean", "sd", "rhat", "ess_bulk", "ess_tail")]))

# Per-chain parameter breakdown to diagnose multimodality
draws_arr <- as_draws_array(fit_wl$diag$draws)
chain_stats <- list()
for (ch in 1:4) {
  ch_draws <- subset_draws(draws_arr, chain = ch)
  ch_sum <- summarize_draws(ch_draws, "mean", "sd", "median")
  ch_sum$chain <- ch
  chain_stats[[ch]] <- ch_sum
}
all_chain_stats <- do.call(rbind, chain_stats)
cat("\n--- Per-Chain Means and Medians ---\n")
for (v in unique(all_chain_stats$variable)) {
  v_sub <- all_chain_stats[all_chain_stats$variable == v, ]
  cat(sprintf("Variable %-18s: Ch1: %.2f (%.2f) | Ch2: %.2f (%.2f) | Ch3: %.2f (%.2f) | Ch4: %.2f (%.2f)\n",
              v,
              v_sub$mean[1], v_sub$median[1],
              v_sub$mean[2], v_sub$median[2],
              v_sub$mean[3], v_sub$median[3],
              v_sub$mean[4], v_sub$median[4]))
}

# Gated variables
gated_vars <- c("phi", "sigma", "n_shift_windows", "n_outliers")
gated <- wl_summary[wl_summary$variable %in% gated_vars, ]
max_rhat_gated <- max(gated$rhat, na.rm = TRUE)
min_ess_gated <- min(gated$ess_bulk, na.rm = TRUE)

cat(sprintf("\nGated Criteria (phi, sigma, n_shift_windows, n_outliers):\n"))
cat(sprintf("  Max R-hat: %.4f (threshold < 1.01)\n", max_rhat_gated))
cat(sprintf("  Min Bulk ESS: %.1f (threshold > 400)\n", min_ess_gated))

# Decision outputs
cat(sprintf("\n--- Detected Features ---\n"))
cat(sprintf("Number of Level Shifts: %d\n", length(fit_wl$cp)))
if (length(fit_wl$cp) > 0) cat("Shift locations:", paste(fit_wl$cp, collapse = ", "), "\n")
cat(sprintf("Number of Point Outliers: %d\n", length(fit_wl$point)))
if (length(fit_wl$point) > 0) cat("Point locations:", paste(fit_wl$point, collapse = ", "), "\n")
cat(sprintf("Number of Collective Segments: %d\n", nrow(fit_wl$collective)))
if (nrow(fit_wl$collective) > 0) {
  print(fit_wl$collective)
}

out_res <- list(
  fit_wl          = fit_wl,
  wl_summary      = wl_summary,
  chain_stats     = all_chain_stats,
  t_wl_total      = t_wl_total,
  t_per_chain_min = t_per_chain_min,
  n_iter          = N_ITER,
  burn_in         = BURN_IN,
  max_rhat_gated  = max_rhat_gated,
  min_ess_gated   = min_ess_gated,
  timestamp       = Sys.time()
)

saveRDS(out_res, file = "output/welllog_long_run.rds")
cat("\nSaved results to output/welllog_long_run.rds\n")
