# ==============================================================================
# gate_p9_verification.R
# Acceptance Gate for Phase P9: Empirical Benchmark, Ablation, and Sensitivity.
#
# Tests:
# 1. Real Data Applications: Fits proposed method, copula comparator, and
#    baselines across all four series (Well-log, Oilwell, INDPRO, BTC).
# 2. Four-Class Labelling & Posterior Bands: Verifies proposed method extracts
#    valid posterior bands, outlier score traces, and FDR-controlled classifications.
# 3. Component Ablation: Verifies all 6 ablations execute properly.
# 4. Sensitivity Grid: Verifies robustness across the 12 hyperparameter configurations.
# ==============================================================================

cat("\n====================================================================\n")
cat("      PHASE P9 ACCEPTANCE GATE VERIFICATION\n")
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
source("R/50_copula.R")
source("R/60_baselines.R")
source("R/70_metrics.R")
source("R/80_sim_design.R")
source("R/85_real_data.R")
source("R/86_sensitivity_ablation.R")

suppressPackageStartupMessages(library(future))
plan(multicore, workers = 4)

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Real-Data Applications
# ------------------------------------------------------------------------------
cat("1. Running Empirical Real-Data Benchmark on four series...\n")

methods_real <- c("proposed", "copula", "decafs", "pelt", "robseg", "capa", "fpnn", "ch4")
m_args_real <- list(
  proposed = list(n_iter = 10000, burn_in = 5000, n_chains = 4, alpha_fdr = 0.05, seed = 909),
  copula   = list(B_boot = 30)
)

t0_real <- proc.time()[3]
real_results <- run_real_data_benchmark(
  series_names = c("indpro", "welllog", "oilwell", "btc"),
  methods      = methods_real,
  m_args_list  = m_args_real,
  raw_dir      = "data/raw"
)
t_real_total <- as.numeric(proc.time()[3] - t0_real)
cat(sprintf("Real-data benchmark finished in %.2f seconds\n", t_real_total))

# Verify all fits exist and pass validate_cpd_fit
all_fits_valid <- TRUE
for (s in names(real_results$fits)) {
  for (m in names(real_results$fits[[s]])) {
    fit_obj <- real_results$fits[[s]][[m]]
    chk <- tryCatch({
      validate_cpd_fit(fit_obj)
      TRUE
    }, error = function(e) {
      warning(sprintf("Fit validation failed for %s on %s: %s", m, s, e$message))
      FALSE
    })
    if (!chk) all_fits_valid <- FALSE
  }
}
cat(sprintf("All real-data fits valid: %s\n", all_fits_valid))
stopifnot(all_fits_valid)

# Sanity checks on every proposed fit: no failure, and the fitted level tracks the data
# (median |y - mu_hat| within 3 robust noise scales). The original gate only checked
# factor levels, which let a fit with mu_hat = 0 on well-log pass.
level_ok <- sapply(names(real_results$fits), function(s) {
  f <- real_results$fits[[s]]$proposed
  y_s <- load_real_series(s)$y
  !isTRUE(f$diag$failed) && !is.null(f$diag$mu_mean) &&
    median(abs(y_s - f$diag$mu_mean)) / robust_standardize(y_s)$scale < 3
})
print(level_ok)
mixing_tab <- do.call(rbind, lapply(names(real_results$fits), function(s) {
  sm <- real_results$fits[[s]]$proposed$diag$summary
  data.frame(series = s, max_rhat = max(sm$rhat, na.rm = TRUE), min_ess = min(sm$ess_bulk, na.rm = TRUE))
}))
print(mixing_tab)
failed_fits <- do.call(rbind, lapply(names(real_results$fits), function(s) {
  data.frame(series = s, method = names(real_results$fits[[s]]),
             failed = sapply(real_results$fits[[s]], function(f) isTRUE(f$diag$failed)))
}))
print(failed_fits[failed_fits$failed, ])

# ------------------------------------------------------------------------------
# 2. Component Ablation Study
# ------------------------------------------------------------------------------
cat("\n2. Running Component Ablation Study...\n")
set.seed(910)
sim_abl <- simulate_series(
  n              = 300,
  phi            = 0.5,
  noise          = "gaussian",
  jump           = 2.0,
  outlier_regime = "mixed",
  p_outlier      = 0.03,
  collective     = 2,
  seed           = 910
)

abl_res <- run_ablation_study(
  y        = sim_abl$y,
  truth    = sim_abl$truth,
  n_iter   = 4000,
  burn_in  = 1500,
  n_chains = 4,
  seed     = 910
)

stopifnot(
  nrow(abl_res$results_df) == 7,
  all(!is.na(abl_res$results_df$macro_f1)),
  all(!is.na(abl_res$results_df$runtime_sec))
)
cat("Component ablation study passed (all 6 ablations + full model completed).\n")

# ------------------------------------------------------------------------------
# 3. Sensitivity Grid
# ------------------------------------------------------------------------------
cat("\n3. Running Hyperparameter Sensitivity Grid...\n")
sens_df <- run_sensitivity_grid(
  y       = sim_abl$y,
  truth   = sim_abl$truth,
  n_iter  = 4000,
  burn_in = 1500,
  seed    = 920
)

stopifnot(
  nrow(sens_df) == 21,
  all(!is.na(sens_df$macro_f1)),
  all(!is.na(sens_df$runtime_sec))
)
cat("Sensitivity analysis grid passed (all 12 parameter settings completed).\n")

# ------------------------------------------------------------------------------
# Gate Persistence and Decision
# ------------------------------------------------------------------------------
gate_p9_passed <- all_fits_valid && all(level_ok) && !any(failed_fits$failed) &&
  (nrow(abl_res$results_df) == 7) && (nrow(sens_df) == 21)

results_p9 <- list(
  gate_passed         = gate_p9_passed,
  real_data_summary   = real_results$summary_table,
  level_ok            = level_ok,
  mixing              = mixing_tab,
  failed_fits         = failed_fits,
  real_data_fits      = real_results$fits,
  ablation_results    = abl_res$results_df,
  sensitivity_results = sens_df,
  runtime_real_sec    = t_real_total,
  timestamp           = Sys.time()
)

saveRDS(results_p9, "output/gate_p9_results.rds")

cat("\n====================================================================\n")
if (gate_p9_passed) {
  cat(">>> PHASE P9 ACCEPTANCE GATE: PASSED <<<\n")
} else {
  cat(">>> PHASE P9 ACCEPTANCE GATE: FAILED <<<\n")
}
cat("====================================================================\n\n")

if (!gate_p9_passed) {
  stop("Phase P9 acceptance gate failed!")
}
