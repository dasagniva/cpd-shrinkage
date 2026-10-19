# ==============================================================================
# gate_p8_verification.R
# Acceptance Gate for Phase P8: Simulation Study Pilot & Budget Verification
#
# Tests:
# 1. Full grid generation: 17 cells (Base + sweeps + stress cells).
# 2. Pilot execution: 5 replicates per cell across benchmark methods.
# 3. Budget check: Total compute estimation, verification of 0 NA fits (no silent failures),
#    and Monte Carlo standard error reporting for every metric.
# ==============================================================================

cat("\n====================================================================\n")
cat("      PHASE P8 ACCEPTANCE GATE VERIFICATION\n")
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

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Verify Grid Generation
# ------------------------------------------------------------------------------
cat("1. Generating simulation grid...\n")
grid <- create_simulation_grid()
cat(sprintf("  Generated %d cells across factors: %s\n",
            nrow(grid), paste(unique(grid$factor_name), collapse = ", ")))
stopifnot(nrow(grid) == 17)

# ------------------------------------------------------------------------------
# 2. Run Pilot Simulation (5 Replicates per Cell)
# ------------------------------------------------------------------------------
cat("\n2. Running Pilot Simulation (5 replicates per cell)...\n")

# Representative methods across families:
# Proposed, Copula comparator, DeCAFS, PELT, Robust FPOP, CAPA, BCP, FPNN
pilot_methods <- c("proposed", "copula", "decafs", "pelt", "robseg", "capa", "fpnn", "ch4")

t_start_pilot <- proc.time()[3]
pilot_output <- run_pilot_simulation(
  grid        = grid,
  n_reps      = 5,
  methods     = pilot_methods,
  seed_base   = 8000,
  m_args_list = list(
    proposed = list(n_iter = 1000, burn_in = 400, n_chains = 1),
    copula   = list(B_boot = 30)
  )
)
pilot_duration <- as.numeric(proc.time()[3] - t_start_pilot)

pilot_df <- pilot_output$data
cat(sprintf("\nPilot simulation finished in %.2f seconds (%d total fits)\n",
            pilot_duration, nrow(pilot_df)))

# ------------------------------------------------------------------------------
# 3. Check for Silent Failures (NA Counts)
# ------------------------------------------------------------------------------
cat("\n3. Checking for silent failures (NA fits)...\n")
na_counts <- tapply(pilot_df$macro_f1, pilot_df$method, function(x) sum(is.na(x)))
cat("NA counts per method:\n")
print(na_counts)

total_na <- sum(na_counts)
pass_no_silent_failures <- (total_na == 0)
cat(sprintf("No silent failures: %s (total NA fits = %d)\n", pass_no_silent_failures, total_na))

# ------------------------------------------------------------------------------
# 4. Monte Carlo SE Reporting & Performance Summary
# ------------------------------------------------------------------------------
cat("\n4. Performance Summary with Monte Carlo Standard Errors:\n")

summary_by_method <- aggregate(
  cbind(cp_f1, macro_f1, runtime_sec) ~ method,
  data = pilot_df,
  FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = sd(x, na.rm = TRUE) / sqrt(length(x)))
)

cat(sprintf("%-12s | %-16s | %-16s | %-16s\n",
            "Method", "CP F1 (SE)", "Macro F1 (SE)", "Runtime (sec)"))
cat("----------------------------------------------------------------------------\n")

for (i in seq_len(nrow(summary_by_method))) {
  m <- summary_by_method$method[i]
  cp_m <- summary_by_method$cp_f1[i, "mean"]
  cp_se <- summary_by_method$cp_f1[i, "se"]
  mf_m <- summary_by_method$macro_f1[i, "mean"]
  mf_se <- summary_by_method$macro_f1[i, "se"]
  rt_m <- summary_by_method$runtime_sec[i, "mean"]

  cat(sprintf("%-12s | %6.3f (%.4f) | %6.3f (%.4f) | %8.3f s\n",
              m, cp_m, cp_se, mf_m, mf_se, rt_m))
}
cat("----------------------------------------------------------------------------\n")

# ------------------------------------------------------------------------------
# 5. Compute Budget Check
# ------------------------------------------------------------------------------
cat("\n5. Compute Budget Check:\n")
total_runtime_pilot <- sum(pilot_df$runtime_sec, na.rm = TRUE)
avg_time_per_sim_rep <- mean(pilot_df$runtime_sec, na.rm = TRUE)

# Estimated cost for 100 replicates of all 17 cells:
total_full_cost_sec <- total_runtime_pilot * (100.0 / 5.0)
cat(sprintf("  Pilot total runtime : %.2f seconds\n", total_runtime_pilot))
cat(sprintf("  Estimated cost for 100 reps across 17 cells: %.1f minutes\n",
            total_full_cost_sec / 60.0))

pass_budget <- (total_full_cost_sec / 3600.0) < 10.0 # Under 10 hours
cat(sprintf("  Budget check passed: %s (estimated %.2f hours < 10 hrs)\n\n",
            pass_budget, total_full_cost_sec / 3600.0))

# ------------------------------------------------------------------------------
# Gate Status and Persistence
# ------------------------------------------------------------------------------
gate_p8_passed <- pass_no_silent_failures && pass_budget

results <- list(
  gate_passed          = gate_p8_passed,
  summary_by_method    = summary_by_method,
  timing_summary       = pilot_output$timing_summary,
  pilot_data           = pilot_df,
  total_runtime_pilot  = total_runtime_pilot,
  estimated_cost_hours = total_full_cost_sec / 3600.0,
  timestamp            = Sys.time()
)

saveRDS(results, "output/gate_p8_results.rds")

cat("====================================================================\n")
if (gate_p8_passed) {
  cat(">>> PHASE P8 ACCEPTANCE GATE: PASSED <<<\n")
} else {
  cat(">>> PHASE P8 ACCEPTANCE GATE: FAILED <<<\n")
}
cat("====================================================================\n\n")

if (!gate_p8_passed) {
  stop("Phase P8 acceptance gate failed!")
}
