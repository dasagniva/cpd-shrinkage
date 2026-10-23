# ==============================================================================
# gate_p9_verification.R
# Acceptance Gate for Phase P9: Real Data, Ablation and Sensitivity.
#
# Inputs (run first):
#   tests/run_real_data_long.R  -> output/real_data_long.rds (all methods on the four
#                                  series; proposed method at 4 x 50,000 iterations)
#   tests/run_p9_replicated.R   -> output/p9_replicated/ (ablation and sensitivity grid,
#                                  20 replicates each of the base cell and the n = 300
#                                  setting)
#
# Gate:
#   1. every real-data fit exists, passes validate_cpd_fit() and did not fail;
#   2. the proposed fit tracks the data (median |y - mu_hat| < 3 robust noise scales);
#   3. the proposed fit mixes on every series (gated quantities of P4: R-hat < 1.01,
#      bulk ESS > 400);
#   4. the ablation and sensitivity studies are complete (all replicates).
# ==============================================================================

cat("\n====================================================================\n")
cat("      PHASE P9 ACCEPTANCE GATE VERIFICATION\n")
cat("====================================================================\n\n")

source("R/80_sim_design.R")
source("R/85_real_data.R")
source("R/90_plots_tables.R")

# ------------------------------------------------------------------------------
# 1-3. Real data
# ------------------------------------------------------------------------------
rd <- readRDS("output/real_data_long.rds")
fits <- rd$fits
all_fits_valid <- all(unlist(lapply(names(fits), function(s) lapply(fits[[s]], function(f) {
  isTRUE(tryCatch({ validate_cpd_fit(f); TRUE }, error = function(e) FALSE))
}))))
failed_fits <- do.call(rbind, lapply(names(fits), function(s) {
  data.frame(series = s, method = names(fits[[s]]), failed = sapply(fits[[s]], function(f) isTRUE(f$diag$failed)))
}))
cat(sprintf("All real-data fits valid: %s; failed fits: %d\n", all_fits_valid, sum(failed_fits$failed)))
if (any(failed_fits$failed)) print(failed_fits[failed_fits$failed, ])

level_ok <- sapply(names(fits), function(s) {
  f <- fits[[s]]$proposed
  y_s <- load_real_series(s)$y
  !isTRUE(f$diag$failed) && median(abs(y_s - f$diag$mu_mean)) / robust_standardize(y_s)$scale < 3
})
print(level_ok)
print(rd$mixing, digits = 4)
mixing_ok <- all(rd$mixing$mixing_ok)
cat(sprintf("Proposed method: %d iterations (%d burn-in) x %d chains\n",
            rd$proposed_args$n_iter, rd$proposed_args$burn_in, rd$proposed_args$n_chains))

# ------------------------------------------------------------------------------
# 4. Replicated ablation and sensitivity
# ------------------------------------------------------------------------------
files <- list.files("output/p9_replicated", pattern = "\\.rds$", full.names = TRUE)
parts <- lapply(files, readRDS)
abl <- do.call(rbind, lapply(parts[vapply(parts, function(d) d$study[1] == "ablation", logical(1))], function(d) {
  d[, c("setting", "replicate", "ablation", "description", "n_cp", "n_point", "n_coll", "cp_prec", "cp_rec", "cp_f1",
        "pt_f1", "pt_fdr", "pt_prauc", "coll_f1", "macro_f1", "phi_hat", "max_rhat", "runtime_sec")]
}))
sens <- do.call(rbind, lapply(parts[vapply(parts, function(d) d$study[1] == "sensitivity", logical(1))], function(d) {
  d[, c("setting", "replicate", "sweep_id", "param", "val", "n_cp", "n_point", "n_coll", "cp_prec", "cp_rec", "cp_f1",
        "pt_f1", "pt_fdr", "pt_prauc", "coll_f1", "macro_f1", "max_rhat", "runtime_sec")]
}))
n_settings <- 2L
abl_complete <- nrow(abl) == n_settings * 20L * 7L
sens_complete <- nrow(sens) == n_settings * 20L * 21L
cat(sprintf("Ablation rows %d (expected %d) | sensitivity rows %d (expected %d)\n",
            nrow(abl), n_settings * 140L, nrow(sens), n_settings * 420L))
abl_order <- c("full", "no_nsp", "no_capa", "indep_outliers", "no_ar1", "standard_hs", "marginal_var")
abl <- abl[order(abl$setting != "base", match(abl$ablation, abl_order), abl$replicate), ]
sens <- sens[order(sens$setting != "base", sens$sweep_id, sens$replicate), ]

gate_p9_passed <- all_fits_valid && !any(failed_fits$failed) && all(level_ok) && mixing_ok && abl_complete && sens_complete

results_p9 <- list(
  gate_passed         = gate_p9_passed,
  real_data_summary   = rd$summary_table,
  level_ok            = level_ok,
  mixing              = rd$mixing,
  failed_fits         = failed_fits,
  real_data_fits      = fits,
  proposed_args       = rd$proposed_args,
  ablation_results    = abl,
  sensitivity_results = sens,
  timestamp           = Sys.time()
)
saveRDS(results_p9, "output/gate_p9_results.rds")

cat("\n====================================================================\n")
cat(sprintf(">>> PHASE P9 ACCEPTANCE GATE: %s <<<\n", if (gate_p9_passed) "PASSED" else "FAILED"))
cat("====================================================================\n\n")
if (!gate_p9_passed) stop("Phase P9 acceptance gate failed!")
