# ==============================================================================
# gate_p8_verification.R
# Acceptance Gate for Phase P8: Full Simulation Study
#
# The study itself is run by tests/run_p8_full.R (resumable, one file per cell x
# replicate in output/p8_full/). The 5-replicate pilot that set the budget is in git
# history (commit 2130b73); its results are not used for tables.
#
# Gate (spec Section 5, P8):
#   1. every cell is complete for every method (planned replicates: 100 for MCMC methods,
#      200 for the others);
#   2. no silent failures: failed fits are counted per method and reported;
#   3. the Monte Carlo SE is reported next to every table entry.
# MCMC convergence of the proposed method is reported per cell (share of fits with gated
# R-hat < 1.01), not gated.
# ==============================================================================

cat("\n====================================================================\n")
cat("      PHASE P8 ACCEPTANCE GATE VERIFICATION (full study)\n")
cat("====================================================================\n\n")

source("R/80_sim_design.R")
source("R/90_plots_tables.R")

full <- collect_p8_full("output/p8_full")
df <- full$data
comp <- full$completeness

# 1. Completeness
incomplete <- comp[!comp$complete, ]
cat(sprintf("Cells x methods complete: %d of %d\n", sum(comp$complete), nrow(comp)))
if (nrow(incomplete) > 0) print(incomplete)
pass_complete <- nrow(incomplete) == 0

# 2. Failures, counted per method (a failed fit has NA metrics and failed = TRUE)
fail_tab <- aggregate(cbind(fits = 1, failed = as.numeric(failed)) ~ method, df, sum)
fail_tab$na_cp_f1 <- tapply(df$cp_f1, df$method, function(x) sum(is.na(x)))[fail_tab$method]
print(fail_tab)
silent <- sum(is.na(df$cp_f1) & !df$failed)
cat(sprintf("Silent failures (NA metric without a failure flag): %d\n", silent))
pass_no_silent <- silent == 0

# 3. Summaries with Monte Carlo SE
metrics <- c("cp_f1", "cp_prec", "cp_rec", "cp_covering", "cp_hausdorff", "pt_f1", "pt_fdr", "pt_prauc",
             "coll_f1", "macro_f1", "runtime_sec")
summ <- do.call(rbind, lapply(intersect(metrics, names(df)), function(m) summarise_p8_metric(df, m)))
pass_se <- all(is.finite(summ$se[summ$n > 1]))

conv <- df[df$method == "proposed", ]
conv_tab <- aggregate(cbind(rhat_ok = as.numeric(max_rhat_gated < 1.01), max_rhat_gated) ~ cell_id, conv, mean)
names(conv_tab) <- c("cell_id", "share_rhat_below_1.01", "mean_max_rhat_gated")
cat("\nProposed method: convergence of the gated quantities by cell\n")
print(conv_tab, digits = 3)

for (m in c("cp_f1", "pt_f1", "macro_f1")) {
  table_sim_metric_latex(
    summ[summ$metric == m, ],
    caption = sprintf("Simulation study: %s by cell and method, mean (Monte Carlo SE); 100 replicates for MCMC methods, 200 otherwise",
                      c(cp_f1 = "change-point F1", pt_f1 = "point-anomaly F1", macro_f1 = "four-class macro-F1")[[m]]),
    output_file = sprintf("output/tab1_sim_%s.tex", m)
  )
}
utils::write.csv(summ, "output/p8_full_summary.csv", row.names = FALSE)

gate_p8_passed <- pass_complete && pass_no_silent && pass_se
saveRDS(list(
  gate_passed  = gate_p8_passed,
  completeness = comp,
  failures     = fail_tab,
  summary      = summ,
  convergence  = conv_tab,
  full_data    = df,
  timestamp    = Sys.time()
), "output/gate_p8_results.rds")

cat("====================================================================\n")
cat(sprintf(">>> PHASE P8 ACCEPTANCE GATE: %s <<<\n", if (gate_p8_passed) "PASSED" else "FAILED"))
cat("====================================================================\n\n")
if (!gate_p8_passed) stop("Phase P8 acceptance gate failed!")
