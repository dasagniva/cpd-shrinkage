# ==============================================================================
# run_real_data_long.R
# Phase P9 real-data benchmark at the MCMC length that converged on well-log
# (4 chains x 50,000 iterations, 25,000 burn-in; tests/run_welllog_long.R). The spec
# length (5,000 + 5,000) does not mix on the n = 4050 series.
#
# The proposed method runs on the four series in parallel (4 series x 4 chains); all
# baselines and the copula comparator run on every series. Mixing is checked on the
# quantities gated in P4 (phi, sigma, n_shift_windows, n_outliers; change C8).
# Output: output/real_data_long.rds, read by tests/gate_p9_verification.R.
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(future)
  library(future.apply)
  library(parallel)
})
source("R/80_sim_design.R")
source("R/85_real_data.R")

SERIES <- c("indpro", "welllog", "oilwell", "btc")
REAL_PROPOSED_ARGS <- list(n_iter = 50000, burn_in = 25000, n_chains = 4, alpha_fdr = 0.05, seed = 909)
BASELINES <- c("copula", "decafs", "abco", "nsp", "capa", "pelt", "robseg", "hsmuce", "bcp", "envcpt", "ch4", "fpnn")
GATED <- c("phi", "sigma", "n_shift_windows", "n_outliers")

data_list <- lapply(setNames(SERIES, SERIES), load_real_series, raw_dir = "data/raw")

t0 <- proc.time()[3]
plan(list(tweak(multicore, workers = length(SERIES)), tweak(multicore, workers = 4)))
log_msg(sprintf("Proposed method on %d series: %d iterations (%d burn-in) x %d chains",
                length(SERIES), REAL_PROPOSED_ARGS$n_iter, REAL_PROPOSED_ARGS$burn_in, REAL_PROPOSED_ARGS$n_chains))
proposed_fits <- future_lapply(SERIES, function(s) {
  run_real_data_method(data_list[[s]]$y, method = "proposed", m_args = REAL_PROPOSED_ARGS)
}, future.seed = TRUE)
names(proposed_fits) <- SERIES
t_proposed <- as.numeric(proc.time()[3] - t0)
plan(sequential)

log_msg("Baselines and copula comparator on every series")
jobs <- expand.grid(method = BASELINES, series = SERIES, stringsAsFactors = FALSE)
baseline_fits <- mclapply(seq_len(nrow(jobs)), function(j) {
  m_args <- if (jobs$method[j] == "copula") list(B_boot = 50) else list()
  run_real_data_method(data_list[[jobs$series[j]]]$y, method = jobs$method[j], m_args = m_args)
}, mc.cores = cpd_workers(nrow(jobs)), mc.preschedule = FALSE)

fits <- lapply(setNames(SERIES, SERIES), function(s) {
  c(list(proposed = proposed_fits[[s]]),
    setNames(baseline_fits[jobs$series == s], jobs$method[jobs$series == s]))
})
summary_table <- do.call(rbind, lapply(SERIES, function(s) {
  do.call(rbind, lapply(names(fits[[s]]), function(m) summarize_real_data_fit(fits[[s]][[m]], series_name = s)))
}))

mixing <- do.call(rbind, lapply(SERIES, function(s) {
  sm <- fits[[s]]$proposed$diag$summary
  g <- sm[sm$variable %in% GATED, ]
  data.frame(series = s, max_rhat_gated = max(g$rhat), min_ess_gated = min(g$ess_bulk),
             max_rhat_all = max(sm$rhat, na.rm = TRUE), min_ess_all = min(sm$ess_bulk, na.rm = TRUE),
             mixing_ok = max(g$rhat) < 1.01 && min(g$ess_bulk) > 400)
}))
print(summary_table, digits = 3)
print(mixing, digits = 4)

saveRDS(list(
  fits            = fits,
  summary_table   = summary_table,
  mixing          = mixing,
  proposed_args   = REAL_PROPOSED_ARGS,
  t_proposed_sec  = t_proposed,
  t_total_sec     = as.numeric(proc.time()[3] - t0),
  timestamp       = Sys.time()
), "output/real_data_long.rds")
log_msg(sprintf("Saved output/real_data_long.rds (%.1f min)", (proc.time()[3] - t0) / 60))
