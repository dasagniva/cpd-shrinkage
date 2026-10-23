# ==============================================================================
# run_real_data_btc_extend.R
# BTC did not meet the P4 mixing criterion at 4 x 50,000 iterations (sigma and
# n_outliers bulk ESS ~190, R-hat 1.02). Refit the proposed method on BTC at
# 4 x 100,000 iterations (50,000 burn-in) and replace it in output/real_data_long.rds.
# ==============================================================================
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
suppressPackageStartupMessages(library(future))
source("R/80_sim_design.R")
source("R/85_real_data.R")
plan(multicore, workers = 4)
GATED <- c("phi", "sigma", "n_shift_windows", "n_outliers")
args_btc <- list(n_iter = 100000, burn_in = 50000, n_chains = 4, alpha_fdr = 0.05, seed = 909)
rd <- readRDS("output/real_data_long.rds")
fit <- run_real_data_method(load_real_series("btc", raw_dir = "data/raw")$y, method = "proposed", m_args = args_btc)
rd$fits$btc$proposed <- fit
rd$summary_table[rd$summary_table$series == "btc" & rd$summary_table$method == "proposed", ] <- summarize_real_data_fit(fit, "btc")
sm <- fit$diag$summary
g <- sm[sm$variable %in% GATED, ]
rd$mixing[rd$mixing$series == "btc", c("max_rhat_gated", "min_ess_gated", "max_rhat_all", "min_ess_all", "mixing_ok")] <-
  list(max(g$rhat), min(g$ess_bulk), max(sm$rhat, na.rm = TRUE), min(sm$ess_bulk, na.rm = TRUE), max(g$rhat) < 1.01 && min(g$ess_bulk) > 400)
rd$mixing$n_iter <- ifelse(rd$mixing$series == "btc", args_btc$n_iter, rd$proposed_args$n_iter)
rd$mixing$burn_in <- ifelse(rd$mixing$series == "btc", args_btc$burn_in, rd$proposed_args$burn_in)
print(as.data.frame(sm)[, c("variable", "mean", "rhat", "ess_bulk")], digits = 3)
print(rd$mixing, digits = 4)
saveRDS(rd, "output/real_data_long.rds")
