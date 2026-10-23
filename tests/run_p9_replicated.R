# ==============================================================================
# run_p9_replicated.R
# Phase P9 ablation (7 configurations) and sensitivity grid (21 configurations), each
# over 20 replicates of two settings (authors' decision, 21 October 2026):
#   base           : the P8 base cell (n = 1000, phi = 0.6, jump 2, isolated 2%, 2 collective)
#   ablation_n300  : the previous single ablation series' setting (n = 300, phi = 0.5,
#                    jump 2, mixed outliers 3%, 2 collective)
# Every configuration is run on the same simulated series (paired comparison). The
# proposed method runs at the spec length (4 chains, 5,000 + 5,000). One file per
# (setting, replicate, study, configuration) in output/p9_replicated/; reruns resume.
#
# Usage: Rscript tests/run_p9_replicated.R [n_workers]
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(parallel)
  library(future)
})
source("R/80_sim_design.R")
source("R/86_sensitivity_ablation.R")
plan(sequential)

args <- commandArgs(trailingOnly = TRUE)
N_WORKERS <- cpd_workers(if (length(args) >= 1) as.integer(args[1]) else NULL) # capped at CPD_MAX_CORES (24)
OUT_DIR <- "output/p9_replicated"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

R_REPS <- 20L
MCMC <- list(n_iter = 10000, burn_in = 5000, n_chains = 4)
base <- create_simulation_grid()[1, ]
SETTINGS <- list(
  base          = list(n = base$n, phi = base$phi, noise = base$noise, jump = base$jump, outlier_regime = base$outlier_regime,
                       p_outlier = base$p_outlier, collective = base$collective, seed0 = 950000L),
  ablation_n300 = list(n = 300L, phi = 0.5, noise = "gaussian", jump = 2.0, outlier_regime = "mixed", p_outlier = 0.03,
                       collective = 2L, seed0 = 960000L)
)
N_ABL <- 7L
N_SENS <- 21L

tasks <- rbind(
  expand.grid(setting = names(SETTINGS), rep = seq_len(R_REPS), study = "ablation", id = seq_len(N_ABL), stringsAsFactors = FALSE),
  expand.grid(setting = names(SETTINGS), rep = seq_len(R_REPS), study = "sensitivity", id = seq_len(N_SENS), stringsAsFactors = FALSE)
)
tasks <- tasks[order(tasks$setting != "base", tasks$rep), ] # n = 1000 first
task_file <- function(i) file.path(OUT_DIR, sprintf("%s__rep%02d__%s%02d.rds", tasks$setting[i], tasks$rep[i], tasks$study[i], tasks$id[i]))

run_task <- function(i) {
  f <- task_file(i)
  if (file.exists(f)) return("cached")
  if (cpd_should_stop(STOP_TIME)) return("stopped") # STOP file or CPD_STOP_AT reached: start nothing new
  st <- SETTINGS[[tasks$setting[i]]]
  rep <- tasks$rep[i]
  sim <- simulate_series(n = st$n, phi = st$phi, noise = st$noise, jump = st$jump, outlier_regime = st$outlier_regime,
                         p_outlier = st$p_outlier, collective = st$collective, seed = st$seed0 + rep)
  invisible(utils::capture.output(res <- if (tasks$study[i] == "ablation") {
    run_ablation_study(sim$y, truth = sim$truth, n_iter = MCMC$n_iter, burn_in = MCMC$burn_in, n_chains = MCMC$n_chains,
                       seed = 910L + 1000L * rep, ablation_ids = tasks$id[i])$results_df
  } else {
    run_sensitivity_grid(sim$y, truth = sim$truth, n_iter = MCMC$n_iter, burn_in = MCMC$burn_in, n_chains = MCMC$n_chains,
                         seed = 920L + 1000L * rep, sweep_ids = tasks$id[i])
  }))
  res$setting <- tasks$setting[i]
  res$replicate <- rep
  res$study <- tasks$study[i]
  tmp <- paste0(f, ".tmp")
  saveRDS(res, tmp)
  file.rename(tmp, f)
  "done"
}

STOP_TIME <- cpd_stop_time()
if (!is.null(STOP_TIME)) log_msg(sprintf("Will stop starting new tasks at %s", format(STOP_TIME)))
todo <- which(!file.exists(vapply(seq_len(nrow(tasks)), task_file, character(1))))
log_msg(sprintf("P9 replicated ablation/sensitivity: %d of %d tasks on %d workers", length(todo), nrow(tasks), N_WORKERS))
t0 <- proc.time()[3]
st <- unlist(mclapply(todo, function(i) tryCatch(run_task(i), error = function(e) paste("error:", conditionMessage(e))),
                      mc.cores = N_WORKERS, mc.preschedule = FALSE))
n_left <- sum(!file.exists(vapply(seq_len(nrow(tasks)), task_file, character(1))))
log_msg(sprintf("Ran %.1f h: %d done, %d errors; %d of %d tasks saved in total, %d left",
                (proc.time()[3] - t0) / 3600, sum(st == "done"), sum(grepl("^error", st)),
                nrow(tasks) - n_left, nrow(tasks), n_left))
if (any(grepl("^error", st))) print(table(st[grepl("^error", st)]))
if (any(st == "stopped")) {
  log_msg("Stopped on request; every finished task is saved. Rerun the same command to resume.")
  quit(status = 75)
}
