# ==============================================================================
# run_p8_full.R
# Phase P8 full simulation study: 17 cells, all methods.
#
# Budget rule of the spec (Section 5, P8): 200 replicates for non-MCMC methods and 100 for
# MCMC methods (proposed, ABCO, BCP, Chapter 4). The proposed method runs at the spec
# length (4 chains, 5,000 burn-in + 5,000 retained).
#
# `targets` is not installed in the project library, so resumability comes from one RDS
# file per (cell, replicate) in output/p8_full/: a rerun skips finished tasks. Each task
# simulates one series and runs every method scheduled for that replicate on it, so the
# methods are compared on the same data.
#
# Usage: Rscript tests/run_p8_full.R [n_workers]
# ==============================================================================

Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(parallel)
  library(future)
})
source("R/80_sim_design.R")
plan(sequential) # chains run sequentially inside each worker; tasks run in parallel

args <- commandArgs(trailingOnly = TRUE)
N_WORKERS <- cpd_workers(if (length(args) >= 1) as.integer(args[1]) else NULL) # capped at CPD_MAX_CORES (24)
OUT_DIR <- "output/p8_full"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

METHODS_MCMC <- c("proposed", "abco", "bcp", "ch4")
METHODS_FAST <- c("decafs", "nsp", "capa", "pelt", "robseg", "hsmuce", "envcpt", "fpnn", "copula")
R_MCMC <- 100L
R_FAST <- 200L
SEED_BASE <- 900000L # distinct from the pilot (8000 + 1000 cell + rep)
PROPOSED_ARGS <- list(n_iter = 10000, burn_in = 5000, n_chains = 4, alpha_fdr = 0.05)

grid <- create_simulation_grid()
tasks <- expand.grid(rep = seq_len(R_FAST), cell = seq_len(nrow(grid)))
# Heaviest tasks first (MCMC replicates, then larger n) for load balancing
tasks <- tasks[order(tasks$rep > R_MCMC, -grid$n[tasks$cell], tasks$rep), ]

task_file <- function(cell, rep) file.path(OUT_DIR, sprintf("%s__rep%03d.rds", grid$cell_id[cell], rep))

# rbind data frames whose metric columns differ by method (missing columns become NA)
bind_fill <- function(lst) {
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) {
    d[setdiff(cols, names(d))] <- NA
    d[cols]
  }))
}

flatten_metrics <- function(ev) {
  v <- unlist(ev[c("cp", "point", "collective", "four_class")])
  v <- v[vapply(v, function(x) is.numeric(x) && length(x) == 1, logical(1))]
  c(v, unlist(ev$fit_info))
}

run_task <- function(i) {
  cell <- tasks$cell[i]
  rep <- tasks$rep[i]
  f <- task_file(cell, rep)
  if (file.exists(f)) return("cached")
  if (cpd_should_stop(STOP_TIME)) return("stopped") # STOP file or CPD_STOP_AT reached: start nothing new
  cl <- grid[cell, ]
  sim_seed <- SEED_BASE + cell * 1000L + rep
  sim <- simulate_series(n = cl$n, phi = cl$phi, noise = cl$noise, jump = cl$jump,
                         outlier_regime = cl$outlier_regime, p_outlier = cl$p_outlier,
                         collective = cl$collective, seed = sim_seed)
  methods <- c(if (rep <= R_MCMC) METHODS_MCMC, METHODS_FAST)
  rows <- lapply(methods, function(m) {
    m_args <- switch(m, proposed = c(PROPOSED_ARGS, list(seed = sim_seed + 1L)), copula = list(B_boot = 50), list())
    t0 <- proc.time()[3]
    ev <- tryCatch(run_method_on_replicate(sim, m, m_args = m_args), error = function(e) NULL)
    wall <- as.numeric(proc.time()[3] - t0)
    met <- if (!is.null(ev)) flatten_metrics(ev) else c()
    row <- data.frame(cell_id = cl$cell_id, factor_name = cl$factor_name, n = cl$n, phi = cl$phi, noise = cl$noise,
                      jump = cl$jump, outlier_reg = cl$outlier_regime, replicate = rep, sim_seed = sim_seed,
                      method = m, wall_sec = wall, failed = is.null(ev) || is.na(met["cp.f1"]),
                      stringsAsFactors = FALSE)
    if (length(met) > 0) row <- cbind(row, as.data.frame(as.list(met), check.names = FALSE))
    row
  })
  out <- bind_fill(rows)
  tmp <- paste0(f, ".tmp")
  saveRDS(out, tmp)
  file.rename(tmp, f)
  "done"
}

STOP_TIME <- cpd_stop_time()
if (!is.null(STOP_TIME)) log_msg(sprintf("Will stop starting new tasks at %s", format(STOP_TIME)))
todo <- which(!file.exists(mapply(task_file, tasks$cell, tasks$rep)))
log_msg(sprintf("P8 full study: %d of %d tasks to run on %d workers", length(todo), nrow(tasks), N_WORKERS))
t0 <- proc.time()[3]
status <- mclapply(todo, function(i) tryCatch(run_task(i), error = function(e) paste("error:", conditionMessage(e))),
                   mc.cores = N_WORKERS, mc.preschedule = FALSE)
st <- unlist(status)
n_left <- sum(!file.exists(mapply(task_file, tasks$cell, tasks$rep)))
log_msg(sprintf("Ran %.1f h: %d done, %d errors; %d of %d tasks saved in total, %d left",
                (proc.time()[3] - t0) / 3600, sum(st == "done"), sum(grepl("^error", st)),
                nrow(tasks) - n_left, nrow(tasks), n_left))
if (any(grepl("^error", st))) print(table(st[grepl("^error", st)]))
if (any(st == "stopped")) {
  log_msg("Stopped on request; every finished task is saved. Rerun the same command to resume.")
  quit(status = 75)
}
