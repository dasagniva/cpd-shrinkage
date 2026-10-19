# ==============================================================================
# gate_p7_verification.R
# Acceptance Gate for Phase P7: Baseline Wrappers and Common Output Schema
#
# Gate: Every wrapper passes validate_cpd_fit() on one simulated series
# and on well-log (n = 4050).
# ==============================================================================

cat("\n====================================================================\n")
cat("      PHASE P7 ACCEPTANCE GATE VERIFICATION\n")
cat("====================================================================\n\n")

source("R/00_utils.R")
source("R/10_dgp.R")
source("R/21_frontend_capa.R")
source("R/32_gibbs_updates.R")
source("R/40_scoring.R")
source("R/41_bfdr.R")
source("R/42_classify.R")
source("R/50_copula.R")
source("R/60_baselines.R")

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Prepare Datasets
# ------------------------------------------------------------------------------
cat("1. Preparing datasets...\n")

# A. Simulated series (n = 200)
set.seed(7001)
sim_data <- simulate_series(
  n              = 200,
  phi            = 0.4,
  sigma          = 1.0,
  noise          = "gaussian",
  k0             = 3,
  jump           = 2.5,
  outlier_regime = "isolated",
  p_outlier      = 0.02,
  collective     = 1,
  seed           = 7001
)
y_sim <- sim_data$y
n_sim <- length(y_sim)
cat(sprintf("  Simulated series ready: n = %d, true CPs = %d, true points = %d\n",
            n_sim, length(sim_data$truth$cp), length(sim_data$truth$point)))

# B. Well-log series (n = 4050)
well_df <- read.csv("data/raw/welllog.csv", header = TRUE)
y_well <- as.numeric(well_df[[1]])
# Remove any NAs if present
y_well <- y_well[!is.na(y_well)]
n_well <- length(y_well)
cat(sprintf("  Well-log series loaded: n = %d (min = %.1f, max = %.1f)\n\n",
            n_well, min(y_well), max(y_well)))

# ------------------------------------------------------------------------------
# 2. Baseline Wrappers to Test
# ------------------------------------------------------------------------------
methods <- c(
  "decafs",
  "pelt",
  "robseg",
  "nsp",
  "capa",
  "hsmuce",
  "bcp",
  "envcpt",
  "ch4",
  "fpnn",
  "copula",
  "abco"
)

sim_results  <- list()
well_results <- list()
all_passed   <- TRUE

cat("2. Running validation across all baseline wrappers...\n")
cat(sprintf("%-12s | %-12s | %-10s | %-12s | %-10s\n",
            "Method", "Sim (sec)", "Sim Valid", "Well (sec)", "Well Valid"))
cat("---------------------------------------------------------------------\n")

for (m in methods) {
  # --- Test on Simulated Series ---
  t0 <- proc.time()[3]
  fit_sim <- tryCatch({
    if (m == "abco") {
      fit_baseline(m, y_sim, nsave = 100, nburn = 100)
    } else if (m == "bcp") {
      fit_baseline(m, y_sim, burnin = 50, mcmc = 100)
    } else if (m == "ch4") {
      fit_baseline(m, y_sim, n_iter = 300, burn_in = 100)
    } else if (m == "copula") {
      fit_baseline(m, y_sim, B_boot = 50)
    } else {
      fit_baseline(m, y_sim)
    }
  }, error = function(e) {
    cat(sprintf("ERROR in %s on sim: %s\n", m, e$message))
    NULL
  })
  t_sim <- as.numeric(proc.time()[3] - t0)

  valid_sim <- FALSE
  if (!is.null(fit_sim)) {
    val_res <- tryCatch({
      validate_cpd_fit(fit_sim, n = n_sim)
      TRUE
    }, error = function(e) FALSE)
    valid_sim <- val_res
  }

  # --- Test on Well-Log Series ---
  t0 <- proc.time()[3]
  fit_well <- tryCatch({
    if (m == "abco") {
      # Scale welllog for ABCO stability
      y_w_scaled <- (y_well - median(y_well)) / mad(y_well)
      fit_baseline(m, y_w_scaled, nsave = 100, nburn = 100)
    } else if (m == "bcp") {
      fit_baseline(m, y_well, burnin = 50, mcmc = 100)
    } else if (m == "ch4") {
      fit_baseline(m, y_well, n_iter = 300, burn_in = 100)
    } else if (m == "copula") {
      fit_baseline(m, y_well, B_boot = 50)
    } else {
      fit_baseline(m, y_well)
    }
  }, error = function(e) {
    cat(sprintf("ERROR in %s on well-log: %s\n", m, e$message))
    NULL
  })
  t_well <- as.numeric(proc.time()[3] - t0)

  valid_well <- FALSE
  if (!is.null(fit_well)) {
    val_res <- tryCatch({
      validate_cpd_fit(fit_well, n = n_well)
      TRUE
    }, error = function(e) FALSE)
    valid_well <- val_res
  }

  cat(sprintf("%-12s | %9.3f s | %-10s | %9.3f s | %-10s\n",
              m, t_sim, ifelse(valid_sim, "PASS", "FAIL"),
              t_well, ifelse(valid_well, "PASS", "FAIL")))

  sim_results[[m]]  <- list(valid = valid_sim, time = t_sim, fit = fit_sim)
  well_results[[m]] <- list(valid = valid_well, time = t_well, fit = fit_well)

  if (!valid_sim || !valid_well) all_passed <- FALSE
}

cat("---------------------------------------------------------------------\n\n")

# Save Gate Results
gate_p7_results <- list(
  gate_passed  = all_passed,
  methods      = methods,
  sim_results  = sim_results,
  well_results = well_results,
  timestamp    = Sys.time()
)

saveRDS(gate_p7_results, "output/gate_p7_results.rds")

cat("====================================================================\n")
if (all_passed) {
  cat(">>> PHASE P7 ACCEPTANCE GATE: PASSED <<<\n")
} else {
  cat(">>> PHASE P7 ACCEPTANCE GATE: FAILED <<<\n")
}
cat("====================================================================\n\n")

if (!all_passed) {
  stop("Phase P7 acceptance gate failed!")
}
