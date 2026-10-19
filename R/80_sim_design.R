# ==============================================================================
# 80_sim_design.R
# Simulation Study Design, Grid Generation, and Execution Engine.
# Implements Section 5 (Phase P8) of CPD2_Code_Workflow.md.
# ==============================================================================

if (file.exists("R/00_utils.R")) {
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
} else if (file.exists("../../R/00_utils.R")) {
  source("../../R/00_utils.R")
  source("../../R/10_dgp.R")
  source("../../R/20_frontend_nsp.R")
  source("../../R/21_frontend_capa.R")
  source("../../R/22_prior_map.R")
  source("../../R/32_gibbs_updates.R")
  source("../../R/40_scoring.R")
  source("../../R/41_bfdr.R")
  source("../../R/42_classify.R")
  source("../../R/31_gibbs_sampler.R")
  source("../../R/50_copula.R")
  source("../../R/60_baselines.R")
  source("../../R/70_metrics.R")
}

#' Create Simulation Study Scenario Grid
#'
#' Generates the 17 experimental cells defined in Section 5 (Phase P8):
#' - Base scenario
#' - One-factor-at-a-time sweeps over n, phi, noise, jump, outlier_regime, collective
#' - 2 stress cells: (phi = 0.9, clustered, t3) and (n = 500, jump = 1, isolated 5%)
#'
#' @return A data frame where each row is a simulation scenario
#' @export
create_simulation_grid <- function() {
  # Base scenario definition
  base_row <- list(
    cell_id        = "base",
    factor_name    = "base",
    n              = 1000L,
    phi            = 0.6,
    noise          = "gaussian",
    jump           = 2.0,
    outlier_regime = "isolated",
    p_outlier      = 0.02,
    collective     = 2L
  )

  grid_list <- list()
  grid_list[[1]] <- base_row

  # 1. Sweep n: 500, 2000
  grid_list[[2]] <- modifyList(base_row, list(cell_id = "n_500", factor_name = "n", n = 500L))
  grid_list[[3]] <- modifyList(base_row, list(cell_id = "n_2000", factor_name = "n", n = 2000L))

  # 2. Sweep phi: 0, 0.3, 0.9
  grid_list[[4]] <- modifyList(base_row, list(cell_id = "phi_0.0", factor_name = "phi", phi = 0.0))
  grid_list[[5]] <- modifyList(base_row, list(cell_id = "phi_0.3", factor_name = "phi", phi = 0.3))
  grid_list[[6]] <- modifyList(base_row, list(cell_id = "phi_0.9", factor_name = "phi", phi = 0.9))

  # 3. Sweep Innovations: t3
  grid_list[[7]] <- modifyList(base_row, list(cell_id = "noise_t3", factor_name = "noise", noise = "t3"))

  # 4. Sweep Jump size: 1, 3
  grid_list[[8]] <- modifyList(base_row, list(cell_id = "jump_1.0", factor_name = "jump", jump = 1.0))
  grid_list[[9]] <- modifyList(base_row, list(cell_id = "jump_3.0", factor_name = "jump", jump = 3.0))

  # 5. Sweep Outlier regime: none, isolated 5%, clustered (patch mean 3), mixed
  grid_list[[10]] <- modifyList(base_row, list(cell_id = "outlier_none", factor_name = "outlier_regime", outlier_regime = "none", p_outlier = 0.0))
  grid_list[[11]] <- modifyList(base_row, list(cell_id = "outlier_iso_5pct", factor_name = "outlier_regime", outlier_regime = "isolated", p_outlier = 0.05))
  grid_list[[12]] <- modifyList(base_row, list(cell_id = "outlier_clustered", factor_name = "outlier_regime", outlier_regime = "clustered", p_outlier = 0.02))
  grid_list[[13]] <- modifyList(base_row, list(cell_id = "outlier_mixed", factor_name = "outlier_regime", outlier_regime = "mixed", p_outlier = 0.04))

  # 6. Sweep Collective anomalies: 0, 4
  grid_list[[14]] <- modifyList(base_row, list(cell_id = "coll_0", factor_name = "collective", collective = 0L))
  grid_list[[15]] <- modifyList(base_row, list(cell_id = "coll_4", factor_name = "collective", collective = 4L))

  # 7. Stress cells:
  # Stress 1: (phi = 0.9, clustered outliers, t3)
  grid_list[[16]] <- modifyList(base_row, list(
    cell_id        = "stress_phi0.9_clust_t3",
    factor_name    = "stress",
    phi            = 0.9,
    outlier_regime = "clustered",
    noise          = "t3"
  ))

  # Stress 2: (n = 500, jump = 1, isolated 5%)
  grid_list[[17]] <- modifyList(base_row, list(
    cell_id        = "stress_n500_jump1_iso5",
    factor_name    = "stress",
    n              = 500L,
    jump           = 1.0,
    outlier_regime = "isolated",
    p_outlier      = 0.05
  ))

  grid_df <- do.call(rbind, lapply(grid_list, as.data.frame, stringsAsFactors = FALSE))
  rownames(grid_df) <- NULL
  grid_df
}

#' Run a Single Simulation Replicate for a Method
#'
#' @param sim Simulated series object returned by `simulate_series`
#' @param method Character identifier of the method
#' @param m_args List of optional method-specific arguments
#'
#' @return A list of evaluated metrics from `evaluate_all_metrics`
#' @export
run_method_on_replicate <- function(sim, method, m_args = list()) {
  y <- sim$y
  truth <- sim$truth
  n <- length(y)

  fit <- tryCatch({
    if (method == "proposed") {
      # Proposed method: prior map + production Gibbs sampler
      # 1. Front-ends
      intervals <- run_frontend_nsp(y, deg = 0, alpha = 0.05)
      capa_res  <- run_frontend_capa(y, prewhiten = TRUE)
      prior_maps <- build_prior_maps(y, intervals = intervals, capa_result = capa_res, kappa_in = 2.0, delta = 1.0)

      n_iter  <- if (!is.null(m_args$n_iter)) m_args$n_iter else 2500
      burn_in <- if (!is.null(m_args$burn_in)) m_args$burn_in else 1000

      gibbs_sampler(
        y               = y,
        kappa           = prior_maps$kappa,
        a_vec           = prior_maps$a_vec,
        n_iter          = n_iter,
        burn_in         = burn_in,
        thin            = 1,
        n_chains        = 2,
        dynamic         = TRUE,
        delta           = 1.0,
        capa_collective = capa_res$collective,
        alpha_fdr       = 0.05
      )
    } else if (method == "copula") {
      fit_copula_comparator(y, alpha = 0.05, B_boot = if (!is.null(m_args$B_boot)) m_args$B_boot else 50)
    } else {
      # Baseline dispatch
      do.call(fit_baseline, c(list(method = method, y = y), m_args))
    }
  }, error = function(e) {
    warning(sprintf("Error running method '%s': %s", method, e$message))
    NULL
  })

  if (is.null(fit)) {
    # Return NA row on failure
    return(list(
      method      = method,
      runtime_sec = NA_real_,
      cp          = list(precision = NA_real_, recall = NA_real_, f1 = NA_real_, hausdorff = NA_real_, covering = NA_real_),
      point       = list(precision = NA_real_, recall = NA_real_, f1 = NA_real_, realised_fdr = NA_real_),
      collective  = list(precision = NA_real_, recall = NA_real_, f1 = NA_real_),
      four_class  = list(macro_f1 = NA_real_)
    ))
  }

  evaluate_all_metrics(truth, fit)
}

#' Run Pilot Simulation Across Scenario Cells
#'
#' Runs 5 replicates per cell across selected methods to benchmark runtime
#' and verify pipeline execution integrity (Phase P8 Pilot).
#'
#' @param grid Optional scenario grid (default `create_simulation_grid()`)
#' @param n_reps Number of pilot replicates per cell (default 5)
#' @param methods Methods to run in pilot (default proposed, copula, decafs, capa, pelt, robseg)
#' @param seed_base Base random seed (default 8000)
#'
#' @return A tidy data frame with all replicate runs and metric summaries
#' @export
run_pilot_simulation <- function(
    grid = NULL,
    n_reps = 5,
    methods = c("proposed", "copula", "decafs", "pelt", "robseg", "capa"),
    seed_base = 8000
) {
  if (is.null(grid)) grid <- create_simulation_grid()

  n_cells <- nrow(grid)
  cat(sprintf("Initiating Pilot Simulation: %d cells x %d replicates x %d methods\n",
              n_cells, n_reps, length(methods)))

  results_rows <- list()
  row_idx <- 1L

  for (c_idx in seq_len(n_cells)) {
    cell <- grid[c_idx, ]
    cat(sprintf("\n--- Cell %2d/%d: %s (n = %d, phi = %.2f, %s, jump = %.1f) ---\n",
                c_idx, n_cells, cell$cell_id, cell$n, cell$phi, cell$outlier_regime, cell$jump))

    for (r in seq_len(n_reps)) {
      sim_seed <- seed_base + c_idx * 1000L + r

      sim <- simulate_series(
        n              = cell$n,
        phi            = cell$phi,
        noise          = cell$noise,
        jump           = cell$jump,
        outlier_regime = cell$outlier_regime,
        p_outlier      = cell$p_outlier,
        collective     = cell$collective,
        seed           = sim_seed
      )

      for (m in methods) {
        t0 <- proc.time()[3]
        eval_res <- run_method_on_replicate(sim, m)
        t_elapsed <- as.numeric(proc.time()[3] - t0)

        results_rows[[row_idx]] <- data.frame(
          cell_id       = cell$cell_id,
          factor_name   = cell$factor_name,
          n             = cell$n,
          phi           = cell$phi,
          noise         = cell$noise,
          jump          = cell$jump,
          outlier_reg   = cell$outlier_regime,
          replicate     = r,
          method        = m,
          runtime_sec   = eval_res$runtime_sec,
          cp_f1         = eval_res$cp$f1,
          cp_prec       = eval_res$cp$precision,
          cp_rec        = eval_res$cp$recall,
          cp_covering   = eval_res$cp$covering,
          pt_f1         = eval_res$point$f1,
          pt_fdr        = eval_res$point$realised_fdr,
          coll_f1       = eval_res$collective$f1,
          macro_f1      = eval_res$four_class$macro_f1,
          stringsAsFactors = FALSE
        )
        row_idx <- row_idx + 1L
      }
    }
  }

  pilot_df <- do.call(rbind, results_rows)
  rownames(pilot_df) <- NULL

  # Compute Budget Summary
  timing_summary <- aggregate(runtime_sec ~ method, data = pilot_df, FUN = function(x) {
    c(mean = mean(x, na.rm = TRUE), median = median(x, na.rm = TRUE), total = sum(x, na.rm = TRUE))
  })

  cat("\n====================================================================\n")
  cat("PILOT EXECUTION COMPLETE\n")
  cat("====================================================================\n")
  print(timing_summary)

  list(
    data           = pilot_df,
    timing_summary = timing_summary,
    grid           = grid
  )
}
