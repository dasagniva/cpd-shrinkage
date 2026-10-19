# ==============================================================================
# 86_sensitivity_ablation.R
# Sensitivity Analysis and Ablation Study (Phase P9).
# Evaluates sensitivity across prior hyperparameters and performs component
# ablations (NSP, CAPA, dynamic outlier coupling, AR(1), HS+, scoring variance).
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
  source("../../R/70_metrics.R")
}

#' Run Comprehensive Component Ablation Study
#'
#' Evaluates the six core ablations defined in Section 5 (Phase P9):
#' 1. full: Full proposed model
#' 2. no_nsp: No significance front-end (kappa_in = 1.0)
#' 3. no_capa: No anomaly front-end (delta = 0, no CAPA collective overlap)
#' 4. indep_outliers: Independent outlier scales (rho = 0, dynamic = FALSE)
#' 5. no_ar1: No AR(1) autoregressive noise (fix phi = 0)
#' 6. standard_hs: Standard Horseshoe instead of Horseshoe+ (eta = 1)
#' 7. marginal_var: Marginal variance instead of innovation variance in o_t (Eq. 8)
#'
#' @param y Numeric vector of observations (simulated or real series)
#' @param truth Optional ground truth list from `simulate_series`
#' @param n_iter Gibbs iterations
#' @param burn_in Gibbs burn-in
#' @param n_chains Number of chains
#' @param seed Master random seed
#'
#' @return A list containing `results_df` and list of `fits`
#' @export
run_ablation_study <- function(
    y,
    truth = NULL,
    n_iter = 2000,
    burn_in = 800,
    n_chains = 2,
    seed = 42
) {
  n <- length(y)
  stopifnot(n >= 20, is.numeric(y))

  ablation_defs <- list(
    list(name = "full", desc = "Full Proposed Model"),
    list(name = "no_nsp", desc = "Ablation 1: No NSP (kappa_in = 1.0)"),
    list(name = "no_capa", desc = "Ablation 2: No CAPA (delta = 0, no coll overlap)"),
    list(name = "indep_outliers", desc = "Ablation 3: Independent outlier scales (rho = 0)"),
    list(name = "no_ar1", desc = "Ablation 4: No AR(1) (phi = 0)"),
    list(name = "standard_hs", desc = "Ablation 5: Horseshoe instead of Horseshoe+"),
    list(name = "marginal_var", desc = "Ablation 6: Marginal variance in scoring")
  )

  cat("\n====================================================================\n")
  cat("     COMPONENT ABLATION STUDY (Section 5, Phase P9)\n")
  cat("====================================================================\n\n")

  # Compute standard front-ends
  intervals <- run_frontend_nsp(y, alpha = 0.05)
  capa_res  <- run_frontend_capa(y, prewhiten = TRUE)

  fits <- list()
  rows <- list()

  for (i in seq_along(ablation_defs)) {
    ab <- ablation_defs[[i]]
    cat(sprintf("Running %s... ", ab$desc))
    t0 <- proc.time()[3]

    # Configure hyperparameters according to ablation
    k_in    <- if (ab$name == "no_nsp") 1.0 else 2.0
    delta   <- if (ab$name == "no_capa") 0.0 else 1.0
    coll    <- if (ab$name == "no_capa") NULL else capa_res$collective
    dynamic <- if (ab$name == "indep_outliers") FALSE else TRUE
    f_phi   <- if (ab$name == "no_ar1") 0.0 else NULL
    hs_plus <- if (ab$name == "standard_hs") FALSE else TRUE
    inn_var <- if (ab$name == "marginal_var") FALSE else TRUE

    pmap <- build_prior_maps(
      n             = n,
      nsp_intervals = if (ab$name == "no_nsp") data.frame(start = integer(0), end = integer(0)) else intervals$intervals,
      capa_results  = if (ab$name == "no_capa") list(point = integer(0), collective = data.frame(start = integer(0), end = integer(0))) else capa_res,
      kappa_in      = k_in
    )

    fit <- gibbs_sampler(
      y               = y,
      kappa           = pmap$kappa,
      a_vec           = pmap$a,
      n_iter          = n_iter,
      burn_in         = burn_in,
      thin            = 1,
      n_chains        = n_chains,
      dynamic         = dynamic,
      delta           = delta,
      capa_collective = coll,
      alpha_fdr       = 0.05,
      seed            = seed + i * 100L,
      fix_phi         = f_phi,
      horseshoe_plus  = hs_plus,
      innovation_var  = inn_var
    )

    t_elapsed <- as.numeric(proc.time()[3] - t0)
    fits[[ab$name]] <- fit

    n_cp <- length(fit$cp)
    n_pt <- length(fit$point)
    n_co <- if (is.data.frame(fit$collective)) nrow(fit$collective) else 0L

    row_data <- list(
      ablation    = ab$name,
      description = ab$desc,
      runtime_sec = fit$runtime_sec,
      n_cp        = n_cp,
      n_point     = n_pt,
      n_coll      = n_co
    )

    # If ground truth is provided, evaluate metrics
    if (!is.null(truth)) {
      ev <- evaluate_all_metrics(fit, truth)
      row_data$cp_f1       <- ev$cp$f1
      row_data$cp_covering <- ev$cp$covering
      row_data$pt_f1       <- ev$point$f1
      row_data$macro_f1    <- ev$four_class$macro_f1
    }

    rows[[i]] <- as.data.frame(row_data, stringsAsFactors = FALSE)
    cat(sprintf("done (%.2fs) | CPs: %d, Pts: %d, Coll: %d\n", t_elapsed, n_cp, n_pt, n_co))
  }

  res_df <- do.call(rbind, rows)
  rownames(res_df) <- NULL

  cat("\n====================================================================\n")
  cat("     ABLATION STUDY COMPLETE\n")
  cat("====================================================================\n\n")
  print(res_df)

  list(
    results_df = res_df,
    fits       = fits,
    timestamp  = Sys.time()
  )
}

#' Run Sensitivity Analysis Grid Across Hyperparameters
#'
#' Evaluates sensitivity across:
#' - alpha (BFDR): {0.05, 0.10}
#' - alpha_nsp: {0.05, 0.10}
#' - kappa_in: {1.0, 2.0, 5.0}
#' - delta: {0.0, 1.0, 2.0}
#' - s_tau: default (1/sqrt(n)) vs p0-based
#' - d_omega, d_psi: {1.0, 1.5, 2.0}
#' - L_max: {10, 20, 30}
#'
#' @param y Numeric vector of observations
#' @param truth Optional ground truth list
#' @param n_iter Gibbs iterations
#' @param burn_in Gibbs burn-in
#' @param seed Master random seed
#'
#' @return A tidy data frame with sensitivity metrics
#' @export
run_sensitivity_grid <- function(
    y,
    truth = NULL,
    n_iter = 1500,
    burn_in = 500,
    seed = 42
) {
  n <- length(y)
  stopifnot(n >= 20, is.numeric(y))

  # Define parameter sweep variations relative to baseline:
  # Baseline: alpha = 0.05, alpha_nsp = 0.05, kappa_in = 2.0, delta = 1.0,
  #           s_tau = "default", d = 1.5, L_max = 20
  base_params <- list(
    param_sweep = "baseline",
    val         = "baseline",
    alpha       = 0.05,
    alpha_nsp   = 0.05,
    kappa_in    = 2.0,
    delta       = 1.0,
    s_tau_type  = "default",
    d_omega     = 1.5,
    d_psi       = 1.5,
    L_max       = 20L
  )

  sweeps <- list()
  sweeps[[1]] <- base_params

  # 1. Sweep alpha: 0.10
  sweeps[[2]] <- modifyList(base_params, list(param_sweep = "alpha", val = "0.10", alpha = 0.10))

  # 2. Sweep alpha_nsp: 0.10
  sweeps[[3]] <- modifyList(base_params, list(param_sweep = "alpha_nsp", val = "0.10", alpha_nsp = 0.10))

  # 3. Sweep kappa_in: 1.0, 5.0
  sweeps[[4]] <- modifyList(base_params, list(param_sweep = "kappa_in", val = "1.0", kappa_in = 1.0))
  sweeps[[5]] <- modifyList(base_params, list(param_sweep = "kappa_in", val = "5.0", kappa_in = 5.0))

  # 4. Sweep delta: 0.0, 2.0
  sweeps[[6]] <- modifyList(base_params, list(param_sweep = "delta", val = "0.0", delta = 0.0))
  sweeps[[7]] <- modifyList(base_params, list(param_sweep = "delta", val = "2.0", delta = 2.0))

  # 5. Sweep s_tau: p0-based (p0 = 0.02)
  sweeps[[8]] <- modifyList(base_params, list(param_sweep = "s_tau", val = "p0_based", s_tau_type = "p0"))

  # 6. Sweep d_omega / d_psi: 1.0, 2.0
  sweeps[[9]]  <- modifyList(base_params, list(param_sweep = "d_tail", val = "1.0", d_omega = 1.0, d_psi = 1.0))
  sweeps[[10]] <- modifyList(base_params, list(param_sweep = "d_tail", val = "2.0", d_omega = 2.0, d_psi = 2.0))

  # 7. Sweep L_max: 10, 30
  sweeps[[11]] <- modifyList(base_params, list(param_sweep = "L_max", val = "10", L_max = 10L))
  sweeps[[12]] <- modifyList(base_params, list(param_sweep = "L_max", val = "30", L_max = 30L))

  cat("\n====================================================================\n")
  cat("     SENSITIVITY ANALYSIS GRID (Section 5, Phase P9)\n")
  cat(sprintf("     Running %d parameter configurations\n", length(sweeps)))
  cat("====================================================================\n\n")

  # Precompute front-ends where possible
  intervals_005 <- run_frontend_nsp(y, alpha = 0.05)
  intervals_010 <- run_frontend_nsp(y, alpha = 0.10)
  capa_res      <- run_frontend_capa(y, prewhiten = TRUE)

  results_list <- list()

  for (s_idx in seq_along(sweeps)) {
    cfg <- sweeps[[s_idx]]
    cat(sprintf("[%2d/%d] Sweep %-12s (val: %-8s)... ",
                s_idx, length(sweeps), cfg$param_sweep, cfg$val))

    t0 <- proc.time()[3]

    # Select NSP intervals based on alpha_nsp
    ints <- if (cfg$alpha_nsp == 0.10) intervals_010 else intervals_005
    pmap <- build_prior_maps(n = n, nsp_intervals = ints$intervals,
                             capa_results = capa_res, kappa_in = cfg$kappa_in)

    # s_tau calculation
    s_tau_val <- if (cfg$s_tau_type == "p0") {
      p0 <- 0.02
      (p0 / (1.0 - p0)) * (1.0 / sqrt(n))
    } else {
      1.0 / sqrt(n)
    }

    fit <- gibbs_sampler(
      y               = y,
      kappa           = pmap$kappa,
      a_vec           = pmap$a,
      n_iter          = n_iter,
      burn_in         = burn_in,
      thin            = 1,
      n_chains        = 1,
      dynamic         = TRUE,
      delta           = cfg$delta,
      alpha_fdr       = cfg$alpha,
      capa_collective = capa_res$collective,
      L_max           = cfg$L_max,
      s_tau           = s_tau_val,
      d_omega         = cfg$d_omega,
      d_psi           = cfg$d_psi,
      seed            = seed + s_idx * 50L
    )

    t_elapsed <- as.numeric(proc.time()[3] - t0)

    n_cp <- length(fit$cp)
    n_pt <- length(fit$point)
    n_co <- if (is.data.frame(fit$collective)) nrow(fit$collective) else 0L

    res_row <- list(
      sweep_id    = s_idx,
      param       = cfg$param_sweep,
      val         = cfg$val,
      n_cp        = n_cp,
      n_point     = n_pt,
      n_coll      = n_co,
      runtime_sec = fit$runtime_sec
    )

    if (!is.null(truth)) {
      ev <- evaluate_all_metrics(fit, truth)
      res_row$cp_f1    <- ev$cp$f1
      res_row$pt_f1    <- ev$point$f1
      res_row$macro_f1 <- ev$four_class$macro_f1
    }

    results_list[[s_idx]] <- as.data.frame(res_row, stringsAsFactors = FALSE)
    cat(sprintf("done (%.2fs) | CPs: %d, Pts: %d, Coll: %d\n",
                t_elapsed, n_cp, n_pt, n_co))
  }

  summary_grid <- do.call(rbind, results_list)
  rownames(summary_grid) <- NULL

  cat("\n====================================================================\n")
  cat("     SENSITIVITY ANALYSIS COMPLETE\n")
  cat("====================================================================\n\n")
  print(summary_grid)

  summary_grid
}
