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
#' 7. marginal_var: Marginal variance instead of innovation variance in o_t (Eq. 8).
#'    o_t only enters the threshold-free point metric (PR-AUC), not the labels (which use
#'    p^psi_t), so this ablation is read off the `pt_prauc` column.
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
      row_data$pt_prauc    <- ev$point$pr_auc
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
#' One-at-a-time sweeps around the default configuration (SPEC_CHANGES.md). Every
#' setting changed from the spec is also swept back to its spec value:
#' - alpha (BFDR): 0.05 (default), 0.10
#' - alpha_nsp: 0.05 (default), 0.10
#' - kappa_in: 2 (default), 1, 5
#' - delta: 1 (default), 0, 2
#' - tau_omega: plug-in with p0 = 0.02 (default), p0 = 0.01, 0.05, sampled with s_tau = 1/sqrt(n) (spec)
#' - d_omega: 1.0 (default), 0.75, 1.5 (spec); d_psi: 1.5 (default = spec), 1.0, 2.0
#' - L_max: 50 (default), 30, 100, 0.1 n (spec)
#' - cp_window h: 5 (default), 0 (spec Eq. 9), 10
#' - min_shift_len m: 4 (default), 1 (spec: no persistence split), 8
#'
#' @param y Numeric vector of observations
#' @param truth Optional ground truth list
#' @param n_iter,burn_in Gibbs iterations and burn-in
#' @param n_chains Number of chains per configuration (default 2)
#' @param seed Master random seed
#'
#' @return A tidy data frame with one row per configuration
#' @export
run_sensitivity_grid <- function(
    y,
    truth = NULL,
    n_iter = 1500,
    burn_in = 500,
    n_chains = 2,
    seed = 42
) {
  n <- length(y)
  stopifnot(n >= 20, is.numeric(y))

  base <- list(param = "default", val = "default", alpha = 0.05, alpha_nsp = 0.05, kappa_in = 2.0,
               delta = 1.0, tau_mode = "plugin", p0 = 0.02, d_omega = 1.0, d_psi = 1.5,
               L_max = 50L, cp_window = 5L, min_shift_len = 4L)
  sw <- function(param, val, ...) modifyList(base, c(list(param = param, val = val), list(...)))
  sweeps <- list(
    base,
    sw("alpha", "0.10", alpha = 0.10),
    sw("alpha_nsp", "0.10", alpha_nsp = 0.10),
    sw("kappa_in", "1", kappa_in = 1.0),
    sw("kappa_in", "5", kappa_in = 5.0),
    sw("delta", "0", delta = 0.0),
    sw("delta", "2", delta = 2.0),
    sw("tau_omega", "plug-in p0=0.01", p0 = 0.01),
    sw("tau_omega", "plug-in p0=0.05", p0 = 0.05),
    sw("tau_omega", "sampled (spec)", tau_mode = "sample"),
    sw("d_omega", "0.75", d_omega = 0.75),
    sw("d_omega", "1.5 (spec)", d_omega = 1.5),
    sw("d_psi", "1.0", d_psi = 1.0),
    sw("d_psi", "2.0", d_psi = 2.0),
    sw("L_max", "30", L_max = 30L),
    sw("L_max", "100", L_max = 100L),
    sw("L_max", "0.1n (spec)", L_max = max(2L, as.integer(floor(0.1 * n)))),
    sw("cp_window", "0 (spec)", cp_window = 0L),
    sw("cp_window", "10", cp_window = 10L),
    sw("min_shift_len", "1 (spec)", min_shift_len = 1L),
    sw("min_shift_len", "8", min_shift_len = 8L)
  )

  cat(sprintf("\n     SENSITIVITY ANALYSIS GRID: %d configurations\n\n", length(sweeps)))

  nsp_by_alpha <- list("0.05" = run_frontend_nsp(y, alpha = 0.05), "0.1" = run_frontend_nsp(y, alpha = 0.10))
  capa_res <- run_frontend_capa(y, prewhiten = TRUE)

  results_list <- vector("list", length(sweeps))
  for (s_idx in seq_along(sweeps)) {
    cfg <- sweeps[[s_idx]]
    ints <- nsp_by_alpha[[as.character(cfg$alpha_nsp)]]
    pmap <- build_prior_maps(n = n, nsp_intervals = ints$intervals, capa_results = capa_res, kappa_in = cfg$kappa_in)

    fit <- gibbs_sampler(
      y               = y,
      kappa           = pmap$kappa,
      a_vec           = pmap$a,
      n_iter          = n_iter,
      burn_in         = burn_in,
      n_chains        = n_chains,
      dynamic         = TRUE,
      delta           = cfg$delta,
      alpha_fdr       = cfg$alpha,
      capa_collective = capa_res$collective,
      L_max           = cfg$L_max,
      tau_omega_mode  = cfg$tau_mode,
      p0              = cfg$p0,
      d_omega         = cfg$d_omega,
      d_psi           = cfg$d_psi,
      cp_window       = cfg$cp_window,
      min_shift_len   = cfg$min_shift_len,
      seed            = seed + s_idx * 50L
    )

    res_row <- list(sweep_id = s_idx, param = cfg$param, val = cfg$val,
                    n_cp = length(fit$cp), n_point = length(fit$point), n_coll = nrow(fit$collective),
                    runtime_sec = fit$runtime_sec, max_rhat = fit$diag$max_rhat)
    if (!is.null(truth)) {
      ev <- evaluate_all_metrics(fit, truth)
      res_row$cp_f1    <- ev$cp$f1
      res_row$pt_f1    <- ev$point$f1
      res_row$pt_prauc <- ev$point$pr_auc
      res_row$coll_f1  <- ev$collective$f1
      res_row$macro_f1 <- ev$four_class$macro_f1
    }
    results_list[[s_idx]] <- as.data.frame(res_row, stringsAsFactors = FALSE)
    cat(sprintf("[%2d/%d] %-14s %-16s | CPs %d, Pts %d, Coll %d\n", s_idx, length(sweeps), cfg$param, cfg$val,
                res_row$n_cp, res_row$n_point, res_row$n_coll))
  }

  summary_grid <- do.call(rbind, results_list)
  rownames(summary_grid) <- NULL
  summary_grid
}
