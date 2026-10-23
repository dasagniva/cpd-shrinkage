# ==============================================================================
# 60_baselines.R
# Baseline Change Point and Anomaly Detection Wrappers.
# Implements Section 5 (Phase P7) and Section 4.7 (Common cpd_fit output schema).
# ==============================================================================

if (file.exists("R/00_utils.R")) {
  source("R/00_utils.R")
  source("R/20_frontend_nsp.R")
  source("R/21_frontend_capa.R")
  source("R/32_gibbs_updates.R")
  source("R/40_scoring.R")
  source("R/41_bfdr.R")
  source("R/42_classify.R")
  source("R/50_copula.R")
  source("R/61_fpnn.R")
} else if (file.exists("../../R/00_utils.R")) {
  source("../../R/00_utils.R")
  source("../../R/20_frontend_nsp.R")
  source("../../R/21_frontend_capa.R")
  source("../../R/32_gibbs_updates.R")
  source("../../R/40_scoring.R")
  source("../../R/41_bfdr.R")
  source("../../R/42_classify.R")
  source("../../R/50_copula.R")
  source("../../R/61_fpnn.R")
}

suppressPackageStartupMessages({
  library(DeCAFS)
  library(dsp)
  library(nsp)
  library(anomaly)
  library(changepoint)
  library(stepR)
  library(bcp)
  library(EnvCpt)
})

# ------------------------------------------------------------------------------
# 1. DeCAFS
# ------------------------------------------------------------------------------
#' DeCAFS Baseline Wrapper
#' @param y Numeric vector of time series observations
#' @param beta Penalty parameter (default 2 * log(length(y)))
#' @param ... Additional arguments passed to `DeCAFS::DeCAFS`
#' @return A validated `cpd_fit` object
#' @export
baseline_decafs <- function(y, beta = NULL, ...) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 5)

  if (is.null(beta)) beta <- 2.0 * log(n)

  fit <- DeCAFS::DeCAFS(y, beta = beta, warningMessage = FALSE, ...)
  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  cps <- as.integer(fit$changepoints)
  cps <- sort(unique(cps[cps >= 2L & cps <= n]))

  class_vec <- rep("nominal", n)
  if (length(cps) > 0) class_vec[cps] <- "shift"

  new_cpd_fit(
    method      = "decafs",
    cp          = cps,
    cp_int      = NULL,
    point       = integer(0),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = NULL,
    runtime_sec = runtime_sec,
    diag        = list(raw = fit, signal = fit$signal, modelParam = fit$modelParameters)
  )
}

# ------------------------------------------------------------------------------
# 2. ABCO (Adaptive Bayesian Changepoint Analysis and Outlier Detection)
# ------------------------------------------------------------------------------
#' ABCO Baseline Wrapper
#' @param y Numeric vector of time series observations
#' @param nsave Number of MCMC samples saved (default 500)
#' @param nburn Number of burn-in iterations (default 500)
#' @param nskip Thinning interval (default 1)
#' @param alpha_fdr FDR cutoff for change points and outliers (default 0.05)
#' @param ... Additional arguments passed to `dsp::dsp_fit`
#' @return A validated `cpd_fit` object
#' @export
baseline_abco <- function(y, nsave = 500, nburn = 500, nskip = 1, alpha_fdr = 0.05, ...) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 5)

  # Model specification: Gaussian changepoint with D = 1 (piecewise constant) and anomaly detection
  spec <- dsp::dsp_spec(family = "gaussian", model = "changepoint", D = 1, useAnom = TRUE)

  fit <- dsp::dsp_fit(
    y          = y,
    model_spec = spec,
    nsave      = nsave,
    nburn      = nburn,
    nskip      = nskip,
    verbose    = FALSE,
    ...
  )
  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  mu_draws <- fit$mcmc_output$mu # (nsave x n)
  mu_hat <- colMeans(mu_draws)

  # Increments omega from ABCO
  omega_draws <- fit$mcmc_output$omega # (nsave x (n-1))
  sigma_est <- stats::mad(y - mu_hat)
  if (sigma_est <= 1e-8) sigma_est <- 1.0
  d_w <- 1.5 * sigma_est
  p_w <- colMeans(abs(omega_draws) > d_w)

  # Observation variance for anomaly scoring
  obs_s2 <- fit$mcmc_output$obs_sigma_t2 # (nsave x n)
  med_s2 <- median(obs_s2)
  if (is.na(med_s2) || med_s2 <= 1e-8) med_s2 <- 1.0
  s_out <- colMeans(obs_s2 / (obs_s2 + med_s2))
  p_psi <- colMeans(obs_s2 > 2.0 * med_s2)

  clf <- classify_and_build_fit(
    y           = y,
    mu_hat      = mu_hat,
    s_out       = s_out,
    p_omega     = p_w,
    p_psi       = p_psi,
    alpha_fdr   = alpha_fdr,
    method      = "abco",
    runtime_sec = runtime_sec,
    diag        = list(raw = fit, mu_hat = mu_hat)
  )

  clf
}

# ------------------------------------------------------------------------------
# 3. NSP (Narrowest Significance Pursuit)
# ------------------------------------------------------------------------------
#' NSP Baseline Wrapper
#' @param y Numeric vector of time series observations
#' @param alpha Significance level (default 0.05)
#' @param ... Additional arguments passed to `run_frontend_nsp` (nsp::nsp_poly_ar)
#' @return A validated `cpd_fit` object
#' @export
baseline_nsp <- function(y, alpha = 0.05, ...) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 5)

  # Same NSP variant as the proposed method's front-end (nsp_poly_ar, spec change C1);
  # nsp_poly_selfnorm is miscalibrated under AR(1) noise (gate P2)
  fit <- run_frontend_nsp(y, alpha = alpha, method = "ar", ...)
  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  ints <- fit$intervals
  if (is.null(ints) || nrow(ints) == 0) {
    cp_ints <- data.frame(start = integer(0), end = integer(0))
    cps <- integer(0)
  } else {
    cp_ints <- data.frame(start = as.integer(ints$start), end = as.integer(ints$end))
    # Change point estimate is midpoint of significant interval
    cps <- as.integer(round((cp_ints$start + cp_ints$end) / 2.0))
    cps <- sort(unique(cps[cps >= 2L & cps <= n]))
  }

  class_vec <- rep("nominal", n)
  if (length(cps) > 0) class_vec[cps] <- "shift"

  new_cpd_fit(
    method      = "nsp",
    cp          = cps,
    cp_int      = cp_ints,
    point       = integer(0),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = NULL,
    runtime_sec = runtime_sec,
    diag        = list(raw = fit, n_intervals = nrow(cp_ints))
  )
}

# ------------------------------------------------------------------------------
# 4. CAPA (Collective and Point Anomaly Detection)
# ------------------------------------------------------------------------------
#' CAPA Baseline Wrapper
#' @param y Numeric vector of time series observations
#' @param prewhiten Logical: prewhiten before CAPA (default FALSE for raw CAPA alone baseline).
#'   Either way the series is robustly centred and scaled, as `anomaly::capa` assumes
#'   standardised data.
#' @param ... Additional arguments passed to `anomaly::capa`
#' @return A validated `cpd_fit` object
#' @export
baseline_capa <- function(y, prewhiten = FALSE, ...) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 10)

  res <- run_frontend_capa(y, prewhiten = isTRUE(prewhiten), ...)
  point_locs <- res$point
  coll_df <- res$collective
  diag_list <- list(phi_hat = res$phi_hat, raw = res$raw)
  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  class_vec <- rep("nominal", n)
  if (nrow(coll_df) > 0) {
    for (r in seq_len(nrow(coll_df))) {
      class_vec[coll_df$start[r]:coll_df$end[r]] <- "collective"
    }
  }
  if (length(point_locs) > 0) {
    pt_only <- point_locs[class_vec[point_locs] == "nominal"]
    if (length(pt_only) > 0) class_vec[pt_only] <- "point"
  }

  new_cpd_fit(
    method      = "capa",
    cp          = integer(0),
    cp_int      = NULL,
    point       = as.integer(point_locs),
    collective  = coll_df,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = NULL,
    runtime_sec = runtime_sec,
    diag        = diag_list
  )
}

# ------------------------------------------------------------------------------
# 5. PELT
# ------------------------------------------------------------------------------
#' PELT Baseline Wrapper
#'
#' `cpt.mean` with the Normal cost assumes independent noise with unit variance, so the
#' series is scaled by its robust marginal noise SD (robust_marginal_scale) first.
#'
#' @param y Numeric vector of time series observations
#' @param penalty Penalty type (default "MBIC")
#' @param ... Additional arguments passed to `changepoint::cpt.mean`
#' @return A validated `cpd_fit` object
#' @export
baseline_pelt <- function(y, penalty = "MBIC", ...) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 5)

  fit <- changepoint::cpt.mean((y - stats::median(y)) / robust_marginal_scale(y), method = "PELT", penalty = penalty, ...)
  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  # In changepoint, cpts gives segment end indices. For a change between t and t+1, end is t.
  # So change point starts at t + 1. Filter out n (last observation).
  raw_cpts <- changepoint::cpts(fit)
  raw_cpts <- raw_cpts[raw_cpts < n]
  cps <- as.integer(raw_cpts + 1L)
  cps <- sort(unique(cps[cps >= 2L & cps <= n]))

  class_vec <- rep("nominal", n)
  if (length(cps) > 0) class_vec[cps] <- "shift"

  new_cpd_fit(
    method      = "pelt",
    cp          = cps,
    cp_int      = NULL,
    point       = integer(0),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = NULL,
    runtime_sec = runtime_sec,
    diag        = list(raw = fit, param_est = fit@param.est)
  )
}

# ------------------------------------------------------------------------------
# 6. Robust FPOP (Biweight)
# ------------------------------------------------------------------------------
#' Robust FPOP Baseline Wrapper
#' @param y Numeric vector of time series observations
#' @param penalty Penalty per change on the standardised scale (default 2 * log(n))
#' @param min_seg_len Minimum segment length for the fallback (default 10)
#' @return A validated `cpd_fit` object
#' @export
baseline_robseg <- function(y, penalty = NULL, min_seg_len = 10) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= min_seg_len)

  fit <- fit_robust_location(y, min_seg_len = min_seg_len, penalty = penalty)
  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  cps <- as.integer(fit$cp)
  cps <- sort(unique(cps[cps >= 2L & cps <= n]))

  class_vec <- rep("nominal", n)
  if (length(cps) > 0) class_vec[cps] <- "shift"

  new_cpd_fit(
    method      = fit$method,
    cp          = cps,
    cp_int      = NULL,
    point       = integer(0),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = NULL,
    runtime_sec = runtime_sec,
    diag        = list(mu_hat = fit$mu_hat, residuals = fit$residuals)
  )
}

# ------------------------------------------------------------------------------
# 7. HSMUCE
# ------------------------------------------------------------------------------
#' HSMUCE Baseline Wrapper
#' @param y Numeric vector of time series observations
#' @param alpha Significance level (default 0.05)
#' @param ... Additional arguments passed to `stepR::stepFit`
#' @return A validated `cpd_fit` object
#' @export
baseline_hsmuce <- function(y, alpha = 0.05, ...) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 5)

  fit <- stepR::stepFit(y, alpha = alpha, family = "hsmuce", ...)
  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  # leftEnd of segments 2..K are the change locations
  cps <- as.integer(fit$leftEnd[-1])
  cps <- sort(unique(cps[cps >= 2L & cps <= n]))

  class_vec <- rep("nominal", n)
  if (length(cps) > 0) class_vec[cps] <- "shift"

  new_cpd_fit(
    method      = "hsmuce",
    cp          = cps,
    cp_int      = NULL,
    point       = integer(0),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = NULL,
    runtime_sec = runtime_sec,
    diag        = list(raw = fit, intervals = fit)
  )
}

# ------------------------------------------------------------------------------
# 8. BCP (Bayesian Product Partition Model)
# ------------------------------------------------------------------------------
#' BCP Baseline Wrapper
#' @param y Numeric vector of time series observations
#' @param p0 Prior probability of change (default 0.05)
#' @param burnin Burn-in iterations (default 200)
#' @param mcmc Retained MCMC iterations (default 500)
#' @param threshold Posterior probability threshold for change points (default 0.5)
#' @param ... Additional arguments passed to `bcp::bcp`
#' @return A validated `cpd_fit` object
#' @export
baseline_bcp <- function(y, p0 = 0.05, burnin = 200, mcmc = 500, threshold = 0.5, ...) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 5)

  fit <- bcp::bcp(y, p0 = p0, burnin = burnin, mcmc = mcmc, ...)
  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  post_p <- as.numeric(fit$posterior.prob)
  cps <- which(post_p >= threshold)
  cps <- sort(unique(cps[cps >= 2L & cps <= n]))

  class_vec <- rep("nominal", n)
  if (length(cps) > 0) class_vec[cps] <- "shift"

  scores_df <- data.frame(
    t     = seq_len(n),
    s_cp  = post_p,
    s_out = rep(NA_real_, n)
  )

  new_cpd_fit(
    method      = "bcp",
    cp          = as.integer(cps),
    cp_int      = NULL,
    point       = integer(0),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = scores_df,
    runtime_sec = runtime_sec,
    diag        = list(raw = fit, posterior_mean = as.numeric(fit$posterior.mean))
  )
}

# ------------------------------------------------------------------------------
# 9. EnvCpt
# ------------------------------------------------------------------------------
#' EnvCpt Baseline Wrapper
#' @param y Numeric vector of time series observations
#' @param ... Additional arguments passed to `EnvCpt::envcpt`
#' @return A validated `cpd_fit` object
#' @export
baseline_envcpt <- function(y, ...) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 15)

  fit <- EnvCpt::envcpt(y, verbose = FALSE, ...)
  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  aic_vals <- stats::AIC(fit)
  best_mod_name <- names(which.min(aic_vals))

  cps <- integer(0)
  if (grepl("cpt", best_mod_name)) {
    cpt_obj <- fit[[best_mod_name]]
    raw_cpts <- changepoint::cpts(cpt_obj)
    raw_cpts <- raw_cpts[raw_cpts < n]
    cps <- as.integer(raw_cpts + 1L)
    cps <- sort(unique(cps[cps >= 2L & cps <= n]))
  }

  class_vec <- rep("nominal", n)
  if (length(cps) > 0) class_vec[cps] <- "shift"

  new_cpd_fit(
    method      = "envcpt",
    cp          = cps,
    cp_int      = NULL,
    point       = integer(0),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = NULL,
    runtime_sec = runtime_sec,
    diag        = list(best_model = best_mod_name, aic = aic_vals)
  )
}

# ------------------------------------------------------------------------------
# 10. Chapter 4 Framework (DeCAFS + HS+ Residual Scoring)
# ------------------------------------------------------------------------------
#' Chapter 4 Framework Baseline Wrapper
#'
#' DeCAFS segmentation followed by horseshoe scoring of the prewhitened, standardised
#' DeCAFS residuals (approximate port of the thesis Chapter 4 framework). Points are
#' flagged by Bayesian FDR on p_t = Pr(|psi_t| > d_psi | residuals), the same
#' exceedance probability used by the proposed method (Eq. 9).
#'
#' @param y Numeric vector of time series observations
#' @param alpha_fdr Bayesian FDR level for residual scoring (default 0.05)
#' @param n_iter Total Gibbs iterations for residual horseshoe (default 1000)
#' @param burn_in Burn-in iterations (default 500)
#' @param d_psi Exceedance threshold in residual SD units (default 1.5)
#' @return A validated `cpd_fit` object
#' @export
baseline_ch4 <- function(y, alpha_fdr = 0.05, n_iter = 1000, burn_in = 500, d_psi = 1.5) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n >= 10)

  # Step 1: DeCAFS segmentation
  decafs_fit <- DeCAFS::DeCAFS(y, warningMessage = FALSE)
  cp_locs <- as.integer(decafs_fit$changepoints)
  cp_locs <- sort(unique(cp_locs[cp_locs >= 2L & cp_locs <= n]))
  mu_hat <- decafs_fit$signal

  # Step 2: prewhiten and standardise residuals
  res <- y - mu_hat
  phi_hat <- as.numeric(decafs_fit$modelParameters$phi)
  if (length(phi_hat) != 1 || is.na(phi_hat) || abs(phi_hat) >= 0.95) phi_hat <- 0.0
  res_pw <- c(sqrt(1 - phi_hat^2) * res[1], res[-1] - phi_hat * res[-n])
  sigma_hat <- stats::mad(res_pw)
  if (sigma_hat <= 1e-12) sigma_hat <- 1.0
  r_std <- res_pw / sigma_hat

  # Step 3: horseshoe on psi in r_std = psi + N(0, 1), global scale fixed at 1/sqrt(n)
  lambda_sq <- rep(1.0, n)
  tau_p <- 1.0 / sqrt(n)
  o_t_accum <- numeric(n)
  exceed <- numeric(n)
  n_saved <- 0L

  for (iter in seq_len(n_iter)) {
    nu2 <- tau_p^2 * lambda_sq
    shrink <- nu2 / (nu2 + 1.0)
    psi_draw <- rnorm(n, mean = shrink * r_std, sd = sqrt(shrink))

    aux <- 1.0 / rgamma(n, shape = 1.0, rate = 1.0 + 1.0 / lambda_sq)
    lambda_sq <- 1.0 / rgamma(n, shape = 1.0, rate = 1.0 / aux + psi_draw^2 / (2.0 * tau_p^2))

    if (iter > burn_in) {
      n_saved <- n_saved + 1L
      o_t_accum <- o_t_accum + shrink
      exceed <- exceed + (abs(psi_draw) > d_psi)
    }
  }

  o_t <- o_t_accum / max(1L, n_saved)
  p_psi <- exceed / max(1L, n_saved)
  flagged_points <- setdiff(bayesian_fdr(p_psi, alpha = alpha_fdr)$flagged, cp_locs)

  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  class_vec <- rep("nominal", n)
  if (length(cp_locs) > 0) class_vec[cp_locs] <- "shift"
  if (length(flagged_points) > 0) class_vec[flagged_points] <- "point"

  new_cpd_fit(
    method      = "ch4",
    cp          = cp_locs,
    cp_int      = NULL,
    point       = as.integer(flagged_points),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = data.frame(t = seq_len(n), s_cp = rep(NA_real_, n), s_out = o_t),
    runtime_sec = runtime_sec,
    diag        = list(decafs = decafs_fit, o_t = o_t, p_psi = p_psi)
  )
}

# ------------------------------------------------------------------------------
# 11. FPNN (JISPS two-phase EV-DeCAFS + Fourier PNN; implementation in 61_fpnn.R)
# ------------------------------------------------------------------------------
#' FPNN Baseline (JISPS EV-DeCAFS + Fourier Probabilistic Neural Network)
#'
#' Port of the ev-decafs-serve training and inference pipeline, run without ground truth:
#' Phase I finds changes with EV-DeCAFS under the flat BIC penalty C log n (C = 2, the
#' untuned default); Phase II labels them with the BOCPD oracle and Algorithm 3, balances
#' the classes with SMOTE, fits the Fourier PNN on the series' own changes and keeps those
#' it classifies as sustained. Recoiled changes are not shifts; the JISPS method gives no
#' point or collective labels, so they are only reported in `diag$recoiled`. With fewer
#' than two Phase-I changes the classifier cannot be trained and the Phase-I changes are
#' returned as they are (`diag$fpnn_trained = FALSE`).
#'
#' @param y Numeric vector of time series observations
#' @param J Number of Fourier harmonics (default 10)
#' @param bic_multiplier BIC multiplier C in alpha_0 = C log n (default 2)
#' @param n_grid Level grid size of the EV-DeCAFS dynamic programme (default 1000)
#' @param window_L Feature window half-width (default 5)
#' @param evi_w,evi_q0 EVI field window half-width and exceedance quantile (50, 0.90)
#' @param kappa_mu_percentile,kappa_S Algorithm-3 thresholds (75, 0.5)
#' @param bocpd_threshold,bocpd_tolerance_fraction BOCPD oracle threshold and match window
#'   as a fraction of n (0.3, 0.02)
#' @param smote_k SMOTE neighbours (default 5)
#' @param seed Seed for SMOTE; the caller's RNG state is restored afterwards
#' @return A validated `cpd_fit` object
#' @export
baseline_fpnn <- function(y, J = 10L, bic_multiplier = 2.0, n_grid = 1000L, window_L = 5L,
                          evi_w = 50L, evi_q0 = 0.90, kappa_mu_percentile = 75, kappa_S = 0.5,
                          bocpd_threshold = 0.3, bocpd_tolerance_fraction = 0.02,
                          smote_k = 5L, seed = 42L) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(is.numeric(y), n > evi_w)

  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
  set.seed(as.integer(seed))

  # Phase I
  ar1 <- fpnn_ar1_params(y)
  lambda_param <- 1.0 / (ar1$sigma_eta_sq + 1e-12)
  gamma <- 1.0 / (ar1$sigma_v_sq + 1e-12)
  xi_field <- fpnn_evi_field(y, w = evi_w, q0 = evi_q0)
  alpha_0 <- bic_multiplier * log(n)
  dp <- ev_decafs_dp_cpp(y, rep(alpha_0, n), lambda_param, gamma, ar1$phi, as.integer(n_grid))
  cps_all <- as.integer(dp$changepoints)

  # Phase II
  trained <- length(cps_all) >= 2
  prob_sustained <- rep(1.0, length(cps_all))
  labels_train <- NULL
  if (trained) {
    X <- fpnn_extract_features(y, cps_all, dp$means, L = window_L, xi_field = xi_field)
    bocpd_cps <- which(fpnn_bocpd(y, phi = ar1$phi, threshold = bocpd_threshold))
    labels_train <- fpnn_make_labels(X, cps_all, bocpd_cps, tol = floor(bocpd_tolerance_fraction * n),
                                     kappa_mu_percentile = kappa_mu_percentile, kappa_S = kappa_S)
    bal <- fpnn_smote(X, labels_train, k = smote_k)
    model <- fourier_pnn_fit(bal$X, bal$y, J = J)
    proba <- fourier_pnn_predict_proba(model, X)
    pred <- model$classes[apply(proba, 1, which.max)]
    prob_sustained <- if ("1" %in% colnames(proba)) proba[, "1"] else rep(0.0, length(cps_all))
    cps <- cps_all[pred == 1L]
  } else {
    cps <- cps_all
  }

  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  class_vec <- rep("nominal", n)
  class_vec[cps] <- "shift"
  s_cp <- numeric(n)
  s_cp[cps_all] <- prob_sustained

  new_cpd_fit(
    method      = "fpnn",
    cp          = cps,
    cp_int      = NULL,
    point       = integer(0),
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = data.frame(t = seq_len(n), s_cp = s_cp, s_out = rep(NA_real_, n)),
    runtime_sec = runtime_sec,
    diag        = list(phase1_cp = cps_all, recoiled = setdiff(cps_all, cps), fpnn_trained = trained,
                       labels_train = labels_train, alpha_0 = alpha_0, ar1 = ar1, J = J)
  )
}

# ------------------------------------------------------------------------------
# 12. Master Dispatcher
# ------------------------------------------------------------------------------
#' Master Baseline Dispatcher
#'
#' Fits a specified baseline method and returns a strictly validated `cpd_fit`.
#'
#' @param method One of "decafs", "abco", "nsp", "capa", "pelt", "robseg",
#'   "hsmuce", "bcp", "envcpt", "ch4", "fpnn", "copula"
#' @param y Numeric vector of time series observations
#' @param ... Arguments forwarded to the specific baseline wrapper
#'
#' @return A validated `cpd_fit` object
#' @export
fit_baseline <- function(method, y, ...) {
  method <- match.arg(tolower(method), c(
    "decafs", "abco", "nsp", "capa", "pelt", "robseg",
    "hsmuce", "bcp", "envcpt", "ch4", "fpnn", "copula"
  ))

  fit <- switch(
    method,
    "decafs" = baseline_decafs(y, ...),
    "abco"   = baseline_abco(y, ...),
    "nsp"    = baseline_nsp(y, ...),
    "capa"   = baseline_capa(y, ...),
    "pelt"   = baseline_pelt(y, ...),
    "robseg" = baseline_robseg(y, ...),
    "hsmuce" = baseline_hsmuce(y, ...),
    "bcp"    = baseline_bcp(y, ...),
    "envcpt" = baseline_envcpt(y, ...),
    "ch4"    = baseline_ch4(y, ...),
    "fpnn"   = baseline_fpnn(y, ...),
    "copula" = fit_copula_comparator(y, ...)
  )

  validate_cpd_fit(fit, n = length(y))
  fit
}
