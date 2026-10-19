# ==============================================================================
# 31_gibbs_sampler.R
# Production Gibbs Sampler for Significance-Guided Global-Local Shrinkage (Phase P4).
# Implements Section 3 (Equations 1-9) and Section 5 (Phase P4).
# ==============================================================================

if (file.exists("R/00_utils.R")) {
  source("R/00_utils.R")
  source("R/32_gibbs_updates.R")
} else if (file.exists("../../R/00_utils.R")) {
  source("../../R/00_utils.R")
  source("../../R/32_gibbs_updates.R")
}

suppressPackageStartupMessages({
  library(spam)
  library(posterior)
  library(future.apply)
})

#' Run a Single Chain of the Production Gibbs Sampler
#'
#' @param chain_id Integer chain identifier
#' @param y Numeric vector of length n
#' @param kappa Numeric vector of length n for local shrinkage multipliers
#' @param a_vec Numeric binary vector of length n for CAPA point offsets
#' @param n_iter Total iterations (burn_in + n_keep)
#' @param burn_in Number of initial burn-in iterations
#' @param thin Thinning interval (default 1)
#' @param dynamic Logical: if TRUE, use dynamic Horseshoe for outliers (rho, delta)
#' @param delta CAPA offset parameter for dynamic outliers (default 0.0)
#' @param s_tau Scale parameter for tau_omega prior (default 1 / sqrt(n))
#' @param seed Optional integer seed for this chain
#'
#' @return A list containing retained draws, traces, and summary metrics
run_gibbs_chain <- function(
    chain_id,
    y,
    kappa = NULL,
    a_vec = NULL,
    n_iter = 10000,
    burn_in = 5000,
    thin = 1,
    dynamic = FALSE,
    delta = 0.0,
    s_tau = NULL,
    seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed) + as.integer(chain_id))
  }

  n <- length(y)
  n_omega <- n - 1
  m <- 2 * n

  if (is.null(kappa)) kappa <- rep(1.0, n)
  kappa_w <- kappa[2:n]
  if (is.null(a_vec)) a_vec <- rep(0.0, n)
  if (is.null(s_tau)) s_tau <- 1.0 / sqrt(n)

  # Precompute banded structure
  banded_struct <- init_banded_structure(n)

  # 1. Warm-start initialization (Section 5, P4)
  # Running median for piecewise constant level
  k_med <- min(11L, 2L * as.integer(floor(n / 10L)) + 1L)
  if (k_med %% 2 == 0) k_med <- k_med + 1L
  mu_init <- as.numeric(runmed(y, k = max(3L, k_med)))
  psi_init <- rep(0.0, n)

  res_init <- y - mu_init
  r1 <- if (n > 5) cor(res_init[-1], res_init[-n]) else 0.5
  phi_curr <- pmax(-0.8, pmin(0.8, ifelse(is.finite(r1), r1, 0.5)))
  sigma_curr <- pmax(1e-3, mad(diff(y)) / sqrt(2.0))

  tau_omega_curr <- s_tau
  lambda_omega_curr <- rep(1.0, n_omega)
  eta_omega_curr <- rep(1.0, n_omega)
  aux_lambda_w <- rep(1.0, n_omega)
  aux_eta_w <- rep(1.0, n_omega)
  aux_tau_w <- 1.0
  c2_tilde <- 1.0
  c2 <- 25.0

  # Compute initial v_omega
  eff_lam_raw <- pmin(kappa_w * eta_omega_curr * lambda_omega_curr, 1e4)
  tau_lam <- pmin(tau_omega_curr * eff_lam_raw, 1e4)
  v_omega_tilde <- (tau_omega_curr^2) * (c2 * (eff_lam_raw^2)) / (c2 + (tau_lam^2) + 1e-12)
  v_omega <- (sigma_curr^2) * v_omega_tilde

  # Outlier scales
  tau_psi_curr <- 1.0 / sqrt(n)
  lambda_psi_curr <- rep(1.0, n)
  eta_psi_curr <- rep(1.0, n)
  aux_lambda_psi <- rep(1.0, n)
  aux_eta_psi <- rep(1.0, n)
  aux_tau_psi <- 1.0
  nu2_curr <- rep(tau_psi_curr^2, n)
  dsp_state <- NULL
  rho_curr <- 0.0
  m_h_curr <- 0.0

  mu_curr  <- mu_init
  psi_curr <- psi_init

  # Retained storage allocation
  n_saved <- floor((n_iter - burn_in) / thin)
  stopifnot(n_saved >= 1)

  # Online accumulator vectors (length n or n_omega)
  mu_sum        <- numeric(n)
  mu2_sum       <- numeric(n)
  psi_sum       <- numeric(n)
  o_t_sum       <- numeric(n)
  p_omega_count <- numeric(n_omega)
  p_psi_count   <- numeric(n)

  # Parameter traces
  trace_phi       <- numeric(n_saved)
  trace_sigma     <- numeric(n_saved)
  trace_tau_omega <- numeric(n_saved)
  trace_tau_psi   <- numeric(n_saved)
  trace_rho       <- if (dynamic) numeric(n_saved) else NULL
  trace_m_h       <- if (dynamic) numeric(n_saved) else NULL

  save_idx <- 0L

  # --- Gibbs Sampling Loop ---
  for (iter in 1:n_iter) {
    # Step 1: Joint (mu, psi) Gaussian draw
    mu_psi_res <- update_mu_psi(
      y             = y,
      phi           = phi_curr,
      sigma         = sigma_curr,
      v_omega       = v_omega,
      nu2           = nu2_curr,
      banded_struct = banded_struct
    )
    mu_curr  <- mu_psi_res$mu
    psi_curr <- mu_psi_res$psi
    omega_curr <- diff(mu_curr)

    # Step 2: Update phi
    e_curr <- y - mu_curr - psi_curr
    phi_res <- update_phi(e = e_curr, phi_curr = phi_curr, sigma = sigma_curr)
    phi_curr <- phi_res$phi

    # Step 3: Compute prewhitened innovations u
    u_curr <- numeric(n)
    u_curr[1] <- sqrt(1.0 - phi_curr^2) * e_curr[1]
    u_curr[2:n] <- e_curr[2:n] - phi_curr * e_curr[1:(n - 1)]

    # Step 4: Update sigma^2
    sigma2_curr <- update_sigma2(
      u             = u_curr,
      omega         = omega_curr,
      v_omega_tilde = v_omega_tilde,
      mu1           = mu_curr[1]
    )
    sigma_curr <- sqrt(pmax(1e-12, sigma2_curr))
    sigma_marg_curr <- sigma_curr / sqrt(1.0 - phi_curr^2)

    # Step 5: Update Horseshoe+ scales on increments omega
    hs_res <- update_hsplus_scales(
      omega         = omega_curr,
      sigma         = sigma_curr,
      tau_omega     = tau_omega_curr,
      lambda_omega  = lambda_omega_curr,
      eta_omega     = eta_omega_curr,
      kappa         = kappa_w,
      aux_lambda    = aux_lambda_w,
      aux_eta       = aux_eta_w,
      aux_tau       = aux_tau_w,
      s_tau         = s_tau,
      c2_tilde      = c2_tilde
    )
    tau_omega_curr    <- hs_res$tau_omega
    lambda_omega_curr <- hs_res$lambda_omega
    eta_omega_curr    <- hs_res$eta_omega
    v_omega           <- hs_res$v_omega
    v_omega_tilde     <- hs_res$v_omega_tilde
    aux_lambda_w      <- hs_res$aux_lambda
    aux_eta_w         <- hs_res$aux_eta
    aux_tau_w         <- hs_res$aux_tau
    c2_tilde          <- hs_res$c2_tilde
    c2                <- hs_res$c2

    # Step 6: Update outlier scales (independent or dynamic)
    out_res <- update_outlier_scales(
      psi            = psi_curr,
      sigma          = sigma_curr,
      tau_psi        = tau_psi_curr,
      lambda_psi     = lambda_psi_curr,
      eta_psi        = eta_psi_curr,
      aux_lambda_psi = aux_lambda_psi,
      aux_eta_psi    = aux_eta_psi,
      aux_tau_psi    = aux_tau_psi,
      dynamic        = dynamic,
      rho            = rho_curr,
      m_h            = m_h_curr,
      delta          = delta,
      a_vec          = a_vec,
      dsp_state      = dsp_state,
      c2             = c2
    )
    nu2_curr <- out_res$nu2
    if (!dynamic) {
      tau_psi_curr    <- out_res$tau_psi
      lambda_psi_curr <- out_res$lambda_psi
      eta_psi_curr    <- out_res$eta_psi
      aux_lambda_psi  <- out_res$aux_lambda_psi
      aux_eta_psi     <- out_res$aux_eta_psi
      aux_tau_psi     <- out_res$aux_tau_psi
    } else {
      dsp_state <- out_res$dsp_state
      rho_curr  <- out_res$rho
      m_h_curr  <- out_res$m_h
    }

    # Step 7: Record retained draws
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      save_idx <- save_idx + 1L

      # Traces
      trace_phi[save_idx]       <- phi_curr
      trace_sigma[save_idx]     <- sigma_curr
      trace_tau_omega[save_idx] <- tau_omega_curr
      trace_tau_psi[save_idx]   <- if (!dynamic) tau_psi_curr else sqrt(mean(nu2_curr))
      if (dynamic) {
        trace_rho[save_idx] <- rho_curr
        trace_m_h[save_idx] <- m_h_curr
      }

      # Accumulate online summaries
      mu_sum  <- mu_sum + mu_curr
      mu2_sum <- mu2_sum + mu_curr^2
      psi_sum <- psi_sum + psi_curr

      # Outlier score o_t (Eq. 8): uses innovation variance sigma^2
      o_t_curr <- nu2_curr / (nu2_curr + sigma2_curr)
      o_t_sum  <- o_t_sum + o_t_curr

      # Indicator counts for p_omega and p_psi (Eq. 9): d = 1.5 * sigma_marg
      d_thresh <- 1.5 * sigma_marg_curr
      p_omega_count <- p_omega_count + as.numeric(abs(omega_curr) > d_thresh)
      p_psi_count   <- p_psi_count + as.numeric(abs(psi_curr) > d_thresh)
    }
  }

  list(
    chain_id      = chain_id,
    n_saved       = n_saved,
    mu_sum        = mu_sum,
    mu2_sum       = mu2_sum,
    psi_sum       = psi_sum,
    o_t_sum       = o_t_sum,
    p_omega_count = p_omega_count,
    p_psi_count   = p_psi_count,
    trace_phi     = trace_phi,
    trace_sigma   = trace_sigma,
    trace_tau_omega = trace_tau_omega,
    trace_tau_psi   = trace_tau_psi,
    trace_rho     = trace_rho,
    trace_m_h     = trace_m_h
  )
}

#' Production Gibbs Sampler for CPD Shrinkage
#'
#' Runs multi-chain Gibbs sampler with warm-start, banded Cholesky updates,
#' Makalic-Schmidt Horseshoe+ scales, and dynamic/independent outlier scales.
#'
#' @param y Numeric vector of time series observations
#' @param kappa Local shrinkage multiplier vector of length n (from build_prior_maps)
#' @param a_vec Binary indicator vector of length n for CAPA point offsets
#' @param n_iter Total iterations per chain (default 10000)
#' @param burn_in Number of burn-in iterations per chain (default 5000)
#' @param thin Thinning interval (default 1)
#' @param n_chains Number of parallel chains (default 4)
#' @param dynamic Logical: if TRUE, use dynamic Horseshoe for outliers (rho, delta)
#' @param delta CAPA offset parameter for dynamic outliers (default 0.0)
#' @param alpha_fdr False discovery rate threshold for candidate selection (default 0.05)
#' @param seed Optional master random seed
#'
#' @return A validated `cpd_fit` object
#' @export
gibbs_sampler <- function(
    y,
    kappa = NULL,
    a_vec = NULL,
    n_iter = 10000,
    burn_in = 5000,
    thin = 1,
    n_chains = 4,
    dynamic = FALSE,
    delta = 0.0,
    alpha_fdr = 0.05,
    seed = NULL
) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(n >= 10, is.numeric(y))

  if (is.null(kappa)) kappa <- rep(1.0, n)
  if (is.null(a_vec)) a_vec <- rep(0.0, n)
  stopifnot(length(kappa) == n, length(a_vec) == n)

  # Run chains
  chain_args <- seq_len(n_chains)
  chains_res <- if (n_chains == 1) {
    list(run_gibbs_chain(
      chain_id = 1,
      y        = y,
      kappa    = kappa,
      a_vec    = a_vec,
      n_iter   = n_iter,
      burn_in  = burn_in,
      thin     = thin,
      dynamic  = dynamic,
      delta    = delta,
      seed     = seed
    ))
  } else {
    future.apply::future_lapply(chain_args, function(cid) {
      run_gibbs_chain(
        chain_id = cid,
        y        = y,
        kappa    = kappa,
        a_vec    = a_vec,
        n_iter   = n_iter,
        burn_in  = burn_in,
        thin     = thin,
        dynamic  = dynamic,
        delta    = delta,
        seed     = seed
      )
    }, future.seed = TRUE)
  }

  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  # Combine results across chains
  total_saved <- sum(sapply(chains_res, function(c) c$n_saved))
  n_saved_per_chain <- chains_res[[1]]$n_saved

  mu_mean  <- Reduce("+", lapply(chains_res, function(c) c$mu_sum)) / total_saved
  mu2_mean <- Reduce("+", lapply(chains_res, function(c) c$mu2_sum)) / total_saved
  mu_sd    <- sqrt(pmax(0.0, mu2_mean - mu_mean^2))

  psi_mean <- Reduce("+", lapply(chains_res, function(c) c$psi_sum)) / total_saved
  o_t_mean <- Reduce("+", lapply(chains_res, function(c) c$o_t_sum)) / total_saved

  p_omega <- Reduce("+", lapply(chains_res, function(c) c$p_omega_count)) / total_saved
  p_psi   <- Reduce("+", lapply(chains_res, function(c) c$p_psi_count)) / total_saved

  # Collect parameter traces into draws_array for posterior package
  mat_phi   <- do.call(cbind, lapply(chains_res, function(c) c$trace_phi))
  mat_sigma <- do.call(cbind, lapply(chains_res, function(c) c$trace_sigma))
  mat_tau_w <- do.call(cbind, lapply(chains_res, function(c) c$trace_tau_omega))
  mat_tau_p <- do.call(cbind, lapply(chains_res, function(c) c$trace_tau_psi))

  # Build 3D array: iterations x chains x variables
  var_names <- c("phi", "sigma", "tau_omega", "tau_psi")
  n_vars <- length(var_names)
  draws_arr <- array(
    NA_real_,
    dim = c(n_saved_per_chain, n_chains, n_vars),
    dimnames = list(
      iteration = seq_len(n_saved_per_chain),
      chain     = seq_len(n_chains),
      variable  = var_names
    )
  )

  draws_arr[, , "phi"]       <- mat_phi
  draws_arr[, , "sigma"]     <- mat_sigma
  draws_arr[, , "tau_omega"] <- mat_tau_w
  draws_arr[, , "tau_psi"]   <- mat_tau_p

  draws_obj <- as_draws_array(draws_arr)
  summ_diag <- summarize_draws(draws_obj)

  # Bayesian FDR selection (Section 3.5)
  fdr_select <- function(p_vec, alpha) {
    k_cand <- integer(0)
    if (length(p_vec) == 0) return(k_cand)
    ord <- order(p_vec, decreasing = TRUE)
    cum_fdr <- cumsum(1.0 - p_vec[ord]) / seq_along(ord)
    valid_k <- which(cum_fdr <= alpha)
    if (length(valid_k) > 0) {
      max_k <- max(valid_k)
      k_cand <- ord[seq_len(max_k)]
    }
    sort(k_cand)
  }

  flagged_cp <- fdr_select(p_omega, alpha_fdr)
  cp_locs <- if (length(flagged_cp) > 0) flagged_cp + 1L else integer(0)

  flagged_point <- fdr_select(p_psi, alpha_fdr)

  # Four-class labelling rule (Section 3.6)
  class_vec <- rep("nominal", n)
  if (length(cp_locs) > 0) {
    class_vec[cp_locs] <- "shift"
  }
  if (length(flagged_point) > 0) {
    pt_only <- setdiff(flagged_point, cp_locs)
    if (length(pt_only) > 0) {
      class_vec[pt_only] <- "point"
    }
  }

  scores_df <- data.frame(
    t     = seq_len(n),
    s_cp  = c(0.0, p_omega),
    s_out = o_t_mean
  )

  new_cpd_fit(
    method      = "proposed",
    cp          = cp_locs,
    cp_int      = NULL,
    point       = flagged_point,
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = scores_df,
    runtime_sec = runtime_sec,
    diag        = list(
      summary   = summ_diag,
      max_rhat  = max(summ_diag$rhat, na.rm = TRUE),
      min_ess   = min(summ_diag$ess_bulk, na.rm = TRUE),
      phi_hat   = mean(mat_phi),
      sigma_hat = mean(mat_sigma),
      mu_mean   = mu_mean,
      mu_sd     = mu_sd,
      psi_mean  = psi_mean,
      draws     = draws_obj
    )
  )
}
