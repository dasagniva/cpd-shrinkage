# ==============================================================================
# 30_stan_reference.R
# Reference Stan model wrapper for Phase P3 validation.
# Wraps stan/hsplus_ar1_ref.stan into the common cpd_fit schema.
# ==============================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

#' Fit Stan Reference Model
#'
#' Fits the reference model with non-centred Horseshoe+ increments and independent
#' outlier scales (rho = 0, delta = 0) under prewhitened AR(1) likelihood.
#' Intended for n <= 300 validation.
#'
#' @param y Numeric vector of time series observations
#' @param kappa Numeric vector of local shrinkage multipliers (default rep(1.0, length(y)))
#' @param alpha_fdr False discovery rate threshold for candidate selection (default 0.05)
#' @param chains Number of MCMC chains (default 2)
#' @param parallel_chains Number of parallel chains (default 2)
#' @param iter_warmup Number of warmup iterations per chain (default 500)
#' @param iter_sampling Number of sampling iterations per chain (default 500)
#' @param adapt_delta NUTS target acceptance rate (default 0.99)
#' @param max_treedepth Maximum NUTS treedepth (default 12)
#' @param init Initialization radius or value (default 0.5)
#' @param seed Random seed for Stan sampling
#' @param stan_file Path to Stan reference model file (default "stan/hsplus_ar1_ref.stan")
#'
#' @return A validated `cpd_fit` object
#' @export
fit_stan_reference <- function(
    y,
    kappa = NULL,
    alpha_fdr = 0.05,
    chains = 2,
    parallel_chains = 2,
    iter_warmup = 500,
    iter_sampling = 500,
    adapt_delta = 0.99,
    max_treedepth = 12,
    init = 0.5,
    seed = NULL,
    stan_file = "stan/hsplus_ar1_ref.stan"
) {
  start_time <- proc.time()[3]
  n <- length(y)
  stopifnot(n >= 10, is.numeric(y))

  if (is.null(kappa)) {
    kappa <- rep(1.0, n)
  }
  stopifnot(length(kappa) == n)

  if (!file.exists(stan_file)) {
    stop("Stan reference model file not found: ", stan_file)
  }

  mod <- cmdstan_model(stan_file)

  stan_data <- list(
    n          = n,
    y          = as.vector(y),
    kappa      = as.vector(kappa),
    s_tau      = 1.0 / sqrt(n),
    a_phi      = 10.0,
    b_phi      = 2.0,
    slab_scale = 5.0,
    slab_df    = 4.0,
    sigma_min  = 1e-4,
    sigma_max  = 1e4
  )

  fit <- mod$sample(
    data            = stan_data,
    seed            = if (!is.null(seed)) as.integer(seed) else sample.int(1e7, 1),
    chains          = chains,
    parallel_chains = parallel_chains,
    iter_warmup     = iter_warmup,
    iter_sampling   = iter_sampling,
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    init            = init,
    show_messages   = FALSE,
    refresh         = 0
  )

  runtime_sec <- proc.time()[3] - start_time

  # Extract posterior draws
  draws <- fit$draws()
  summ <- fit$summary(c("phi", "sigma", "tau_omega", "tau_psi"))

  # Posterior means
  phi_hat <- as.numeric(summ$mean[summ$variable == "phi"])
  sigma_hat <- as.numeric(summ$mean[summ$variable == "sigma"])
  sigma_marg_hat <- sigma_hat / sqrt(1.0 - phi_hat^2)

  # Outlier scores o_t from generated quantities
  o_t_draws <- fit$draws("o_t")
  o_t_mat <- as_draws_matrix(o_t_draws)
  o_t_mean <- colMeans(o_t_mat)

  # Increment scores p_omega (Pr(|omega_t| > 1.5 * sigma_marg))
  omega_draws <- fit$draws("omega")
  omega_mat <- as_draws_matrix(omega_draws)
  d_omega <- 1.5 * sigma_marg_hat
  p_omega <- colMeans(abs(omega_mat) > d_omega)

  # Outlier draws p_psi (Pr(|psi_t| > 1.5 * sigma_marg))
  psi_draws <- fit$draws("psi")
  psi_mat <- as_draws_matrix(psi_draws)
  d_psi <- 1.5 * sigma_marg_hat
  p_psi <- colMeans(abs(psi_mat) > d_psi)

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
  # omega_t corresponds to change between t and t+1 (1-indexed index t+1 in 1..n)
  cp_locs <- if (length(flagged_cp) > 0) flagged_cp + 1L else integer(0)

  flagged_point <- fdr_select(p_psi, alpha_fdr)

  # Four-class labelling rule (Section 3.6)
  class_vec <- rep("nominal", n)
  if (length(cp_locs) > 0) {
    class_vec[cp_locs] <- "shift"
  }
  if (length(flagged_point) > 0) {
    # Outliers that are not sustained shifts
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

  mu_draws <- fit$draws("mu")
  mu_mat <- as_draws_matrix(mu_draws)
  mu_mean <- colMeans(mu_mat)

  new_cpd_fit(
    method      = "stan_ref",
    cp          = cp_locs,
    cp_int      = NULL,
    point       = flagged_point,
    collective  = NULL,
    class       = factor(class_vec, levels = CPD_CLASS_LEVELS),
    scores      = scores_df,
    runtime_sec = as.numeric(runtime_sec),
    diag        = list(
      summary   = summ,
      max_rhat  = max(summ$rhat, na.rm = TRUE),
      min_ess   = min(summ$ess_bulk, na.rm = TRUE),
      phi_hat   = phi_hat,
      sigma_hat = sigma_hat,
      mu_mean   = mu_mean
    )
  )
}
