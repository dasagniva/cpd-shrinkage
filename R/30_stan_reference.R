# ==============================================================================
# 30_stan_reference.R
# Reference Stan model wrapper for Phase P3 validation.
# Wraps stan/hsplus_ar1_ref.stan into the common cpd_fit schema.
# ==============================================================================

if (file.exists("R/00_utils.R")) {
  source("R/00_utils.R")
  source("R/32_gibbs_updates.R")
  source("R/40_scoring.R")
  source("R/41_bfdr.R")
  source("R/42_classify.R")
} else if (file.exists("../../R/00_utils.R")) {
  source("../../R/00_utils.R")
  source("../../R/32_gibbs_updates.R")
  source("../../R/40_scoring.R")
  source("../../R/41_bfdr.R")
  source("../../R/42_classify.R")
}

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

#' Fit Stan Reference Model
#'
#' Fits the reference model with non-centred horseshoe+ increments and independent
#' outlier scales (rho = 0, delta = 0) under the prewhitened AR(1) likelihood.
#' Like the Gibbs sampler, the model is fitted to robust_standardize(y) and levels
#' are back-transformed. Intended for n <= 300 validation.
#'
#' @param y Numeric vector of time series observations
#' @param kappa Numeric vector of local shrinkage multipliers (default rep(1.0, length(y)))
#' @param alpha_fdr Bayesian FDR level for candidate selection (default 0.05)
#' @param capa_collective,L_max,r Collective-anomaly rule inputs (Section 3.6)
#' @param d_omega,d_psi Threshold multipliers of the marginal noise SD (Eq. 9)
#' @param use_slab Logical; TRUE uses the regularised slab (not the Gibbs model)
#' @param chains,parallel_chains,iter_warmup,iter_sampling,adapt_delta,max_treedepth,init NUTS settings
#' @param seed Random seed for Stan sampling
#' @param stan_file Path to Stan reference model file
#'
#' @return A validated `cpd_fit` object; `diag` holds posterior summaries and
#'   90% interval endpoints for mu_t and o_t
#' @export
fit_stan_reference <- function(
    y,
    kappa = NULL,
    alpha_fdr = 0.05,
    capa_collective = NULL,
    L_max = NULL,
    r = 0.3,
    d_omega = 1.5,
    d_psi = 1.5,
    use_slab = FALSE,
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

  if (is.null(kappa)) kappa <- rep(1.0, n)
  stopifnot(length(kappa) == n)
  if (!file.exists(stan_file)) stop("Stan reference model file not found: ", stan_file)

  std <- robust_standardize(y)
  mod <- cmdstan_model(stan_file)

  stan_data <- list(
    n          = n,
    y          = as.vector(std$z),
    kappa      = as.vector(kappa),
    s_tau      = 1.0 / sqrt(n),
    a_phi      = 10.0,
    b_phi      = 2.0,
    mu1_sd     = sqrt(MU1_PRIOR_VAR),
    use_slab   = as.integer(isTRUE(use_slab)),
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

  # Global draws on the same scale as the Gibbs sampler (sigma, tau_psi in data units)
  gdraws <- fit$draws(c("phi", "sigma", "tau_omega", "tau_psi"))
  gdraws[, , "sigma"] <- gdraws[, , "sigma"] * std$scale
  gdraws[, , "tau_psi"] <- gdraws[, , "tau_psi"] * std$scale
  summ <- summarize_draws(gdraws)
  draws_mat <- function(v) as_draws_matrix(fit$draws(v))
  phi_d   <- as.vector(draws_mat("phi"))
  sigma_d <- as.vector(draws_mat("sigma"))
  sigma_marg_d <- sigma_d / sqrt(1.0 - phi_d^2)

  omega_mat <- draws_mat("omega")
  psi_mat   <- draws_mat("psi")
  mu_mat    <- draws_mat("mu")
  o_t_mat   <- draws_mat("o_t")

  # Eq. 9 with per-draw thresholds, as in the Gibbs sampler
  p_omega <- colMeans(abs(omega_mat) > d_omega * sigma_marg_d)
  p_psi   <- colMeans(abs(psi_mat) > d_psi * sigma_marg_d)

  q90 <- function(m) apply(m, 2, stats::quantile, probs = c(0.05, 0.95), names = FALSE)
  mu_q <- q90(mu_mat) * std$scale + std$center
  o_q  <- q90(o_t_mat)

  sampler_diag <- fit$diagnostic_summary(quiet = TRUE)

  diag_list <- list(
    summary       = summ,
    max_rhat      = max(summ$rhat, na.rm = TRUE),
    min_ess       = min(summ$ess_bulk, na.rm = TRUE),
    phi_hat       = mean(phi_d),
    sigma_hat     = mean(sigma_d) * std$scale,
    draws         = gdraws,
    tau_omega_hat = mean(as.vector(draws_mat("tau_omega"))),
    mu_mean       = colMeans(mu_mat) * std$scale + std$center,
    mu_sd         = apply(mu_mat, 2, stats::sd) * std$scale,
    mu_q05        = mu_q[1, ],
    mu_q95        = mu_q[2, ],
    o_t_q05       = o_q[1, ],
    o_t_q95       = o_q[2, ],
    num_divergent = sum(sampler_diag$num_divergent),
    standardization = std[c("center", "scale")]
  )

  classify_and_build_fit(
    y               = y,
    mu_hat          = diag_list$mu_mean,
    s_out           = colMeans(o_t_mat),
    p_omega         = p_omega,
    p_psi           = p_psi,
    alpha_fdr       = alpha_fdr,
    capa_collective = capa_collective,
    method          = "stan_ref",
    runtime_sec     = as.numeric(runtime_sec),
    diag            = diag_list,
    L_max           = L_max,
    r               = r
  )
}
