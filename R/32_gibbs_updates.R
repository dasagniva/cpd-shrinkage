# ==============================================================================
# 32_gibbs_updates.R
# Conditional updates for the production Gibbs sampler (Phase P4).
# Implements:
#   1. update_mu_psi: Joint Gaussian draw via O(n) banded Cholesky (spam)
#   2. update_phi: Truncated Gaussian with MH correction for stationary eps_1
#   3. update_sigma2: Conjugate Inverse-Gamma update
#   4. update_hsplus_scales: Makalic-Schmidt (2016) inverse-gamma hierarchy
#   5. update_outlier_scales: Makalic-Schmidt (independent) or dsp::sampleDSP (dynamic)
# ==============================================================================

suppressPackageStartupMessages({
  library(spam)
})

find_project_root <- function() {
  cur <- getwd()
  while (cur != "/" && !file.exists(file.path(cur, "DESCRIPTION")) && !file.exists(file.path(cur, "renv.lock"))) {
    cur <- dirname(cur)
  }
  cur
}
proj_root <- find_project_root()
local_gsl <- file.path(proj_root, "inst/libs/usr/lib/x86_64-linux-gnu/libgsl.so.27")
local_cblas <- file.path(proj_root, "inst/libs/usr/lib/x86_64-linux-gnu/libgslcblas.so.0")
if (file.exists(local_gsl) && file.exists(local_cblas)) {
  tryCatch({
    dyn.load(local_cblas)
    dyn.load(local_gsl)
  }, error = function(e) NULL)
}

#' Initialize Banded Matrix Index Structure for spam
#'
#' @param n Length of the time series
#' @return A list with precomputed spam indices and empty template matrix
#' @export
init_banded_structure <- function(n) {
  m <- 2 * n
  rows <- c(1:m, 1:(m - 1), 2:m, 1:(m - 2), 3:m, 1:(m - 3), 4:m)
  cols <- c(1:m, 2:m, 1:(m - 1), 3:m, 1:(m - 2), 4:m, 1:(m - 3))
  ord <- order(rows, cols)
  rows <- rows[ord]
  cols <- cols[ord]

  Q_tmpl <- spam(list(i = rows, j = cols, values = rep(1.0, length(rows))), nrow = m, ncol = m)
  list(n = n, m = m, rows = rows, cols = cols, Q_tmpl = Q_tmpl)
}

#' Joint (mu, psi) Conditional Gaussian Draw
#'
#' Samples the interleaved state vector x = (mu_1, psi_1, mu_2, psi_2, ..., mu_n, psi_n)
#' from its exact conditional Gaussian distribution using an O(n) banded Cholesky.
#'
#' @param y Observed time series vector of length n
#' @param phi Autoregressive parameter in (-0.95, 0.95)
#' @param sigma Innovation standard deviation > 0
#' @param v_omega Vector of length n - 1 with prior variances of omega_t (scaled by sigma^2)
#' @param nu2 Vector of length n with prior variances of psi_t
#' @param banded_struct Optional precomputed banded structure from init_banded_structure
#'
#' @return A list with elements `mu` (length n) and `psi` (length n)
#' @export
update_mu_psi <- function(y, phi, sigma, v_omega, nu2, banded_struct = NULL) {
  n <- length(y)
  m <- 2 * n

  if (is.null(banded_struct) || banded_struct$n != n) {
    banded_struct <- init_banded_structure(n)
  }

  sigma2 <- sigma^2
  inv_sigma2 <- 1.0 / sigma2

  # Band diagonals: d0 (main), d1 (dist 1), d2 (dist 2), d3 (dist 3)
  d0 <- numeric(m)
  d1 <- numeric(m - 1)
  d2 <- numeric(m - 2)
  d3 <- numeric(m - 3)

  phi2 <- phi^2

  # 1. Likelihood contribution to Q
  for (t in 1:n) {
    idx_mu  <- 2 * t - 1
    idx_psi <- 2 * t

    c_diag <- if (t == 1) 1.0 else if (t < n) (1.0 + phi2) else 1.0
    c_diag_s <- c_diag * inv_sigma2

    d0[idx_mu]  <- d0[idx_mu]  + c_diag_s
    d0[idx_psi] <- d0[idx_psi] + c_diag_s

    # Superdiagonal 1: between mu_t and psi_t
    d1[idx_mu]  <- d1[idx_mu]  + c_diag_s

    if (t > 1) {
      m_phi_s <- -phi * inv_sigma2
      # Superdiagonal 1: between psi_{t-1} and mu_t
      d1[2 * t - 2] <- d1[2 * t - 2] + m_phi_s
      # Superdiagonal 2: between mu_{t-1} and mu_t
      d2[2 * t - 3] <- d2[2 * t - 3] + m_phi_s
      # Superdiagonal 2: between psi_{t-1} and psi_t
      d2[2 * t - 2] <- d2[2 * t - 2] + m_phi_s
      # Superdiagonal 3: between mu_{t-1} and psi_t
      d3[2 * t - 3] <- d3[2 * t - 3] + m_phi_s
    }
  }

  # 2. Prior contributions to Q
  # Base mean level prior: mu_1 ~ N(0, 100 * sigma^2)
  d0[1] <- d0[1] + 1.0 / (100.0 * sigma2)

  # Increments: mu_t - mu_{t-1} ~ N(0, v_omega[t-1])
  for (t in 2:n) {
    w <- 1.0 / v_omega[t - 1]
    idx_prev <- 2 * (t - 1) - 1
    idx_curr <- 2 * t - 1
    d0[idx_prev] <- d0[idx_prev] + w
    d0[idx_curr] <- d0[idx_curr] + w
    d2[idx_prev] <- d2[idx_prev] - w
  }

  # Outliers: psi_t ~ N(0, nu2[t])
  for (t in 1:n) {
    d0[2 * t] <- d0[2 * t] + 1.0 / nu2[t]
  }

  # 3. Assemble spam matrix Q
  rows <- c(1:m, 1:(m - 1), 2:m, 1:(m - 2), 3:m, 1:(m - 3), 4:m)
  cols <- c(1:m, 2:m, 1:(m - 1), 3:m, 1:(m - 2), 4:m, 1:(m - 3))
  vals <- c(d0, d1, d1, d2, d2, d3, d3)

  # Fast spam construction
  Q <- spam(list(i = rows, j = cols, values = vals), nrow = m, ncol = m)
  diag.spam(Q) <- diag.spam(Q) + 1e-8

  # 4. Compute RHS vector b = M^T * \tilde{y} / sigma^2
  u_y <- numeric(n)
  u_y[1] <- sqrt(1.0 - phi2) * y[1]
  u_y[2:n] <- y[2:n] - phi * y[1:(n - 1)]

  Pt_u <- numeric(n)
  Pt_u[1] <- sqrt(1.0 - phi2) * u_y[1] - phi * u_y[2]
  if (n > 2) {
    for (t in 2:(n - 1)) {
      Pt_u[t] <- u_y[t] - phi * u_y[t + 1]
    }
  }
  Pt_u[n] <- u_y[n]

  b <- numeric(m)
  for (t in 1:n) {
    val_t <- Pt_u[t] * inv_sigma2
    b[2 * t - 1] <- val_t
    b[2 * t]     <- val_t
  }

  # 5. O(n) Banded Cholesky draw: R^T R = Q
  R <- tryCatch(chol(Q), error = function(e) {
    diag.spam(Q) <- diag.spam(Q) + 1e-6
    chol(Q)
  })
  z <- rnorm(m)
  w <- forwardsolve(R, b)
  x <- backsolve(R, w + z)

  # 6. De-interleave state vector
  mu  <- x[seq(1, m, by = 2)]
  psi <- x[seq(2, m, by = 2)]

  list(mu = mu, psi = psi)
}

#' Update Autoregressive Persistence phi
#'
#' Conditional draw of phi using truncated Gaussian proposal with Metropolis-Hastings
#' correction for stationary distribution of eps_1 and Beta(a_phi, b_phi) prior.
#'
#' @param e Residual vector e = y - mu - psi of length n
#' @param phi_curr Current value of phi
#' @param sigma Innovation standard deviation > 0
#' @param a_phi Beta prior parameter 1 (default 10.0)
#' @param b_phi Beta prior parameter 2 (default 2.0)
#'
#' @return A list with `phi` (scalar numeric) and `accepted` (logical)
#' @export
update_phi <- function(e, phi_curr, sigma, a_phi = 10.0, b_phi = 2.0) {
  n <- length(e)
  stopifnot(n >= 2)

  e_lag <- e[1:(n - 1)]
  e_curr <- e[2:n]

  S_xx <- sum(e_lag^2)
  S_xy <- sum(e_curr * e_lag)

  # Truncated Gaussian proposal
  V_prop <- sigma^2 / (S_xx + 1e-12)
  m_prop <- S_xy / (S_xx + 1e-12)

  # Propose candidate phi_cand
  phi_cand <- rnorm(1, mean = m_prop, sd = sqrt(V_prop))

  # Truncate to (-0.95, 0.95)
  if (phi_cand <= -0.95 || phi_cand >= 0.95) {
    return(list(phi = phi_curr, accepted = FALSE))
  }

  # Log-acceptance probability:
  # Proposal cancels the sum_{t=2}^n (e_t - phi e_{t-1})^2 part.
  # Remaining terms:
  # 1. Stationary term: 0.5 * log(1 - phi^2) - (1 - phi^2) * e[1]^2 / (2 * sigma^2)
  # 2. Beta prior: beta_lpdf((phi + 1)/2 | a_phi, b_phi)
  log_stat_cand <- 0.5 * log(1.0 - phi_cand^2) - (1.0 - phi_cand^2) * e[1]^2 / (2.0 * sigma^2)
  log_stat_curr <- 0.5 * log(1.0 - phi_curr^2) - (1.0 - phi_curr^2) * e[1]^2 / (2.0 * sigma^2)

  log_prior_cand <- dbeta((phi_cand + 1.0) / 2.0, a_phi, b_phi, log = TRUE)
  log_prior_curr <- dbeta((phi_curr + 1.0) / 2.0, a_phi, b_phi, log = TRUE)

  log_alpha <- (log_stat_cand - log_stat_curr) + (log_prior_cand - log_prior_curr)

  if (log(runif(1)) < log_alpha) {
    list(phi = phi_cand, accepted = TRUE)
  } else {
    list(phi = phi_curr, accepted = FALSE)
  }
}

#' Update Innovation Variance sigma2
#'
#' Inverse-Gamma conditional draw combining prewhitened likelihood innovations
#' and scaled increments omega_t.
#'
#' @param u Prewhitened innovations vector of length n
#' @param omega Trend increments vector of length n - 1
#' @param v_omega_tilde Prior scale multipliers for omega (without sigma^2)
#' @param mu1 Base level scalar mu_1
#'
#' @return Numeric scalar sigma2 > 0
#' @export
update_sigma2 <- function(u, omega, v_omega_tilde, mu1 = 0.0) {
  n <- length(u)
  n_omega <- length(omega)

  # Sum of squared innovations
  ss_u <- sum(u^2)

  # Sum of standardized squared increments
  ss_omega <- sum(omega^2 / (v_omega_tilde + 1e-12))

  # Base level contribution: mu_1 ~ N(0, 100 * sigma^2)
  ss_mu1 <- (mu1^2) / 100.0

  # Jeffreys prior pi(sigma^2) proportional to 1 / sigma^2
  shape <- (n + n_omega + 1.0) / 2.0
  rate  <- (ss_u + ss_omega + ss_mu1) / 2.0

  1.0 / rgamma(1, shape = shape, rate = rate)
}

#' 1D Slice Sampler for Global Scale tau ~ C+(0, s_tau)
#'
#' Implements slice sampling on eta = log(tau) (Neal 2003, CPD2_Code_Workflow.md line 251)
#'
#' @param tau_curr Current value of tau > 0
#' @param S_sum Sum of squared standardized primitives
#' @param K Number of items (n - 1 for increments, n for outliers)
#' @param s_tau Prior scale parameter
#' @param w Initial step-out window width (default 1.0)
#'
#' @return Updated scalar value tau > 0
#' @export
slice_sample_tau <- function(tau_curr, S_sum, K, s_tau, w = 1.0) {
  log_post <- function(eta) {
    -K * eta - 0.5 * S_sum * exp(-2.0 * eta) - log1p(exp(2.0 * (eta - log(s_tau))))
  }

  eta0 <- log(max(1e-12, tau_curr))
  y_level <- log_post(eta0) - rexp(1)

  # Step out
  L <- eta0 - runif(1) * w
  R <- L + w
  while (log_post(L) > y_level) L <- L - w
  while (log_post(R) > y_level) R <- R + w

  # Shrinkage
  for (iter in 1:100) {
    eta_prop <- runif(1, L, R)
    if (log_post(eta_prop) >= y_level) {
      return(exp(eta_prop))
    }
    if (eta_prop < eta0) L <- eta_prop else R <- eta_prop
  }
  exp(eta_prop)
}

#' Update Horseshoe+ Scales via Makalic-Schmidt (2016) Scheme with Slice Sampling for tau
#'
#' Updates local scales lambda_omega, secondary scales eta_omega, global scale tau_omega,
#' and their auxiliary parameters using exact conjugate Inverse-Gamma draws and slice sampling.
#'
#' @param omega Trend increments of length n - 1
#' @param sigma Innovation standard deviation
#' @param tau_omega Current global scale scalar
#' @param lambda_omega Current local scales vector of length n - 1
#' @param eta_omega Current secondary scales vector of length n - 1
#' @param kappa Prior inflation multiplier vector of length n - 1
#' @param aux_lambda Auxiliary vector for lambda of length n - 1
#' @param aux_eta Auxiliary vector for eta of length n - 1
#' @param aux_tau Auxiliary scalar for tau
#' @param s_tau Scale parameter for global scale Cauchy prior (default 1 / sqrt(n))
#' @param slab_scale Slab scale (default 5.0)
#' @param slab_df Slab degrees of freedom (default 4.0)
#'
#' @return A list with updated `tau_omega`, `lambda_omega`, `eta_omega`, `v_omega`, `v_omega_tilde`,
#'   and updated auxiliary parameters `aux_lambda`, `aux_eta`, `aux_tau`, `c2`
#' @export
update_hsplus_scales <- function(
    omega,
    sigma,
    tau_omega,
    lambda_omega,
    eta_omega,
    kappa,
    aux_lambda,
    aux_eta,
    aux_tau,
    s_tau = 0.1,
    c2_tilde = 1.0,
    slab_scale = 5.0,
    slab_df = 4.0
) {
  n_omega <- length(omega)
  sigma2 <- sigma^2
  tau2 <- tau_omega^2

  # 1. Update regularised slab parameter c2
  # c2_tilde ~ Inv-Gamma(slab_df/2, slab_df/2)
  c2_tilde <- 1.0 / rgamma(1, shape = slab_df / 2.0, rate = slab_df / 2.0)
  c2 <- (slab_scale^2) * c2_tilde

  # 2. Update local scales lambda_omega
  # theta_t = omega_t / (sigma * tau_omega * kappa_t * eta_omega_t) ~ N(0, lambda_t^2)
  eff_denom_lam <- sigma * tau_omega * kappa * eta_omega + 1e-12
  theta_lam <- omega / eff_denom_lam

  rate_lam <- (1.0 / aux_lambda) + 0.5 * (theta_lam^2)
  lambda2 <- 1.0 / rgamma(n_omega, shape = 1.0, rate = rate_lam)
  lambda_omega <- sqrt(pmax(1e-12, lambda2))

  # Auxiliary aux_lambda
  aux_lambda <- 1.0 / rgamma(n_omega, shape = 1.0, rate = 1.0 + 1.0 / lambda2)

  # 3. Update secondary scales eta_omega (Horseshoe+)
  # theta_eta_t = omega_t / (sigma * tau_omega * kappa_t * lambda_t) ~ N(0, eta_t^2)
  eff_denom_eta <- sigma * tau_omega * kappa * lambda_omega + 1e-12
  theta_eta <- omega / eff_denom_eta

  rate_eta <- (1.0 / aux_eta) + 0.5 * (theta_eta^2)
  eta2 <- 1.0 / rgamma(n_omega, shape = 1.0, rate = rate_eta)
  eta_omega <- sqrt(pmax(1e-12, eta2))

  # Auxiliary aux_eta
  aux_eta <- 1.0 / rgamma(n_omega, shape = 1.0, rate = 1.0 + 1.0 / eta2)

  # 4. Update global scale tau_omega via slice sampler (Section 5, Phase P4)
  eff_lam_raw <- pmin(kappa * eta_omega * lambda_omega, 1e4)
  zeta_tau <- omega / (sigma * eff_lam_raw + 1e-12)
  tau_omega <- slice_sample_tau(tau_omega, S_sum = sum(zeta_tau^2), K = n_omega, s_tau = s_tau)

  # 5. Compute regularised variances v_omega and v_omega_tilde
  tau_lam <- pmin(tau_omega * eff_lam_raw, 1e4)
  eff_lam_sq <- (c2 * (eff_lam_raw^2)) / (c2 + (tau_lam^2) + 1e-12)

  # v_omega_tilde: scale without sigma^2
  v_omega_tilde <- pmax(1e-12, pmin(1e4, (tau_omega^2) * eff_lam_sq))
  # v_omega: actual variance including sigma^2
  v_omega <- pmax(1e-12, pmin(1e4, sigma2 * v_omega_tilde))

  list(
    tau_omega     = tau_omega,
    lambda_omega  = lambda_omega,
    eta_omega     = eta_omega,
    v_omega       = v_omega,
    v_omega_tilde = v_omega_tilde,
    aux_lambda    = aux_lambda,
    aux_eta       = aux_eta,
    aux_tau       = aux_tau,
    c2_tilde      = c2_tilde,
    c2            = c2
  )
}

#' Update Outlier Scales (Independent or Dynamic)
#'
#' @param psi Outlier vector of length n
#' @param sigma Innovation standard deviation
#' @param tau_psi Global outlier scale scalar (for independent model)
#' @param lambda_psi Local outlier scales of length n
#' @param eta_psi Secondary outlier scales of length n
#' @param aux_lambda_psi Auxiliary parameters
#' @param aux_eta_psi Auxiliary parameters
#' @param aux_tau_psi Auxiliary parameter
#' @param dynamic Logical: if TRUE, use dsp::sampleDSP; if FALSE, use independent Makalic-Schmidt
#' @param rho Persistence parameter for dynamic outlier log-scales
#' @param m_h Base mean for dynamic outlier log-scales
#' @param delta CAPA offset parameter
#' @param a_vec Binary indicator vector of CAPA point anomalies
#' @param dsp_state Internal state list for dsp::sampleDSP
#'
#' @return A list with updated outlier variances `nu2`, `tau_psi`, and scales
#' @export
update_outlier_scales <- function(
    psi,
    sigma,
    tau_psi = 0.1,
    lambda_psi = NULL,
    eta_psi = NULL,
    aux_lambda_psi = NULL,
    aux_eta_psi = NULL,
    aux_tau_psi = NULL,
    dynamic = FALSE,
    rho = 0.0,
    m_h = 0.0,
    delta = 0.0,
    a_vec = NULL,
    dsp_state = NULL,
    c2 = 25.0
) {
  n <- length(psi)
  sigma2 <- sigma^2

  if (!dynamic) {
    # Independent Horseshoe+ outlier scales (Section 3.3 with rho = 0, delta = 0)
    if (is.null(lambda_psi)) lambda_psi <- rep(1.0, n)
    if (is.null(eta_psi)) eta_psi <- rep(1.0, n)
    if (is.null(aux_lambda_psi)) aux_lambda_psi <- rep(1.0, n)
    if (is.null(aux_eta_psi)) aux_eta_psi <- rep(1.0, n)
    if (is.null(aux_tau_psi)) aux_tau_psi <- 1.0

    # 1. Update lambda_psi
    denom_lam <- tau_psi * eta_psi + 1e-12
    theta_lam <- psi / denom_lam
    rate_lam <- (1.0 / aux_lambda_psi) + 0.5 * (theta_lam^2)
    lam2 <- 1.0 / rgamma(n, shape = 1.0, rate = rate_lam)
    lambda_psi <- sqrt(pmax(1e-12, lam2))
    aux_lambda_psi <- 1.0 / rgamma(n, shape = 1.0, rate = 1.0 + 1.0 / lam2)

    # 2. Update eta_psi
    denom_eta <- tau_psi * lambda_psi + 1e-12
    theta_eta <- psi / denom_eta
    rate_eta <- (1.0 / aux_eta_psi) + 0.5 * (theta_eta^2)
    eta2 <- 1.0 / rgamma(n, shape = 1.0, rate = rate_eta)
    eta_psi <- sqrt(pmax(1e-12, eta2))
    aux_eta_psi <- 1.0 / rgamma(n, shape = 1.0, rate = 1.0 + 1.0 / eta2)

    # 3. Update tau_psi via slice sampler
    eff_lam_raw <- pmin(eta_psi * lambda_psi, 1e4)
    zeta_tau <- psi / (eff_lam_raw + 1e-12)
    s_tau_p <- 1.0 / sqrt(n)
    tau_psi <- slice_sample_tau(tau_psi, S_sum = sum(zeta_tau^2), K = n, s_tau = s_tau_p)

    # Regularised outlier variance
    tau_lam_p <- pmin(tau_psi * eff_lam_raw, 1e4)
    eff_lam_sq <- (c2 * (eff_lam_raw^2)) / (c2 + (tau_lam_p^2) + 1e-12)
    nu2 <- pmax(1e-12, (tau_psi^2) * eff_lam_sq)

    list(
      nu2            = nu2,
      tau_psi        = tau_psi,
      lambda_psi     = lambda_psi,
      eta_psi        = eta_psi,
      aux_lambda_psi = aux_lambda_psi,
      aux_eta_psi    = aux_eta_psi,
      aux_tau_psi    = aux_tau_psi
    )
  } else {
    # Dynamic outlier log-scales via dsp::sampleDSP (Section 3.3)
    suppressPackageStartupMessages(library(dsp))
    ns_dsp <- asNamespace("dsp")

    if (is.null(dsp_state)) {
      dsp_state <- list(
        ht          = matrix(log(pmax(1e-6, psi^2)), ncol = 1),
        dhs_mean    = m_h,
        dhs_phi     = rho,
        sigma_eta_t = matrix(1.0, nrow = n - 1, ncol = 1),
        sigma_eta_0 = 1.0,
        dhs_mean0   = m_h
      )
    }

    # Subtract mean offset delta * a_t
    mean_offset <- if (!is.null(a_vec)) delta * a_vec else rep(0.0, n)
    ht_in <- dsp_state$ht - matrix(mean_offset, ncol = 1)

    evol_params <- list(
      ht          = ht_in,
      dhs_mean    = dsp_state$dhs_mean,
      dhs_phi     = dsp_state$dhs_phi,
      sigma_eta_t = dsp_state$sigma_eta_t,
      sigma_eta_0 = dsp_state$sigma_eta_0,
      dhs_mean0   = dsp_state$dhs_mean0
    )

    dsp_out <- ns_dsp$sampleDSP(
      omega         = matrix(psi, ncol = 1),
      evolParams    = evol_params,
      loc           = NULL,
      alphaPlusBeta = 1.0
    )

    # Add back mean offset
    ht_out <- dsp_out$ht + matrix(mean_offset, ncol = 1)
    nu2 <- as.vector(exp(ht_out))
    nu2 <- pmax(1e-12, nu2)

    list(
      nu2       = nu2,
      dsp_state = list(
        ht          = ht_out,
        dhs_mean    = dsp_out$dhs_mean,
        dhs_phi     = dsp_out$dhs_phi,
        sigma_eta_t = dsp_out$sigma_eta_t,
        sigma_eta_0 = dsp_out$sigma_eta_0,
        dhs_mean0   = dsp_out$dhs_mean0
      ),
      rho = as.numeric(dsp_out$dhs_phi),
      m_h = as.numeric(dsp_out$dhs_mean)
    )
  }
}
