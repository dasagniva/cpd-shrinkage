# ==============================================================================
# 32_gibbs_updates.R
# Conditional updates for the production Gibbs sampler (Phase P4).
# All updates assume the series has been standardised with robust_standardize()
# (00_utils.R), so every prior scale below is in units of the robust noise scale.
# Implements:
#   1. update_mu_psi: Joint Gaussian draw via O(n) banded Cholesky (spam)
#   2. update_phi: Exact truncated-Gaussian proposal + MH correction for eps_1
#   3. update_sigma2: Conjugate Inverse-Gamma update
#   4. update_hsplus_scales: Makalic-Schmidt (2016) horseshoe+ hierarchy (Eq. 3-5),
#      with a partially collapsed (lambda-marginal) update of the global scale
#   5. update_outlier_scales: Makalic-Schmidt (independent) or dsp::sampleDSP (Eq. 6-7)
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

#' Prior variance of mu_1 on the standardised scale (sd = 100 robust noise units)
MU1_PRIOR_VAR <- 1e4

#' Initialize Banded Matrix Index Structure for spam
#'
#' Stores the fixed sparsity pattern of the interleaved precision Q (bandwidth 3)
#' and a cache for the symbolic Cholesky factorisation, so that every iteration
#' only refills numeric entries and reuses the symbolic structure.
#'
#' @param n Length of the time series
#' @return A list with precomputed spam indices, template matrix and chol cache
#' @export
init_banded_structure <- function(n) {
  m <- 2 * n
  rows <- c(1:m, 1:(m - 1), 2:m, 1:(m - 2), 3:m, 1:(m - 3), 4:m)
  cols <- c(1:m, 2:m, 1:(m - 1), 3:m, 1:(m - 2), 4:m, 1:(m - 3))
  ord <- order(rows, cols)

  Q_tmpl <- spam(list(i = rows[ord], j = cols[ord], values = rep(1.0, length(rows))), nrow = m, ncol = m)
  cache <- new.env(parent = emptyenv())
  cache$chol <- NULL
  list(n = n, m = m, rows = rows[ord], cols = cols[ord], ord = ord, Q_tmpl = Q_tmpl, cache = cache)
}

#' Assemble the interleaved precision matrix Q and right-hand side b
#'
#' @inheritParams update_mu_psi
#' @return A list with spam matrix `Q` and numeric vector `b`
#' @keywords internal
build_mu_psi_system <- function(y, phi, sigma, v_omega, nu2, banded_struct) {
  n <- length(y)
  m <- 2 * n
  inv_s2 <- 1.0 / sigma^2
  idx_mu <- seq(1, m, by = 2)
  idx_psi <- seq(2, m, by = 2)
  w <- 1.0 / v_omega

  # Likelihood: A^T P^T P A / sigma^2, with (P^T P) tridiagonal in time
  c_diag <- c(1.0, rep(1.0 + phi^2, max(0, n - 2)), 1.0)[seq_len(n)]
  off <- -phi * inv_s2

  d0 <- numeric(m)
  d0[idx_mu] <- c_diag * inv_s2
  d0[idx_psi] <- c_diag * inv_s2 + 1.0 / nu2
  d0[1] <- d0[1] + 1.0 / MU1_PRIOR_VAR
  d0[idx_mu[-n]] <- d0[idx_mu[-n]] + w # mu_{t-1} side of omega_t
  d0[idx_mu[-1]] <- d0[idx_mu[-1]] + w # mu_t side of omega_t

  d1 <- numeric(m - 1)
  d1[seq(1, m - 1, by = 2)] <- c_diag * inv_s2 # (mu_t, psi_t)
  if (n > 1) d1[seq(2, m - 2, by = 2)] <- off # (psi_{t-1}, mu_t)

  d2 <- numeric(m - 2)
  if (n > 1) {
    d2[seq(1, m - 3, by = 2)] <- off - w # (mu_{t-1}, mu_t)
    d2[seq(2, m - 2, by = 2)] <- off # (psi_{t-1}, psi_t)
  }

  d3 <- numeric(m - 3)
  if (n > 1) d3[seq(1, m - 3, by = 2)] <- off # (mu_{t-1}, psi_t); even slots stay 0

  vals <- c(d0, d1, d1, d2, d2, d3, d3)
  Q <- banded_struct$Q_tmpl
  Q@entries <- vals[banded_struct$ord] + ifelse(banded_struct$rows == banded_struct$cols, 1e-8, 0)

  # RHS b = A^T P^T P y / sigma^2
  u <- c(sqrt(1.0 - phi^2) * y[1], y[-1] - phi * y[-n])
  Ptu <- u
  Ptu[1] <- sqrt(1.0 - phi^2) * u[1]
  if (n > 1) Ptu[-n] <- Ptu[-n] - phi * u[-1]
  b <- rep(Ptu * inv_s2, each = 2)

  list(Q = Q, b = b)
}

#' Factorise the (mu, psi) conditional and evaluate the scale-dependent marginal likelihood
#'
#' Computes the banded Cholesky factor R of Q (reusing the cached symbolic
#' factorisation) and, with (mu, psi) integrated out,
#'   log p(y | v_omega, nu2, phi, sigma) = -1/2 sum log v_omega - 1/2 sum log nu2
#'                                         - 1/2 log|Q| + 1/2 b' Q^{-1} b + const,
#' where the constant does not depend on v_omega or nu2.
#'
#' @inheritParams update_mu_psi
#' @return A list with factor `R`, `w = R^{-T} b` and `logml`
#' @export
factor_mu_psi <- function(y, phi, sigma, v_omega, nu2, banded_struct) {
  sys <- build_mu_psi_system(y, phi, sigma, v_omega, nu2, banded_struct)

  cache <- banded_struct$cache
  R <- NULL
  if (!is.null(cache$chol)) {
    R <- tryCatch(update(cache$chol, sys$Q), error = function(e) NULL)
  }
  if (is.null(R)) {
    R <- tryCatch(chol(sys$Q), error = function(e) {
      Q2 <- sys$Q
      diag.spam(Q2) <- diag.spam(Q2) + 1e-6
      chol(Q2)
    })
  }
  cache$chol <- R

  w <- forwardsolve(R, sys$b)
  logml <- -0.5 * sum(log(v_omega)) - 0.5 * sum(log(nu2)) - sum(log(diag(R))) + 0.5 * sum(w^2)
  list(R = R, w = w, logml = logml)
}

#' Joint (mu, psi) Conditional Gaussian Draw
#'
#' Samples the interleaved state vector x = (mu_1, psi_1, mu_2, psi_2, ..., mu_n, psi_n)
#' from its exact conditional Gaussian distribution using an O(n) banded Cholesky.
#'
#' @param y Observed (standardised) time series vector of length n
#' @param phi Autoregressive parameter in [-0.95, 0.95]
#' @param sigma Innovation standard deviation > 0
#' @param v_omega Vector of length n - 1 with prior variances of omega_t (including sigma^2)
#' @param nu2 Vector of length n with prior variances of psi_t
#' @param banded_struct Optional precomputed banded structure from init_banded_structure
#' @param fac Optional factorisation from factor_mu_psi for the same inputs
#'
#' @return A list with elements `mu` (length n) and `psi` (length n)
#' @export
update_mu_psi <- function(y, phi, sigma, v_omega, nu2, banded_struct = NULL, fac = NULL) {
  n <- length(y)
  m <- 2 * n
  if (is.null(banded_struct) || banded_struct$n != n) {
    banded_struct <- init_banded_structure(n)
  }
  if (is.null(fac)) fac <- factor_mu_psi(y, phi, sigma, v_omega, nu2, banded_struct)

  # x = Q^{-1} b + R^{-1} z, with pivoting handled by spam's solvers
  x <- backsolve(fac$R, fac$w + rnorm(m))
  list(mu = x[seq(1, m, by = 2)], psi = x[seq(2, m, by = 2)])
}

#' Draw from a univariate truncated normal N(mean, sd^2) restricted to [lower, upper]
#'
#' Exact inverse-CDF sampler that works in log space, so it stays accurate when
#' [lower, upper] lies far in the tail of the untruncated normal.
#'
#' @param mean,sd Parameters of the untruncated normal
#' @param lower,upper Truncation bounds (lower < upper)
#' @return A single draw in [lower, upper]
#' @export
rtnorm1 <- function(mean, sd, lower, upper) {
  stopifnot(sd > 0, lower < upper)
  a <- (lower - mean) / sd
  b <- (upper - mean) / sd

  lower_tail_draw <- function(a, b) {
    # Region [a, b] with b <= 0: sample uniformly between Phi(a) and Phi(b) in log space
    la <- pnorm(a, log.p = TRUE)
    lb <- pnorm(b, log.p = TRUE)
    lp <- lb + log1p(-runif(1) * (-expm1(la - lb)))
    qnorm(lp, log.p = TRUE)
  }

  z <- if (b <= 0) {
    lower_tail_draw(a, b)
  } else if (a >= 0) {
    -lower_tail_draw(-b, -a)
  } else {
    qnorm(runif(1, pnorm(a), pnorm(b)))
  }
  min(upper, max(lower, mean + sd * z))
}

#' Update Autoregressive Persistence phi
#'
#' Independence Metropolis-Hastings step. The proposal is the Gaussian part of the
#' conditional, N(S_xy / S_xx, sigma^2 / S_xx), truncated exactly to |phi| <= phi_max.
#' The acceptance ratio therefore contains only the stationary eps_1 term and the
#' Beta prior on (phi + 1) / 2 (Section 3.4).
#'
#' @param e Residual vector e = y - mu - psi of length n
#' @param phi_curr Current value of phi
#' @param sigma Innovation standard deviation > 0
#' @param a_phi Beta prior parameter 1 (default 10.0)
#' @param b_phi Beta prior parameter 2 (default 2.0)
#' @param phi_max Truncation bound (default 0.95)
#'
#' @return A list with `phi` (scalar numeric) and `accepted` (logical)
#' @export
update_phi <- function(e, phi_curr, sigma, a_phi = 10.0, b_phi = 2.0, phi_max = 0.95) {
  n <- length(e)
  stopifnot(n >= 2)

  e_lag <- e[-n]
  e_curr <- e[-1]
  S_xx <- sum(e_lag^2) + 1e-12
  S_xy <- sum(e_curr * e_lag)

  phi_cand <- rtnorm1(S_xy / S_xx, sigma / sqrt(S_xx), -phi_max, phi_max)

  log_rest <- function(phi) {
    0.5 * log1p(-phi^2) - (1.0 - phi^2) * e[1]^2 / (2.0 * sigma^2) +
      dbeta((phi + 1.0) / 2.0, a_phi, b_phi, log = TRUE)
  }

  if (log(runif(1)) < log_rest(phi_cand) - log_rest(phi_curr)) {
    list(phi = phi_cand, accepted = TRUE)
  } else {
    list(phi = phi_curr, accepted = FALSE)
  }
}

#' Update Innovation Variance sigma2
#'
#' Inverse-Gamma conditional combining the prewhitened likelihood (n terms) and the
#' sigma-scaled increments omega_t (n - 1 terms), under the Jeffreys prior
#' pi(sigma^2) proportional to 1 / sigma^2. The mu_1 prior does not involve sigma.
#'
#' @param u Prewhitened innovations vector of length n
#' @param omega Trend increments vector of length n - 1
#' @param v_omega_tilde Prior variances of omega without the sigma^2 factor
#'
#' @return Numeric scalar sigma2 > 0
#' @export
update_sigma2 <- function(u, omega, v_omega_tilde) {
  shape <- (length(u) + length(omega)) / 2.0
  rate <- (sum(u^2) + sum(omega^2 / v_omega_tilde)) / 2.0
  1.0 / rgamma(1, shape = shape, rate = rate)
}

#' Generic univariate slice sampler (stepping out + shrinkage; Neal 2003)
#'
#' @param x0 Current value
#' @param log_f Log density (up to a constant)
#' @param w Initial window width
#' @return A new draw leaving log_f invariant
#' @keywords internal
slice_sample_1d <- function(x0, log_f, w = 1.0) {
  y_level <- log_f(x0) - rexp(1)
  L <- x0 - runif(1) * w
  R <- L + w
  while (log_f(L) > y_level) L <- L - w
  while (log_f(R) > y_level) R <- R + w
  repeat {
    x1 <- runif(1, L, R)
    if (log_f(x1) >= y_level) return(x1)
    if (x1 < x0) L <- x1 else R <- x1
  }
}

#' Partially collapsed update of a global horseshoe scale
#'
#' With the Makalic-Schmidt representation lambda_t^2 | nu_t ~ IG(1/2, 1/nu_t),
#' integrating lambda_t out gives x_t | tau, nu_t ~ Cauchy(0, tau a_t sqrt(2 / nu_t)),
#' where a_t is the remaining (non-lambda) scale. tau is drawn from
#'   p(tau | x, nu) propto C+(tau; 0, s_tau) prod_t Cauchy(x_t; 0, tau a_t sqrt(2 / nu_t))
#' by slice sampling on log tau; lambda must then be drawn from its full conditional
#' given the new tau (partially collapsed Gibbs, van Dyk & Park 2008).
#'
#' @param tau_curr Current global scale
#' @param x Coefficients (omega_t / sigma or psi_t)
#' @param a Remaining scale multipliers (kappa_t * eta_t, or eta_t)
#' @param nu Makalic-Schmidt auxiliaries of lambda
#' @param s_tau Half-Cauchy prior scale
#' @return Updated tau
#' @export
collapsed_tau_update <- function(tau_curr, x, a, nu, s_tau) {
  r2 <- (x / (a * sqrt(2.0 / nu)))^2
  K <- length(x)
  log_f <- function(e) {
    -(K - 1) * e - sum(log1p(r2 * exp(-2.0 * e))) - log1p(exp(2.0 * (e - log(s_tau))))
  }
  exp(slice_sample_1d(log(max(1e-300, tau_curr)), log_f))
}

#' 1D Slice Sampler for a Global Scale tau ~ C+(0, s_tau)
#'
#' Samples eta = log(tau) (Neal 2003) from the conditional
#'   p(tau | .) propto tau^{-K} exp(-S / (2 tau^2)) / (1 + (tau / s_tau)^2),
#' including the Jacobian d tau / d eta = tau of the log transform.
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
    -(K - 1) * eta - 0.5 * S_sum * exp(-2.0 * eta) - log1p(exp(2.0 * (eta - log(s_tau))))
  }

  eta0 <- log(max(1e-300, tau_curr))
  y_level <- log_post(eta0) - rexp(1)

  # Step out
  L <- eta0 - runif(1) * w
  R <- L + w
  while (log_post(L) > y_level) L <- L - w
  while (log_post(R) > y_level) R <- R + w

  # Shrinkage
  repeat {
    eta_prop <- runif(1, L, R)
    if (log_post(eta_prop) >= y_level) return(exp(eta_prop))
    if (eta_prop < eta0) L <- eta_prop else R <- eta_prop
  }
}

#' Update Horseshoe+ Scales on the Increments (Makalic-Schmidt 2016)
#'
#' Model (Eq. 3-5) in product form: omega_t ~ N(0, sigma^2 tau^2 kappa_t^2 eta_t^2 lambda_t^2),
#' lambda_t, eta_t ~ C+(0, 1), tau ~ C+(0, s_tau). Each half-Cauchy uses the
#' inverse-gamma auxiliary representation. With `horseshoe_plus = FALSE`, eta_t is
#' fixed at 1 (ordinary horseshoe, ablation 5).
#'
#' @param omega Trend increments of length n - 1
#' @param sigma Innovation standard deviation
#' @param tau_omega Current global scale scalar
#' @param lambda_omega Current local scales vector of length n - 1
#' @param eta_omega Current secondary scales vector of length n - 1
#' @param kappa Prior inflation multiplier vector of length n - 1
#' @param aux_lambda Auxiliary vector for lambda of length n - 1
#' @param aux_eta Auxiliary vector for eta of length n - 1
#' @param s_tau Scale parameter for global scale Cauchy prior
#' @param horseshoe_plus Logical; FALSE fixes eta_t = 1
#' @param update_tau Logical; FALSE keeps tau_omega fixed (plug-in value)
#'
#' @return A list with updated `tau_omega`, `lambda_omega`, `eta_omega`, `v_omega`,
#'   `v_omega_tilde`, `aux_lambda`, `aux_eta`
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
    s_tau,
    horseshoe_plus = TRUE,
    update_tau = TRUE
) {
  n_omega <- length(omega)
  w2 <- (omega / sigma)^2

  # 1. tau | omega, eta, nu_lambda with lambda integrated out (partially collapsed)
  if (isTRUE(update_tau)) {
    tau_omega <- collapsed_tau_update(tau_omega, omega / sigma, kappa * eta_omega, aux_lambda, s_tau)
  }

  # 2. lambda_t^2 | . ~ IG(1, 1/aux + w2 / (2 tau^2 kappa^2 eta^2)), then its auxiliary
  lambda2 <- 1.0 / rgamma(n_omega, shape = 1.0,
                          rate = 1.0 / aux_lambda + 0.5 * w2 / (tau_omega * kappa * eta_omega)^2)
  aux_lambda <- 1.0 / rgamma(n_omega, shape = 1.0, rate = 1.0 + 1.0 / lambda2)
  lambda_omega <- sqrt(lambda2)

  # 3. eta_t^2 (horseshoe+ only)
  if (isTRUE(horseshoe_plus)) {
    eta2 <- 1.0 / rgamma(n_omega, shape = 1.0,
                         rate = 1.0 / aux_eta + 0.5 * w2 / (tau_omega * kappa * lambda_omega)^2)
    aux_eta <- 1.0 / rgamma(n_omega, shape = 1.0, rate = 1.0 + 1.0 / eta2)
    eta_omega <- sqrt(eta2)
  } else {
    eta_omega <- rep(1.0, n_omega)
    aux_eta <- rep(1.0, n_omega)
  }

  local_sd <- kappa * eta_omega * lambda_omega
  v_omega_tilde <- pmin(1e12, pmax(1e-12, (tau_omega * local_sd)^2))

  list(
    tau_omega     = tau_omega,
    lambda_omega  = lambda_omega,
    eta_omega     = eta_omega,
    v_omega       = sigma^2 * v_omega_tilde,
    v_omega_tilde = v_omega_tilde,
    aux_lambda    = aux_lambda,
    aux_eta       = aux_eta
  )
}

#' Non-centred Metropolis update of the log-volatility persistence rho
#'
#' In dsp's dynamic horseshoe, h_t = m_h + h~_t with h~_1 = e_1 and
#' h~_t = rho h~_{t-1} + e_t (t >= 2), where the innovations e_t have rho-free
#' Polya-Gamma scale-mixture priors. Holding (e_1, ..., e_n) fixed and proposing rho'
#' by a random walk, h~(rho') follows from the AR recursion and the acceptance ratio
#' contains only the likelihood psi_t ~ N(0, exp(h_t)) and the prior
#' (rho + 1) / 2 ~ Beta(a, b). Several steps are taken per call.
#'
#' @param psi_tilde Offset-adjusted outliers (length n)
#' @param ht Current log-volatilities h_t = m_h + h~_t (length n)
#' @param m_h Current mean of the log-volatilities
#' @param rho Current persistence
#' @param prior Beta prior parameters for (rho + 1) / 2 (default c(10, 2), as in dsp)
#' @param n_steps Number of random-walk steps (default 5)
#' @param step Random-walk standard deviation (default 0.1)
#' @return A list with updated `rho` and `ht`
#' @export
rho_noncentred_update <- function(psi_tilde, ht, m_h, rho, prior = c(10, 2), n_steps = 5, step = 0.1) {
  n <- length(ht)
  h_til <- ht - m_h
  e <- c(h_til[1], h_til[-1] - rho * h_til[-n])
  loglik <- function(hh) sum(stats::dnorm(psi_tilde, 0, exp((m_h + hh) / 2), log = TRUE))
  logprior <- function(r) stats::dbeta((r + 1) / 2, prior[1], prior[2], log = TRUE)
  ll_curr <- loglik(h_til)
  for (k in seq_len(n_steps)) {
    rho_prop <- rho + step * rnorm(1)
    if (abs(rho_prop) >= 1) next
    h_prop <- as.numeric(stats::filter(e, rho_prop, method = "recursive"))
    ll_prop <- loglik(h_prop)
    if (log(runif(1)) < ll_prop - ll_curr + logprior(rho_prop) - logprior(rho)) {
      rho <- rho_prop
      h_til <- h_prop
      ll_curr <- ll_prop
    }
  }
  list(rho = rho, ht = m_h + h_til)
}

#' Update Outlier Scales (Independent or Dynamic)
#'
#' Independent (rho = 0, delta = 0): psi_t ~ N(0, tau_psi^2 eta_t^2 lambda_t^2) with
#' horseshoe+ scales and tau_psi ~ C+(0, 1/sqrt(n)).
#'
#' Dynamic (Eq. 6-7): h_t = log nu_t^2 = h~_t + delta a_t, where h~_t follows the
#' dynamic horseshoe AR(1) with mean m_h and persistence rho. Because
#' psi_t exp(-delta a_t / 2) ~ N(0, exp(h~_t)), dsp::sampleDSP is applied to the
#' offset-adjusted psi, so the CAPA offset enters the likelihood for h~.
#' The global-scale prior is centred at log(1/n), matching tau_psi ~ C+(0, 1/sqrt(n)).
#'
#' @param psi Outlier vector of length n
#' @param tau_psi Global outlier scale scalar (independent model)
#' @param lambda_psi,eta_psi Local / secondary outlier scales of length n
#' @param aux_lambda_psi,aux_eta_psi Auxiliary parameters
#' @param dynamic Logical: TRUE uses dsp::sampleDSP; FALSE uses independent horseshoe+
#' @param delta CAPA offset parameter
#' @param a_vec Binary indicator vector of CAPA point anomalies
#' @param dsp_state Internal state list for dsp::sampleDSP (NULL to initialise)
#'
#' @return A list with updated outlier variances `nu2` and the sampler state
#' @export
update_outlier_scales <- function(
    psi,
    tau_psi = NULL,
    lambda_psi = NULL,
    eta_psi = NULL,
    aux_lambda_psi = NULL,
    aux_eta_psi = NULL,
    dynamic = FALSE,
    delta = 0.0,
    a_vec = NULL,
    dsp_state = NULL
) {
  n <- length(psi)
  s_tau_psi <- 1.0 / sqrt(n)

  if (!dynamic) {
    if (is.null(tau_psi)) tau_psi <- s_tau_psi
    if (is.null(lambda_psi)) lambda_psi <- rep(1.0, n)
    if (is.null(eta_psi)) eta_psi <- rep(1.0, n)
    if (is.null(aux_lambda_psi)) aux_lambda_psi <- rep(1.0, n)
    if (is.null(aux_eta_psi)) aux_eta_psi <- rep(1.0, n)

    p2 <- psi^2
    tau_psi <- collapsed_tau_update(tau_psi, psi, eta_psi, aux_lambda_psi, s_tau_psi)

    lam2 <- 1.0 / rgamma(n, shape = 1.0, rate = 1.0 / aux_lambda_psi + 0.5 * p2 / (tau_psi * eta_psi)^2)
    aux_lambda_psi <- 1.0 / rgamma(n, shape = 1.0, rate = 1.0 + 1.0 / lam2)
    lambda_psi <- sqrt(lam2)

    eta2 <- 1.0 / rgamma(n, shape = 1.0, rate = 1.0 / aux_eta_psi + 0.5 * p2 / (tau_psi * lambda_psi)^2)
    aux_eta_psi <- 1.0 / rgamma(n, shape = 1.0, rate = 1.0 + 1.0 / eta2)
    eta_psi <- sqrt(eta2)

    local_sd <- eta_psi * lambda_psi

    list(
      nu2            = pmin(1e12, pmax(1e-12, (tau_psi * local_sd)^2)),
      tau_psi        = tau_psi,
      lambda_psi     = lambda_psi,
      eta_psi        = eta_psi,
      aux_lambda_psi = aux_lambda_psi,
      aux_eta_psi    = aux_eta_psi
    )
  } else {
    ns_dsp <- asNamespace("dsp")
    offset <- if (!is.null(a_vec)) delta * a_vec else rep(0.0, n)
    psi_tilde <- matrix(psi * exp(-0.5 * offset), ncol = 1)

    if (is.null(dsp_state)) {
      dsp_state <- ns_dsp$initDHS(psi_tilde)
      dsp_state$sigma_wt <- NULL
      # sampleDSP skips the rho update when dhs_phi is exactly 0 or 1
      if (!is.finite(dsp_state$dhs_phi) || abs(dsp_state$dhs_phi) < 0.05 || abs(dsp_state$dhs_phi) > 0.95) {
        dsp_state$dhs_phi <- 0.5
      }
    }

    dsp_out <- ns_dsp$sampleDSP(
      omega         = psi_tilde,
      evolParams    = dsp_state,
      sigma_e       = s_tau_psi,
      loc           = NULL,
      prior_dhs_phi = c(10, 2),
      alphaPlusBeta = 1.0
    )
    dsp_out$sigma_wt <- NULL

    # Non-centred update of rho: h~ and rho are strongly coupled in dsp's centred
    # update, so rho is also moved holding the log-volatility innovations fixed
    nc <- rho_noncentred_update(as.vector(psi_tilde), as.vector(dsp_out$ht), dsp_out$dhs_mean, dsp_out$dhs_phi)
    dsp_out$ht <- matrix(nc$ht, ncol = 1)
    dsp_out$dhs_phi <- nc$rho

    list(
      nu2       = pmin(1e12, pmax(1e-12, exp(as.vector(dsp_out$ht) + offset))),
      dsp_state = dsp_out,
      rho       = as.numeric(dsp_out$dhs_phi),
      m_h       = as.numeric(dsp_out$dhs_mean)
    )
  }
}
