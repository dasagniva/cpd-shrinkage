# ==============================================================================
# 40_scoring.R
# Posterior Outlier Scoring and Tail Exceedance Probabilities.
# Implements Section 3.5, Equations (8) and (9).
# ==============================================================================

#' Compute Outlier Score o_t (Equation 8)
#'
#' Computes the posterior mean outlier score o_t across M retained draws:
#'   o_t = (1 / M) * sum_{i=1}^M [ nu_t^{2(i)} / (nu_t^{2(i)} + sigma_{eff}^{2(i)}) ]
#'
#' By default, sigma_{eff}^{2(i)} = sigma^{2(i)} (the innovation variance).
#' If `use_marginal_variance = TRUE`, sigma_{eff}^{2(i)} = sigma^{2(i)} / (1 - phi^{2(i)}),
#' as specified for the sensitivity ablation.
#'
#' @param nu2_draws Matrix or array of posterior draws for nu_t^2.
#'   Can be of dimensions (M, n) where M is iterations and n is series length,
#'   or (n, M).
#' @param sigma2_draws Numeric vector of length M containing posterior draws of sigma^2.
#' @param phi_draws Optional numeric vector of length M containing posterior draws of phi.
#'   Required if `use_marginal_variance = TRUE`.
#' @param use_marginal_variance Logical; if TRUE, uses marginal variance sigma^2 / (1 - phi^2).
#'   Default is FALSE (uses innovation variance sigma^2).
#'
#' @return Numeric vector of length n containing posterior outlier scores o_t in [0, 1].
#' @export
compute_outlier_score <- function(
    nu2_draws,
    sigma2_draws,
    phi_draws = NULL,
    use_marginal_variance = FALSE
) {
  # Coerce to matrix
  if (is.vector(nu2_draws)) {
    nu2_mat <- matrix(nu2_draws, nrow = 1)
  } else {
    nu2_mat <- as.matrix(nu2_draws)
  }

  M <- length(sigma2_draws)
  stopifnot(M >= 1, all(sigma2_draws > 0))

  # Ensure rows are iterations (M) and columns are time (n)
  if (nrow(nu2_mat) != M) {
    if (ncol(nu2_mat) == M) {
      nu2_mat <- t(nu2_mat)
    } else {
      stop(sprintf(
        "Dimension mismatch: nu2_draws has dim (%d, %d), but sigma2_draws has length %d",
        nrow(nu2_mat), ncol(nu2_mat), M
      ))
    }
  }

  n <- ncol(nu2_mat)

  # Compute denominator scale sigma_{eff}^{2(i)}
  if (isTRUE(use_marginal_variance)) {
    if (is.null(phi_draws) || length(phi_draws) != M) {
      stop("phi_draws of length M must be supplied when use_marginal_variance = TRUE")
    }
    denom_scale <- sigma2_draws / pmax(0.01, 1.0 - phi_draws^2)
  } else {
    denom_scale <- sigma2_draws
  }

  # Vectorized calculation over M x n matrix
  # nu2 / (nu2 + denom_scale)
  # denom_scale is of length M, matching nrow of nu2_mat
  ratios <- nu2_mat / (nu2_mat + denom_scale)
  # Replace any NaN (e.g. 0/0) with 0
  ratios[is.na(ratios)] <- 0.0

  o_t <- colMeans(ratios)
  # Guarantee strictly bounded in [0, 1]
  pmax(0.0, pmin(1.0, o_t))
}

#' Compute Tail Exceedance Probabilities for Trend Increments p_omega (Equation 9)
#'
#' Evaluates p_t^omega = Pr(|omega_t| > d_omega | y).
#' Threshold d_omega = c_omega * sigma / sqrt(1 - phi^2), default c_omega = 1.5.
#'
#' @param omega_draws Matrix of posterior draws for omega_t of dimensions (M, n - 1) or (n - 1, M),
#'   or (M, n) if t=1 is included.
#' @param sigma_draws Numeric vector of length M of sigma draws.
#' @param phi_draws Numeric vector of length M of phi draws.
#' @param c_omega Multiplier for marginal noise standard deviation (default 1.5).
#' @param d_omega Optional fixed numeric scalar or vector of length M specifying threshold directly.
#'
#' @return Numeric vector of length matching columns of omega_draws containing probabilities in [0, 1].
#' @export
compute_tail_prob_omega <- function(
    omega_draws,
    sigma_draws,
    phi_draws,
    c_omega = 1.5,
    d_omega = NULL
) {
  if (is.vector(omega_draws)) {
    w_mat <- matrix(omega_draws, nrow = 1)
  } else {
    w_mat <- as.matrix(omega_draws)
  }

  M <- length(sigma_draws)
  stopifnot(M >= 1, length(phi_draws) == M)

  if (nrow(w_mat) != M) {
    if (ncol(w_mat) == M) {
      w_mat <- t(w_mat)
    } else {
      stop(sprintf(
        "Dimension mismatch: omega_draws has dim (%d, %d), but sigma_draws has length %d",
        nrow(w_mat), ncol(w_mat), M
      ))
    }
  }

  if (is.null(d_omega)) {
    sd_marginal <- sigma_draws / sqrt(pmax(0.01, 1.0 - phi_draws^2))
    thresh <- c_omega * sd_marginal
  } else {
    thresh <- if (length(d_omega) == 1) rep(d_omega, M) else d_omega
    stopifnot(length(thresh) == M)
  }

  # Compare absolute draws against threshold for each iteration
  abs_w <- abs(w_mat)
  exceed <- abs_w > thresh
  p_omega <- colMeans(exceed)
  pmax(0.0, pmin(1.0, p_omega))
}

#' Compute Tail Exceedance Probabilities for Outliers p_psi (Equation 9)
#'
#' Evaluates p_t^psi = Pr(|psi_t| > d_psi | y).
#' Threshold d_psi = c_psi * sigma / sqrt(1 - phi^2), default c_psi = 1.5.
#'
#' @param psi_draws Matrix of posterior draws for psi_t of dimensions (M, n) or (n, M).
#' @param sigma_draws Numeric vector of length M of sigma draws.
#' @param phi_draws Numeric vector of length M of phi draws.
#' @param c_psi Multiplier for marginal noise standard deviation (default 1.5).
#' @param d_psi Optional fixed numeric scalar or vector of length M specifying threshold directly.
#'
#' @return Numeric vector of length n containing outlier probabilities in [0, 1].
#' @export
compute_tail_prob_psi <- function(
    psi_draws,
    sigma_draws,
    phi_draws,
    c_psi = 1.5,
    d_psi = NULL
) {
  if (is.vector(psi_draws)) {
    psi_mat <- matrix(psi_draws, nrow = 1)
  } else {
    psi_mat <- as.matrix(psi_draws)
  }

  M <- length(sigma_draws)
  stopifnot(M >= 1, length(phi_draws) == M)

  if (nrow(psi_mat) != M) {
    if (ncol(psi_mat) == M) {
      psi_mat <- t(psi_mat)
    } else {
      stop(sprintf(
        "Dimension mismatch: psi_draws has dim (%d, %d), but sigma_draws has length %d",
        nrow(psi_mat), ncol(psi_mat), M
      ))
    }
  }

  if (is.null(d_psi)) {
    sd_marginal <- sigma_draws / sqrt(pmax(0.01, 1.0 - phi_draws^2))
    thresh <- c_psi * sd_marginal
  } else {
    thresh <- if (length(d_psi) == 1) rep(d_psi, M) else d_psi
    stopifnot(length(thresh) == M)
  }

  abs_psi <- abs(psi_mat)
  exceed <- abs_psi > thresh
  p_psi <- colMeans(exceed)
  pmax(0.0, pmin(1.0, p_psi))
}

#' Compute Full Set of Posterior Decision Scores
#'
#' Comprehensive convenience wrapper to compute o_t, p^omega_t, and p^psi_t from draws.
#'
#' @param nu2_draws Matrix (M x n) of nu^2 draws
#' @param psi_draws Matrix (M x n) of psi draws
#' @param omega_draws Matrix (M x (n-1)) of omega draws
#' @param sigma_draws Vector (M) of sigma draws
#' @param phi_draws Vector (M) of phi draws
#' @param c_omega Multiplier for omega threshold (default 1.5)
#' @param c_psi Multiplier for psi threshold (default 1.5)
#' @param use_marginal_variance Logical (default FALSE)
#'
#' @return A list with elements:
#'   \item{o_t}{Numeric vector (length n) of outlier scores}
#'   \item{p_omega}{Numeric vector (length n-1) of increment probabilities}
#'   \item{p_psi}{Numeric vector (length n) of outlier probabilities}
#' @export
compute_all_scores <- function(
    nu2_draws,
    psi_draws,
    omega_draws,
    sigma_draws,
    phi_draws,
    c_omega = 1.5,
    c_psi = 1.5,
    use_marginal_variance = FALSE
) {
  sigma2_draws <- sigma_draws^2

  o_t <- compute_outlier_score(
    nu2_draws              = nu2_draws,
    sigma2_draws           = sigma2_draws,
    phi_draws              = phi_draws,
    use_marginal_variance  = use_marginal_variance
  )

  p_omega <- compute_tail_prob_omega(
    omega_draws = omega_draws,
    sigma_draws = sigma_draws,
    phi_draws   = phi_draws,
    c_omega     = c_omega
  )

  p_psi <- compute_tail_prob_psi(
    psi_draws   = psi_draws,
    sigma_draws = sigma_draws,
    phi_draws   = phi_draws,
    c_psi       = c_psi
  )

  list(
    o_t     = o_t,
    p_omega = p_omega,
    p_psi   = p_psi
  )
}
