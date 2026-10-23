# ==============================================================================
# 31_gibbs_sampler.R
# Production Gibbs Sampler for Significance-Guided Global-Local Shrinkage (Phase P4).
# Implements Section 3 (Equations 1-9) and Section 5 (Phase P4).
#
# The series is standardised with robust_standardize() before sampling; all
# reported levels (mu, psi, sigma, tau_psi) are back-transformed to data units.
# Scores (o_t, p^omega_t, p^psi_t) are scale-invariant.
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
  library(spam)
  library(posterior)
  library(future.apply)
})

#' Warm-start level for the Gibbs sampler (Section 5, P4)
#'
#' Uses the DeCAFS fit when available, otherwise a running median.
#'
#' @param z Standardised series
#' @return Numeric vector of initial mu values
#' @keywords internal
warm_start_mu <- function(z) {
  n <- length(z)
  mu0 <- NULL
  if (requireNamespace("DeCAFS", quietly = TRUE) && n >= 20) {
    mu0 <- tryCatch(as.numeric(DeCAFS::DeCAFS(z, warningMessage = FALSE)$signal), error = function(e) NULL)
  }
  if (is.null(mu0) || length(mu0) != n || any(!is.finite(mu0))) {
    k_med <- min(11L, 2L * as.integer(floor(n / 10L)) + 1L)
    mu0 <- as.numeric(stats::runmed(z, k = max(3L, k_med)))
  }
  mu0
}

#' Minimum-persistence split of one posterior draw (spec change C5)
#'
#' The split is applied to the total signal s = mu + psi, so it does not depend on how
#' the sampler divided a feature between the level and the outlier process (the two
#' are not separately identified for short features). A running median of width
#' 2m - 1 keeps level changes and anything lasting at least m points; excursions
#' lasting at most m - 1 points are removed:
#'   mu* = runmed(mu + psi, 2m - 1),  psi* = (mu + psi) - mu*.
#' Short spikes therefore end up in psi* even if the draw held them in mu, and a
#' persistent block held in psi (possible in the dynamic outlier model when rho > 0)
#' ends up in the level mu*. mu* + psi* = mu + psi, so the fitted signal is unchanged.
#'
#' @param mu,psi One draw of the level and outlier processes
#' @param m Minimum number of points a level shift must persist (m = 1: no split)
#' @return A list with the persistent level `mu` and the transient outliers `psi`
#' @export
persistence_split <- function(mu, psi, m) {
  k <- 2L * as.integer(m) - 1L
  if (k <= 1L || length(mu) < k) return(list(mu = mu, psi = psi))
  total <- mu + psi
  mu_star <- as.numeric(stats::runmed(total, k = k, endrule = "median", algorithm = "Turlach"))
  list(mu = mu_star, psi = total - mu_star)
}

#' Run a Single Chain of the Production Gibbs Sampler
#'
#' @param chain_id Integer chain identifier
#' @param y Numeric vector of length n (already standardised)
#' @param kappa Numeric vector of length n for local shrinkage multipliers
#' @param a_vec Numeric binary vector of length n for CAPA point offsets
#' @param n_iter Total iterations (burn_in + retained)
#' @param burn_in Number of initial burn-in iterations
#' @param thin Thinning interval (default 1)
#' @param dynamic Logical: if TRUE, use dynamic horseshoe for outliers (rho, delta)
#' @param delta CAPA offset parameter for dynamic outliers (default 0.0)
#' @param s_tau Scale parameter for tau_omega prior (default 1 / sqrt(n))
#' @param seed Optional integer seed (only for direct calls; gibbs_sampler uses RNG streams)
#' @param fix_phi Optional fixed value of phi (ablation 4)
#' @param horseshoe_plus Logical; FALSE gives the ordinary horseshoe (ablation 5)
#' @param innovation_var Logical; FALSE uses the marginal variance in Eq. 8 (ablation 6)
#' @param d_omega,d_psi Threshold multipliers of the marginal noise SD (Eq. 9)
#' @param cp_window Half-width h of the windowed shift probability (0 = Eq. 9 exactly)
#' @param tau_omega_fixed Optional fixed value of tau_omega (plug-in); NULL samples it
#' @param min_shift_len Minimum persistence m of a level shift. In every retained draw,
#'   the total signal mu + psi is split into a persistent part mu* = runmed(mu + psi, 2m - 1)
#'   and a transient part psi* = mu + psi - mu* (persistence_split). Shift and outlier
#'   probabilities are computed from (mu*, psi*). m = 1 disables the split.
#' @param marginal_scales Logical; add Metropolis updates of log tau_omega and of the
#'   outlier global scale (tau_psi, or exp(m_h / 2) in the dynamic model) with
#'   (mu, psi) integrated out (Johndrow, Orenstein & Bhattacharya 2020). Step sizes
#'   adapt during burn-in only.
#' @param phi_prior Shape parameters (a_phi, b_phi) of the Beta prior on (phi + 1) / 2
#' @param sigma2_prior (shape, rate) of an inverse-gamma prior on sigma^2; c(0, 0) is the
#'   Jeffreys prior of the spec
#' @param geweke_fun Optional function of the current state. When given, the chain becomes
#'   the successive-conditional simulator of Geweke (2004): after every sweep, y is redrawn
#'   from p(y | mu, psi, phi, sigma), and after burn-in geweke_fun(state) is recorded in
#'   `geweke_draws` (one row per iteration). Scores and traces are then not computed.
#' @param init_state Optional complete initial state (mu, psi, phi, sigma, tau_omega,
#'   lambda_omega, eta_omega, aux_lambda_omega, aux_eta_omega, nu2, and either tau_psi,
#'   lambda_psi, eta_psi, aux_lambda_psi, aux_eta_psi or the dsp state list `dsp_state`)
#'   replacing the warm start; used to start the Geweke test from exact joint draws
#'
#' @return A list containing online summaries and global-parameter traces
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
    seed = NULL,
    fix_phi = NULL,
    horseshoe_plus = TRUE,
    innovation_var = TRUE,
    d_omega = 1.0,
    d_psi = 1.5,
    cp_window = 5L,
    tau_omega_fixed = NULL,
    min_shift_len = 4L,
    marginal_scales = TRUE,
    phi_prior = c(10, 2),
    sigma2_prior = c(0, 0),
    geweke_fun = NULL,
    init_state = NULL
) {
  if (!is.null(seed)) set.seed(as.integer(seed))
  # Closures created inside the updates (e.g. dsp's slice-sampler target) would otherwise
  # be byte-compiled afresh on every iteration; the heavy work is vectorised or in C++
  old_jit <- compiler::enableJIT(0)
  on.exit(compiler::enableJIT(old_jit), add = TRUE)

  n <- length(y)
  n_omega <- n - 1
  h_win <- as.integer(cp_window)

  if (is.null(kappa)) kappa <- rep(1.0, n)
  kappa_w <- kappa[2:n]
  if (is.null(a_vec)) a_vec <- rep(0.0, n)
  if (is.null(s_tau)) s_tau <- 1.0 / sqrt(n)

  banded_struct <- init_banded_structure(n)

  # 1. Warm start: robust segmentation for mu, psi = 0
  mu_curr <- warm_start_mu(y)
  psi_curr <- rep(0.0, n)
  res_init <- y - mu_curr
  r1 <- suppressWarnings(stats::cor(res_init[-1], res_init[-n]))
  phi_curr <- if (!is.null(fix_phi)) as.numeric(fix_phi) else pmax(-0.9, pmin(0.9, ifelse(is.finite(r1), r1, 0.5)))
  sigma_curr <- max(1e-3, stats::mad(diff(res_init)) / sqrt(2.0))

  # Local scales start consistent with the warm-start increments. Under horseshoe
  # shrinkage a closed site (small lambda * eta) almost never reopens, so sites where
  # the robust segmentation already shows a jump must start open:
  # lambda_t * eta_t = max(1, |omega0_t| / (sigma tau kappa_t)).
  tau_omega_curr <- if (!is.null(tau_omega_fixed)) tau_omega_fixed else s_tau
  open_ratio <- pmax(1.0, abs(diff(mu_curr)) / (sigma_curr * tau_omega_curr * kappa_w))
  lambda_omega_curr <- sqrt(open_ratio)
  eta_omega_curr <- if (isTRUE(horseshoe_plus)) sqrt(open_ratio) else rep(1.0, n_omega)
  if (!isTRUE(horseshoe_plus)) lambda_omega_curr <- open_ratio
  aux_lambda_w <- rep(1.0, n_omega)
  aux_eta_w <- rep(1.0, n_omega)
  v_omega_tilde <- (tau_omega_curr * kappa_w * lambda_omega_curr * eta_omega_curr)^2
  v_omega <- sigma_curr^2 * v_omega_tilde

  tau_psi_curr <- 1.0 / sqrt(n)
  lambda_psi_curr <- rep(1.0, n)
  eta_psi_curr <- rep(1.0, n)
  aux_lambda_psi <- rep(1.0, n)
  aux_eta_psi <- rep(1.0, n)
  nu2_curr <- rep(tau_psi_curr^2, n)
  dsp_state <- NULL
  rho_curr <- NA_real_
  m_h_curr <- NA_real_

  if (!is.null(init_state)) {
    s0 <- init_state
    mu_curr <- s0$mu
    psi_curr <- s0$psi
    phi_curr <- s0$phi
    sigma_curr <- s0$sigma
    tau_omega_curr <- s0$tau_omega
    lambda_omega_curr <- s0$lambda_omega
    eta_omega_curr <- s0$eta_omega
    aux_lambda_w <- s0$aux_lambda_omega
    aux_eta_w <- s0$aux_eta_omega
    v_omega_tilde <- pmin(1e12, pmax(1e-12, (tau_omega_curr * kappa_w * lambda_omega_curr * eta_omega_curr)^2))
    v_omega <- sigma_curr^2 * v_omega_tilde
    nu2_curr <- s0$nu2
    if (!dynamic) {
      tau_psi_curr <- s0$tau_psi
      lambda_psi_curr <- s0$lambda_psi
      eta_psi_curr <- s0$eta_psi
      aux_lambda_psi <- s0$aux_lambda_psi
      aux_eta_psi <- s0$aux_eta_psi
    } else {
      dsp_state <- s0$dsp_state
      rho_curr <- dsp_state$dhs_phi
      m_h_curr <- dsp_state$dhs_mean
    }
  }

  n_saved <- floor((n_iter - burn_in) / thin)
  stopifnot(n_saved >= 1)

  mu_sum        <- numeric(n)
  mu2_sum       <- numeric(n)
  psi_sum       <- numeric(n)
  o_t_sum       <- numeric(n)
  abs_omega_sum <- numeric(n_omega)
  p_omega_count <- numeric(n_omega)
  p_win_count   <- numeric(n_omega)
  p_psi_count   <- numeric(n)

  trace_phi       <- numeric(n_saved)
  trace_sigma     <- numeric(n_saved)
  trace_tau_omega <- numeric(n_saved)
  trace_tau_psi   <- numeric(n_saved)
  trace_n_shift   <- numeric(n_saved)
  trace_n_out     <- numeric(n_saved)
  trace_rho       <- if (dynamic) numeric(n_saved) else NULL
  trace_m_h       <- if (dynamic) numeric(n_saved) else NULL
  phi_accept      <- 0L

  # Index pairs for the windowed shift probability: mu_{t+h} - mu_{t-h-1}, t = 2..n
  if (h_win > 0) {
    t_idx <- 2:n
    hi <- pmin(n, t_idx + h_win)
    lo <- pmax(1L, t_idx - h_win - 1L)
  }

  save_idx <- 0L
  geweke_draws <- NULL

  # Adaptive random-walk step sizes (log scale) for the marginal scale updates
  log_half_cauchy <- function(x, s) -log1p((x / s)^2)
  s_psi <- 1.0 / sqrt(n)
  step_tau <- 0.1
  step_c <- 0.1
  acc_tau <- 0L
  acc_c <- 0L
  scale_accept <- c(tau_omega = 0L, outlier = 0L)

  for (iter in seq_len(n_iter)) {
    fac <- NULL

    # Step 0: global scales with (mu, psi) integrated out (blocked with Step 1)
    if (isTRUE(marginal_scales)) {
      fac <- factor_mu_psi(y, phi_curr, sigma_curr, v_omega, nu2_curr, banded_struct)

      # 0a. log tau_omega (skipped when tau_omega is fixed at its plug-in value)
      if (is.null(tau_omega_fixed)) {
      tau_prop <- tau_omega_curr * exp(step_tau * rnorm(1))
      ratio2 <- (tau_prop / tau_omega_curr)^2
      fac_prop <- factor_mu_psi(y, phi_curr, sigma_curr, v_omega * ratio2, nu2_curr, banded_struct)
      log_acc <- fac_prop$logml - fac$logml +
        log_half_cauchy(tau_prop, s_tau) - log_half_cauchy(tau_omega_curr, s_tau) + log(tau_prop / tau_omega_curr)
      if (log(runif(1)) < log_acc) {
        tau_omega_curr <- tau_prop
        v_omega <- v_omega * ratio2
        v_omega_tilde <- v_omega_tilde * ratio2
        fac <- fac_prop
        acc_tau <- acc_tau + 1L
        if (iter > burn_in) scale_accept[1] <- scale_accept[1] + 1L
      }
      }

      # 0b. log outlier global scale (tau_psi, or exp(m_h / 2) holding the centred log-vols fixed)
      if (!dynamic || !is.null(dsp_state)) {
        c_cur <- if (!dynamic) tau_psi_curr else exp(m_h_curr / 2)
        c_prop <- c_cur * exp(step_c * rnorm(1))
        ratio2 <- (c_prop / c_cur)^2
        fac_prop <- factor_mu_psi(y, phi_curr, sigma_curr, v_omega, nu2_curr * ratio2, banded_struct)
        log_acc <- fac_prop$logml - fac$logml +
          log_half_cauchy(c_prop, s_psi) - log_half_cauchy(c_cur, s_psi) + log(c_prop / c_cur)
        if (log(runif(1)) < log_acc) {
          nu2_curr <- nu2_curr * ratio2
          if (!dynamic) {
            tau_psi_curr <- c_prop
          } else {
            shift <- log(ratio2)
            dsp_state$ht <- dsp_state$ht + shift
            dsp_state$dhs_mean <- dsp_state$dhs_mean + shift
            dsp_state$dhs_mean0 <- dsp_state$dhs_mean0 + shift
            m_h_curr <- m_h_curr + shift
          }
          fac <- fac_prop
          acc_c <- acc_c + 1L
          if (iter > burn_in) scale_accept[2] <- scale_accept[2] + 1L
        }
      }

      if (iter <= burn_in && iter %% 50 == 0) {
        step_tau <- step_tau * exp(acc_tau / 50 - 0.35)
        step_c <- step_c * exp(acc_c / 50 - 0.35)
        acc_tau <- 0L
        acc_c <- 0L
      }
    }

    # Step 1: Joint (mu, psi) Gaussian draw
    mp <- update_mu_psi(y, phi_curr, sigma_curr, v_omega, nu2_curr, banded_struct, fac = fac)
    mu_curr <- mp$mu
    psi_curr <- mp$psi
    omega_curr <- diff(mu_curr)

    # Step 2: phi
    e_curr <- y - mu_curr - psi_curr
    if (is.null(fix_phi)) {
      phi_res <- update_phi(e = e_curr, phi_curr = phi_curr, sigma = sigma_curr,
                            a_phi = phi_prior[1], b_phi = phi_prior[2])
      phi_curr <- phi_res$phi
      if (iter > burn_in) phi_accept <- phi_accept + as.integer(phi_res$accepted)
    }

    # Step 3: sigma^2 from prewhitened innovations and scaled increments
    u_curr <- c(sqrt(1.0 - phi_curr^2) * e_curr[1], e_curr[-1] - phi_curr * e_curr[-n])
    sigma2_curr <- update_sigma2(u_curr, omega_curr, v_omega_tilde, prior = sigma2_prior)
    sigma_curr <- sqrt(sigma2_curr)
    sigma_marg_curr <- sigma_curr / sqrt(1.0 - phi_curr^2)

    # Step 4: horseshoe(+) scales on increments
    hs_res <- update_hsplus_scales(
      omega          = omega_curr,
      sigma          = sigma_curr,
      tau_omega      = tau_omega_curr,
      lambda_omega   = lambda_omega_curr,
      eta_omega      = eta_omega_curr,
      kappa          = kappa_w,
      aux_lambda     = aux_lambda_w,
      aux_eta        = aux_eta_w,
      s_tau          = s_tau,
      horseshoe_plus = horseshoe_plus,
      update_tau     = is.null(tau_omega_fixed)
    )
    tau_omega_curr    <- hs_res$tau_omega
    lambda_omega_curr <- hs_res$lambda_omega
    eta_omega_curr    <- hs_res$eta_omega
    aux_lambda_w      <- hs_res$aux_lambda
    aux_eta_w         <- hs_res$aux_eta
    v_omega           <- hs_res$v_omega
    v_omega_tilde     <- hs_res$v_omega_tilde

    # Step 5: outlier scales (independent or dynamic)
    out_res <- update_outlier_scales(
      psi            = psi_curr,
      tau_psi        = tau_psi_curr,
      lambda_psi     = lambda_psi_curr,
      eta_psi        = eta_psi_curr,
      aux_lambda_psi = aux_lambda_psi,
      aux_eta_psi    = aux_eta_psi,
      dynamic        = dynamic,
      delta          = delta,
      a_vec          = a_vec,
      dsp_state      = dsp_state
    )
    nu2_curr <- out_res$nu2
    if (!dynamic) {
      tau_psi_curr    <- out_res$tau_psi
      lambda_psi_curr <- out_res$lambda_psi
      eta_psi_curr    <- out_res$eta_psi
      aux_lambda_psi  <- out_res$aux_lambda_psi
      aux_eta_psi     <- out_res$aux_eta_psi
    } else {
      dsp_state <- out_res$dsp_state
      rho_curr  <- out_res$rho
      m_h_curr  <- out_res$m_h
    }

    # Geweke (2004) successive-conditional simulator: redraw y | theta, record test functions
    if (!is.null(geweke_fun)) {
      u_new <- c(rnorm(1, 0, sigma_curr / sqrt(1 - phi_curr^2)), sigma_curr * rnorm(n - 1))
      y <- mu_curr + psi_curr + as.numeric(stats::filter(u_new, phi_curr, method = "recursive"))
      if (iter > burn_in) {
        g <- geweke_fun(list(y = y, mu = mu_curr, psi = psi_curr, phi = phi_curr, sigma = sigma_curr,
                             tau_omega = tau_omega_curr, lambda_omega = lambda_omega_curr,
                             eta_omega = eta_omega_curr, nu2 = nu2_curr, tau_psi = tau_psi_curr,
                             rho = rho_curr, m_h = m_h_curr))
        if (is.null(geweke_draws)) geweke_draws <- matrix(NA_real_, n_iter - burn_in, length(g),
                                                          dimnames = list(NULL, names(g)))
        geweke_draws[iter - burn_in, ] <- g
      }
      next
    }

    # Step 6: retained draws
    if (iter > burn_in && ((iter - burn_in) %% thin == 0)) {
      save_idx <- save_idx + 1L

      trace_phi[save_idx]       <- phi_curr
      trace_sigma[save_idx]     <- sigma_curr
      trace_tau_omega[save_idx] <- tau_omega_curr
      trace_tau_psi[save_idx]   <- if (!dynamic) tau_psi_curr else exp(m_h_curr / 2)
      if (dynamic) {
        trace_rho[save_idx] <- rho_curr
        trace_m_h[save_idx] <- m_h_curr
      }

      # Minimum-persistence split of this draw: excursions of mu shorter than
      # min_shift_len points are reassigned from the level to the outlier process
      split <- persistence_split(mu_curr, psi_curr, min_shift_len)
      mu_star <- split$mu
      psi_star <- split$psi
      omega_star <- diff(mu_star)

      mu_sum  <- mu_sum + mu_star
      mu2_sum <- mu2_sum + mu_star^2
      psi_sum <- psi_sum + psi_star

      # Eq. 8: innovation variance (default) or marginal variance (ablation 6)
      denom_var <- if (isTRUE(innovation_var)) sigma2_curr else sigma_marg_curr^2
      o_t_sum <- o_t_sum + nu2_curr / (nu2_curr + denom_var)

      # Eq. 9 with thresholds proportional to the marginal noise SD
      d_w <- d_omega * sigma_marg_curr
      abs_omega_sum <- abs_omega_sum + abs(omega_star)
      p_omega_count <- p_omega_count + (abs(omega_star) > d_w)
      p_psi_count   <- p_psi_count + (abs(psi_star) > d_psi * sigma_marg_curr)

      exceed_psi <- abs(psi_star) > d_psi * sigma_marg_curr
      exceed_w <- if (h_win > 0) abs(mu_star[hi] - mu_star[lo]) > d_w else abs(omega_star) > d_w
      if (h_win > 0) p_win_count <- p_win_count + exceed_w
      # Decision functionals (invariant to how a spike is split between mu and psi)
      trace_n_shift[save_idx] <- sum(exceed_w)
      trace_n_out[save_idx]   <- sum(exceed_psi)
    }
  }

  list(
    chain_id        = chain_id,
    n_saved         = n_saved,
    mu_sum          = mu_sum,
    mu2_sum         = mu2_sum,
    psi_sum         = psi_sum,
    o_t_sum         = o_t_sum,
    abs_omega_sum   = abs_omega_sum,
    p_omega_count   = p_omega_count,
    p_win_count     = p_win_count,
    p_psi_count     = p_psi_count,
    phi_accept      = phi_accept,
    scale_accept    = scale_accept,
    trace_phi       = trace_phi,
    trace_sigma     = trace_sigma,
    trace_tau_omega = trace_tau_omega,
    trace_tau_psi   = trace_tau_psi,
    trace_n_shift   = trace_n_shift,
    trace_n_out     = trace_n_out,
    trace_rho       = trace_rho,
    trace_m_h       = trace_m_h,
    geweke_draws    = geweke_draws
  )
}

#' Production Gibbs Sampler for CPD Shrinkage
#'
#' Runs a multi-chain Gibbs sampler on the robustly standardised series with
#' warm-start, banded Cholesky (mu, psi) updates, horseshoe+ increment scales and
#' dynamic/independent outlier scales. Chains use independent L'Ecuyer-CMRG
#' streams derived from `seed`; set a `future::plan()` to run them in parallel.
#'
#' @param y Numeric vector of time series observations (any location/scale)
#' @param kappa Local shrinkage multiplier vector of length n (from build_prior_maps)
#' @param a_vec Binary indicator vector of length n for CAPA point offsets
#' @param n_iter Total iterations per chain (default 10000)
#' @param burn_in Number of burn-in iterations per chain (default 5000)
#' @param thin Thinning interval (default 1)
#' @param n_chains Number of chains (default 4)
#' @param dynamic Logical: if TRUE, use dynamic horseshoe for outliers (rho, delta)
#' @param delta CAPA offset parameter for dynamic outliers (default 0.0)
#' @param alpha_fdr Bayesian FDR level for candidate selection (default 0.05)
#' @param capa_collective Optional CAPA collective segments for rule 3.6(1)
#' @param L_max,r Collective-anomaly rule parameters (Section 3.6)
#' @param seed Optional master random seed
#' @param s_tau Scale of the tau_omega prior (default 1 / sqrt(n))
#' @param fix_phi,horseshoe_plus,innovation_var,d_omega,d_psi,cp_window,min_shift_len,marginal_scales
#'   See run_gibbs_chain. Defaults follow SPEC_CHANGES.md: d_omega = 1.0, d_psi = 1.5,
#'   cp_window = 5, min_shift_len = 4 (the spec's Eq. 9 is d_omega = d_psi = 1.5,
#'   cp_window = 0, min_shift_len = 1).
#' @param tau_omega_mode "plugin" (default) fixes tau_omega at (p0 / (1 - p0)) / sqrt(n);
#'   "sample" gives it a C+(0, s_tau) prior and samples it
#' @param p0 Expected fraction of time points carrying a change (default 0.02)
#' @param phi_prior Shape parameters (a_phi, b_phi) of the Beta prior on (phi + 1) / 2
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
    capa_collective = NULL,
    L_max          = NULL,
    r              = 0.3,
    seed           = NULL,
    s_tau          = NULL,
    fix_phi        = NULL,
    horseshoe_plus = TRUE,
    innovation_var = TRUE,
    d_omega        = 1.0,
    d_psi          = 1.5,
    cp_window      = 5L,
    tau_omega_mode = c("plugin", "sample"),
    p0             = 0.02,
    min_shift_len  = 4L,
    marginal_scales = TRUE,
    phi_prior      = c(10, 2)
) {
  start_time <- proc.time()[3]
  tau_omega_mode <- match.arg(tau_omega_mode)
  n <- length(y)
  stopifnot(n >= 10, is.numeric(y), all(is.finite(y)))

  if (is.null(kappa)) kappa <- rep(1.0, n)
  if (is.null(a_vec)) a_vec <- rep(0.0, n)
  stopifnot(length(kappa) == n, length(a_vec) == n)

  std <- robust_standardize(y)
  z <- std$z

  # Plug-in global scale (Piironen & Vehtari 2017 form, with p0 the expected fraction
  # of time points that carry a change); spec Section 9 fallback, change C3
  tau_omega_fixed <- if (tau_omega_mode == "plugin") (p0 / (1 - p0)) / sqrt(n) else NULL

  if (!is.null(seed)) set.seed(as.integer(seed))

  chain_fun <- function(cid) {
    run_gibbs_chain(
      chain_id       = cid,
      y              = z,
      kappa          = kappa,
      a_vec          = a_vec,
      n_iter         = n_iter,
      burn_in        = burn_in,
      thin           = thin,
      dynamic        = dynamic,
      delta          = delta,
      s_tau          = s_tau,
      fix_phi        = fix_phi,
      horseshoe_plus = horseshoe_plus,
      innovation_var = innovation_var,
      d_omega        = d_omega,
      d_psi          = d_psi,
      cp_window      = cp_window,
      tau_omega_fixed = tau_omega_fixed,
      min_shift_len  = min_shift_len,
      marginal_scales = marginal_scales,
      phi_prior      = phi_prior
    )
  }

  chains_res <- future.apply::future_lapply(seq_len(n_chains), chain_fun, future.seed = TRUE)

  runtime_sec <- as.numeric(proc.time()[3] - start_time)

  total_saved <- sum(vapply(chains_res, function(c) c$n_saved, numeric(1)))
  n_saved_per_chain <- chains_res[[1]]$n_saved
  pool <- function(field) Reduce("+", lapply(chains_res, function(c) c[[field]])) / total_saved

  mu_mean_z  <- pool("mu_sum")
  mu_sd_z    <- sqrt(pmax(0.0, pool("mu2_sum") - mu_mean_z^2))
  psi_mean_z <- pool("psi_sum")
  o_t_mean   <- pool("o_t_sum")
  p_omega    <- pool("p_omega_count")
  p_psi      <- pool("p_psi_count")
  abs_omega  <- pool("abs_omega_sum")
  p_win      <- if (cp_window > 0) pool("p_win_count") else NULL

  # Global-parameter traces (sigma and tau_psi back-transformed to data units)
  trace_mat <- function(field, mult = 1.0) do.call(cbind, lapply(chains_res, function(c) c[[field]])) * mult
  var_list <- list(
    phi       = trace_mat("trace_phi"),
    sigma     = trace_mat("trace_sigma", std$scale),
    tau_omega = trace_mat("trace_tau_omega"),
    tau_psi   = trace_mat("trace_tau_psi", std$scale)
  )
  if (dynamic) {
    var_list$rho <- trace_mat("trace_rho")
    var_list$m_h <- trace_mat("trace_m_h")
  }
  # Decision functionals: per-draw number of exceeding shift windows and outliers
  var_list$n_shift_windows <- trace_mat("trace_n_shift")
  var_list$n_outliers <- trace_mat("trace_n_out")
  if (!is.null(fix_phi)) var_list$phi <- NULL
  if (!is.null(tau_omega_fixed)) var_list$tau_omega <- NULL

  draws_arr <- array(
    NA_real_,
    dim = c(n_saved_per_chain, n_chains, length(var_list)),
    dimnames = list(iteration = NULL, chain = NULL, variable = names(var_list))
  )
  for (v in names(var_list)) draws_arr[, , v] <- var_list[[v]]
  draws_obj <- as_draws_array(draws_arr)
  summ_diag <- summarize_draws(draws_obj)

  phi_draws <- trace_mat("trace_phi")
  diag_list <- list(
    summary        = summ_diag,
    max_rhat       = if (n_chains > 1) max(summ_diag$rhat, na.rm = TRUE) else NA_real_,
    min_ess        = min(summ_diag$ess_bulk, na.rm = TRUE),
    phi_hat        = mean(phi_draws),
    phi_accept     = sum(vapply(chains_res, function(c) c$phi_accept, numeric(1))) / total_saved,
    scale_accept   = Reduce("+", lapply(chains_res, function(c) c$scale_accept)) / total_saved,
    sigma_hat      = mean(var_list$sigma),
    mu_mean        = mu_mean_z * std$scale + std$center,
    mu_sd          = mu_sd_z * std$scale,
    psi_mean       = psi_mean_z * std$scale,
    abs_omega_mean = abs_omega * std$scale,
    p_omega_point  = p_omega,
    p_omega_window = p_win,
    cp_window      = as.integer(cp_window),
    min_shift_len  = as.integer(min_shift_len),
    tau_omega      = if (!is.null(tau_omega_fixed)) c(mode = "plugin", value = tau_omega_fixed, p0 = p0) else c(mode = "sample"),
    standardization = std[c("center", "scale")],
    n_chains       = n_chains,
    n_saved        = n_saved_per_chain,
    draws          = draws_obj
  )

  classify_and_build_fit(
    y               = y,
    mu_hat          = diag_list$mu_mean,
    s_out           = o_t_mean,
    p_omega         = if (cp_window > 0) p_win else p_omega,
    p_psi           = p_psi,
    alpha_fdr       = alpha_fdr,
    capa_collective = capa_collective,
    method          = "proposed",
    runtime_sec     = runtime_sec,
    diag            = diag_list,
    L_max           = L_max,
    r               = r,
    cp_window       = cp_window,
    omega_strength  = abs_omega
  )
}

#' Fit the Full Proposed Method (front-ends + prior map + Gibbs sampler)
#'
#' Single entry point used by the simulation study, the real-data benchmark and the
#' acceptance gates, so all of them run the same configuration. Requires
#' 20_frontend_nsp.R, 21_frontend_capa.R and 22_prior_map.R to be sourced.
#'
#' @param y Numeric series
#' @param alpha_nsp NSP significance level (default 0.05)
#' @param kappa_in Local-scale inflation on the candidate set (default 2)
#' @param delta CAPA offset in the dynamic outlier model (default 1)
#' @param dynamic Logical; dynamic horseshoe for outliers (default TRUE)
#' @param ... Further arguments passed to gibbs_sampler (n_iter, burn_in, n_chains,
#'   alpha_fdr, d_omega, d_psi, cp_window, seed, ...)
#'
#' @return A validated `cpd_fit`; `diag$frontends` holds the NSP intervals and CAPA output
#' @export
fit_proposed <- function(y, alpha_nsp = 0.05, kappa_in = 2.0, delta = 1.0, dynamic = TRUE, ...) {
  n <- length(y)
  nsp_res <- run_frontend_nsp(y, alpha = alpha_nsp)
  capa_res <- run_frontend_capa(y, prewhiten = TRUE)
  pm <- build_prior_maps(n = n, nsp_intervals = nsp_res$intervals, capa_results = capa_res, kappa_in = kappa_in)

  fit <- gibbs_sampler(
    y               = y,
    kappa           = pm$kappa,
    a_vec           = pm$a,
    dynamic         = dynamic,
    delta           = if (dynamic) delta else 0.0,
    capa_collective = capa_res$collective,
    ...
  )
  fit$diag$frontends <- list(nsp_intervals = nsp_res$intervals, capa_point = capa_res$point,
                             capa_collective = capa_res$collective, capa_phi_hat = capa_res$phi_hat)
  fit
}
