# ==============================================================================
# test-32-gibbs-updates.R
# Unit tests for conditional Gibbs sampler updates (Phase P4).
# ==============================================================================

source("../../R/00_utils.R")
source("../../R/10_dgp.R")
source("../../R/32_gibbs_updates.R")

test_that("init_banded_structure initializes correct dimensions", {
  n <- 50
  bs <- init_banded_structure(n)
  expect_equal(bs$n, 50)
  expect_equal(bs$m, 100)
  expect_equal(dim(bs$Q_tmpl), c(100, 100))
  expect_s4_class(bs$Q_tmpl, "spam")
})

test_that("update_mu_psi produces valid dimensions and finite values", {
  n <- 40
  y <- rnorm(n)
  phi <- 0.5
  sigma <- 1.0
  v_omega <- rep(0.1, n - 1)
  nu2 <- rep(0.05, n)

  res <- update_mu_psi(y, phi, sigma, v_omega, nu2)
  expect_equal(length(res$mu), n)
  expect_equal(length(res$psi), n)
  expect_true(all(is.finite(res$mu)))
  expect_true(all(is.finite(res$psi)))
})

test_that("update_phi produces values in (-0.95, 0.95)", {
  n <- 50
  set.seed(42)
  e <- arima.sim(list(ar = 0.6), n = n)
  phi_curr <- 0.5
  sigma <- 1.0

  for (i in 1:10) {
    res <- update_phi(e, phi_curr, sigma)
    expect_true(res$phi > -0.95 && res$phi < 0.95)
    expect_true(is.logical(res$accepted))
    phi_curr <- res$phi
  }
})

test_that("update_sigma2 returns positive finite scalar", {
  n <- 50
  u <- rnorm(n, 0, 1)
  omega <- rnorm(n - 1, 0, 0.1)
  v_omega_tilde <- rep(0.01, n - 1)

  sig2 <- update_sigma2(u, omega, v_omega_tilde)
  expect_true(is.numeric(sig2))
  expect_equal(length(sig2), 1)
  expect_true(sig2 > 0)
  expect_true(is.finite(sig2))
})

test_that("update_hsplus_scales updates all Horseshoe+ components", {
  n <- 40
  n_omega <- n - 1
  omega <- rnorm(n_omega, 0, 0.05)
  sigma <- 1.0
  tau_omega <- 0.1
  lambda_omega <- rep(1.0, n_omega)
  eta_omega <- rep(1.0, n_omega)
  kappa <- rep(1.0, n_omega)
  aux_lambda <- rep(1.0, n_omega)
  aux_eta <- rep(1.0, n_omega)

  res <- update_hsplus_scales(
    omega = omega,
    sigma = sigma,
    tau_omega = tau_omega,
    lambda_omega = lambda_omega,
    eta_omega = eta_omega,
    kappa = kappa,
    aux_lambda = aux_lambda,
    aux_eta = aux_eta,
    s_tau = 1.0 / sqrt(n)
  )

  expect_true(res$tau_omega > 0)
  expect_equal(length(res$lambda_omega), n_omega)
  expect_equal(length(res$eta_omega), n_omega)
  expect_equal(length(res$v_omega), n_omega)
  expect_true(all(res$v_omega > 0))
  expect_true(all(is.finite(res$v_omega)))
})

test_that("update_outlier_scales works for independent model", {
  n <- 30
  psi <- rnorm(n, 0, 0.1)

  res <- update_outlier_scales(psi, dynamic = FALSE)
  expect_equal(length(res$nu2), n)
  expect_true(all(res$nu2 > 0))
  expect_true(all(is.finite(res$nu2)))
  expect_true(res$tau_psi > 0)
})

# ------------------------------------------------------------------------------
# Correctness checks against dense / numerical references
# ------------------------------------------------------------------------------

dense_mu_psi_system <- function(y, phi, sigma, v_omega, nu2) {
  n <- length(y)
  A <- matrix(0, n, 2 * n)
  A[cbind(1:n, 2 * (1:n) - 1)] <- 1
  A[cbind(1:n, 2 * (1:n))] <- 1
  P <- diag(n)
  P[1, 1] <- sqrt(1 - phi^2)
  P[cbind(2:n, 1:(n - 1))] <- -phi
  D <- matrix(0, n - 1, 2 * n)
  D[cbind(1:(n - 1), 2 * (2:n) - 1)] <- 1
  D[cbind(1:(n - 1), 2 * (1:(n - 1)) - 1)] <- -1
  Q <- t(A) %*% t(P) %*% P %*% A / sigma^2 + t(D) %*% diag(1 / v_omega, n - 1) %*% D
  Q[cbind(2 * (1:n), 2 * (1:n))] <- Q[cbind(2 * (1:n), 2 * (1:n))] + 1 / nu2
  Q[1, 1] <- Q[1, 1] + 1 / MU1_PRIOR_VAR
  Q <- Q + diag(1e-8, 2 * n)
  b <- as.vector(t(A) %*% t(P) %*% P %*% y) / sigma^2
  list(Q = Q, b = b)
}

test_that("banded precision and RHS match the dense construction", {
  set.seed(1)
  n <- 7
  y <- rnorm(n)
  v <- runif(n - 1, 0.05, 2)
  nu2 <- runif(n, 0.01, 3)
  bs <- init_banded_structure(n)
  sys <- build_mu_psi_system(y, 0.7, 1.3, v, nu2, bs)
  ref <- dense_mu_psi_system(y, 0.7, 1.3, v, nu2)
  expect_equal(as.matrix(sys$Q), ref$Q, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(sys$b, ref$b, tolerance = 1e-10)
})

test_that("update_mu_psi draws have the exact conditional mean and covariance", {
  set.seed(2)
  n <- 5
  y <- rnorm(n)
  v <- runif(n - 1, 0.1, 1)
  nu2 <- runif(n, 0.1, 1)
  ref <- dense_mu_psi_system(y, 0.4, 0.8, v, nu2)
  Sigma <- solve(ref$Q)
  m <- as.vector(Sigma %*% ref$b)

  bs <- init_banded_structure(n)
  N <- 20000
  draws <- t(replicate(N, {
    r <- update_mu_psi(y, 0.4, 0.8, v, nu2, bs)
    as.vector(rbind(r$mu, r$psi))
  }))
  se <- sqrt(diag(Sigma) / N)
  expect_true(all(abs(colMeans(draws) - m) < 5 * se))
  expect_equal(cov(draws), Sigma, tolerance = 0.05, ignore_attr = TRUE)
})

test_that("rtnorm1 is exact in the far tail and in the body", {
  set.seed(3)
  # Far tail: N(2, 0.05^2) truncated to [-0.95, 0.95] -> approx exponential below 0.95
  x <- replicate(5000, rtnorm1(2, 0.05, -0.95, 0.95))
  expect_true(all(x >= -0.95 & x <= 0.95))
  expect_equal(mean(x), 0.95 - 0.05^2 / 1.05, tolerance = 2e-4)

  # Body: N(0.3, 0.4^2) truncated to [-0.95, 0.95]
  a <- (-0.95 - 0.3) / 0.4
  b <- (0.95 - 0.3) / 0.4
  mean_ref <- 0.3 + 0.4 * (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
  x2 <- replicate(20000, rtnorm1(0.3, 0.4, -0.95, 0.95))
  expect_equal(mean(x2), mean_ref, tolerance = 0.01)
})

test_that("slice_sample_tau targets the half-Cauchy conditional", {
  set.seed(4)
  K <- 20
  S <- 3
  s_tau <- 0.3
  dens <- function(tau) tau^(-K) * exp(-S / (2 * tau^2)) / (1 + (tau / s_tau)^2)
  grid <- seq(1e-3, 5, length.out = 20000)
  w <- dens(grid)
  mean_ref <- sum(grid * w) / sum(w)

  tau <- 0.5
  draws <- numeric(20000)
  for (i in seq_along(draws)) {
    tau <- slice_sample_tau(tau, S, K, s_tau)
    draws[i] <- tau
  }
  expect_equal(mean(draws), mean_ref, tolerance = 0.02)
})

test_that("update_phi leaves the exact phi conditional invariant", {
  set.seed(5)
  e <- as.numeric(arima.sim(list(ar = 0.6), n = 60))
  sigma <- 1.0
  log_cond <- function(phi) {
    u <- c(sqrt(1 - phi^2) * e[1], e[-1] - phi * e[-60])
    0.5 * log(1 - phi^2) - sum(u^2) / (2 * sigma^2) + dbeta((phi + 1) / 2, 10, 2, log = TRUE)
  }
  grid <- seq(-0.95, 0.95, length.out = 5000)
  lw <- sapply(grid, log_cond)
  w <- exp(lw - max(lw))
  mean_ref <- sum(grid * w) / sum(w)

  phi <- 0
  draws <- numeric(20000)
  for (i in seq_along(draws)) {
    phi <- update_phi(e, phi, sigma)$phi
    draws[i] <- phi
  }
  expect_equal(mean(draws), mean_ref, tolerance = 0.01)
  expect_true(all(abs(draws) <= 0.95))
})

test_that("phi update moves when the conditional mass sits at the 0.95 bound", {
  set.seed(6)
  e <- cumsum(rnorm(200)) # random-walk residuals push phi to the bound
  phis <- numeric(200)
  phi <- 0.8
  for (i in seq_along(phis)) {
    phi <- update_phi(e, phi, 1.0)$phi
    phis[i] <- phi
  }
  expect_gt(sd(phis), 0)
  expect_true(all(phis <= 0.95))
})

test_that("ordinary horseshoe option keeps eta fixed at 1 in the variance", {
  set.seed(7)
  n_omega <- 30
  omega <- c(rnorm(n_omega - 1, 0, 0.01), 3)
  res <- update_hsplus_scales(
    omega = omega, sigma = 1, tau_omega = 0.1,
    lambda_omega = rep(1, n_omega), eta_omega = rep(5, n_omega), kappa = rep(1, n_omega),
    aux_lambda = rep(1, n_omega), aux_eta = rep(1, n_omega), s_tau = 0.2,
    horseshoe_plus = FALSE
  )
  expect_equal(res$eta_omega, rep(1, n_omega))
  expect_equal(res$v_omega_tilde, (res$tau_omega * res$lambda_omega)^2, tolerance = 1e-12)
})

test_that("dynamic outlier update samples rho and applies the CAPA offset", {
  set.seed(8)
  n <- 200
  psi <- rnorm(n, 0, 0.05)
  psi[c(50, 120)] <- c(6, -5)
  a <- numeric(n)
  a[c(50, 120)] <- 1
  st <- NULL
  rhos <- numeric(50)
  for (i in 1:50) {
    res <- update_outlier_scales(psi, dynamic = TRUE, delta = 2, a_vec = a, dsp_state = st)
    st <- res$dsp_state
    rhos[i] <- res$rho
  }
  expect_gt(sd(rhos), 0)
  expect_equal(length(res$nu2), n)
  expect_true(all(is.finite(res$nu2) & res$nu2 > 0))
})

test_that("collapsed_tau_update targets the lambda-marginal conditional", {
  set.seed(9)
  x <- c(rnorm(15, 0, 0.02), 2.5)
  a <- runif(16, 0.5, 2)
  nu <- runif(16, 0.2, 3)
  s_tau <- 0.1
  dens <- function(tau) {
    sc <- tau * a * sqrt(2 / nu)
    exp(sum(dcauchy(x, 0, sc, log = TRUE))) / (1 + (tau / s_tau)^2)
  }
  grid <- exp(seq(log(1e-4), log(10), length.out = 20000))
  w <- sapply(grid, dens) * c(diff(grid), 0)
  mean_ref <- sum(grid * w) / sum(w)

  tau <- 0.05
  draws <- numeric(20000)
  for (i in seq_along(draws)) {
    tau <- collapsed_tau_update(tau, x, a, nu, s_tau)
    draws[i] <- tau
  }
  expect_equal(mean(draws), mean_ref, tolerance = 0.02)
})

test_that("factor_mu_psi log marginal likelihood matches the dense Gaussian marginal", {
  set.seed(11)
  n <- 6
  y <- rnorm(n)
  bs <- init_banded_structure(n)
  logml_dense <- function(v, nu2, phi = 0.5, sigma = 0.9) {
    # y = mu + psi + eps; mu = mu_1 + cumsum(omega); Cov(y) built explicitly
    L <- lower.tri(diag(n), diag = TRUE) * 1
    C_mu <- MU1_PRIOR_VAR * matrix(1, n, n) + L[, -1, drop = FALSE] %*% diag(v, n - 1) %*% t(L[, -1, drop = FALSE])
    C_eps <- sigma^2 / (1 - phi^2) * phi^abs(outer(1:n, 1:n, "-"))
    S <- C_mu + diag(nu2) + C_eps
    -0.5 * determinant(S)$modulus - 0.5 * sum(y * solve(S, y))
  }
  v1 <- runif(n - 1, 0.1, 1); nu1 <- runif(n, 0.1, 1)
  v2 <- v1 * 3; nu2 <- nu1 / 2
  d_dense <- logml_dense(v2, nu2) - logml_dense(v1, nu1)
  d_band <- factor_mu_psi(y, 0.5, 0.9, v2, nu2, bs)$logml - factor_mu_psi(y, 0.5, 0.9, v1, nu1, bs)$logml
  expect_equal(as.numeric(d_band), as.numeric(d_dense), tolerance = 1e-5)
})

test_that("rho_noncentred_update preserves innovations and targets the right conditional", {
  set.seed(12)
  n <- 30
  e <- rnorm(n, 0, 0.8)
  m_h <- -1
  # Conditional of rho given fixed innovations e: p(rho) * prod N(psi; 0, exp(m_h + h(rho)))
  psi <- rnorm(n, 0, exp((m_h + as.numeric(stats::filter(e, 0.6, method = "recursive"))) / 2))
  log_cond <- function(r) {
    h <- as.numeric(stats::filter(e, r, method = "recursive"))
    sum(dnorm(psi, 0, exp((m_h + h) / 2), log = TRUE)) + dbeta((r + 1) / 2, 10, 2, log = TRUE)
  }
  grid <- seq(-0.999, 0.999, length.out = 4000)
  lw <- sapply(grid, log_cond)
  w <- exp(lw - max(lw))
  mean_ref <- sum(grid * w) / sum(w)

  rho <- 0.2
  ht <- m_h + as.numeric(stats::filter(e, rho, method = "recursive"))
  draws <- numeric(15000)
  for (i in seq_along(draws)) {
    out <- rho_noncentred_update(psi, ht, m_h, rho)
    rho <- out$rho
    ht <- out$ht
    draws[i] <- rho
  }
  h_til <- ht - m_h
  expect_equal(c(h_til[1], h_til[-1] - rho * h_til[-n]), e, tolerance = 1e-8) # innovations unchanged
  expect_equal(mean(draws), mean_ref, tolerance = 0.02)
})

test_that("C++ banded factorisation matches the dense precision and the spam reference", {
  set.seed(321)
  n <- 30
  y <- cumsum(rnorm(n))
  v <- exp(rnorm(n - 1, -2, 2))
  nu2 <- exp(rnorm(n, -3, 2))
  bs <- init_banded_structure(n)
  fac <- factor_mu_psi(y, 0.6, 1.1, v, nu2, bs)
  m <- 2 * n
  L <- matrix(0, m, m)
  for (i in seq_len(m)) for (d in 0:3) if (i - d >= 1) L[i, i - d] <- fac$L[i, d + 1]
  Q_dense <- as.matrix(build_mu_psi_system(y, 0.6, 1.1, v, nu2, bs)$Q)
  expect_equal(L %*% t(L), Q_dense, tolerance = 1e-10)

  ref <- factor_mu_psi_spam(y, 0.6, 1.1, v, nu2, bs)
  expect_equal(fac$logml, ref$logml, tolerance = 1e-8)
  mean_cpp <- band_backsolve_cpp(fac$L, fac$w)
  expect_equal(mean_cpp, as.numeric(solve(Q_dense, build_mu_psi_system(y, 0.6, 1.1, v, nu2, bs)$b)), tolerance = 1e-8)
})

test_that("tridiagonal sampler and AR(1) filter match dense / stats::filter references", {
  set.seed(5)
  n <- 25
  qd <- runif(n, 2, 4)
  qo <- runif(n - 1, -0.9, 0.9)
  b <- rnorm(n)
  Q <- diag(qd)
  Q[cbind(1:(n - 1), 2:n)] <- qo
  Q[cbind(2:n, 1:(n - 1))] <- qo
  expect_equal(tridiag_sample_cpp(qd, qo, b, rep(0, n)), as.numeric(solve(Q, b)), tolerance = 1e-10)
  z <- rnorm(n)
  U <- chol(Q)
  expect_equal(tridiag_sample_cpp(qd, qo, b, z), as.numeric(solve(U, solve(t(U), b) + z)), tolerance = 1e-10)
  e <- rnorm(40)
  expect_equal(ar1_filter_cpp(e, 0.7), as.numeric(stats::filter(e, 0.7, method = "recursive")), tolerance = 1e-12)
})

test_that("sample_dsp_fast reproduces dsp::sampleDSP draw for draw", {
  ns <- asNamespace("dsp")
  set.seed(11)
  psi <- matrix(c(rnorm(60, 0, 0.05), 3, -4, rnorm(38, 0, 0.05)), ncol = 1)
  st <- ns$initDHS(psi)
  st$sigma_wt <- NULL
  st$dhs_phi <- 0.5
  for (k in 1:5) {
    set.seed(100 + k)
    a <- ns$sampleDSP(omega = psi, evolParams = st, sigma_e = 0.1, loc = NULL, prior_dhs_phi = c(10, 2), alphaPlusBeta = 1)
    set.seed(100 + k)
    b <- sample_dsp_fast(omega = psi, evolParams = st, sigma_e = 0.1, prior_dhs_phi = c(10, 2), alphaPlusBeta = 1)
    expect_equal(as.numeric(b$ht), as.numeric(a$ht), tolerance = 1e-8)
    expect_equal(b$dhs_phi, a$dhs_phi, tolerance = 1e-8)
    expect_equal(b$dhs_mean, a$dhs_mean, tolerance = 1e-8)
    expect_equal(as.numeric(b$sigma_eta_t), as.numeric(a$sigma_eta_t), tolerance = 1e-8)
    a$sigma_wt <- NULL
    st <- a
  }
})
