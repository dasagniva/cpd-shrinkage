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
  aux_tau <- 1.0

  res <- update_hsplus_scales(
    omega = omega,
    sigma = sigma,
    tau_omega = tau_omega,
    lambda_omega = lambda_omega,
    eta_omega = eta_omega,
    kappa = kappa,
    aux_lambda = aux_lambda,
    aux_eta = aux_eta,
    aux_tau = aux_tau,
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
  sigma <- 1.0

  res <- update_outlier_scales(psi, sigma, dynamic = FALSE)
  expect_equal(length(res$nu2), n)
  expect_true(all(res$nu2 > 0))
  expect_true(all(is.finite(res$nu2)))
  expect_true(res$tau_psi > 0)
})
