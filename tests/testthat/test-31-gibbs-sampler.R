# ==============================================================================
# test-31-gibbs-sampler.R
# Unit tests for the production Gibbs sampler (Phase P4).
# ==============================================================================

source("../../R/00_utils.R")
source("../../R/10_dgp.R")
source("../../R/40_scoring.R")
source("../../R/41_bfdr.R")
source("../../R/42_classify.R")
source("../../R/32_gibbs_updates.R")
source("../../R/31_gibbs_sampler.R")

test_that("gibbs_sampler runs and returns a validated cpd_fit object", {
  n <- 50
  sim <- simulate_series(
    n              = n,
    phi            = 0.5,
    sigma          = 1.0,
    noise          = "gaussian",
    k0             = 1,
    jump           = 3.0,
    outlier_regime = "none",
    collective     = 0,
    seed           = 12345
  )

  fit <- gibbs_sampler(
    y         = sim$y,
    n_iter    = 300,
    burn_in   = 100,
    thin      = 1,
    n_chains  = 2,
    dynamic   = FALSE,
    alpha_fdr = 0.05,
    seed      = 12345
  )

  # Check cpd_fit schema
  expect_silent(validate_cpd_fit(fit, n = n))
  expect_s3_class(fit, "cpd_fit")
  expect_equal(fit$method, "proposed")

  # Check scores
  expect_equal(nrow(fit$scores), n)
  expect_true(all(c("t", "s_cp", "s_out") %in% names(fit$scores)))
  expect_true(all(fit$scores$s_cp >= 0 & fit$scores$s_cp <= 1))
  expect_true(all(fit$scores$s_out >= 0 & fit$scores$s_out <= 1))

  # Check diagnostics
  expect_true(!is.null(fit$diag$summary))
  expect_true(fit$diag$phi_hat > -0.95 && fit$diag$phi_hat < 0.95)
  expect_true(fit$diag$sigma_hat > 0)
  expect_true(is.finite(fit$diag$max_rhat))
  expect_equal(length(fit$diag$mu_mean), n)
})

test_that("gibbs_sampler runs with dynamic outliers enabled", {
  n <- 50
  sim <- simulate_series(
    n              = n,
    phi            = 0.4,
    sigma          = 1.0,
    noise          = "gaussian",
    k0             = 0,
    jump           = 0,
    outlier_regime = "isolated",
    p_outlier      = 0.04,
    collective     = 0,
    seed           = 54321
  )

  fit <- gibbs_sampler(
    y         = sim$y,
    n_iter    = 150,
    burn_in   = 50,
    n_chains  = 1,
    dynamic   = TRUE,
    seed      = 54321
  )

  expect_silent(validate_cpd_fit(fit, n = n))
  expect_equal(nrow(fit$scores), n)
})

test_that("persistence_split moves short excursions of mu to psi and keeps shifts", {
  mu <- c(rep(0, 20), rep(3, 20))
  mu[8:9] <- 5      # 2-point excursion: shorter than m = 4
  mu[30] <- -2      # 1-point excursion inside the second level
  psi <- rep(0, 40)
  sp <- persistence_split(mu, psi, m = 4)
  expect_equal(sp$mu, c(rep(0, 20), rep(3, 20)))
  expect_equal(sp$mu + sp$psi, mu + psi)
  expect_equal(sp$psi[8:9], c(5, 5))
  expect_equal(sp$psi[30], -5)
  # Excursions lasting m points or more stay in the level
  mu2 <- c(rep(0, 20), rep(4, 4), rep(0, 16))
  expect_equal(persistence_split(mu2, rep(0, 40), m = 4)$mu, mu2)
  expect_equal(persistence_split(mu, psi, m = 1)$mu, mu)
})

test_that("persistence_split is applied to mu + psi jointly", {
  n <- 60
  # A spike held in mu and a 12-point block held in psi
  mu <- rep(0, n)
  mu[10] <- 6
  psi <- rep(0, n)
  psi[30:41] <- 3
  sp <- persistence_split(mu, psi, m = 4)
  expect_equal(sp$mu + sp$psi, mu + psi)
  expect_equal(sp$mu[10], 0)       # spike moved to the outlier process
  expect_equal(sp$psi[10], 6)
  expect_equal(sp$mu[30:41], rep(3, 12)) # persistent block moved to the level
  expect_equal(sp$psi[30:41], rep(0, 12))
})
