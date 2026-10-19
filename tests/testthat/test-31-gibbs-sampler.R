# ==============================================================================
# test-31-gibbs-sampler.R
# Unit tests for the production Gibbs sampler (Phase P4).
# ==============================================================================

source("../../R/00_utils.R")
source("../../R/10_dgp.R")
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
