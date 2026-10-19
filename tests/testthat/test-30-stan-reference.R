# ==============================================================================
# test-30-stan-reference.R
# Unit tests for the Stan reference model and R interface.
# ==============================================================================

source("../../R/00_utils.R")
source("../../R/10_dgp.R")
source("../../R/40_scoring.R")
source("../../R/41_bfdr.R")
source("../../R/42_classify.R")
source("../../R/30_stan_reference.R")

test_that("Stan reference model file exists and compiles", {
  stan_file <- "../../stan/hsplus_ar1_ref.stan"
  expect_true(file.exists(stan_file))

  mod <- cmdstan_model(stan_file)
  expect_s3_class(mod, "CmdStanModel")
})

test_that("fit_stan_reference returns a valid cpd_fit object on short series", {
  skip_if_not(dir.exists("../../.cmdstan/cmdstan-2.40.0"))

  # Small test series (n = 50) with 1 change point
  sim <- simulate_series(
    n              = 50,
    phi            = 0.5,
    sigma          = 1.0,
    noise          = "gaussian",
    k0             = 1,
    jump           = 3.0,
    outlier_regime = "none",
    collective     = 0,
    seed           = 12345
  )

  fit <- fit_stan_reference(
    y             = sim$y,
    chains        = 1,
    iter_warmup   = 150,
    iter_sampling = 150,
    adapt_delta   = 0.99,
    seed          = 12345,
    stan_file     = "../../stan/hsplus_ar1_ref.stan"
  )

  # Validate cpd_fit schema
  expect_silent(validate_cpd_fit(fit, n = 50))
  expect_s3_class(fit, "cpd_fit")
  expect_equal(fit$method, "stan_ref")

  # Check scores
  expect_equal(nrow(fit$scores), 50)
  expect_true(all(c("t", "s_cp", "s_out") %in% names(fit$scores)))
  expect_true(all(fit$scores$s_cp >= 0 & fit$scores$s_cp <= 1))
  expect_true(all(fit$scores$s_out >= 0 & fit$scores$s_out <= 1))

  # Check diagnostics
  expect_true(!is.null(fit$diag$phi_hat))
  expect_true(fit$diag$phi_hat >= -0.95 && fit$diag$phi_hat <= 0.95)
  expect_true(fit$diag$sigma_hat > 0)
  expect_true(is.finite(fit$diag$max_rhat))
})

test_that("fit_stan_reference handles input validation errors cleanly", {
  expect_error(fit_stan_reference(y = 1:5), "n >= 10")
  expect_error(fit_stan_reference(y = rnorm(20), kappa = rep(1, 10)), "length\\(kappa\\) == n")
})
