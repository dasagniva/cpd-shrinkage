# ==============================================================================
# test-50-copula.R
# Unit tests for Phase P6: Copula Markov Chain Comparator
# ==============================================================================

library(testthat)

if (!exists("fit_copula_comparator")) {
  source("../../R/00_utils.R")
  source("../../R/50_copula.R")
}

test_that("biweight_loss computes expected Tukey bisquare loss", {
  c_val <- 4.685
  e_zero <- 0.0
  expect_equal(biweight_loss(e_zero, c = c_val), 0.0)

  # Inside tuning boundary: e = 1.0
  u <- 1.0 / c_val
  expected_loss <- (c_val^2 / 6.0) * (1.0 - (1.0 - u^2)^3)
  expect_equal(biweight_loss(1.0, c = c_val), expected_loss)

  # Outside tuning boundary: e = 10.0 > 4.685
  expect_equal(biweight_loss(10.0, c = c_val), c_val^2 / 6.0)
})

test_that("fit_robust_location recovers shift and returns residuals", {
  n <- 60
  # Jump at t = 31 from 0 to 5, plus an isolated outlier at t = 15
  y <- c(rep(0.0, 30), rep(5.0, 30))
  y[15] <- 12.0 # Outlier

  fit_loc <- fit_robust_location(y, min_seg_len = 10)
  expect_true(is.list(fit_loc))
  expect_equal(length(fit_loc$mu_hat), n)
  expect_equal(length(fit_loc$residuals), n)

  # The location estimate around 1..10 should be close to 0
  expect_lt(abs(fit_loc$mu_hat[5]), 0.5)
  # The location estimate around 40..50 should be close to 5
  expect_lt(abs(fit_loc$mu_hat[45] - 5.0), 0.5)
  # Outlier residual at 15 should be large
  expect_gt(abs(fit_loc$residuals[15]), 8.0)
})

test_that("simulate_copula_markov generates stationary Markov uniform sequences", {
  cop_g <- rvinecopulib::bicop_dist("gaussian", 0, 0.6)
  n <- 200
  u <- simulate_copula_markov(cop_g, n, seed = 123)

  expect_equal(length(u), n)
  expect_true(all(u > 0 & u < 1))

  # Rank correlation (Kendall tau) should be positive and close to theoretical
  tau_emp <- cor(u[-n], u[-1], method = "kendall")
  # Gaussian with rho = 0.6 has tau = (2/pi) * asin(0.6) approx 0.4097
  expect_gt(tau_emp, 0.25)
  expect_lt(tau_emp, 0.55)
})

test_that("calibrate_copula_threshold returns quantile approximately equal to alpha", {
  cop <- rvinecopulib::bicop_dist("clayton", 0, 2.0)
  thresh <- calibrate_copula_threshold(cop, n = 50, alpha = 0.05, B_boot = 50, seed = 42)

  expect_true(is.numeric(thresh))
  expect_length(thresh, 1)
  expect_gt(thresh, 0.01)
  expect_lt(thresh, 0.15)
})

test_that("fit_copula_comparator returns a validated cpd_fit object", {
  set.seed(999)
  n <- 50
  y <- rnorm(n)
  y[25] <- 6.0 # Add outlier

  fit <- fit_copula_comparator(y, alpha = 0.05, B_boot = 30, seed = 42)

  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n))
  expect_equal(fit$method, "copula_markov")
  expect_true("lambda_L" %in% names(fit$diag))
  expect_true("lambda_U" %in% names(fit$diag))
  expect_true("family" %in% names(fit$diag))
  expect_true("threshold" %in% names(fit$diag))

  # Outlier at 25 should be detected in point anomalies
  expect_true(25L %in% fit$point)
  expect_equal(as.character(fit$class[25]), "point")
})

test_that("copula Markov comparator recovers family and parameters (Gate requirement)", {
  n <- 500

  # 1. Gaussian recovery
  cop_g <- rvinecopulib::bicop_dist("gaussian", 0, 0.6)
  u_g <- simulate_copula_markov(cop_g, n, seed = 101)
  fit_g <- rvinecopulib::bicop(cbind(u_g[-n], u_g[-1]), family_set = c("gaussian", "t", "clayton", "gumbel", "frank"))
  expect_true(fit_g$family %in% c("gaussian", "t"))
  expect_lt(abs(fit_g$parameters[1] - 0.6), 0.15)

  # 2. Gumbel recovery
  cop_gu <- rvinecopulib::bicop_dist("gumbel", 0, 2.0)
  u_gu <- simulate_copula_markov(cop_gu, n, seed = 102)
  fit_gu <- rvinecopulib::bicop(cbind(u_gu[-n], u_gu[-1]), family_set = c("gaussian", "t", "clayton", "gumbel", "frank"))
  expect_equal(fit_gu$family, "gumbel")
  expect_lt(abs(fit_gu$parameters[1] - 2.0), 0.35)
  # Gumbel has upper tail dependence
  td_gu <- rvinecopulib::tail_dep(fit_gu)
  expect_gt(td_gu["upper", "upper"], 0.3)
})

test_that("null false-flag rate is approximately alpha (Gate requirement)", {
  # Under clean AR(1) or null copula series, false-flag rate should be near nominal alpha
  set.seed(2026)
  n <- 100
  alpha <- 0.05
  B <- 10
  flag_counts <- numeric(B)

  for (b in 1:B) {
    # Generate null stationary series: AR(1) phi = 0.5
    eps <- numeric(n)
    eps[1] <- rnorm(1)
    for (t in 2:n) eps[t] <- 0.5 * eps[t - 1] + rnorm(1)

    fit <- fit_copula_comparator(eps, alpha = alpha, B_boot = 50, seed = 100 + b)
    flag_counts[b] <- length(fit$point) / n
  }

  mean_flag_rate <- mean(flag_counts)
  # Empirical false-flag rate should be around 0.05 (+- 0.04 over 10 small series)
  expect_gt(mean_flag_rate, 0.01)
  expect_lt(mean_flag_rate, 0.10)
})
