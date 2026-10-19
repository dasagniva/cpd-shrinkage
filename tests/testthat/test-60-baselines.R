# ==============================================================================
# test-60-baselines.R
# Unit tests for Phase P7: Baseline Wrappers and Common Output Schema
# ==============================================================================

library(testthat)

if (!exists("fit_baseline")) {
  source("../../R/00_utils.R")
  source("../../R/10_dgp.R")
  source("../../R/21_frontend_capa.R")
  source("../../R/32_gibbs_updates.R")
  source("../../R/40_scoring.R")
  source("../../R/41_bfdr.R")
  source("../../R/42_classify.R")
  source("../../R/50_copula.R")
  source("../../R/60_baselines.R")
}

# Test series with a change point at 31 and an outlier at 15
set.seed(42)
n_test <- 60
y_test <- c(rnorm(30, 0, 1), rnorm(30, 3, 1))
y_test[15] <- 8.0

test_that("baseline_decafs returns a validated cpd_fit", {
  fit <- baseline_decafs(y_test)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "decafs")
  expect_true(length(fit$cp) >= 1)
})

test_that("baseline_nsp returns a validated cpd_fit with intervals", {
  fit <- baseline_nsp(y_test)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "nsp")
  expect_true(is.data.frame(fit$cp_int))
})

test_that("baseline_capa returns a validated cpd_fit with point anomalies", {
  fit <- baseline_capa(y_test)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "capa")
  expect_true(15L %in% fit$point)
})

test_that("baseline_pelt returns a validated cpd_fit", {
  fit <- baseline_pelt(y_test)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "pelt")
})

test_that("baseline_robseg returns a validated cpd_fit", {
  fit <- baseline_robseg(y_test)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "robseg")
})

test_that("baseline_hsmuce returns a validated cpd_fit", {
  fit <- baseline_hsmuce(y_test)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "hsmuce")
})

test_that("baseline_bcp returns a validated cpd_fit", {
  fit <- baseline_bcp(y_test, burnin = 50, mcmc = 100)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "bcp")
})

test_that("baseline_envcpt returns a validated cpd_fit", {
  fit <- baseline_envcpt(y_test)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "envcpt")
})

test_that("baseline_ch4 returns a validated cpd_fit with CP and point anomalies", {
  fit <- baseline_ch4(y_test, n_iter = 200, burn_in = 100)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "ch4")
  expect_true(length(fit$cp) >= 1)
  expect_equal(length(fit$scores$s_out), n_test)
})

test_that("baseline_fpnn returns a validated cpd_fit", {
  fit <- baseline_fpnn(y_test, window_len = 10, J = 6)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "fpnn_placeholder")
})

test_that("baseline_abco returns a validated cpd_fit", {
  fit <- baseline_abco(y_test, nsave = 100, nburn = 100)
  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n_test))
  expect_equal(fit$method, "abco")
})

test_that("fit_baseline master dispatcher works across methods", {
  methods <- c("decafs", "nsp", "capa", "pelt", "robseg", "hsmuce", "fpnn")
  labels <- c(fpnn = "fpnn_placeholder")
  for (m in methods) {
    fit_m <- fit_baseline(m, y_test)
    expect_s3_class(fit_m, "cpd_fit")
    expect_equal(fit_m$method, if (m %in% names(labels)) labels[[m]] else m)
    expect_silent(validate_cpd_fit(fit_m, n = n_test))
  }
})
