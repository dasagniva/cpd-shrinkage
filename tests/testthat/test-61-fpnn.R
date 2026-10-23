# ==============================================================================
# test-61-fpnn.R
# The R port of the JISPS EV-DeCAFS + FPNN pipeline against the Python implementation
# (ev-decafs-serve), stage by stage, on a stored reference series.
# ==============================================================================

library(testthat)

if (!exists("fit_baseline")) {
  source("../../R/00_utils.R")
  source("../../R/60_baselines.R")
}

ref <- jsonlite::fromJSON(test_path("fixtures", "fpnn_python_reference.json"))
y <- ref$y
n <- length(y)

test_that("AR(1) estimates and EVI field match the Python implementation", {
  ar1 <- fpnn_ar1_params(y)
  expect_equal(ar1$phi, ref$ar1$phi, tolerance = 1e-10)
  expect_equal(ar1$sigma_v_sq, ref$ar1$sigma_v_sq, tolerance = 1e-10)
  expect_equal(ar1$sigma_eta_sq, ref$ar1$sigma_eta_sq, tolerance = 1e-10)
  expect_equal(fpnn_evi_field(y), ref$xi, tolerance = 1e-10)
})

test_that("EV-DeCAFS dynamic programme, features, BOCPD, labels and FPNN match Python", {
  ar1 <- fpnn_ar1_params(y)
  dp <- ev_decafs_dp_cpp(y, rep(2 * log(n), n), 1 / (ar1$sigma_eta_sq + 1e-12), 1 / (ar1$sigma_v_sq + 1e-12),
                         ar1$phi, 1000L)
  expect_equal(dp$changepoints, as.integer(ref$cps))
  expect_equal(dp$means, ref$means, tolerance = 1e-10)
  expect_equal(dp$cost, ref$cost, tolerance = 1e-10)

  X <- fpnn_extract_features(y, dp$changepoints, dp$means, L = 5, xi_field = fpnn_evi_field(y))
  expect_equal(unname(X), ref$X, tolerance = 1e-10)

  bocpd <- which(fpnn_bocpd(y, ar1$phi, threshold = 0.3))
  expect_equal(bocpd, as.integer(ref$bocpd))

  labels <- fpnn_make_labels(X, dp$changepoints, bocpd, tol = floor(0.02 * n))
  expect_equal(labels, as.integer(ref$labels))

  model <- fourier_pnn_fit(X, labels, J = 10)
  expect_equal(unname(fourier_pnn_predict_proba(model, X)), ref$proba, tolerance = 1e-10)
})

test_that("SMOTE balances the classes by interpolating within the minority class", {
  set.seed(1)
  X <- rbind(matrix(rnorm(20), 10), matrix(rnorm(6, 5), 3))
  lab <- c(rep(0L, 10), rep(1L, 3))
  bal <- fpnn_smote(X, lab, k = 5)
  expect_equal(as.vector(table(bal$y)), c(10, 10))
  new_pts <- bal$X[-(1:13), , drop = FALSE]
  minority <- X[lab == 1L, ]
  expect_true(all(new_pts[, 1] >= min(minority[, 1]) - 1e-12 & new_pts[, 1] <= max(minority[, 1]) + 1e-12))
})

test_that("baseline_fpnn keeps only sustained Phase-I changes and restores the caller's RNG", {
  set.seed(99)
  before <- .Random.seed
  fit <- baseline_fpnn(y)
  expect_identical(.Random.seed, before)
  expect_equal(fit$method, "fpnn")
  expect_true(all(fit$cp %in% fit$diag$phase1_cp))
  expect_setequal(c(fit$cp, fit$diag$recoiled), fit$diag$phase1_cp)
})
