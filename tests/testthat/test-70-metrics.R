# ==============================================================================
# test-70-metrics.R
# Unit tests for Section 6 Evaluation Metrics (70_metrics.R)
# ==============================================================================

library(testthat)

if (!exists("evaluate_all_metrics")) {
  source("../../R/00_utils.R")
  source("../../R/10_dgp.R")
  source("../../R/70_metrics.R")
}

test_that("evaluate_change_points calculates exact precision, recall, F1, Hausdorff, and covering", {
  n <- 100
  truth_cp <- c(25L, 60L)
  # Estimated: 26 (within margin of 25), 80 (false positive)
  est_cp   <- c(26L, 80L)

  eval_res <- evaluate_change_points(truth_cp, est_cp, n, base_margin = 5)

  # Margin for n <= 1000 is 5 (spec Section 6); 2000 points give 10
  expect_equal(eval_res$margin, 5L)
  expect_equal(evaluate_change_points(truth_cp, est_cp, 2000, base_margin = 5)$margin, 10L)
  # Est 26 matches 25 (|26-25|=1 <= 5). Est 80 matches none. Prec = 1/2 = 0.5
  expect_equal(eval_res$precision, 0.5)
  # True 25 matches 26. True 60 matches none. Rec = 1/2 = 0.5
  expect_equal(eval_res$recall, 0.5)
  expect_equal(eval_res$f1, 0.5)
  expect_equal(eval_res$count_error, 0L)

  # Hausdorff distance:
  # min dist for 25 is 1 (to 26). min dist for 60 is 20 (to 80). max is 20.
  # min dist for 26 is 1. min dist for 80 is 20 (to 60). max is 20.
  expect_equal(eval_res$hausdorff, 20.0)
  # Covering metric should be in [0, 1]
  expect_gt(eval_res$covering, 0.5)
  expect_lt(eval_res$covering, 1.0)
})

test_that("evaluate_nsp_intervals computes interval coverage correctly", {
  truth_cp <- c(30L, 70L)
  # 2 intervals: [25, 35] (covers 30), [80, 90] (covers nothing)
  ints <- data.frame(start = c(25L, 80L), end = c(35L, 90L))

  eval_nsp <- evaluate_nsp_intervals(truth_cp, ints, alpha_nsp = 0.05)
  expect_equal(eval_nsp$n_intervals, 2L)
  expect_equal(eval_nsp$coverage, 0.5)
  expect_equal(eval_nsp$median_width, 11.0)
})

test_that("evaluate_point_anomalies computes exact metrics and PR-AUC", {
  truth_pt <- c(10L, 25L, 40L)
  est_pt   <- c(10L, 25L, 50L) # 2 TP, 1 FP, 1 FN

  n <- 60
  o_t <- rep(0.01, n)
  o_t[truth_pt] <- 0.95
  o_t[50] <- 0.80

  eval_pt <- evaluate_point_anomalies(truth_pt, est_pt, o_t = o_t, n = n)
  expect_equal(eval_pt$precision, 2 / 3)
  expect_equal(eval_pt$recall, 2 / 3)
  expect_equal(eval_pt$f1, 2 / 3)
  expect_equal(eval_pt$realised_fdr, 1 / 3)
  expect_gt(eval_pt$pr_auc, 0.7)
})

test_that("evaluate_collective_anomalies matches segments based on Jaccard overlap", {
  # True segment: [20, 35]
  truth_coll <- data.frame(start = 20L, end = 35L)
  # Est segment: [21, 34] (inter = 14, union = 16, jacc = 14/16 = 0.875 >= 0.5)
  est_coll   <- data.frame(start = 21L, end = 34L)

  eval_coll <- evaluate_collective_anomalies(truth_coll, est_coll, min_overlap = 0.5)
  expect_equal(eval_coll$precision, 1.0)
  expect_equal(eval_coll$recall, 1.0)
  expect_equal(eval_coll$f1, 1.0)
  expect_equal(eval_coll$start_error, 1.0)
  expect_equal(eval_coll$end_error, 1.0)
})

test_that("evaluate_four_class computes confusion matrix and Macro-F1", {
  truth_class <- factor(c("nominal", "shift", "collective", "point"), levels = CPD_CLASS_LEVELS)
  est_class   <- factor(c("nominal", "shift", "collective", "nominal"), levels = CPD_CLASS_LEVELS)

  eval_clf <- evaluate_four_class(truth_class, est_class)
  expect_equal(nrow(eval_clf$confusion_matrix), 4L)
  expect_equal(ncol(eval_clf$confusion_matrix), 4L)

  # Class "point": TP = 0, FP = 0, FN = 1 -> F1 = 0
  expect_equal(eval_clf$per_class$point$f1, 0.0)
  # Macro-F1 is average of 4 classes
  expect_gt(eval_clf$macro_f1, 0.5)
  expect_lt(eval_clf$macro_f1, 1.0)
})

test_that("evaluate_all_metrics runs end-to-end on simulate_series and a cpd_fit", {
  sim <- simulate_series(
    n              = 60,
    phi            = 0.3,
    sigma          = 1.0,
    k0             = 1,
    jump           = 3.0,
    outlier_regime = "isolated",
    p_outlier      = 0.02,
    collective     = 0,
    seed           = 123
  )

  # Create a dummy valid cpd_fit matching truth
  dummy_fit <- new_cpd_fit(
    method = "test_method",
    cp     = sim$truth$cp,
    point  = sim$truth$point,
    class  = sim$truth$class
  )

  metrics <- evaluate_all_metrics(sim$truth, dummy_fit)
  expect_equal(metrics$method, "test_method")
  expect_equal(metrics$cp$f1, 1.0)
  expect_equal(metrics$four_class$macro_f1, 1.0)
})
