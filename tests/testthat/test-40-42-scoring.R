# ==============================================================================
# test-40-42-scoring.R
# Unit tests for Phase P5: Scoring, Bayesian FDR, Classification
# ==============================================================================

library(testthat)

# Source modules if running interactively
if (!exists("compute_outlier_score")) {
  source("../../R/00_utils.R")
  source("../../R/40_scoring.R")
  source("../../R/41_bfdr.R")
  source("../../R/42_classify.R")
}

# ------------------------------------------------------------------------------
# 1. Tests for 40_scoring.R
# ------------------------------------------------------------------------------
test_that("compute_outlier_score calculates exact ratios and respects innovation vs marginal variance", {
  M <- 4
  n <- 3

  # nu^2 values: row = iter, col = time
  nu2_mat <- matrix(c(
    0.0, 1.0, 4.0,
    0.0, 2.0, 8.0,
    0.0, 1.0, 9.0,
    0.0, 0.0, 3.0
  ), nrow = M, ncol = n, byrow = TRUE)

  sigma2_draws <- c(1.0, 2.0, 1.0, 1.0)
  phi_draws <- c(0.6, 0.6, 0.6, 0.6) # 1 - phi^2 = 0.64

  # Default: innovation variance
  # Iter 1: [0/(0+1)=0, 1/(1+1)=0.5, 4/(4+1)=0.8]
  # Iter 2: [0/(0+2)=0, 2/(2+2)=0.5, 8/(8+2)=0.8]
  # Iter 3: [0/(0+1)=0, 1/(1+1)=0.5, 9/(9+1)=0.9]
  # Iter 4: [0/(0+1)=0, 0/(0+1)=0.0, 3/(3+1)=0.75]
  # Means: col 1 = 0, col 2 = (0.5 + 0.5 + 0.5 + 0) / 4 = 0.375, col 3 = (0.8 + 0.8 + 0.9 + 0.75) / 4 = 0.8125
  o_t <- compute_outlier_score(nu2_mat, sigma2_draws, use_marginal_variance = FALSE)

  expect_equal(length(o_t), n)
  expect_equal(o_t[1], 0.0)
  expect_equal(o_t[2], 0.375)
  expect_equal(o_t[3], 0.8125)

  # Marginal variance ablation: denom = sigma^2 / 0.64
  # Iter 1: denom = 1 / 0.64 = 1.5625 -> col 2: 1 / (1 + 1.5625) = 1 / 2.5625 approx 0.39024
  o_t_marg <- compute_outlier_score(nu2_mat, sigma2_draws, phi_draws = phi_draws, use_marginal_variance = TRUE)
  expect_equal(length(o_t_marg), n)
  expect_true(all(o_t_marg >= 0 & o_t_marg <= 1))
  # With denom larger (since 1 - phi^2 < 1), outlier score using marginal variance is strictly smaller
  expect_true(o_t_marg[2] < o_t[2])
  expect_true(o_t_marg[3] < o_t[3])
})

test_that("compute_tail_prob_omega and compute_tail_prob_psi evaluate exceedance correctly", {
  M <- 5
  n <- 4

  # omega draws for n-1=3 increments
  w_mat <- matrix(c(
    0.1, 3.5, 0.2,
    0.2, 4.0, 0.1,
    0.0, 2.5, 0.5,
    0.1, 1.0, 0.1,
    0.3, 0.5, 0.2
  ), nrow = M, ncol = 3, byrow = TRUE)

  sigma_draws <- rep(1.0, M)
  phi_draws <- rep(0.0, M) # sd_marg = 1.0, d_omega = 1.5 * 1.0 = 1.5

  # Col 1: all <= 1.5 -> prob 0.0
  # Col 2: 3.5 (>1.5), 4.0 (>1.5), 2.5 (>1.5), 1.0 (<=1.5), 0.5 (<=1.5) -> 3/5 = 0.60
  # Col 3: all <= 1.5 -> prob 0.0
  p_w <- compute_tail_prob_omega(w_mat, sigma_draws, phi_draws, c_omega = 1.5)

  expect_equal(length(p_w), 3)
  expect_equal(p_w[1], 0.0)
  expect_equal(p_w[2], 0.6)
  expect_equal(p_w[3], 0.0)

  # Check psi tail exceedance
  psi_mat <- matrix(c(
    0.0, 0.0, 5.0, 0.0,
    0.1, 0.0, 4.5, 0.2,
    0.0, 0.0, 6.0, 0.1,
    0.0, 0.1, 5.5, 0.0,
    0.2, 0.0, 7.0, 0.0
  ), nrow = M, ncol = 4, byrow = TRUE)

  p_psi <- compute_tail_prob_psi(psi_mat, sigma_draws, phi_draws, c_psi = 1.5)
  expect_equal(length(p_psi), 4)
  expect_equal(p_psi[3], 1.0)
  expect_equal(p_psi[1], 0.0)
})

# ------------------------------------------------------------------------------
# 2. Tests for 41_bfdr.R
# ------------------------------------------------------------------------------
test_that("bayesian_fdr matches known hand-calculated threshold and discoveries", {
  # Scenario: 5 hypotheses
  p <- c(0.70, 0.99, 0.91, 0.98, 0.94) # unsorted
  # Sorted descending:
  # 1: 0.99 -> cum FDR = (1 - 0.99) / 1 = 0.010 <= 0.05
  # 2: 0.98 -> cum FDR = (0.01 + 0.02) / 2 = 0.015 <= 0.05
  # 3: 0.94 -> cum FDR = (0.03 + 0.06) / 3 = 0.030 <= 0.05
  # 4: 0.91 -> cum FDR = (0.09 + 0.09) / 4 = 0.045 <= 0.05
  # 5: 0.70 -> cum FDR = (0.18 + 0.30) / 5 = 0.096 > 0.05

  res <- bayesian_fdr(p, alpha = 0.05)
  expect_s3_class(res, "bfdr_result")
  expect_equal(res$k, 4L)
  expect_equal(res$threshold, 0.91)
  expect_equal(res$est_fdr, 0.045)
  # Preserves time order: indices of 0.99 (2), 0.98 (4), 0.94 (5), 0.91 (3)
  expect_equal(res$flagged, c(2L, 3L, 4L, 5L))

  # If alpha is 0.02:
  # k=1: 0.01 <= 0.02
  # k=2: 0.015 <= 0.02
  # k=3: 0.03 > 0.02
  res_strict <- bayesian_fdr(p, alpha = 0.02)
  expect_equal(res_strict$k, 2L)
  expect_equal(res_strict$flagged, c(2L, 4L))
})

test_that("bayesian_fdr handles boundary cases: empty, null, and certainty", {
  # Empty input
  res_empty <- bayesian_fdr(numeric(0), alpha = 0.05)
  expect_equal(res_empty$k, 0L)
  expect_equal(length(res_empty$flagged), 0L)

  # All low probabilities (no discovery)
  p_low <- c(0.1, 0.2, 0.05, 0.15) # max p = 0.2, 1 - 0.2 = 0.8 > 0.05
  res_low <- bayesian_fdr(p_low, alpha = 0.05)
  expect_equal(res_low$k, 0L)
  expect_equal(length(res_low$flagged), 0L)
  expect_equal(res_low$est_fdr, 0.0)

  # All 1.0 (certain discoveries)
  p_high <- rep(1.0, 6)
  res_high <- bayesian_fdr(p_high, alpha = 0.05)
  expect_equal(res_high$k, 6L)
  expect_equal(res_high$flagged, 1:6)
  expect_equal(res_high$est_fdr, 0.0)
})

# ------------------------------------------------------------------------------
# 3. Tests for 42_classify.R
# ------------------------------------------------------------------------------
test_that("classify_four_class labels sustained shifts correctly", {
  n <- 50
  mu_hat <- c(rep(0.0, 20), rep(3.0, 30)) # Step jump at t = 21

  flagged_w <- c(21L)
  flagged_p <- integer(0)

  res <- classify_four_class(n, flagged_w, flagged_p, mu_hat)

  expect_equal(res$cp, 21L)
  expect_equal(nrow(res$collective), 0L)
  expect_equal(length(res$point), 0L)
  expect_equal(as.character(res$class[21]), "shift")
  expect_equal(sum(res$class == "shift"), 1L)
  expect_equal(sum(res$class == "nominal"), 49L)
})

test_that("classify_four_class labels collective anomalies correctly via return condition", {
  n <- 60
  # Mean shifts up at t=15 from 0 to 4, then returns at t=25 back to 0.1
  mu_hat <- numeric(n)
  mu_hat[1:14] <- 0.0
  mu_hat[15:24] <- 4.0
  mu_hat[25:60] <- 0.1

  # Increments: t1=15 (jump +4), t2=25 (jump -3.9)
  flagged_w <- c(15L, 25L)
  flagged_p <- integer(0)

  # sep = 25 - 15 = 10 <= L_max (floor(0.1*60)=6, so set L_max = 15)
  res <- classify_four_class(n, flagged_w, flagged_p, mu_hat, L_max = 15, r = 0.3)

  expect_equal(nrow(res$collective), 1L)
  expect_equal(res$collective$start[1], 15L)
  expect_equal(res$collective$end[1], 24L) # omega_25 is the return step; excursion is 15..24
  expect_equal(length(res$cp), 0L) # Neither 15 nor 25 should be sustained shift!
  expect_equal(sum(res$class == "collective"), 10L) # 15..24
  expect_equal(as.character(res$class[15]), "collective")
  expect_equal(as.character(res$class[24]), "collective")
  expect_equal(as.character(res$class[25]), "nominal")
  expect_equal(as.character(res$class[14]), "nominal")
})

test_that("classify_four_class uses CAPA overlap alternative when return condition fails", {
  n <- 80
  # Jump at 20 (+3), return at 35 to 1.5 (return diff = 1.5 > 0.3*3=0.9, so return cond fails)
  mu_hat <- numeric(n)
  mu_hat[1:19] <- 0.0
  mu_hat[20:34] <- 3.0
  mu_hat[35:80] <- 1.5

  flagged_w <- c(20L, 35L)
  flagged_p <- integer(0)

  # Without CAPA: return condition fails, so both become shifts
  res_no_capa <- classify_four_class(n, flagged_w, flagged_p, mu_hat, L_max = 20, r = 0.3)
  expect_equal(nrow(res_no_capa$collective), 0L)
  expect_equal(res_no_capa$cp, c(20L, 35L))

  # With CAPA collective segment [20, 35] having 100% overlap
  capa_df <- data.frame(start = 20L, end = 35L)
  res_capa <- classify_four_class(n, flagged_w, flagged_p, mu_hat, capa_collective = capa_df, L_max = 20, r = 0.3)
  expect_equal(nrow(res_capa$collective), 1L)
  expect_equal(res_capa$collective$start[1], 20L)
  expect_equal(res_capa$collective$end[1], 34L)
  expect_equal(length(res_capa$cp), 0L)
  expect_equal(as.character(res_capa$class[20]), "collective")
})

test_that("classify_four_class enforces strict sequential precedence (first match wins)", {
  n <- 40
  mu_hat <- rep(0.0, n)
  mu_hat[10:18] <- 2.5
  mu_hat[19:40] <- 0.05

  # Increments at 10 and 19 (collective)
  flagged_w <- c(10L, 19L, 30L) # 30 is a shift
  # Outliers: one at 14 (inside collective), one at 30 (coincides with shift), one at 5 (isolated)
  flagged_p <- c(5L, 14L, 30L)

  res <- classify_four_class(n, flagged_w, flagged_p, mu_hat, L_max = 12)

  # 1. Collective: 10..18 (omega_19 is the return step)
  expect_equal(nrow(res$collective), 1L)
  expect_equal(res$collective$start[1], 10L)
  expect_equal(res$collective$end[1], 18L)
  expect_equal(as.character(res$class[14]), "collective") # Outlier 14 absorbed into collective

  # 2. Shift: 30
  expect_equal(res$cp, 30L)
  expect_equal(as.character(res$class[30]), "shift") # Outlier 30 absorbed into shift (rule 2 beats rule 3)

  # 3. Point: 5
  expect_equal(res$point, 5L)
  expect_equal(as.character(res$class[5]), "point")

  # 4. Nominal: all other points
  expect_equal(sum(res$class == "nominal"), 40 - 9 - 1 - 1)
})

# ------------------------------------------------------------------------------
# 4. Dedicated Identifiability Test: Lone spike labelled 'point', NOT shifts
# ------------------------------------------------------------------------------
test_that("lone spike is classified as point anomaly and NOT two opposite-signed shifts", {
  n <- 50
  mu_hat <- rep(0.0, n) # Clean baseline

  # Lone spike at t = 25
  # Under an identifiability failure, increments at 25 and 26 would be flagged
  # But here, only the outlier at 25 is flagged:
  flagged_w <- integer(0)
  flagged_p <- c(25L)

  res <- classify_four_class(n, flagged_w, flagged_p, mu_hat)

  expect_equal(res$point, 25L)
  expect_equal(length(res$cp), 0L)
  expect_equal(nrow(res$collective), 0L)
  expect_equal(as.character(res$class[25]), "point")

  # Even IF consecutive increments t=25 and t=26 were flagged by mistake:
  # sep = 26 - 25 = 1 < 2, so it CANNOT be a collective anomaly!
  flagged_w_mistake <- c(25L, 26L)
  res_mistake <- classify_four_class(n, flagged_w_mistake, flagged_p, mu_hat)
  expect_equal(nrow(res_mistake$collective), 0L) # sep = 1 cannot form collective!
})

test_that("classify_and_build_fit builds fully valid cpd_fit object", {
  n <- 30
  y <- rnorm(n)
  mu_hat <- rep(0.0, n)
  s_out <- runif(n, 0, 0.1)
  s_out[12] <- 0.95

  p_omega <- rep(0.01, n - 1L)
  p_psi <- rep(0.01, n)
  p_psi[12] <- 0.99

  fit <- classify_and_build_fit(
    y = y,
    mu_hat = mu_hat,
    s_out = s_out,
    p_omega = p_omega,
    p_psi = p_psi,
    alpha_fdr = 0.05,
    method = "test_method"
  )

  expect_s3_class(fit, "cpd_fit")
  expect_silent(validate_cpd_fit(fit, n = n))
  expect_equal(fit$point, 12L)
  expect_equal(length(fit$cp), 0L)
  expect_equal(as.character(fit$class[12]), "point")
})

test_that("windowed shift probabilities collapse each flagged run to one change", {
  n <- 100
  y <- c(rep(0, 50), rep(3, 50))
  mu_hat <- y
  p_win <- rep(0.0, n - 1)
  p_win[46:54] <- 0.99 # run of flagged times around the change at t = 51
  strength <- rep(0.01, n - 1)
  strength[50] <- 3 # E|omega| peaks at omega_51 = mu_51 - mu_50
  fit <- classify_and_build_fit(y, mu_hat, s_out = rep(0, n), p_omega = p_win, p_psi = rep(0, n),
                                alpha_fdr = 0.05, cp_window = 5L, omega_strength = strength)
  expect_equal(fit$cp, 51L)
  expect_equal(sum(fit$class == "shift"), 1L)
})

test_that("flagged runs with short gaps collapse to a single change", {
  n <- 100
  y <- c(rep(0, 50), rep(3, 50))
  p_win <- rep(0.0, n - 1)
  p_win[c(46:49, 51:54)] <- 0.99 # one unflagged time inside the run
  strength <- rep(0.01, n - 1)
  strength[50] <- 3
  fit <- classify_and_build_fit(y, y, s_out = rep(0, n), p_omega = p_win, p_psi = rep(0, n),
                                alpha_fdr = 0.05, cp_window = 5L, omega_strength = strength)
  expect_length(fit$cp, 1L)
})

test_that("nested or interleaved collective pairs give one merged segment", {
  n <- 60
  # Levels: 0 -> 2 at t = 10 -> 4 at t = 20 -> 2 at t = 30 -> 0 at t = 40.
  # Rule 1 pairs (10, 40) and then (20, 30), whose segments 10-39 and 20-29 overlap.
  mu_hat <- rep(0, n)
  mu_hat[10:19] <- 2
  mu_hat[20:29] <- 4
  mu_hat[30:39] <- 2
  res <- classify_four_class(n, c(10L, 20L, 30L, 40L), integer(0), mu_hat, L_max = 50, r = 0.3)
  expect_equal(nrow(res$collective), 1L)
  expect_equal(res$collective$start, 10L)
  expect_equal(res$collective$end, 39L)
  expect_equal(sum(res$class == "collective"), 30L)
})
