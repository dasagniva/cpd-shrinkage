# ==============================================================================
# test-80-sim-design.R
# Unit tests for Phase P8 Simulation Design (80_sim_design.R)
# ==============================================================================

library(testthat)

if (!exists("create_simulation_grid")) {
  source("../../R/00_utils.R")
  source("../../R/10_dgp.R")
  source("../../R/20_frontend_nsp.R")
  source("../../R/21_frontend_capa.R")
  source("../../R/22_prior_map.R")
  source("../../R/32_gibbs_updates.R")
  source("../../R/40_scoring.R")
  source("../../R/41_bfdr.R")
  source("../../R/42_classify.R")
  source("../../R/31_gibbs_sampler.R")
  source("../../R/50_copula.R")
  source("../../R/60_baselines.R")
  source("../../R/70_metrics.R")
  source("../../R/80_sim_design.R")
}

test_that("create_simulation_grid creates 17 valid scenario cells", {
  grid <- create_simulation_grid()
  expect_true(is.data.frame(grid))
  expect_equal(nrow(grid), 17L)

  # Check required columns
  expected_cols <- c("cell_id", "factor_name", "n", "phi", "noise", "jump", "outlier_regime", "p_outlier", "collective")
  expect_true(all(expected_cols %in% names(grid)))

  # Base cell
  expect_equal(grid$cell_id[1], "base")
  expect_equal(grid$n[1], 1000L)
  expect_equal(grid$phi[1], 0.6)
  expect_equal(grid$jump[1], 2.0)

  # Stress cells
  expect_true("stress_phi0.9_clust_t3" %in% grid$cell_id)
  expect_true("stress_n500_jump1_iso5" %in% grid$cell_id)
})

test_that("run_method_on_replicate executes and evaluates metrics on simulated series", {
  sim <- simulate_series(
    n              = 60,
    phi            = 0.3,
    sigma          = 1.0,
    k0             = 1,
    jump           = 3.0,
    outlier_regime = "isolated",
    p_outlier      = 0.02,
    collective     = 0,
    seed           = 42
  )

  # Run pelt
  res_pelt <- run_method_on_replicate(sim, "pelt")
  expect_equal(res_pelt$method, "pelt")
  expect_true(is.numeric(res_pelt$runtime_sec))
  expect_true("f1" %in% names(res_pelt$cp))

  # Run decafs
  res_decafs <- run_method_on_replicate(sim, "decafs")
  expect_equal(res_decafs$method, "decafs")
  expect_true("macro_f1" %in% names(res_decafs$four_class))
})

test_that("run_pilot_simulation runs on miniature grid and returns tidy results", {
  # Miniature 2-cell grid for testing
  mini_grid <- data.frame(
    cell_id        = c("mini_1", "mini_2"),
    factor_name    = c("mini", "mini"),
    n              = c(50L, 50L),
    phi            = c(0.2, 0.5),
    noise          = c("gaussian", "gaussian"),
    jump           = c(2.0, 3.0),
    outlier_regime = c("none", "isolated"),
    p_outlier      = c(0.0, 0.02),
    collective     = c(0L, 0L),
    stringsAsFactors = FALSE
  )

  pilot_res <- run_pilot_simulation(
    grid      = mini_grid,
    n_reps    = 2,
    methods   = c("pelt", "decafs"),
    seed_base = 9000
  )

  expect_true(is.list(pilot_res))
  expect_true(is.data.frame(pilot_res$data))
  expect_equal(nrow(pilot_res$data), 2 * 2 * 2) # 2 cells x 2 reps x 2 methods = 8 rows
  expect_true(all(c("cell_id", "method", "macro_f1", "runtime_sec") %in% names(pilot_res$data)))
  expect_false(any(is.na(pilot_res$data$runtime_sec)))
})
