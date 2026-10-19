# ==============================================================================
# test-85-86-real-data-ablation.R
# Unit tests for Phase P9: Real Data Loaders, Ablations, and Sensitivity Grid.
# ==============================================================================

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
source("../../R/85_real_data.R")
source("../../R/86_sensitivity_ablation.R")

test_that("load_real_series accurately loads all four benchmark series", {
  # 1. Well-log: 4,050 points
  wl <- load_real_series("welllog", raw_dir = "../../data/raw")
  expect_equal(wl$name, "welllog")
  expect_equal(wl$n, 4050)
  expect_false(anyNA(wl$y))
  expect_true(is.numeric(wl$y))

  # 2. Oilwell: 4,000 points
  ow <- load_real_series("oilwell", raw_dir = "../../data/raw")
  expect_equal(ow$name, "oilwell")
  expect_equal(ow$n, 4000)
  expect_false(anyNA(ow$y))

  # 3. INDPRO: 312 points (monthly log growth)
  ip <- load_real_series("indpro", raw_dir = "../../data/raw")
  expect_equal(ip$name, "indpro")
  expect_equal(ip$n, 312)
  expect_false(anyNA(ip$y))
  expect_equal(length(ip$dates), 312)

  # 4. BTC: 2,191 points (2018-01-01 to 2023-12-31)
  btc <- load_real_series("btc", raw_dir = "../../data/raw")
  expect_equal(btc$name, "btc")
  expect_equal(btc$n, 2191)
  expect_false(anyNA(btc$y))
  expect_equal(btc$dates[1], "2018-01-01")
  expect_equal(btc$dates[length(btc$dates)], "2023-12-31")
})

test_that("run_real_data_method and summarize_real_data_fit produce valid cpd_fit and summaries", {
  # Mini slice of INDPRO for fast smoke test
  ip <- load_real_series("indpro", raw_dir = "../../data/raw")
  y_sub <- ip$y[1:50]

  fit_pelt <- run_real_data_method(y_sub, method = "pelt")
  expect_true(inherits(fit_pelt, "cpd_fit"))
  expect_false(is.null(fit_pelt$class))

  sum_df <- summarize_real_data_fit(fit_pelt, series_name = "indpro_sub")
  expect_equal(nrow(sum_df), 1)
  expect_equal(sum_df$series, "indpro_sub")
  expect_equal(sum_df$method, "pelt")
  expect_true(is.numeric(sum_df$runtime_sec))
})

test_that("run_ablation_study executes all six ablations plus full model", {
  # Fast synthetic series
  set.seed(901)
  sim <- simulate_series(n = 60, phi = 0.2, noise = "gaussian", jump = 2.0,
                         outlier_regime = "isolated", p_outlier = 0.03, collective = 1, seed = 901)

  abl_res <- run_ablation_study(
    y        = sim$y,
    truth    = sim$truth,
    n_iter   = 100,
    burn_in  = 40,
    n_chains = 1,
    seed     = 901
  )

  expect_true(is.list(abl_res))
  expect_true(is.data.frame(abl_res$results_df))
  # 7 rows: full + 6 ablations
  expect_equal(nrow(abl_res$results_df), 7)
  expect_equal(abl_res$results_df$ablation,
               c("full", "no_nsp", "no_capa", "indep_outliers", "no_ar1", "standard_hs", "marginal_var"))
  expect_true(all(!is.na(abl_res$results_df$runtime_sec)))
  expect_true(all(!is.na(abl_res$results_df$macro_f1)))
})

test_that("run_sensitivity_grid executes all parameter sweeps", {
  set.seed(902)
  sim <- simulate_series(n = 60, phi = 0.2, noise = "gaussian", jump = 2.0,
                         outlier_regime = "isolated", p_outlier = 0.03, collective = 1, seed = 902)

  sens_df <- run_sensitivity_grid(
    y       = sim$y,
    truth   = sim$truth,
    n_iter  = 100,
    burn_in = 40,
    seed    = 902
  )

  expect_true(is.data.frame(sens_df))
  # 12 parameter configurations
  expect_equal(nrow(sens_df), 12)
  expect_true(all(!is.na(sens_df$runtime_sec)))
  expect_true(all(!is.na(sens_df$macro_f1)))
})
