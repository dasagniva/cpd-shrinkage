# ==============================================================================
# test-90-plots-tables.R
# Unit tests for Phase P10: Plotting and LaTeX table generators.
# ==============================================================================

source("../../R/00_utils.R")
source("../../R/90_plots_tables.R")

test_that("plot_fit_panel generates valid patchwork ggplot object", {
  n <- 50
  y <- rnorm(n)
  fit <- new_cpd_fit(
    method = "test",
    cp = c(15L, 35L),
    point = c(8L, 42L),
    class = factor(sample(c("nominal", "shift", "collective", "point"), n, replace = TRUE),
                   levels = c("nominal", "shift", "collective", "point")),
    scores = data.frame(t = 1:n, s_cp = runif(n), s_out = runif(n)),
    diag = list(mu_mean = rep(0, n), mu_sd = rep(0.1, n))
  )

  p <- plot_fit_panel(fit, y, title = "Test Fit")
  expect_true(inherits(p, "patchwork") || inherits(p, "ggplot"))
})

test_that("plot_identifiability_quadrants generates valid ggplot object", {
  n <- 50
  fit <- new_cpd_fit(
    method = "test",
    class = factor(sample(c("nominal", "shift", "collective", "point"), n, replace = TRUE),
                   levels = c("nominal", "shift", "collective", "point")),
    scores = data.frame(t = 1:n, s_cp = runif(n), s_out = runif(n))
  )

  p <- plot_identifiability_quadrants(fit, alpha_fdr = 0.05, output_file = NULL)
  expect_true(inherits(p, "ggplot"))
})

test_that("table generation functions output valid LaTeX code", {
  # Mock simulation summary df
  mock_pilot_df <- data.frame(
    method      = c("proposed", "decafs", "pelt"),
    cp_f1       = c(0.85, 0.70, 0.60),
    cp_covering = c(0.90, 0.75, 0.65),
    pt_f1       = c(0.80, 0.40, 0.10),
    coll_f1     = c(0.75, 0.20, 0.05),
    macro_f1    = c(0.82, 0.55, 0.40),
    runtime_sec = c(5.2, 0.1, 0.05),
    stringsAsFactors = FALSE
  )

  tex1 <- table_sim_summary_latex(mock_pilot_df, output_file = NULL)
  expect_true(is.character(tex1))
  expect_true(grepl("tabular|table", tex1))

  # Mock real data summary
  mock_real_df <- data.frame(
    series       = c("welllog", "oilwell"),
    method       = c("proposed", "proposed"),
    n_cp         = c(10, 8),
    n_point      = c(5, 3),
    n_coll_segs  = c(2, 1),
    n_coll_pts   = c(25, 12),
    phi_est      = c(0.65, 0.45),
    runtime_sec  = c(42.5, 38.2),
    stringsAsFactors = FALSE
  )

  tex2 <- table_real_data_latex(mock_real_df, output_file = NULL)
  expect_true(is.character(tex2))
  expect_true(grepl("tabular|table", tex2))

  # Mock ablation summary
  mock_abl_df <- data.frame(
    ablation    = c("full", "no_nsp"),
    description = c("Full Proposed", "No NSP"),
    n_cp        = c(10, 8),
    n_point     = c(5, 5),
    n_coll      = c(2, 2),
    cp_f1       = c(0.88, 0.75),
    macro_f1    = c(0.84, 0.71),
    runtime_sec = c(12.0, 11.5),
    stringsAsFactors = FALSE
  )

  tex3 <- table_ablation_latex(mock_abl_df, output_file = NULL)
  expect_true(is.character(tex3))
  expect_true(grepl("tabular|table", tex3))

  # Mock sensitivity summary
  mock_sens_df <- data.frame(
    param       = c("baseline", "alpha"),
    val         = c("baseline", "0.10"),
    n_cp        = c(10, 11),
    n_point     = c(5, 6),
    n_coll      = c(2, 2),
    cp_f1       = c(0.88, 0.86),
    macro_f1    = c(0.84, 0.83),
    runtime_sec = c(12.0, 12.1),
    stringsAsFactors = FALSE
  )

  tex4 <- table_sensitivity_latex(mock_sens_df, output_file = NULL)
  expect_true(is.character(tex4))
  expect_true(grepl("tabular|table", tex4))
})
