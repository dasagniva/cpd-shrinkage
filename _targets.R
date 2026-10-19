# ==============================================================================
# _targets.R
# Reproducible Targets Workflow Pipeline (Phase P10).
# Defines the end-to-end target graph for simulation benchmark, empirical
# applications, ablation analysis, sensitivity grid, tables, and figures.
# ==============================================================================

if (requireNamespace("targets", quietly = TRUE)) {
  library(targets)

  # Source all R functions
  lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source)

  tar_option_set(
    packages = c(
      "spam", "rvinecopulib", "anomaly", "nsp",
      "changepoint", "robseg", "stepR", "bcp", "EnvCpt",
      "ggplot2", "patchwork", "knitr", "kableExtra"
    ),
    format = "rds"
  )

  list(
    # 1. Experimental Grid
    tar_target(sim_grid, create_simulation_grid()),

    # 2. Real Datasets
    tar_target(wl_data, load_real_series("welllog", raw_dir = "data/raw")),
    tar_target(ow_data, load_real_series("oilwell", raw_dir = "data/raw")),
    tar_target(ip_data, load_real_series("indpro", raw_dir = "data/raw")),
    tar_target(bc_data, load_real_series("btc", raw_dir = "data/raw")),

    # 3. Real Data Application Benchmark
    tar_target(real_benchmark_results, {
      run_real_data_benchmark(
        series_names = c("welllog", "oilwell", "indpro", "btc"),
        methods      = c("proposed", "copula", "decafs", "pelt", "robseg", "capa", "fpnn", "ch4"),
        raw_dir      = "data/raw"
      )
    }),

    # 4. Component Ablation Study
    tar_target(ablation_results, {
      set.seed(42)
      sim <- simulate_series(n = 500, phi = 0.5, noise = "gaussian", jump = 2.0,
                             outlier_regime = "mixed", p_outlier = 0.03, collective = 2, seed = 42)
      run_ablation_study(sim$y, truth = sim$truth, n_iter = 2000, burn_in = 800)
    }),

    # 5. Sensitivity Analysis Grid
    tar_target(sensitivity_results, {
      set.seed(42)
      sim <- simulate_series(n = 500, phi = 0.5, noise = "gaussian", jump = 2.0,
                             outlier_regime = "mixed", p_outlier = 0.03, collective = 2, seed = 42)
      run_sensitivity_grid(sim$y, truth = sim$truth, n_iter = 1500, burn_in = 500)
    }),

    # 6. Tables (LaTeX)
    tar_target(tab1_file, {
      if (file.exists("output/gate_p8_results.rds")) {
        p8_res <- readRDS("output/gate_p8_results.rds")
        table_sim_summary_latex(p8_res$pilot_data, output_file = "output/tab1_sim_summary.tex")
      }
      "output/tab1_sim_summary.tex"
    }, format = "file"),

    tar_target(tab2_file, {
      table_real_data_latex(real_benchmark_results$summary_table, output_file = "output/tab2_real_data.tex")
      "output/tab2_real_data.tex"
    }, format = "file"),

    tar_target(tab3_file, {
      table_ablation_latex(ablation_results$results_df, output_file = "output/tab3_ablation.tex")
      "output/tab3_ablation.tex"
    }, format = "file"),

    tar_target(tab4_file, {
      table_sensitivity_latex(sensitivity_results, output_file = "output/tab4_sensitivity.tex")
      "output/tab4_sensitivity.tex"
    }, format = "file"),

    # 7. Figures (PNG)
    tar_target(fig1_file, {
      data_list <- list(welllog = wl_data, oilwell = ow_data, indpro = ip_data, btc = bc_data)
      fits_list <- list(
        welllog = real_benchmark_results$fits$welllog$proposed,
        oilwell = real_benchmark_results$fits$oilwell$proposed,
        indpro  = real_benchmark_results$fits$indpro$proposed,
        btc     = real_benchmark_results$fits$btc$proposed
      )
      plot_four_panel_real_data(fits_list, data_list, output_file = "output/fig1_real_data_four_panel.png")
      "output/fig1_real_data_four_panel.png"
    }, format = "file"),

    tar_target(fig3_file, {
      fit_wl <- real_benchmark_results$fits$welllog$proposed
      plot_identifiability_quadrants(fit_wl, alpha_fdr = 0.05, output_file = "output/fig3_identifiability_quadrants.png")
      "output/fig3_identifiability_quadrants.png"
    }, format = "file")
  )
}
