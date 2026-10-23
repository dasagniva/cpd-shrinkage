# ==============================================================================
# run_pipeline.R
# Master Pipeline Execution and Publication Artifact Generator (Phase P10).
# Generates all publication-grade figures (Figures 1-3) and LaTeX tables (Tables 1-4)
# from persisted simulation, empirical benchmark, ablation, and sensitivity results.
# ==============================================================================

cat("\n====================================================================\n")
cat("   MASTER PUBLICATION PIPELINE: FIGURES & LATEX TABLES (Phase P10)\n")
cat("====================================================================\n\n")

source("R/00_utils.R")
source("R/85_real_data.R")
source("R/90_plots_tables.R")

dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1. Generate Table 1 & Figure 2 from the Phase P8 full simulation study
# ------------------------------------------------------------------------------
if (file.exists("output/gate_p8_results.rds")) {
  cat("1. Processing Phase P8 full simulation study results...\n")
  source("R/80_sim_design.R")
  p8_res <- readRDS("output/gate_p8_results.rds")
  stopifnot(!is.null(p8_res$full_data)) # the 5-replicate pilot is not used for tables

  # Table 1 (a-c): cells x methods, mean (Monte Carlo SE)
  labels <- c(cp_f1 = "change-point F1", pt_f1 = "point-anomaly F1", macro_f1 = "four-class macro-F1")
  for (m in names(labels)) {
    table_sim_metric_latex(
      p8_res$summary[p8_res$summary$metric == m, ],
      caption = sprintf("Simulation study: %s by cell and method, mean (Monte Carlo SE); 100 replicates for MCMC methods, 200 otherwise", labels[[m]]),
      output_file = sprintf("output/tab1_sim_%s.tex", m)
    )
  }

  # Figure 2: Simulation Factor Effects Across Benchmark Methods
  cat("   Generating Figure 2 (output/fig2_sim_factor_effects.png)... ")
  plot_sim_factor_effects(p8_res$full_data, output_file = "output/fig2_sim_factor_effects.png")
} else {
  warning("output/gate_p8_results.rds not found. Skipping Table 1 and Figure 2.")
}

# ------------------------------------------------------------------------------
# 2. Generate Tables 2-4 and Figures 1 & 3 from Phase P9 Empirical Results
# ------------------------------------------------------------------------------
if (file.exists("output/gate_p9_results.rds")) {
  cat("\n2. Processing Phase P9 Empirical & Ablation Results...\n")
  p9_res <- readRDS("output/gate_p9_results.rds")

  # Table 2: Empirical Benchmark Summary
  cat("   Generating Table 2 (output/tab2_real_data.tex)... ")
  table_real_data_latex(p9_res$real_data_summary, output_file = "output/tab2_real_data.tex")

  # Tables 3-4: replicated ablation and sensitivity, mean (MC SE) over 20 replicates per setting
  metrics <- c("cp_prec", "cp_rec", "cp_f1", "pt_f1", "pt_fdr", "coll_f1", "macro_f1")
  cat("   Generating Table 3 (output/tab3_ablation.tex)... ")
  table_replicated_latex(p9_res$ablation_results, by = c("setting", "ablation"), metrics = c(metrics, "phi_hat"),
                         caption = "Component ablation: mean (Monte Carlo SE) over 20 replicates of each setting",
                         output_file = "output/tab3_ablation.tex")
  cat("   Generating Table 4 (output/tab4_sensitivity.tex)... ")
  table_replicated_latex(p9_res$sensitivity_results, by = c("setting", "param", "val"), metrics = metrics,
                         caption = "Sensitivity grid: mean (Monte Carlo SE) over 20 replicates of each setting",
                         output_file = "output/tab4_sensitivity.tex")

  # Load real datasets for plotting
  cat("\n   Loading empirical series for visualization...\n")
  data_list <- list(
    welllog = load_real_series("welllog", raw_dir = "data/raw"),
    oilwell = load_real_series("oilwell", raw_dir = "data/raw"),
    indpro  = load_real_series("indpro", raw_dir = "data/raw"),
    btc     = load_real_series("btc", raw_dir = "data/raw")
  )

  fits_list <- list(
    welllog = p9_res$real_data_fits$welllog$proposed,
    oilwell = p9_res$real_data_fits$oilwell$proposed,
    indpro  = p9_res$real_data_fits$indpro$proposed,
    btc     = p9_res$real_data_fits$btc$proposed
  )

  # Figure 1: Four-Panel Real Data Benchmark
  cat("   Generating Figure 1 (output/fig1_real_data_four_panel.png)... ")
  plot_four_panel_real_data(fits_list, data_list, output_file = "output/fig1_real_data_four_panel.png")

  # Figure 3: Identifiability Separation in Exceedance Space (Well-log)
  cat("   Generating Figure 3 (output/fig3_identifiability_quadrants.png)... ")
  plot_identifiability_quadrants(
    fit         = p9_res$real_data_fits$welllog$proposed,
    alpha_fdr   = 0.05,
    output_file = "output/fig3_identifiability_quadrants.png"
  )
} else {
  warning("output/gate_p9_results.rds not found. Skipping Tables 2-4 and Figures 1, 3.")
}

cat("\n====================================================================\n")
cat("   PUBLICATION ARTIFACTS GENERATION COMPLETE\n")
cat("====================================================================\n\n")
