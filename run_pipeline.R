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
# 1. Generate Table 1 & Figure 2 from Phase P8 Simulation Results
# ------------------------------------------------------------------------------
if (file.exists("output/gate_p8_results.rds")) {
  cat("1. Processing Phase P8 Simulation Study Results...\n")
  p8_res <- readRDS("output/gate_p8_results.rds")
  pilot_df <- p8_res$pilot_data

  # Table 1: Simulation Benchmark Summary
  cat("   Generating Table 1 (output/tab1_sim_summary.tex)... ")
  table_sim_summary_latex(pilot_df, output_file = "output/tab1_sim_summary.tex")

  # Figure 2: Simulation Factor Effects Across Benchmark Methods
  cat("   Generating Figure 2 (output/fig2_sim_factor_effects.png)... ")
  plot_sim_factor_effects(pilot_df, output_file = "output/fig2_sim_factor_effects.png")
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

  # Table 3: Component Ablation Study
  cat("   Generating Table 3 (output/tab3_ablation.tex)... ")
  table_ablation_latex(p9_res$ablation_results, output_file = "output/tab3_ablation.tex")

  # Table 4: Hyperparameter Sensitivity Analysis
  cat("   Generating Table 4 (output/tab4_sensitivity.tex)... ")
  table_sensitivity_latex(p9_res$sensitivity_results, output_file = "output/tab4_sensitivity.tex")

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
