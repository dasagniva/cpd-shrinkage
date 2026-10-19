# ==============================================================================
# 90_plots_tables.R
# Publication-Grade Visualizations and LaTeX Table Generators (Phase P10).
# Generates Figures 1-4 and Tables 1-4 as specified in Section 5 (Phase P10)
# and Section 6 (Reporting Standards).
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(knitr)
  library(kableExtra)
})

# Universal 4-class palette (colorblind-safe, high contrast)
CPD_PALETTE <- c(
  nominal    = "#757575", # Slate grey
  shift      = "#1565C0", # Deep cobalt blue
  collective = "#F57C00", # Vivid amber / orange
  point      = "#D32F2F"  # Crimson red
)

#' Plot a Single Time Series Fit Panel (Top: mu + points; Bottom: Outlier Score)
#'
#' @param fit A validated `cpd_fit` object
#' @param y Numeric vector of observed data
#' @param dates Optional vector of dates/timestamps
#' @param title Character plot title
#'
#' @return A `patchwork` ggplot object
#' @export
plot_fit_panel <- function(fit, y, dates = NULL, title = "") {
  n <- length(y)
  t_vec <- if (!is.null(dates) && length(dates) == n) {
    if (inherits(dates, "Date") || inherits(dates, "POSIXt")) dates else as.Date(dates)
  } else {
    seq_len(n)
  }

  cls <- if (!is.null(fit$class) && length(fit$class) == n) {
    as.character(fit$class)
  } else {
    rep("nominal", n)
  }
  cls <- factor(cls, levels = names(CPD_PALETTE))

  # Extract mu and uncertainty bands
  mu_hat <- if (!is.null(fit$diag$mu_mean) && length(fit$diag$mu_mean) == n) {
    fit$diag$mu_mean
  } else {
    NULL
  }
  mu_sd <- if (!is.null(fit$diag$mu_sd) && length(fit$diag$mu_sd) == n) {
    fit$diag$mu_sd
  } else {
    NULL
  }

  df_top <- data.frame(
    t     = t_vec,
    y     = y,
    class = cls
  )

  if (!is.null(mu_hat)) {
    df_top$mu <- mu_hat
    if (!is.null(mu_sd)) {
      df_top$mu_lower <- mu_hat - 1.645 * mu_sd # 90% credible band
      df_top$mu_upper <- mu_hat + 1.645 * mu_sd
    }
  }

  # Build top plot
  p_top <- ggplot(df_top, aes(x = t))

  # Shaded 90% credible band
  if (!is.null(mu_hat) && !is.null(mu_sd)) {
    p_top <- p_top + geom_ribbon(aes(ymin = mu_lower, ymax = mu_upper),
                                 fill = "#BBDEFB", alpha = 0.5)
  }

  # Points colored by 4-class classification
  p_top <- p_top +
    geom_point(aes(y = y, color = class), size = 1.2, alpha = 0.8) +
    scale_color_manual(values = CPD_PALETTE, drop = FALSE, name = "Classification")

  # Posterior mean mu line
  if (!is.null(mu_hat)) {
    p_top <- p_top + geom_line(aes(y = mu), color = "#0D47A1", linewidth = 0.9)
  }

  # Vertical lines at detected change points
  if (length(fit$cp) > 0) {
    cp_t <- t_vec[fit$cp]
    p_top <- p_top + geom_vline(xintercept = as.numeric(cp_t), linetype = "dashed",
                                color = "#1565C0", alpha = 0.7, linewidth = 0.5)
  }

  p_top <- p_top +
    labs(title = title, x = NULL, y = "Value") +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 12)
    )

  # Build bottom plot: outlier scores o_t
  s_out <- if (!is.null(fit$scores$s_out) && length(fit$scores$s_out) == n) {
    fit$scores$s_out
  } else {
    rep(0.0, n)
  }

  df_bot <- data.frame(t = t_vec, o_t = s_out)

  p_bot <- ggplot(df_bot, aes(x = t, y = o_t)) +
    geom_col(fill = "#D32F2F", width = 1.0, alpha = 0.75) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "black", linewidth = 0.5) +
    ylim(0, 1) +
    labs(x = if (!is.null(dates)) "Time" else "Index t", y = expression(o[t])) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())

  # Combine using patchwork
  (p_top / p_bot) + plot_layout(heights = c(3.2, 1))
}

#' Generate Figure 1: Four-Panel Real-Data Benchmark
#'
#' @param fits_list Named list of fits: `list(welllog = ..., oilwell = ..., indpro = ..., btc = ...)`
#' @param data_list Named list of datasets loaded from `load_real_series`
#' @param output_file Path to save PNG
#'
#' @return Combined ggplot object
#' @export
plot_four_panel_real_data <- function(fits_list, data_list, output_file = "output/fig1_real_data_four_panel.png") {
  p_wl <- plot_fit_panel(fits_list$welllog, data_list$welllog$y, title = "A. Nuclear Magnetic Resonance Well-Log (n = 4,050)")
  p_ow <- plot_fit_panel(fits_list$oilwell, data_list$oilwell$y, title = "B. Oilwell Downhole Pressure (n = 4,000)")
  p_ip <- plot_fit_panel(fits_list$indpro, data_list$indpro$y, dates = data_list$indpro$dates,
                         title = "C. US Industrial Production Monthly Log Growth (n = 312)")
  p_bc <- plot_fit_panel(fits_list$btc, data_list$btc$y, dates = data_list$btc$dates,
                         title = "D. Bitcoin Daily Log-Price 2018-2023 (n = 2,191)")

  combined <- (p_wl | p_ow) / (p_ip | p_bc) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggsave(output_file, combined, width = 16, height = 11, dpi = 300)
    cat(sprintf("Figure 1 saved to %s\n", output_file))
  }

  combined
}

#' Generate Figure 2: Simulation Factor Effects Across Benchmark Methods
#'
#' @param pilot_df Data frame from `run_pilot_simulation`
#' @param output_file Path to save PNG
#'
#' @return ggplot object
#' @export
plot_sim_factor_effects <- function(pilot_df, output_file = "output/fig2_sim_factor_effects.png") {
  # Aggregate Macro-F1 across factors and methods
  agg_df <- aggregate(
    macro_f1 ~ factor_name + cell_id + method,
    data = pilot_df,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = sd(x, na.rm = TRUE) / sqrt(length(x)))
  )

  df_plot <- data.frame(
    factor_name = agg_df$factor_name,
    cell_id     = agg_df$cell_id,
    method      = agg_df$method,
    macro_f1    = agg_df$macro_f1[, "mean"],
    se          = agg_df$macro_f1[, "se"]
  )

  p <- ggplot(df_plot, aes(x = cell_id, y = macro_f1, color = method, group = method)) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = macro_f1 - se, ymax = macro_f1 + se), width = 0.2) +
    facet_wrap(~factor_name, scales = "free_x") +
    labs(
      title = "Figure 2: Four-Class Macro-F1 Across Simulation Design Factors",
      subtitle = "Evaluates Proposed Global-Local Shrinkage against Benchmark Comparators",
      x = "Experimental Cell ID",
      y = "Macro-F1 (with MC Standard Error)",
      color = "Method"
    ) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggsave(output_file, p, width = 14, height = 9, dpi = 300)
    cat(sprintf("Figure 2 saved to %s\n", output_file))
  }

  p
}

#' Generate Figure 3: Identifiability Separation in Exceedance Space (p_omega vs p_psi)
#'
#' @param fit Proposed `cpd_fit` object with scores
#' @param alpha_fdr FDR threshold line
#' @param output_file Path to save PNG
#'
#' @return ggplot object
#' @export
plot_identifiability_quadrants <- function(fit, alpha_fdr = 0.05, output_file = "output/fig3_identifiability_quadrants.png") {
  stopifnot(!is.null(fit$scores$s_cp), !is.null(fit$scores$s_out), !is.null(fit$class))

  df_quad <- data.frame(
    p_omega = fit$scores$s_cp,
    p_psi   = fit$scores$s_out,
    class   = fit$class
  )

  p <- ggplot(df_quad, aes(x = p_omega, y = p_psi, color = class)) +
    geom_point(alpha = 0.75, size = 1.8) +
    geom_vline(xintercept = alpha_fdr, linetype = "dashed", color = "black", alpha = 0.6) +
    geom_hline(yintercept = alpha_fdr, linetype = "dashed", color = "black", alpha = 0.6) +
    scale_color_manual(values = CPD_PALETTE, name = "Class") +
    annotate("text", x = 0.8, y = 0.1, label = "Sustained Shifts\n(High p_omega, Low p_psi)", size = 3.5, fontface = "italic") +
    annotate("text", x = 0.1, y = 0.8, label = "Point Anomalies\n(Low p_omega, High p_psi)", size = 3.5, fontface = "italic") +
    annotate("text", x = 0.1, y = 0.05, label = "Nominal", size = 3.5, fontface = "italic") +
    labs(
      title = "Figure 3: Identifiability Separation in Exceedance Space",
      subtitle = "Decoupling Sustained Mean Shifts from Point Outliers via Tail Probabilities",
      x = expression(Tail~Exceedance~Probability~p[t]^omega~(Sustained~Shift)),
      y = expression(Outlier~Posterior~Score~o[t]~or~p[t]^psi~(Point~Outlier))
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right", panel.grid.minor = element_blank())

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    ggsave(output_file, p, width = 8, height = 6.5, dpi = 300)
    cat(sprintf("Figure 3 saved to %s\n", output_file))
  }

  p
}

#' Format and Output LaTeX Table 1: Simulation Benchmark Summary
#'
#' @param pilot_df Data frame from `run_pilot_simulation`
#' @param output_file Path to save .tex file
#'
#' @return Formatted character string of LaTeX code
#' @export
table_sim_summary_latex <- function(pilot_df, output_file = "output/tab1_sim_summary.tex") {
  summary_by_method <- aggregate(
    cbind(cp_f1, cp_covering, pt_f1, coll_f1, macro_f1, runtime_sec) ~ method,
    data = pilot_df,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), se = sd(x, na.rm = TRUE) / sqrt(length(x)))
  )

  df_tab <- data.frame(
    Method      = summary_by_method$method,
    CP_F1       = sprintf("%.3f (%.3f)", summary_by_method$cp_f1[, "mean"], summary_by_method$cp_f1[, "se"]),
    Covering    = sprintf("%.3f (%.3f)", summary_by_method$cp_covering[, "mean"], summary_by_method$cp_covering[, "se"]),
    Point_F1    = sprintf("%.3f (%.3f)", summary_by_method$pt_f1[, "mean"], summary_by_method$pt_f1[, "se"]),
    Coll_F1     = sprintf("%.3f (%.3f)", summary_by_method$coll_f1[, "mean"], summary_by_method$coll_f1[, "se"]),
    Macro_F1    = sprintf("%.3f (%.3f)", summary_by_method$macro_f1[, "mean"], summary_by_method$macro_f1[, "se"]),
    Runtime_sec = sprintf("%.2f s", summary_by_method$runtime_sec[, "mean"]),
    stringsAsFactors = FALSE
  )

  colnames(df_tab) <- c("Method", "CP F1 (SE)", "Covering (SE)", "Point F1 (SE)", "Coll F1 (SE)", "Macro-F1 (SE)", "Mean Runtime")

  tex_code <- kable(df_tab, format = "latex", booktabs = TRUE, caption = "Simulation Benchmark Performance Across 17 Design Cells (Monte Carlo SE in parentheses)") %>%
    kable_styling(latex_options = c("striped", "hold_position"))

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    writeLines(as.character(tex_code), output_file)
    cat(sprintf("Table 1 saved to %s\n", output_file))
  }

  tex_code
}

#' Format and Output LaTeX Table 2: Empirical Real-Data Benchmark Summary
#'
#' @param summary_df Data frame from `run_real_data_benchmark`
#' @param output_file Path to save .tex file
#'
#' @return Formatted character string of LaTeX code
#' @export
table_real_data_latex <- function(summary_df, output_file = "output/tab2_real_data.tex") {
  df_tab <- summary_df[, c("series", "method", "n_cp", "n_point", "n_coll_segs", "n_coll_pts", "phi_est", "runtime_sec")]
  colnames(df_tab) <- c("Dataset", "Method", "CPs", "Points", "Coll. Segs", "Coll. Pts", "phi_hat", "Runtime (s)")

  df_tab$`phi_hat` <- ifelse(is.na(df_tab$`phi_hat`), "--", sprintf("%.3f", df_tab$`phi_hat`))
  df_tab$`Runtime (s)` <- sprintf("%.2f", df_tab$`Runtime (s)`)

  tex_code <- kable(df_tab, format = "latex", booktabs = TRUE, caption = "Empirical Application Feature Detections and Runtime Across Four Benchmark Series") %>%
    kable_styling(latex_options = c("striped", "hold_position"))

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    writeLines(as.character(tex_code), output_file)
    cat(sprintf("Table 2 saved to %s\n", output_file))
  }

  tex_code
}

#' Format and Output LaTeX Table 3: Component Ablation Study Results
#'
#' @param ablation_df Data frame from `run_ablation_study`
#' @param output_file Path to save .tex file
#'
#' @return Formatted character string of LaTeX code
#' @export
table_ablation_latex <- function(ablation_df, output_file = "output/tab3_ablation.tex") {
  cols_keep <- intersect(names(ablation_df), c("ablation", "description", "n_cp", "n_point", "n_coll", "cp_f1", "macro_f1", "runtime_sec"))
  df_tab <- ablation_df[, cols_keep]

  if ("cp_f1" %in% names(df_tab)) df_tab$cp_f1 <- sprintf("%.3f", df_tab$cp_f1)
  if ("macro_f1" %in% names(df_tab)) df_tab$macro_f1 <- sprintf("%.3f", df_tab$macro_f1)
  if ("runtime_sec" %in% names(df_tab)) df_tab$runtime_sec <- sprintf("%.2f s", df_tab$runtime_sec)

  tex_code <- kable(df_tab, format = "latex", booktabs = TRUE, caption = "Component Ablation Study: Dissecting Significance Guidance, Anomaly Coupling, and Innovation Scaling") %>%
    kable_styling(latex_options = c("striped", "hold_position"))

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    writeLines(as.character(tex_code), output_file)
    cat(sprintf("Table 3 saved to %s\n", output_file))
  }

  tex_code
}

#' Format and Output LaTeX Table 4: Sensitivity Analysis Grid
#'
#' @param sens_df Data frame from `run_sensitivity_grid`
#' @param output_file Path to save .tex file
#'
#' @return Formatted character string of LaTeX code
#' @export
table_sensitivity_latex <- function(sens_df, output_file = "output/tab4_sensitivity.tex") {
  cols_keep <- intersect(names(sens_df), c("param", "val", "n_cp", "n_point", "n_coll", "cp_f1", "macro_f1", "runtime_sec"))
  df_tab <- sens_df[, cols_keep]

  if ("cp_f1" %in% names(df_tab)) df_tab$cp_f1 <- sprintf("%.3f", df_tab$cp_f1)
  if ("macro_f1" %in% names(df_tab)) df_tab$macro_f1 <- sprintf("%.3f", df_tab$macro_f1)
  if ("runtime_sec" %in% names(df_tab)) df_tab$runtime_sec <- sprintf("%.2f s", df_tab$runtime_sec)

  tex_code <- kable(df_tab, format = "latex", booktabs = TRUE, caption = "Hyperparameter Sensitivity Grid: Robustness Across Significance Levels, Priors, and Truncation Bounds") %>%
    kable_styling(latex_options = c("striped", "hold_position"))

  if (!is.null(output_file)) {
    dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
    writeLines(as.character(tex_code), output_file)
    cat(sprintf("Table 4 saved to %s\n", output_file))
  }

  tex_code
}
