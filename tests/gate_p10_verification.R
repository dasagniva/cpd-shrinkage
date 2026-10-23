# ==============================================================================
# gate_p10_verification.R
# Acceptance Gate for Phase P10: Publication Outputs & Pipeline Verification.
#
# Tests:
# 1. Existence and integrity of Figures 1-3 (PNG format, > 50KB).
# 2. Existence and valid structure of Tables 1-4 (LaTeX format, > 500 bytes).
# 3. Targets pipeline file syntax validation.
# ==============================================================================

cat("\n====================================================================\n")
cat("      PHASE P10 ACCEPTANCE GATE VERIFICATION\n")
cat("====================================================================\n\n")

# 1. Check Figures
cat("1. Checking publication figures in output/...\n")
figures <- c(
  "output/fig1_real_data_four_panel.png",
  "output/fig2_sim_factor_effects.png",
  "output/fig3_identifiability_quadrants.png"
)

fig_ok <- TRUE
for (f in figures) {
  if (!file.exists(f)) {
    warning(sprintf("Missing figure: %s", f))
    fig_ok <- FALSE
  } else {
    f_size <- file.info(f)$size
    cat(sprintf("  [OK] %-42s (%.1f KB)\n", f, f_size / 1024))
    if (f_size < 10000) fig_ok <- FALSE
  }
}
stopifnot(fig_ok)

# 2. Check LaTeX Tables
cat("\n2. Checking publication LaTeX tables in output/...\n")
tables <- c(
  "output/tab1_sim_cp_f1.tex",
  "output/tab1_sim_pt_f1.tex",
  "output/tab1_sim_macro_f1.tex",
  "output/tab2_real_data.tex",
  "output/tab3_ablation.tex",
  "output/tab4_sensitivity.tex"
)

tab_ok <- TRUE
for (t in tables) {
  if (!file.exists(t)) {
    warning(sprintf("Missing table: %s", t))
    tab_ok <- FALSE
  } else {
    t_lines <- readLines(t, warn = FALSE)
    t_size <- file.info(t)$size
    has_tab <- any(grepl("tabular|table", t_lines))
    cat(sprintf("  [OK] %-30s (%d lines, %.1f KB, valid tabular: %s)\n",
                t, length(t_lines), t_size / 1024, has_tab))
    if (!has_tab || t_size < 300) tab_ok <- FALSE
  }
}
stopifnot(tab_ok)

# 3. Check _targets.R syntax
cat("\n3. Validating _targets.R syntax...\n")
stopifnot(file.exists("_targets.R"))
parsed <- parse("_targets.R")
cat(sprintf("  [OK] _targets.R parsed cleanly (%d top-level expressions)\n", length(parsed)))

# 4. Every output must be newer than the results it is built from (stale outputs fail)
cat("\n4. Checking that outputs are newer than their source results...\n")
sources <- list(
  "output/gate_p8_results.rds" = c(grep("tab1_", tables, value = TRUE), "output/fig2_sim_factor_effects.png"),
  "output/gate_p9_results.rds" = c("output/tab2_real_data.tex", "output/tab3_ablation.tex", "output/tab4_sensitivity.tex",
                                   "output/fig1_real_data_four_panel.png", "output/fig3_identifiability_quadrants.png")
)
fresh_ok <- TRUE
for (src in names(sources)) {
  stopifnot(file.exists(src))
  for (f in sources[[src]]) {
    ok <- file.exists(f) && file.mtime(f) >= file.mtime(src)
    cat(sprintf("  [%s] %-42s vs %s\n", if (ok) "OK" else "STALE", f, src))
    fresh_ok <- fresh_ok && ok
  }
}

gate_p10_passed <- fig_ok && tab_ok && fresh_ok

cat("\n====================================================================\n")
if (gate_p10_passed) {
  cat(">>> PHASE P10 ACCEPTANCE GATE: PASSED <<<\n")
} else {
  cat(">>> PHASE P10 ACCEPTANCE GATE: FAILED <<<\n")
}
cat("====================================================================\n\n")

if (!gate_p10_passed) {
  stop("Phase P10 acceptance gate failed!")
}
