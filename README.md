# Significance-Guided Global–Local Shrinkage for Autocorrelated Time Series

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

> **Status (23 October 2026):** gates P2-P10 pass on the current code, and the full
> simulation study (P8) and the replicated ablation and sensitivity study (P9) are complete.
> Read Section 4 before citing any result: at the spec MCMC length only 17% of the
> simulation fits of the proposed method reach R-hat < 1.01, and the ablation does not
> isolate a contribution of each component.

Official implementation and reproducibility suite for:
> **Significance-Guided Global–Local Shrinkage for Distinguishing Change Points, Collective Anomalies and Point Outliers in Autocorrelated Time Series**  
> *Agniva Das and Muralidharan Kunnummal (2026).*

---

## 1. Overview and Key Contributions

When time series exhibit serial autocorrelation alongside sparse anomalies, standard segmentation methods face an acute **identifiability problem**:
- Classical change-point algorithms (e.g. PELT, DeCAFS) frequently break single isolated spikes into pairs of spurious opposite-signed shifts.
- Standard anomaly detectors (e.g. CAPA) can mistake gradual or sustained level shifts for bursts of anomalies.
- Heavy-tailed innovations and dynamic volatility clustering can degrade boundary localization.

This package implements the **Significance-Guided Global–Local Shrinkage Framework (SG-GLS)** to solve this problem:
1. **Significance-Informed Prior Modulation:** Narrowest Significance Pursuit (NSP) intervals and robustly prewhitened Collective and Point Anomaly (CAPA) detections map into local shrinkage scale multipliers $\kappa_t$ (Eq. 5) and dynamic offset indicators $a_t$ (Eq. 7).
2. **Joint Banded Gaussian Updates:** An interleaved parameterization stacks $x = (\mu_1, \psi_1, \mu_2, \psi_2, \dots)^\top$, yielding a penta-diagonal precision matrix $Q$ drawn in $O(n)$ time via a banded Cholesky factorization (C++, `src/gibbs_fast.cpp`; the `spam` version is kept as a reference).
3. **Dynamic Outlier Coupling:** Persistent contamination segments are decoupled from isolated spikes through Pólya–Gamma augmented dynamic scale processes (`dsp`).
4. **Bayesian FDR Classification:** Outlier scores $o_t = \nu_t^2 / (\nu_t^2 + \sigma^2)$ and tail exceedance probabilities $p_t^\omega, p_t^\psi$ partition time series observations into four mutually exclusive classes: **Nominal**, **Sustained Shift**, **Collective Anomaly**, and **Point Outlier**.

---

## 2. Mathematical Formulation

The generative model for $t = 1, \dots, n$ is:
  $$y_t = \mu_t + \psi_t + \varepsilon_t,$$
  $$\varepsilon_t = \phi \varepsilon_{t-1} + u_t, \quad u_t \stackrel{\text{iid}}{\sim} \mathcal{N}(0, \sigma^2), \quad |\phi| < 1,$$
  $$\omega_t = \mu_{t+1} - \mu_t \mid \tau_\omega, \lambda_t, \eta_t, \sigma \sim \mathcal{N}(0, \sigma^2 \tau_\omega^2 \kappa_t^2 \lambda_t^2 \eta_t^2), $$
  $$\psi_t \mid \nu_t \sim \mathcal{N}(0, \nu_t^2), \quad \log \nu_t^2 = m_h + h_t + \delta a_t,$$
where $\kappa_t \ge 1$ inflates the local Horseshoe+ scale on candidate significance intervals $\mathcal{I}$, and $a_t \in \{0, 1\}$ provides dynamic offset guidance from robust prewhitened anomaly screening.

---

## 3. Repository Architecture

```text
cpd-shrinkage/
├── R/
│   ├── 00_utils.R              # Common cpd_fit S3 schema, validation, RNG streams
│   ├── 10_dgp.R                # Multi-regime synthetic data generating process (DGP)
│   ├── 20_frontend_nsp.R       # Narrowest Significance Pursuit (NSP) front-end
│   ├── 21_frontend_capa.R      # Robustly prewhitened CAPA front-end
│   ├── 22_prior_map.R          # Candidate set I construction and (kappa_t, a_t) maps
│   ├── 30_stan_reference.R     # Stan reference model wrapper and CmdStan compilation
│   ├── 31_gibbs_sampler.R      # Multi-chain Production Gibbs Sampler engine
│   ├── 32_gibbs_updates.R      # Full conditional updates (banded Cholesky, slice, sample_dsp_fast)
│   ├── 40_scoring.R            # Outlier score o_t (Eq. 8) and tail exceedance p_t (Eq. 9)
│   ├── 41_bfdr.R               # Bayesian False Discovery Rate (BFDR) selection rule
│   ├── 42_classify.R           # Sequential four-class classification and cpd_fit builder
│   ├── 50_copula.R             # Copula Markov comparator (AIC family selection, bootstrap)
│   ├── 60_baselines.R          # Baseline wrappers under the unified cpd_fit schema
│   ├── 61_fpnn.R               # JISPS / FPNN baseline (port of ev-decafs-serve)
│   ├── 70_metrics.R            # Change point F1, covering, Hausdorff, PR-AUC, macro-F1
│   ├── 80_sim_design.R         # 17-cell simulation grid, replicate runner, MC SE reporting
│   ├── 85_real_data.R          # Loaders and runners for 4 empirical series
│   ├── 86_sensitivity_ablation.R # 7 ablation configurations and 21-row sensitivity grid
│   └── 90_plots_tables.R       # Publication ggplot2 figures and LaTeX kable tables
├── src/
│   ├── gibbs_fast.cpp          # Banded Cholesky, tridiagonal sampler, AR(1) filter (Rcpp)
│   └── ev_decafs_dp.cpp        # EV-DeCAFS grid dynamic programme for FPNN (Rcpp)
├── stan/
│   └── hsplus_ar1_ref.stan     # Reference Stan implementation with non-centred priors
├── data/
│   ├── raw/                    # Well-log, Oilwell, US Industrial Production, Bitcoin
│   └── processed/              # Formatted benchmarks
├── tests/
│   ├── testthat/               # Unit tests, incl. exactness checks of each Gibbs update
│   ├── gate_p2_verification.R  # NSP front-end null calibration gate
│   ├── gate_p3_verification.R  # Reference Stan vs Gibbs verification gate
│   ├── gate_p4_verification.R  # Production Gibbs convergence and speed gate
│   ├── gate_p5_verification.R  # BFDR error control and spike identifiability gate
│   ├── gate_p6_verification.R  # Copula Markov comparator calibration gate
│   ├── gate_p7_verification.R  # 12 Baseline wrappers validation gate
│   ├── gate_p8_verification.R  # Full simulation study gate (completeness, failures, convergence)
│   ├── gate_p9_verification.R  # Empirical applications, ablation, sensitivity gate
│   ├── gate_p10_verification.R # Publication outputs and LaTeX table gate
│   ├── test_geweke_getting_it_right.R # Geweke (2004) joint-distribution test (part of P4)
│   ├── run_welllog_long.R      # Well-log run at the converged length (P4 mixing and speed)
│   ├── run_real_data_long.R    # Four real series, all methods (P9); + run_real_data_btc_extend.R
│   ├── run_p8_full.R           # Full simulation study, resumable (output/p8_full/)
│   └── run_p9_replicated.R     # Replicated ablation and sensitivity, resumable (output/p9_replicated/)
├── output/                     # Gate results (.rds), figures (PNG) and LaTeX tables (.tex)
├── run_all.sh                  # Full rerun in dependency order, with stop/resume
├── run_pipeline.R              # Figures and tables from the saved gate results
└── DESCRIPTION                 # R package metadata
```

---

## 4. Verification Gates: Current Status

Each phase in `CPD2_Code_Workflow.md` has an acceptance gate. Gate scripts save their
results to `output/gate_pX_results.rds` and stop with an error if a criterion fails.

| Phase | Gate script | Status (last run) |
|---|---|---|
| P0 | `renv::status()` | `robseg` locked; `SBC` is still recorded from an unknown source |
| P1 | `tests/testthat/test-10-dgp.R` | Unit tests pass (all tests: 468 checks, 0 failed, 23 Oct) |
| P2 | `tests/gate_p2_verification.R` | Passed for the pipeline front-end `nsp_poly_ar` (false-interval rate <= 0.015 at phi = 0, 0.6, 0.9); the spec's `nsp_poly_selfnorm` fails (0.38, 1.00), see SPEC_CHANGES.md C1 (19 Oct) |
| P3 | `tests/gate_p3_verification.R` | Passed: parameter recovery and SBC (KS p: phi 0.010, sigma 0.051, tau_omega 0.60) (21 Oct; Stan model unchanged since) |
| P4 | `tests/gate_p4_verification.R` | Passed: Stan agreement, Getting-it-Right test (both configurations pass, planted error detected), well-log mixing at 4 x 50,000 iterations (gated R-hat < 1.01, ESS > 400), 8.0 min per chain (23 Oct) |
| P5 | `tests/gate_p5_verification.R` | Passed: realised shift and point FDR 0; power 0.27 (shift) and 0.93 (point), against a required 0.2; spike identifiability 1.00 (21 Oct) |
| P6 | `tests/gate_p6_verification.R` | Passed: copula family recovery and null calibration (21 Oct) |
| P7 | `tests/gate_p7_verification.R` | Passed after the baseline fixes (robseg, PELT/CAPA scaling, ch4) (21 Oct) |
| P8 | `tests/gate_p8_verification.R` | Passed: full study complete, 17 cells x 13 methods, 0 failed fits, 0 silent NA metrics (23 Oct). The gate checks completeness; convergence is reported, see below |
| P9 | `tests/gate_p9_verification.R` | Passed: all real-data fits valid and tracking the data, mixing on all four series, replicated ablation (280 rows) and sensitivity (840 rows) complete (23 Oct) |
| P10 | `tests/gate_p10_verification.R` | Passed: all figures and tables present and newer than their gate results (23 Oct) |

**Deviations from the specification** are recorded in [`SPEC_CHANGES.md`](SPEC_CHANGES.md)
(NSP front-end, windowed shift probability, plug-in tau_omega, L_max = 50, minimum shift
persistence, Bitcoin log returns). Each is swept back to its spec value in the sensitivity grid.

**Main results** (details and Monte Carlo SEs in `output/p8_full_summary.csv` and `output/tab*.tex`)

- *Simulation (P8), unweighted mean over the 17 cells.* The proposed method has the highest
  point-outlier F1 (0.93) and macro-F1 (0.62; next best: copula comparator 0.45). Its
  change-point F1 (0.50) is second to H-SMUCE (0.54).
- *Real data (P9).* Well-log: 12 change points, 21 point outliers. Oilwell: 8 and 12.
  INDPRO: 0 change points, 7 point outliers, 1 collective anomaly. Bitcoin (log returns):
  1 change point, 155 point outliers, 3 collective anomalies.
- *Ablation (P9, 20 replicates each of the base cell and an n = 300 setting).* Removing the
  NSP or CAPA front-end, the dynamic outlier model, horseshoe+ or innovation-variance scoring changes
  change-point and point F1 by no more than about 0.08 (base). Only removing AR(1) has a large
  effect: change-point recall rises (base 0.47 to 0.66; n = 300 0.05 to 0.52), but precision
  falls (base 0.91 to 0.49) and false change points rise from about 2.6 to 6.7 per series.
- *Sensitivity (P9).* The spec values of the changed settings (d_omega = 1.5, per-point
  shift probability, sampled tau_omega) lower change-point recall on the base cell (F1 0.21,
  0.06 and 0.37, against 0.57 for the defaults). Sampled tau_omega is the only setting that
  raises recall at n = 300 (F1 0.28 against 0.06). Most other settings move F1 by less than 0.05.

**Known limitations**

1. *Convergence in the simulation study.* At the spec length (4 chains x 5,000 + 5,000) only
   17% of the 1,700 proposed fits reach gated R-hat < 1.01 (72% < 1.05, 85% < 1.1). The share
   ranges from 0-1% (phi = 0.9, t3 noise, the phi = 0.9 stress cell) to about 50% (jump 1.0,
   no collective anomalies). It is reported per cell in `output/gate_p8_results.rds`
   (`$convergence`). The real-data fits use longer chains and all mix.
2. *Conservative change-point detection.* Change-point precision is 0.8-1.0 and recall 0.2-0.5
   at alpha = 0.05. On short series dense in structure (n = 300) the AR(1) term absorbs level
   shifts (phi_hat about 0.75 for a true 0.5) and recall is about 0.05. This is a property of the
   posterior, not of mixing (SPEC_CHANGES.md, decisions of 21 October).
3. *Weakly identified outlier globals.* tau_psi, m_h and rho mix slowly on well-log (rho ESS
   about 55 per 200,000 draws); they are reported, not gated (C8).
4. *Collective anomalies* are rarely recovered by any method (proposed: collective F1 0.22 in P8).

---

## 5. Quickstart and Reproduction

### Running Unit Tests
To run all unit test files:
```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Full rerun
`targets` is not installed in the project library, so `run_all.sh` runs every step in
dependency order. The long steps are resumable: each P8 replicate and each P9 fit is saved as
its own file when it finishes, and a rerun skips finished work.
```bash
./run_all.sh                       # all steps 1-10
./run_all.sh 8 3 4 9 10            # exactly these steps, in this order
CPD_STOP_AT=18:30 ./run_all.sh 6   # stop starting new tasks at 18:30; rerun to resume
./run_all.sh stop                  # finish the tasks in progress, save and exit
./run_all.sh status                # what is saved and what is left
```
| Step | Content | Time (16 workers) |
|---|---|---|
| 1 | P3 gate (Stan reference, SBC) | ~40 min |
| 2 | Geweke Getting-it-Right test | ~6 min |
| 3 | Long well-log run (P4 mixing and speed) | ~10 min |
| 4 | P4 gate | ~20 min |
| 5 | P5, P6, P7 gates | ~20 min |
| 6 | Full simulation study (P8, 3,400 tasks) | 7.7 h |
| 7 | Real data at the converged length (skipped if `output/real_data_long.rds` exists; `REDO_REAL=1` reruns it) | ~45 min |
| 8 | Replicated ablation and sensitivity (P9, 1,120 fits) | 3.3 h |
| 9 | P8 and P9 gates | ~1 min |
| 10 | Tables and figures, P10 gate, unit tests | ~1 min |

Core use is capped by `CPD_MAX_CORES` (default 24) and the worker count by `N_WORKERS`
(default 16; throughput on the development machine does not improve beyond 16-20). BLAS is
single-threaded. Progress is logged to `output/run_all.log`.

`output/p8_full/` (3,400 files, read-only, backup `output/p8_full_backup_2026-10-22.tar.gz`)
and `output/p9_replicated/` (1,120 files) hold the completed studies; delete or move them to
rerun P8 or P9 from scratch.

### Figures and tables
`Rscript run_pipeline.R` (step 10) regenerates these from the saved gate results:
- **`output/fig1_real_data_four_panel.png`**: Four-panel empirical application plot (Well-log, Oilwell, US INDPRO, Bitcoin) with shaded 90% posterior credible bands for $\mu$ and bottom strip for $o_t$.
- **`output/fig2_sim_factor_effects.png`**: Multi-panel simulation performance comparison across sample size, autocorrelation, innovation distribution, jump size, and outlier regimes.
- **`output/fig3_identifiability_quadrants.png`**: Scatter plot in $(p_t^\omega, p_t^\psi)$ exceedance space separating sustained shifts from isolated point outliers.
- **`output/tab1_sim_cp_f1.tex`, `tab1_sim_pt_f1.tex`, `tab1_sim_macro_f1.tex`**: change-point, point-outlier and macro-F1 by cell and method, mean (Monte Carlo SE) over the full P8 replicates (written by the P8 gate).
- **`output/tab2_real_data.tex`**: feature detections and runtime across the four empirical series.
- **`output/tab3_ablation.tex`**: the seven ablation configurations, mean (MC SE) over 20 replicates of each setting.
- **`output/tab4_sensitivity.tex`**: the 21 sensitivity configurations, mean (MC SE) over 20 replicates of each setting.