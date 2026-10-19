# Significance-Guided Global–Local Shrinkage for Autocorrelated Time Series

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

> **Status (October 2026): work in progress.** The sampler was corrected after an audit
> (branch `fix/audit-2026-10`); acceptance gates P3-P10 must be re-run and several design
> decisions are open (see Section 4). No results in `output/` should be cited.

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
2. **Joint Banded Gaussian Updates:** An interleaved parameterization stacks $x = (\mu_1, \psi_1, \mu_2, \psi_2, \dots)^\top$, yielding a penta-diagonal precision matrix $Q$ drawn in $O(n)$ time via sparse banded Cholesky factorizations (`spam`).
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
│   ├── 32_gibbs_updates.R      # Conditional full conditional updates (Banded Cholesky, slice)
│   ├── 40_scoring.R            # Outlier score o_t (Eq. 8) and tail exceedance p_t (Eq. 9)
│   ├── 41_bfdr.R               # Bayesian False Discovery Rate (BFDR) selection rule
│   ├── 42_classify.R           # Sequential four-class classification and cpd_fit builder
│   ├── 50_copula.R             # Copula Markov comparator (AIC family selection, bootstrap)
│   ├── 60_baselines.R          # Baseline wrappers under the unified cpd_fit schema (FPNN is a placeholder)
│   ├── 70_metrics.R            # Change point F1, covering, Hausdorff, PR-AUC, macro-F1
│   ├── 80_sim_design.R         # 17-cell simulation grid, pilot runner, MC SE reporting
│   ├── 85_real_data.R          # Loaders and runners for 4 empirical series
│   ├── 86_sensitivity_ablation.R # 6 component ablations and 12-cell sensitivity grid
│   └── 90_plots_tables.R       # Publication ggplot2 figures and LaTeX kable tables
├── stan/
│   └── hsplus_ar1_ref.stan     # Reference Stan implementation with non-centred priors
├── data/
│   ├── raw/                    # Well-log, Oilwell, US Industrial Production, Bitcoin
│   └── processed/              # Formatted benchmarks
├── tests/
│   ├── testthat/               # Unit tests, incl. exactness checks of each Gibbs update
│   ├── gate_p3_verification.R  # Reference Stan vs Gibbs verification gate
│   ├── gate_p4_verification.R  # Production Gibbs convergence and speed gate
│   ├── gate_p5_verification.R  # BFDR error control and spike identifiability gate
│   ├── gate_p6_verification.R  # Copula Markov comparator calibration gate
│   ├── gate_p7_verification.R  # 12 Baseline wrappers validation gate
│   ├── gate_p8_verification.R  # 17-cell Simulation pilot and budget gate
│   ├── gate_p9_verification.R  # Empirical applications, ablation, sensitivity gate
│   └── gate_p10_verification.R # Publication outputs and LaTeX table gate
├── output/                     # Generated figures (PNG) and LaTeX tables (.tex)
├── _targets.R                  # Full reproducible workflow pipeline
├── run_pipeline.R              # Standalone master execution script
└── DESCRIPTION                 # R package metadata
```

---

## 4. Verification Gates: Current Status

Each phase in `CPD2_Code_Workflow.md` has an acceptance gate. Gate scripts save their
results to `output/gate_pX_results.rds` and stop with an error if a criterion fails.

| Phase | Gate script | Status |
|---|---|---|
| P0 | `renv::status()` | `robseg` now locked; `SBC` is still recorded from an unknown source |
| P1 | `tests/testthat/test-10-dgp.R` | Unit tests pass |
| P2 | `tests/gate_p2_verification.R` | Passed for the pipeline front-end `nsp_poly_ar` (false-interval rate <= 0.015 at phi = 0, 0.6, 0.9); the spec's `nsp_poly_selfnorm` fails (0.38, 1.00) - see SPEC_CHANGES.md C1 |
| P3 | `tests/gate_p3_verification.R` | SBC now simulates from the fitted model; never run to completion before (no saved results); re-run needed |
| P4 | `tests/gate_p4_verification.R` | Spec thresholds restored (R-hat < 1.01, ESS > 400, pipeline configuration). Global scales still mix slowly (see below). Getting-it-Right test not implemented |
| P5 | `tests/gate_p5_verification.R` | Now requires power > 0.2 so that FDR = 0 cannot pass vacuously; re-run needed |
| P6 | `tests/gate_p6_verification.R` | Passed (copula family recovery and null calibration) |
| P7 | `tests/gate_p7_verification.R` | Re-run needed after baseline fixes (robseg, PELT/CAPA scaling, ch4) |
| P8 | `tests/gate_p8_verification.R` | Only a 5-replicate pilot exists; full study not run |
| P9 | `tests/gate_p9_verification.R` | Now checks that fitted levels track the data; re-run needed |
| P10 | `tests/gate_p10_verification.R` | Checks files only; `tar_make()` does not yet run the simulation study |

**Deviations from the specification** are recorded in [`SPEC_CHANGES.md`](SPEC_CHANGES.md)
(NSP front-end, windowed shift probability, plug-in tau_omega, L_max = 50, minimum shift
persistence, Bitcoin log returns). Each is swept back to its spec value in the sensitivity grid.

**Known open issues**

1. *Mixing* of the outlier global scale and rho in the dynamic outlier model.
2. *FPNN* is a placeholder (`fpnn_placeholder`); the JISPS implementation is not yet ported.
3. The Geweke (2004) joint-distribution test for P4 is not implemented.

---

## 5. Quickstart and Reproduction

### Running Unit Tests
To run all unit test files:
```bash
Rscript -e 'testthat::test_dir("tests/testthat")'
```

### Reproducing All Publication Artifacts
To execute the master pipeline and generate all publication-grade figures and LaTeX tables:
```bash
Rscript run_pipeline.R
```

This regenerates the figures and tables below from the saved gate results. They are
**not** current: they were produced before the audit fixes and must be regenerated after
the gates are re-run.
- **`output/fig1_real_data_four_panel.png`**: Four-panel empirical application plot (Well-log, Oilwell, US INDPRO, Bitcoin) with shaded 90% posterior credible bands for $\mu$ and bottom strip for $o_t$.
- **`output/fig2_sim_factor_effects.png`**: Multi-panel simulation performance comparison across sample size, autocorrelation, innovation distribution, jump size, and outlier regimes.
- **`output/fig3_identifiability_quadrants.png`**: Scatter plot in $(p_t^\omega, p_t^\psi)$ exceedance space demonstrating complete decoupling of sustained shifts from isolated point outliers.
- **`output/tab1_sim_summary.tex`**: LaTeX table of simulation benchmark metrics (CP F1, covering, Point F1, Collective F1, Macro-F1, runtime) with Monte Carlo standard errors.
- **`output/tab2_real_data.tex`**: LaTeX table of feature detections and runtime across all four empirical series.
- **`output/tab3_ablation.tex`**: LaTeX table for the six component ablations.
- **`output/tab4_sensitivity.tex`**: LaTeX table for the twelve hyperparameter sensitivity configurations.
