# Data Provenance and Environment Specifications

This directory contains raw datasets and records package version specifications for the project:
*Significance-Guided Global–Local Shrinkage for Distinguishing Change Points, Collective Anomalies and Point Outliers in Autocorrelated Time Series* (Das & Kunnummal, 2026).

---

## 1. Raw Datasets (`data/raw/`)

| File | Series | Source / Provenance | Size / Observations | SHA-256 Checksum |
|---|---|---|---|---|
| `welllog.csv` | Nuclear magnetic resonance well-log | Standard benchmark (`changepoint.influence::welldata` raw, uncleaned) | 4,050 points | `d53e53bf227f0504e119a51b411c5847a331055f25b89f0739708db11f03293f` |
| `oilwell.csv` | Oilwell pressure sensor | JISPS benchmark (Das & Kunnummal, 2026) | 4,000 points | `eaffa35b1c625ca082429913c37c6f90d8b202b9223bebec680a6d8574e34769` |
| `us_indpro.csv` | US Industrial Production Index (`INDPRO`) | Federal Reserve Economic Data (FRED) monthly log growth | 312 observations (2000-01 to 2026) | `3cc36f6f026b54c7139f9b8c295e35a4d3874590bc0251036cc1dd11c0ab2cb6` |
| `btc_usd.csv` | Bitcoin / USD daily log price | CoinGecko / Thesis benchmark (1 Jan 2018 – 31 Dec 2023) | 2,191 points | `ec16ea2cdb5eea8cffcb49e41373a2bd6e3498078e7950f7142798956fd1c157` |

All files in `data/raw/` are frozen and read-only.

---

## 2. Package Versions and Function Signatures (Audited)

As required in Section 2, the function signatures of key dynamic shrinkage and front-end packages were audited against installed versions:

### `nsp` (Version 1.0.0)
- **Signature:** `nsp_poly_selfnorm(y, M = 1000, thresh.val = NULL, power = 1/2, min.size = 20, alpha = 0.1, deg = 0, eps = 0.03, c = exp(1 + 2 * eps), overlap = FALSE)`
- **Usage:** Set `deg = 0` for piecewise-constant mean change point detection under autocorrelated noise.

### `anomaly` (Version 4.3.3)
- **Signature:** `capa(x, beta, beta_tilde, type = "meanvar", min_seg_len = 10, max_seg_len = Inf, max_lag = 0)`
- **Accessors:** `point_anomalies(res)`, `collective_anomalies(res)`.
- **Usage:** Requires prewhitening for autocorrelated series ($\tilde y_t = y_t - \hat\phi y_{t-1}$).

### `dsp` (Version 1.7.0)
- **Functions:** `sampleDSP`, `sampleFastGaussian`, `initChol_spam`, `sampleAR1`, `initDHS`, `abco`.
