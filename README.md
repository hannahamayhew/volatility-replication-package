# Replication Package
## "Has the Great Moderation Survived? Volatility Break Evidence from FRED-MD, 1959–2026"

---

## Overview

This package reproduces all tables and figures in the paper. The analysis applies the Sensier and van Dijk (2004) volatility-break methodology to 122–123 FRED-MD series across three sample windows: Historical (1959:03–1999:12), Pre-GFC (1959:03–2007:12), and Full Sample (1959:03–2026:01).

**Pre-computed outputs are included** in `output/tables/` and `output/figures/`. You do not need to re-run the code to inspect the results. Instructions for re-running from scratch are in Section 3 below.

---

## 1. Requirements

- **MATLAB** (any version from R2019b onward; no toolboxes required beyond base MATLAB)
- The FRED-MD February 2026 vintage, `2026-02-MD.csv`, included in `data_raw/`

---

## 2. Package Structure

```
volatility_replicationpackage_final/
├── README.md                       ← this file
├── master_run.m                    ← runs the full pipeline end-to-end
├── config.m                        ← paths and estimation parameters
│
├── data_raw/
│   └── 2026-02-MD.csv              ← FRED-MD February 2026 vintage (included)
│
├── helpers/                        ← shared econometric functions
│   ├── uv_supw_test.m
│   ├── uv_onebreak_search.m
│   ├── uv_recursive_breaks.m
│   ├── cm_supw_test.m
│   ├── cm_onebreak_search.m
│   ├── cv_supw_test.m
│   ├── cv_recursive_breaks.m
│   ├── select_ar_lag_bic.m
│   ├── hansen_supf_pval.m
│   ├── hansen_supf_pval_general.m
│   ├── smooth_break_search_mean.m
│   ├── logistic_transition.m
│   ├── prepare_missing.m
│   ├── remove_outliers.m
│   ├── factors_em.m
│   └── [other utilities]
│
├── src/
│   ├── data_prep/
│   │   ├── run_data_prep.m         ← loads, transforms, and saves the panel
│   │   └── run_stationarity_tests.m
│   ├── uv/
│   │   ├── build_uv_object.m       ← constructs z_t = sqrt(pi/2)|y_t - mean(y)|
│   │   ├── run_uv_historical.m     ← UV one-break test, Historical window
│   │   ├── run_uv_2007.m           ← UV one-break test, Pre-GFC window
│   │   ├── run_uv_fullsample.m     ← UV one-break test, Full sample
│   │   └── run_uv_smooth.m         ← UV smooth-break (logistic) estimation
│   ├── cm_cv/
│   │   ├── run_cm_*.m              ← AR(p) conditional-mean break models
│   │   ├── run_cv_*.m              ← CV break tests (uses CM residuals)
│   │   ├── run_mb_fixed74.m        ← recursive multi-break estimation (main spec)
│   │   ├── run_mb_fixed74_smooth_char.m  ← smooth-break characterization for MB
│   │   ├── run_mb_robustness_samples.m   ← balanced panel + no-flag robustness
│   │   └── run_mb_lm_stationarity_robust.m  ← Labor Market stationarity robustness
│   └── tables_figures/
│       ├── run_tables_figures.m    ← merges results, produces summary CSVs
│       └── produce_paper_figures.m ← generates all paper figures (Figs 1–4b, A1–A7)
│
├── data_intermediate/              ← .mat files saved between stages (auto-generated)
│
└── output/
    ├── tables/                     ← all results as .csv files (pre-computed)
    └── figures/                    ← all paper figures as .png and .pdf (pre-computed)
```

---

## 3. How to Reproduce

### Option A — Inspect pre-computed results (no re-run needed)

All tables and figures from the paper are already in:
- `output/tables/` — CSV files for every analysis stage
- `output/figures/` — PNG and PDF versions of all paper figures

### Option B — Re-run the full pipeline from scratch

1. Open MATLAB
2. Navigate to this folder (the one containing `master_run.m`)
3. Run:

```matlab
master_run
```

That is the only command needed. The script runs all stages in order and writes results to `output/tables/` and `output/figures/`. Total runtime is approximately 30–60 minutes depending on hardware.

### Option C — Re-run individual stages

```matlab
run('config.m')                                          % must run first every session
run('src/data_prep/run_data_prep.m')                     % Stage 1
run('src/uv/run_uv_historical.m')                        % Stage 2a
run('src/uv/run_uv_2007.m')                              % Stage 2b
run('src/uv/run_uv_fullsample.m')                        % Stage 2c
run('src/cm_cv/run_cm_historical.m')                     % Stage 3a
run('src/cm_cv/run_cv_historical.m')                     % Stage 3b (needs 3a)
run('src/cm_cv/run_cm_2007.m')                           % Stage 3d
run('src/cm_cv/run_cv_2007.m')                           % Stage 3e (needs 3d)
run('src/cm_cv/run_cm_fullsample.m')                     % Stage 3g
run('src/cm_cv/run_cv_fullsample.m')                     % Stage 3h (needs 3g)
run('src/cm_cv/run_mb_fixed74.m')                        % Stage 3j — main multi-break (needs all CM/CV)
run('src/uv/run_uv_smooth.m')                            % Stage 3k
run('src/cm_cv/run_cv_smooth.m')                         % Stage 3l
run('src/tables_figures/run_tables_figures.m')           % Stage 4a
run('src/tables_figures/produce_paper_figures.m')        % Stage 4b — paper figures
```

---

## 4. Paper–Output Correspondence

### Tables

| Paper table | Source file(s) |
|-------------|---------------|
| Table 1 (Sample Windows) | Descriptive; no computation |
| Table 2 (One-Break Results Summary) | `output/tables/uv_hist_results.csv`, `uv_2007_results.csv`, `uv_full_results.csv`, `cv_hist_results.csv`, `cv_2007_results.csv`, `cv_full_results.csv`, `cm_hist_results.csv` |
| Table 3 (Single-Break Abruptness by Sector) | `output/tables/mb_fixed74_smooth_char_results.csv` |
| Table 4 (Recursive MB Distribution) | `output/tables/mb_fixed74_results.csv` |
| Table 5 (Multiple-Break Abruptness by Sector) | `output/tables/mb_fixed74_smooth_char_results.csv` |
| Table 6 (Sector-Level Recursive Breaks) | `output/tables/mb_fixed74_group_summary.csv` |

### Figures

| Paper figure | Output file |
|-------------|-------------|
| Figure 1a (UV break magnitude, Real vs. Nominal) | `output/figures/figure4a_magnitude_uv_real_nominal.*` |
| Figure 1b (CV break magnitude, Real vs. Nominal) | `output/figures/figure4b_magnitude_cv_real_nominal.*` |
| Figure A1 (Output and Income) | `output/figures/figureA1_output_income_break_magnitude.*` |
| Figure A2 (Labor Market) | `output/figures/figureA2_labor_market_break_magnitude.*` |
| Figure A3 (Orders and Inventories) | `output/figures/figureA3_orders_inventories_break_magnitude.*` |
| Figure A4 (Housing/Construction) | `output/figures/figureA4_housing_construction_break_magnitude.*` |
| Figure A5 (Prices) | `output/figures/figureA5_prices.*` |
| Figure A6 (Money and Credit) | `output/figures/figureA6_money_credit_break_magnitude.*` |
| Figure A7 (Interest and Exchange Rates) | `output/figures/figureA7_interest_fx_break_magnitude.*` |

Each figure is saved in both `.png` (300 dpi) and `.pdf` (vector) formats.

---

## 5. Key Estimation Parameters

All parameters are set in `config.m`. The values below match the paper exactly:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `TRIM_FRAC` | 0.15 | 15% trimming for break-date candidates (Sensier & van Dijk 2004) |
| `MAX_BREAKS` | 5 | Maximum recursive breaks per series |
| `MIN_SEG_LEN` | 60 | Hard floor on segment length (months) |
| `MIN_SEG_EXT` | 74 | Fixed min-segment for Pre-GFC and Full Sample windows |
| `PMAX` | 12 | Maximum AR lag for BIC selection |
| Significance level | 5% | Hansen (1997) asymptotic p-values throughout |

The historical window uses a sample-relative 15% × T minimum segment rule, which equals 74 months for T = 490. The fixed 74-month rule is carried over to the Pre-GFC and Full Sample windows for comparability.

---

## 6. Data

**Source:** McCracken, M.W. and Ng, S. (2016). "FRED-MD: A Monthly Database for Macroeconomic Research." *Journal of Business & Economic Statistics*, 34(4), 574–589.

**Vintage:** February 2026 (`2026-02-MD.csv`), included in `data_raw/`. The same vintage can be downloaded from: https://research.stlouisfed.org/econ/mccracken/fred-databases/

**Coverage:** 126 raw series, January 1959 – January 2026. After McCracken-Ng stationarity transformations, the panel begins 1959:03. Series with fewer than 80% non-missing observations are excluded, yielding 122–123 series depending on window.

---

## 7. References

- Sensier, M. and van Dijk, D. (2004). "Testing for Volatility Changes in U.S. Macroeconomic Time Series." *Review of Economics and Statistics*, 86(3), 833–839.
- Bai, J. (1997). "Estimating Multiple Breaks One at a Time." *Econometric Theory*, 13(3), 315–352.
- Bai, J. and Perron, P. (1998). "Estimating and Testing Linear Models with Multiple Structural Changes." *Econometrica*, 66(1), 47–78.
- Hansen, B.E. (1997). "Approximate Asymptotic P Values for Structural-Change Tests." *Journal of Business & Economic Statistics*, 15(1), 60–67.
- McCracken, M.W. and Ng, S. (2016). "FRED-MD: A Monthly Database for Macroeconomic Research." *Journal of Business & Economic Statistics*, 34(4), 574–589.
