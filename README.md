# FluCast: A Regime-Conditional BMA Ensemble for ILI Forecasting

**Karl Reger** • INF 599 Final Project • Northern Arizona University

This repository contains the FluCast forecasting project: a
regime-conditional Bayesian Model Averaging (BMA) ensemble for
weighted influenza-like-illness (wILI) forecasts on the FluSight ILI
Sandbox Hub. Submitted as `KReger-bma_ensemble`.

## Headline result

On a 4-cell validation comparison (BMA vs equal-weight × primary
8-model pool vs expanded 10-model pool), the BMA-expanded ensemble
won by every margin that mattered:

| Cell | Pool | Method | LOO total WIS | Mean WIS |
|---|---|---|---:|---:|
| **BMA-expanded**    | expanded (10) | BMA   | **614.0** | **0.245** |
| BMA-primary         | primary (8)   | BMA   | 630.3 | 0.251 |
| Equal-primary       | primary (8)   | Equal | 732.3 | 0.292 |
| Equal-expanded      | expanded (10) | Equal | 750.6 | 0.299 |

*Validation: leave-one-date-out cross-validation across 2,508 forecast
cells (57 origin dates × 11 locations × 4 horizons + quantile
aggregation). Test-period evaluation reserved for the instructor.*

The BMA-expanded ensemble beats equal-weight aggregation by 14% on
the primary pool and 18% on the expanded pool, leads in every regime
(growing, declining, peaking), and improves on its best individual
component (KReger-nnetar_bc_bs at validation mean WIS 0.284) by
13.7%.

The full report is `docs/final_report.pdf` (3 pages + references).

## Repository structure

```
.
├── README.md                       (this file)
├── docs/
│   ├── final_report.tex            LaTeX source for the final report
│   ├── final_report.pdf            Compiled report
│   └── weights_heatmap.pdf         Figure 1 of the report
├── scripts/                        Pipeline scripts (see scripts/README.md)
│   ├── reproduce_full.R            Driver: full from-scratch reproduction
│   ├── reproduce_smart.R           Driver: skip phases with existing outputs
│   ├── reproduce_final.R           Driver: final assembly only
│   ├── forecast_helpers.R          Phase 1 helpers
│   ├── run_validation.R            Phase 1: baseline validation
│   ├── run_test.R                  Phase 1: baseline test forecasts
│   ├── run_all_candidates.R        Phase 2: candidate generation
│   ├── run_bsts_validation.R       Phase 2: bsts (long-running)
│   ├── score_candidates_expanded.R Phase 4: score 20-model pool
│   ├── twin_diagnostics_expanded.R Phase 4: diversity pruning
│   ├── compute_bma_weights_expanded.R   Phase 4: BMA weight estimation
│   ├── generate_ensembles_matrix.R Phase 4: 4-cell comparison
│   ├── promote_ensemble.R          Phase 5: promotion to model-output/
│   ├── score_validation_coverage.R Phase 5: validation PI coverage
│   ├── regime_helpers.R            Phase 2 helper: regime classifier
│   └── archive/                    Superseded earlier-phase scripts
├── analysis/                       Intermediate artifacts
│   ├── phase2/decisions.md         Phase 2 decision log
│   ├── phase4/matrix_loo_comparison.csv     Validation 4-cell matrix
│   ├── phase4/weights_log_score_bma_expanded.csv  BMA weights (frozen)
│   ├── phase4/validation_coverage.csv       Validation PI coverage
│   └── phase4/staging/             Per-ensemble staged forecasts
├── model-output/                   Submitted forecasts (mirrors hub)
│   ├── KReger-bma_ensemble/        ← The designated submission (134 CSVs)
│   ├── KReger-snaive_bc_bs/        Phase 1 baseline
│   ├── KReger-{arima,ets,nnetar,...}/   Component models
│   └── (other classmate models for the expanded pool)
├── model-metadata/
│   └── KReger-bma_ensemble.yml     Submission metadata
├── target-data/
│   ├── oracle-output.csv           Truth data (FluSight hub format)
│   └── time-series.csv             wILI time series
└── hub-config/
    ├── tasks.json                  Hub-defined origin dates and quantile levels
    └── (other hub config)
```

## How to reproduce

Three reproduction paths are supported, each calling the same underlying
phase scripts but at different scope and runtime:

### 1. Full reproduction (`reproduce_full.R`)

Runs the entire pipeline from raw data to the final submission. Use
this to verify reproducibility from inputs alone with no shortcuts.

```bash
Rscript scripts/reproduce_full.R
```

Runtime: **6–12 hours** (dominated by candidate scoring and bsts MCMC).
Phases run unconditionally regardless of whether outputs already exist.

### 2. Smart reproduction (`reproduce_smart.R`) — recommended

Same pipeline, but skips phases whose sentinel output files already
exist. If you've cloned this repo (which has all intermediate artifacts
checked in), this completes in minutes.

```bash
Rscript scripts/reproduce_smart.R
```

Runtime: **minutes if intermediates exist; same as full if cold**.

To force a single phase to re-run, delete its sentinel output (the
script reports each phase's sentinel) and re-run the driver.

### 3. Final assembly (`reproduce_final.R`)

Skips directly to the final promotion step using existing intermediate
artifacts. Useful for spot-checking that intermediate artifacts produce
the expected submission CSVs.

```bash
Rscript scripts/reproduce_final.R
```

Runtime: **~30 seconds**.

## Requirements

- **R ≥ 4.3** with packages: `tidyverse`, `fpp3`, `scoringutils ≥ 2.x`,
  `bsts`, `hubValidations`, `fable`, `fabletools`, `tsibble`,
  `purrr`, `readr`
- Working directory: project root

To install:

```r
install.packages(c("tidyverse", "fpp3", "scoringutils", "bsts",
                   "hubValidations", "fable", "fabletools", "tsibble"))
```

## Pipeline overview

Five phases produce the submission:

1. **Baseline** — Generate the assigned-model baseline (KReger-snaive_bc_bs).
2. **Candidates** — Generate 11 KReger candidate models, score, prune for diversity into a primary 8-model pool.
3. **Primary BMA** — Compute regime-conditional BMA weights and ensemble.
4. **Expanded pool** — Add 9 hub-submitted models, re-prune to 10, generate all four ensembles (BMA × {primary, expanded} and Equal × {primary, expanded}) for comparison.
5. **Submission** — Promote the winning BMA-expanded ensemble to `model-output/` and validate.

Detailed per-script documentation is in [`scripts/README.md`](scripts/README.md).
The Phase 2 decision log is at [`analysis/phase2/decisions.md`](analysis/phase2/decisions.md).

## Key methodological choices

- **Diversity pruning via hierarchical clustering** of the validation-period point-forecast error correlation matrix (complete linkage at threshold 1−0.93)
- **Regime-conditional BMA weighting** with three regimes (growing, declining, peaking) classified via a 4-step priority cascade
- **Softmax of negative WIS at τ=0.050** (LOO-tuned), preferred over inverse-WIS and stacking
- **Validation-only frozen weights** to prevent test-period leakage
- **No-floor weight allocation** — the softmax discounts weak components without zeroing them out

The full methodology is in `docs/final_report.pdf`.

## Author

**Karl Reger** • [kcr28@nau.edu](mailto:kcr28@nau.edu)
INF 599, Spring 2026
Northern Arizona University

Submission PR: [sjfox/ATSF2026#18](https://github.com/sjfox/ATSF2026/pull/18)

## License

Coursework submitted for INF 599 at Northern Arizona University.
