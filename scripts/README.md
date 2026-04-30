# Forecasting Scripts

This directory contains the FluCast forecasting pipeline scripts. The
high-level project overview is in [`../README.md`](../README.md).

## Quick navigation

If you want to **reproduce the project**, start with one of the driver scripts:

| Driver | Runtime | Purpose |
|---|---|---|
| `reproduce_smart.R` | minutes (recommended) | Skips phases whose outputs already exist |
| `reproduce_full.R` | 6–12 hours | Top-to-bottom from raw data |
| `reproduce_final.R` | ~30 seconds | Final assembly only |

```bash
Rscript scripts/reproduce_smart.R   # run from project root
```

If you want to **add or run an individual model**, see the engine
architecture section below.

If you're **looking for a specific phase's script**, see the phase
index below.

---

## Phase index

The work proceeds through five phases. Scripts are organized below
by phase; each script's own header block explains its inputs, outputs,
and how to run it.

### Phase 1 — Baseline submission

Generates the assigned-model baseline (`KReger-snaive_bc_bs`).

| Script | Purpose |
|---|---|
| `forecast_engine.R` | Generic, config-driven forecast runner (used by all KReger models) |
| `forecast_helpers.R` | Shared utilities: data loading, hub formatting, validation |
| `run_validation.R` | Run a config's validation forecasts (Seasons 1–2) |
| `run_test.R` | Run a config's test forecasts (Seasons 3–5) |
| `cdf_to_quantiles.R` | Helper: instructor-supplied CDF→quantile conversion |
| `get_step_ahead_model_output.R` | Helper: instructor-supplied horizon-N extraction |

### Phase 2 — Candidate generation, scoring, and pruning

Generates 11 KReger candidate models and prunes for diversity.

| Script | Purpose |
|---|---|
| `run_all_candidates.R` | Batch-run all 11 candidate model configs |
| `run_bsts_validation.R` | bsts validation forecasts (separate due to ~4.6h MCMC runtime) |
| `analyze_candidates.R` | Compute correlations, regimes, model summaries |
| `regime_helpers.R` | Helper: 4-step regime classification cascade |
| `archive/score_candidates.R` | Score candidates with scoringutils (Phase 2 → 8-model pool) |
| `archive/twin_diagnostics.R` | Diversity-prune to primary 8-model pool |

Outputs: `analysis/phase2/regimes.csv`, `analysis/phase2/wis_summary.csv`,
`analysis/phase2/decisions.md` (decision log).

### Phase 3 — Primary pool BMA ensemble

| Script | Purpose |
|---|---|
| `archive/compute_bma_weights.R` | Three weighting methods on primary pool (BMA, inverse-WIS, stacking) |
| `archive/generate_ensemble.R` | Linear-pool ensemble using primary BMA weights |

Outputs: `analysis/phase3/weights_log_score_bma.csv`,
`analysis/phase3/staging/ensemble_bma_primary/`.

### Phase 4 — Expanded pool (with classmate models + delphi-epicast)

Adds 9 hub-submitted models, re-prunes to 10, generates all four
ensembles for the comparison matrix that drives the report's central
finding.

| Script | Purpose |
|---|---|
| `score_candidates_expanded.R` | Score the full 20-candidate expanded pool |
| `twin_diagnostics_expanded.R` | Diversity-prune to expanded 10-model pool |
| `compute_bma_weights_expanded.R` | BMA weights on expanded pool |
| `generate_ensembles_matrix.R` | Generate all four ensembles + LOO-CV comparison |

Outputs: `analysis/phase4/shortlist_expanded.csv`,
`analysis/phase4/weights_log_score_bma_expanded.csv`,
`analysis/phase4/matrix_loo_comparison.csv` (the headline 4-cell result),
`analysis/phase4/staging/ensemble_*/` (four directories, 134 CSVs each).

### Phase 5 — Submission and validation diagnostics

| Script | Purpose |
|---|---|
| `promote_ensemble.R` | Copy BMA-expanded staging to model-output/, write metadata YAML |
| `score_validation_coverage.R` | Compute 50% / 95% PI coverage on validation (for report §3 PIC) |

Outputs: `model-output/KReger-bma_ensemble/` (134 CSVs),
`model-metadata/KReger-bma_ensemble.yml`,
`analysis/phase4/validation_coverage.csv`.

### Auxiliary

| Script | Purpose |
|---|---|
| `presentation_figures.R` | Figures for the in-class presentation (not part of the submission pipeline) |

### Archive

The `archive/` subdirectory contains earlier-phase scripts that were
superseded by their Phase 4 expanded equivalents. They remain in the
pipeline because Phase 3 produces the primary pool BMA ensemble used
in the four-cell comparison; archiving signals that they are stable
historical scripts rather than active development.

---

## Engine architecture

The KReger forecast pipeline is **model-agnostic**. Model-specific
configuration lives in `../configs/` — these scripts never need to be
edited to add a new model.

```
configs/snaive_bc_bs.R ──┐
configs/ets_log.R ───────┤
configs/your_model.R ────┴──> forecast_engine.R ──> forecast_helpers.R
                                    │
                                    ├── model-output/{team}-{model}/*.csv
                                    └── model-metadata/{team}-{model}.yml
```

### `forecast_engine.R`

The main entry point. Provides two functions:

- **`run_forecast(config)`** — loads data, preprocesses (e.g., lambda
  estimation), writes metadata YAML, and optionally runs the full forecast
  loop across all origin dates. Returns a result object with an ad-hoc
  `result$forecast(date)` function for single-date forecasts.

- **`run_diagnostics(result)`** — generates visual QA plots from a completed
  forecast run: fan charts, sanity dashboard (median vs actuals + PI width),
  and transform parameter summaries.

Key features handled by the engine:

- Per-date deterministic seeding for reproducibility
- Skip-if-exists for resuming interrupted runs
- Overwrite protection (configurable)
- Per-date error handling (failures logged, loop continues)
- Progress reporting with ETA
- Post-loop spot-check validation
- Hub naming constraint enforcement on startup
- Metadata YAML generation from config fields

### `forecast_helpers.R`

Shared utility functions and constants sourced automatically by the
engine. These implement the hub submission contract — formatting,
validation, and file I/O that every model needs regardless of its
forecasting method.

**Constants:**
`DATE_COL`, `QUANTILE_LEVELS`, `HUB_LOCATIONS`, `HUB_HORIZONS`,
`EXPECTED_ROWS`, `EXPECTED_COLS`

**Functions:**

| Function | Purpose |
|---|---|
| `load_and_prepare_tsibble()` | Load target CSV, convert to tsibble |
| `get_origin_dates()` | Parse origin dates from tasks.json |
| `compute_quantiles_from_sims()` | Empirical quantiles from simulation paths |
| `format_for_hub()` | Shape data into the 8-column hub CSV schema |
| `validate_forecast_df()` | 15 programmatic quality checks per forecast |
| `write_hub_csv()` | Write CSV with hub filename convention |
| `plot_quantile_fan()` | Single-location fan chart for spot-checking |

### Usage (single model)

```r
# Run a model end-to-end
source("scripts/forecast_engine.R")
source("configs/snaive_bc_bs.R")
config$run_loop <- TRUE
result <- run_forecast(config)

# Visual QA
run_diagnostics(result)

# Ad-hoc single-date forecast (uses deterministic seed)
df <- result$forecast(as.Date("2017-01-14"), write_file = FALSE)
```

---

## Common conventions

- All scripts assume the working directory is the project root
- Phase artifacts go to `analysis/phaseN/`, never to `scripts/` itself
- Roxygen-style comments document non-trivial functions
- Heavy comments explain non-obvious design choices inline
- Scripts use `getwd()`-based path resolution; no hardcoded absolute
  paths (so the project clones cleanly to any location)

## Dependencies

Loaded automatically by the engine:

- `fpp3` (includes tsibble, fable, fabletools, feasts)
- `tidyverse`

Phase-specific scripts additionally use:

- `scoringutils` (>=2.0) — Phase 2/4 WIS scoring
- `bsts` — Phase 2 bsts candidate
- `hubValidations` — submission validation
- `patchwork` — diagnostic plots
