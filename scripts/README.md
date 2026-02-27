# Forecasting Scripts

This directory contains the model-agnostic forecasting pipeline for the
FluSight ILI Sandbox Hub (ATSF2026). Model-specific configuration lives
in `configs/` — these scripts never need to be edited to add a new model.

## Architecture

```
configs/snaive_bc_bs.R ──┐
configs/ets_log.R ───────┤
configs/your_model.R ────┴──> forecast_engine.R ──> forecast_helpers.R
                                    │
                                    ├── model-output/{team}-{model}/*.csv
                                    └── model-metadata/{team}-{model}.yml
```

## Files

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

Shared utility functions and constants sourced automatically by the engine.
These implement the hub submission contract — formatting, validation, and
file I/O that every model needs regardless of its forecasting method.

**Constants:**
- `DATE_COL`, `QUANTILE_LEVELS`, `HUB_LOCATIONS`, `HUB_HORIZONS`
- `EXPECTED_ROWS`, `EXPECTED_COLS`

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

## Usage

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

## Other scripts in this directory

The following scripts are part of the original hub repository and are not
modified by our forecasting pipeline:

- `get_step_ahead_model_output.R` — instructor's CDF-to-quantile pipeline
  for the hist-avg and delphi-epicast models
- `cdf_to_quantiles.R` — helper function used by the above script

## Dependencies

Loaded automatically by the engine:
- `fpp3` (includes tsibble, fable, fabletools, feasts)
- `tidyverse`

Required only for `run_diagnostics()`:
- `patchwork`

Required only for hub validation (optional, not called by the engine):
- `hubValidations`