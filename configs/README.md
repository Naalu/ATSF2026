# Model Configurations

Each file in this directory defines a complete model configuration for the
FluSight ILI Sandbox Hub forecast engine (`scripts/forecast_engine.R`).

## Quick start

```r
source("scripts/forecast_engine.R")
source("configs/snaive_bc_bs.R")   # loads `config`
result <- run_forecast(config)      # setup + optional forecast loop
run_diagnostics(result)             # visual QA plots
```

## Current models

| Config file | Model ID | Description |
|---|---|---|
| `snaive_bc_bs.R` | `KReger-snaive_bc_bs` | Seasonal Naive + Box-Cox (guerrero) + bootstrapped residuals |
| `ets_log.R` | `KReger-ets_log` | ETS (auto-selected) + log(x+0.001) + bootstrapped residuals |

## Adding a new model

1. Copy an existing config file (e.g., `cp snaive_bc_bs.R my_new_model.R`)
2. Change these fields:
   - `model_abbr` — unique identifier, ≤15 chars, alphanumeric + underscore
   - `model_fn` — the fable model formula (or `NULL` for non-fable models)
   - `transform` — `"box_cox_guerrero"` or `"none"`
   - `simulate_method` — `"bootstrap"` or `"parametric"`
   - `metadata` — update `model_name`, `methods`, and `methods_long`
3. Run it:
   ```r
   source("scripts/forecast_engine.R")
   source("configs/my_new_model.R")
   config$run_loop <- TRUE
   result <- run_forecast(config)
   ```

The engine automatically creates `model-output/KReger-my_new_model/` and
`model-metadata/KReger-my_new_model.yml`.

## Config fields reference

### Identity
| Field | Type | Description |
|---|---|---|
| `team_abbr` | string | Team name (≤15 chars, `^[a-zA-Z0-9_]+$`) |
| `model_abbr` | string | Model name (same constraints) |

### Paths
| Field | Type | Description |
|---|---|---|
| `hub_path` | string | Absolute path to local ATSF2026 repo clone |

### Simulation
| Field | Type | Description |
|---|---|---|
| `n_sims` | integer | Number of bootstrap/parametric replicates (1000 for production, 100 for dev) |
| `min_train_weeks` | integer | Minimum training data required (52 for seasonal models) |

### Run control
| Field | Type | Description |
|---|---|---|
| `dev_mode` | logical | `TRUE` = process first 5 dates only |
| `run_loop` | logical | `TRUE` = run full forecast loop on `run_forecast()` |
| `overwrite` | logical | `TRUE` = overwrite existing CSVs and metadata |
| `base_seed` | integer or NULL | Per-date seed: `set.seed(base_seed + as.integer(date))`. NULL = stochastic. |

### Model specification
| Field | Type | Description |
|---|---|---|
| `transform` | string | `"box_cox_guerrero"` or `"none"` |
| `model_fn` | function or NULL | `function(lambda) { ... }` returning a fable model spec. NULL for custom models. |
| `simulate_method` | string | `"bootstrap"`, `"parametric"`, or `"custom"` |
| `custom_simulate_fn` | function or NULL | Escape hatch for non-fable models. Must return tibble with `location`, `target_end_date`, `.sim`. |

### Metadata
| Field | Type | Description |
|---|---|---|
| `metadata$team_name` | string | Human-readable team name |
| `metadata$model_name` | string | Human-readable model name (used in plot titles) |
| `metadata$model_contributors` | list of lists | Each with `name`, `affiliation`, `email`, optional `orcid` |
| `metadata$license` | string | One of: CC0-1.0, CC-BY-4.0, CC-BY_SA-4.0, PPDL, ODC-by, ODbL, OGL-3.0 |
| `metadata$designated_model` | logical | Eligible for hub ensemble? (max 2 per team) |
| `metadata$data_inputs` | string | Data sources used |
| `metadata$methods` | string | Brief description (≤200 characters) |
| `metadata$methods_long` | string | Full methodology description |
| `metadata$ensemble_of_models` | logical | Is this an ensemble? |
| `metadata$ensemble_of_hub_models` | logical | Ensemble of other hub models? |

## Non-fable models

For models outside the fable ecosystem, set `model_fn = NULL` and provide
`custom_simulate_fn`:

```r
config$model_fn <- NULL
config$simulate_method <- "custom"
config$custom_simulate_fn <- function(train_data, h, n_sims, preprocess_result) {
  # train_data: tsibble filtered to <= origin_date
  # h: number of horizons (4)
  # n_sims: number of simulation paths
  # preprocess_result: list (e.g., lambdas if Box-Cox)
  #
  # Must return a tibble with columns:
  #   location (chr), target_end_date (Date), .sim (dbl)
}
```

The engine handles everything downstream: quantile computation, hub formatting,
validation, CSV writing, metadata generation, and diagnostics.