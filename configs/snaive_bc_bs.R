###############################################################################
# CONFIG: configs/snaive_bc_bs.R
#
# Model configuration for: SNAIVE + Box-Cox + Bootstrapped Residuals
# Engine:  scripts/forecast_engine.R
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/snaive_bc_bs.R")
#   result <- run_forecast(config)
#
# To create a NEW fable model, copy this file and change:
#   - model_abbr, model_name, methods, methods_long (identity + description)
#   - model_fn   (the fable formula)
#   - transform  (if not using Box-Cox)
#   See the examples at the bottom of this file.
#
# To create a NON-fable model (e.g., GAM, external API), set:
#   - model_fn = NULL
#   - custom_simulate_fn = your function
#   See the "Custom model" example at the bottom.
#
###############################################################################

config <- list(
  
  # ---- Identity ----
  # Hub naming rules: each must be <=15 chars, alphanumeric + underscore only.
  # The engine validates these on startup.
  team_abbr  = "KReger",
  model_abbr = "snaive_bc_bs",
  
  # ---- Paths ----
  # Root of the local ATSF2026 repo clone. All other paths derived from this.
  hub_path = "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026",
  
  # ---- Simulation settings ----
  # 1000 gives >=10 obs in each tail for the 0.01/0.99 quantiles.
  # Use 100 for fast dev/testing.
  n_sims = 1000,
  
  # Minimum weeks of training data required. SNAIVE with lag("year") needs
  # at least 52 weeks (one full seasonal cycle).
  min_train_weeks = 52,
  
  # ---- Run control ----
  dev_mode  = FALSE,  # TRUE = process first 5 dates only
  run_loop  = FALSE,  # TRUE = run full forecast loop on source()
  overwrite = FALSE,  # TRUE = overwrite existing CSVs; FALSE = skip them
  
  # ---- Reproducibility ----
  # Per-date deterministic seed: set.seed(base_seed + as.integer(origin_date))
  # Each date gets a unique, repeatable seed. Rerunning a single date always
  # produces the same simulations regardless of what other dates were processed.
  # Set to NULL to disable seeding (stochastic variation between runs).
  base_seed = 12345,
  
  # ---- Transform ----
  # How the engine preprocesses data before model fitting.
  #   "box_cox_guerrero" -- estimate per-location lambda via guerrero()
  #   "none"             -- no preprocessing; model_fn receives NA for lambda
  transform = "box_cox_guerrero",
  
  # ---- Model specification (fable path) ----
  # A function that receives a lambda value and returns a fable model
  # definition (the thing that goes inside model()). The engine calls this
  # once per location, passing the location-specific lambda.
  #
  # For models without Box-Cox, the lambda argument is NA -- just ignore it.
  model_fn = function(lambda) {
    SNAIVE(box_cox(observation, lambda) ~ lag("year"))
  },
  
  # ---- Simulation method ----
  # "bootstrap"  -- generate(bootstrap = TRUE): resample fitted residuals
  # "parametric" -- generate(bootstrap = FALSE): assume normal residuals
  simulate_method = "bootstrap",
  
  # ---- Custom simulation function (non-fable escape hatch) ----
  # Set model_fn = NULL and provide this instead. Must accept:
  #   (train_data, h, n_sims, preprocess_result)
  # Must return a tibble with columns: location, target_end_date, .sim
  custom_simulate_fn = NULL,
  
  # ---- Metadata ----
  # These fields are written to model-metadata/{team_abbr}-{model_abbr}.yml
  # by the engine. See model-metadata/README.md for field descriptions.
  # The engine calls write_metadata() automatically during run_forecast().
  metadata = list(
    team_name    = "Reger",
    model_name   = "Seasonal Naive with Box-Cox and Bootstrap",
    model_contributors = list(
      list(
        name        = "Karl Reger",
        affiliation = "Northern Arizona University",
        email       = "kcr28@nau.edu"
      )
    ),
    license             = "CC-BY-4.0",
    designated_model    = TRUE,
    data_inputs         = "ILI%",
    methods             = paste0(
      "Seasonal naive forecasts with Box-Cox variance stabilization and ",
      "bootstrapped residual prediction intervals via fable generate()."
    ),
    methods_long        = paste0(
      "Uses SNAIVE() with a yearly lag and Box-Cox transformation (lambda ",
      "estimated via guerrero method per location) to forecast weekly ",
      "weighted ILI percentages. Prediction intervals are generated from ",
      "1000 bootstrapped residual paths using fable generate(bootstrap=TRUE), ",
      "with empirical quantiles extracted at 23 levels. Forecasted values ",
      "are floored at zero to ensure non-negativity."
    ),
    ensemble_of_models      = FALSE,
    ensemble_of_hub_models  = FALSE
  )
)


###############################################################################
# EXAMPLES: How other models would look as config files
#
# Copy this file, change the fields below, save as e.g. configs/ets_log.R.
# The engine handles metadata YAML generation, output directory creation,
# and all formatting/validation automatically.
#
# ---- ETS with log transform ----
# config$model_abbr <- "ets_log"
# config$transform  <- "none"
# config$model_fn   <- function(lambda) ETS(log(observation))
# config$metadata$model_name <- "ETS with Log Transform"
# config$metadata$methods    <- "ETS state space model with log transform..."
# config$metadata$methods_long <- "Uses fable's ETS() with log transform..."
#
# ---- ARIMA (no transform) ----
# config$model_abbr <- "arima_raw"
# config$transform  <- "none"
# config$model_fn   <- function(lambda) ARIMA(observation)
# config$metadata$model_name <- "ARIMA"
# config$metadata$methods    <- "Auto ARIMA with parametric intervals..."
#
# ---- Custom non-fable model ----
# config$model_fn   <- NULL
# config$simulate_method <- "custom"
# config$custom_simulate_fn <- function(train_data, h, n_sims, preprocess_result) {
#   # Must return tibble with: location, target_end_date, .sim
# }
###############################################################################