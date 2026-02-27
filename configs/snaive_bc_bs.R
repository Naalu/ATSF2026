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
#   - model_abbr (and team_abbr if different team)
#   - model_fn   (the fable formula)
#   - transform  (if not using Box-Cox)
#   See the "Other fable models" examples at the bottom of this file.
#
# To create a NON-fable model (e.g., GAM, external API), set:
#   - model_fn = NULL
#   - custom_simulate_fn = your function
#   See the "Custom model" example at the bottom.
#
###############################################################################

config <- list(
  
  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "snaive_bc_bs",
  
  # ---- Paths ----
  # Root of the local ATSF2026 repo clone. All other paths derived from this.
  hub_path = "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026",
  
  # ---- Simulation settings ----
  # 1000 gives ≥10 obs in each tail for the 0.01/0.99 quantiles.
  # Use 100 for fast dev/testing.
  n_sims = 1000,
  
  # Minimum weeks of training data required. SNAIVE with lag("year") needs
  # at least 52 weeks (one full seasonal cycle).
  min_train_weeks = 52,
  
  # ---- Run control ----
  dev_mode = FALSE,   # TRUE = process first 5 dates only
  run_loop = FALSE,   # TRUE = run full forecast loop on source()
  
  # ---- Transform ----
  # How the engine preprocesses data before model fitting.
  #   "box_cox_guerrero" — estimate per-location lambda via guerrero(), pass
  #                        to model_fn as the `lambda` argument
  #   "none"             — no preprocessing; model_fn receives NA for lambda
  transform = "box_cox_guerrero",
  
  # ---- Model specification (fable path) ----
  # A function that receives a lambda value and returns a fable model
  # definition (the thing that goes inside model()). The engine calls this
  # once per location, passing the location-specific lambda.
  #
  # For models without Box-Cox, the lambda argument is NA — just ignore it.
  model_fn = function(lambda) {
    SNAIVE(box_cox(observation, lambda) ~ lag("year"))
  },
  
  # ---- Simulation method ----
  # "bootstrap"  — generate(bootstrap = TRUE): resample fitted residuals
  # "parametric" — generate(bootstrap = FALSE): assume normal residuals
  simulate_method = "bootstrap",
  
  # ---- Custom simulation function (non-fable escape hatch) ----
  # Set model_fn = NULL and provide this instead. Must accept:
  #   (train_data, h, n_sims, preprocess_result)
  # Must return a tibble with columns: location, target_end_date, .sim
  custom_simulate_fn = NULL
)


###############################################################################
# EXAMPLES: How other models would look as config files
#
# These are NOT executed — they're here as reference for when you create
# new model configs. Copy this file, change the fields below, and save as
# e.g. configs/ets_log.R
#
# ---- ETS with log transform (no Box-Cox lambda needed) ----
# config <- list(
#   team_abbr  = "KReger",
#   model_abbr = "ets_log",
#   hub_path   = "/Users/chrisreger/.../ATSF2026",
#   n_sims     = 1000,
#   min_train_weeks = 52,
#   dev_mode   = FALSE,
#   run_loop   = FALSE,
#   transform  = "none",
#   model_fn   = function(lambda) ETS(log(observation)),
#   simulate_method = "bootstrap",
#   custom_simulate_fn = NULL
# )
#
# ---- ARIMA (no transform) ----
# config <- list(
#   ...same boilerplate...
#   model_abbr = "arima_raw",
#   transform  = "none",
#   model_fn   = function(lambda) ARIMA(observation),
#   simulate_method = "parametric",
# )
#
# ---- SNAIVE + Box-Cox but parametric residuals ----
# config <- list(
#   ...same boilerplate...
#   model_abbr = "snaive_bc_pr",
#   transform  = "box_cox_guerrero",
#   model_fn   = function(lambda) SNAIVE(box_cox(observation, lambda) ~ lag("year")),
#   simulate_method = "parametric",
# )
#
# ---- Custom non-fable model (e.g. a GAM) ----
# config <- list(
#   ...same boilerplate...
#   model_abbr = "custom_gam",
#   transform  = "none",
#   model_fn   = NULL,  # not using fable
#   simulate_method = "custom",
#   custom_simulate_fn = function(train_data, h, n_sims, preprocess_result) {
#     # Your code here. Must return tibble with:
#     #   location, target_end_date, .sim
#     # The engine handles quantile computation, formatting, validation.
#   }
# )
###############################################################################