###############################################################################
# CONFIG: configs/ets_log.R
#
# Model configuration for: ETS + Log Transform + Bootstrapped Residuals
# Engine:  scripts/forecast_engine.R
#
# ETS (Error-Trend-Seasonal) is a state space model. fable's ETS()
# automatically selects the best error, trend, and seasonal components
# via AICc. The log transform stabilizes variance (equivalent to
# Box-Cox with lambda = 0) and is applied inside the model formula.
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/ets_log.R")
#   result <- run_forecast(config)
#
###############################################################################

config <- list(
  
  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "ets_log",
  
  # ---- Paths ----
  hub_path = "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026",
  
  # ---- Simulation settings ----
  n_sims = 1000,
  
  # ETS can fit with less data than SNAIVE, but benefits from at least one
  # full seasonal cycle for the automatic seasonal component selection.
  min_train_weeks = 52,
  
  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,
  
  # ---- Reproducibility ----
  # Same base_seed as snaive_bc_bs so that for any given origin date, the
  # RNG state is determined solely by the date. Different models produce
  # different simulations anyway because the fitted model differs.
  base_seed = 12345,
  
  # ---- Transform ----
  # "none" because the log transform is embedded in the model formula itself
  # (fable handles the back-transformation automatically). No external lambda
  # estimation is needed.
  transform = "none",
  
  # ---- Model specification ----
  # ETS() with log(observation + offset) lets fable auto-select
  # error/trend/seasonal components. The small offset (0.001) prevents
  # log(0) = -Inf for the ~10 zero-valued summer observations in the ILI
  # data. Fable parses the expression and inverts it (exp(y) - 0.001)
  # automatically during back-transformation.
  # The lambda argument is ignored (NA) since transform = "none".
  model_fn = function(lambda) {
    ETS(log(observation + 0.001))
  },
  
  # ---- Simulation method ----
  simulate_method = "bootstrap",
  
  # ---- Custom simulation function ----
  custom_simulate_fn = NULL,
  
  # ---- Metadata ----
  metadata = list(
    team_name    = "Reger",
    model_name   = "ETS with Log Transform and Bootstrap",
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
      "ETS state space model with log variance stabilization and ",
      "bootstrapped residual prediction intervals via fable generate()."
    ),
    methods_long        = paste0(
      "Uses fable's ETS() with automatic error/trend/seasonal component ",
      "selection (minimizing AICc) to forecast weekly weighted ILI ",
      "percentages. A log(x + 0.001) transformation stabilizes variance ",
      "across flu seasons; the small offset prevents log(0) for rare ",
      "zero-valued summer observations. The model is fit independently ",
      "per location for each origin date. Prediction intervals are ",
      "generated via bootstrap resampling of fitted residuals (1000 ",
      "replicates) using fable::generate(bootstrap = TRUE). Empirical ",
      "quantiles are extracted at 23 levels (0.01 to 0.99). All forecast ",
      "values are floored at zero to satisfy the hub's non-negativity ",
      "constraint."
    ),
    ensemble_of_models      = FALSE,
    ensemble_of_hub_models  = FALSE
  )
)