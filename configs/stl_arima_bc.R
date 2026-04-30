###############################################################################
# CONFIG: configs/stl_arima_bc.R
#
# Model configuration for: STL Decomposition + ARIMA on Remainder + Box-Cox
# Engine:  scripts/forecast_engine.R
#
# This model uses a two-stage approach:
#   Stage 1: STL (Seasonal-Trend decomposition using Loess) separates the
#            series into seasonal + trend + remainder components. STL uses
#            local regression (loess), making it robust to outliers and
#            able to capture a seasonal shape that changes gradually.
#   Stage 2: ARIMA is fit to the seasonally adjusted remainder, capturing
#            short-term dynamics that STL's smooth decomposition misses.
#
# WHY THIS DIFFERS FROM PLAIN ARIMA:
# Plain ARIMA handles seasonality via seasonal differencing (removing it
# algebraically) or seasonal AR/MA terms (modeling it parametrically).
# STL+ARIMA handles seasonality NON-PARAMETRICALLY via loess smoothing,
# then ARIMA only models the leftover short-term patterns. The seasonal
# component is re-added when forecasting.
#
# The practical difference: STL's seasonal component is estimated from
# the entire training history using robust local regression, while
# ARIMA's seasonal terms are parametric and estimated jointly with the
# non-seasonal terms. STL is more flexible for irregular seasonal
# shapes, while ARIMA's seasonal terms are more constrained.
#
# This model is implemented using fabletools::decomposition_model(),
# which handles the decompose -> model remainder -> recombine pipeline
# automatically, including proper forecast generation.
#
# POTENTIAL ISSUE: decomposition_model() with generate(bootstrap=TRUE)
# resamples residuals from the ARIMA fit on the remainder. This is
# conceptually clean but may have edge cases. If this config fails at
# runtime, switch to the custom_simulate_fn path (see hist_week.R for
# the pattern) and implement the STL+ARIMA pipeline manually.
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/stl_arima_bc.R")
#   result <- run_forecast(config)
#
###############################################################################

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "stl_arima_bc",   # 12 chars, within 15-char limit

  # ---- Paths ----
  hub_path = normalizePath(".", mustWork = FALSE),

  # ---- Simulation settings ----
  n_sims = 1000,

  # STL needs at least 2 full cycles for a stable seasonal decomposition,
  # plus ARIMA needs enough remainder observations. 104 weeks is safe.
  min_train_weeks = 104,

  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,

  # ---- Reproducibility ----
  base_seed = 12345,

  # ---- Transform ----
  # Box-Cox applied inside the STL formula. The decomposition operates
  # in the transformed space, which makes the additive decomposition
  # assumption (observation = seasonal + trend + remainder) more
  # appropriate when variance scales with level.
  transform = "box_cox_guerrero",

  # ---- Model specification ----
  # decomposition_model() from fabletools:
  #   1st arg: STL decomposition formula with season(period = 52)
  #            The period must be explicit because our tsibble index
  #            is Date (not yearweek), so auto-detection may not work.
  #   2nd arg: ARIMA model applied to the season_adjust component
  #            (the trend + remainder after removing seasonality).
  #            Auto-ARIMA selects the best order for this smoother series.
  #
  # When forecasting, fabletools automatically re-adds the seasonal
  # component (using the last seasonal cycle from STL, like SNAIVE)
  # to the ARIMA forecast of the seasonally adjusted series.
  model_fn = function(lambda) {
    decomposition_model(
      STL(box_cox(observation, lambda) ~ season(period = 52)),
      ARIMA(season_adjust)
    )
  },

  # ---- Simulation method ----
  # Bootstrap residuals from the ARIMA fit on the remainder component.
  # The seasonal component is treated as deterministic (same pattern
  # as the most recent cycle), so all forecast uncertainty comes from
  # the ARIMA model's residuals.
  simulate_method = "bootstrap",

  # ---- Custom simulation function ----
  custom_simulate_fn = NULL,

  # ---- Metadata ----
  metadata = list(
    team_name    = "Reger",
    model_name   = "STL Decomposition with ARIMA and Box-Cox",
    model_contributors = list(
      list(
        name        = "Karl Reger",
        affiliation = "Northern Arizona University",
        email       = "kcr28@nau.edu"
      )
    ),
    license             = "CC-BY-4.0",
    designated_model    = FALSE,
    data_inputs         = "ILI%",
    methods             = paste0(
      "STL seasonal decomposition with ARIMA on remainder, Box-Cox ",
      "stabilization, bootstrapped intervals via fable generate()."
    ),
    methods_long        = paste0(
      "Two-stage model using fabletools decomposition_model(). ",
      "Stage 1: STL decomposes Box-Cox transformed ILI% into seasonal ",
      "(period 52, loess-based) and seasonally adjusted components. ",
      "Stage 2: Auto-ARIMA models the seasonally adjusted remainder. ",
      "Forecasts recombine the ARIMA predictions with the last ",
      "seasonal cycle. Lambda estimated via guerrero per location. ",
      "Prediction intervals from 1000 bootstrapped residual paths. ",
      "Empirical quantiles at 23 levels, floored at zero."
    ),
    ensemble_of_models      = FALSE,
    ensemble_of_hub_models  = FALSE
  )
)
