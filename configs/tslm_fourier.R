###############################################################################
# CONFIG: configs/tslm_fourier.R
#
# Model configuration for: TSLM with Fourier Seasonal Terms + Bootstrap
# Engine:  scripts/forecast_engine.R
#
# TSLM (Time Series Linear Model) is a regression model where the response
# is a function of time-based predictors. With Fourier terms, the seasonal
# pattern is captured by K pairs of sine/cosine functions at different
# frequencies, creating a smooth approximation of the annual cycle.
#
# The key conceptual difference from other models:
#   - SNAIVE: Seasonality is a step function (copy last year's exact value)
#   - ETS:    Seasonality evolves via state updates (can adapt over time)
#   - ARIMA:  Seasonality via differencing (removes it, models remainder)
#   - TSLM:   Seasonality is a FIXED smooth curve (same shape every year)
#
# Fourier(K=5) uses 10 terms (5 sine + 5 cosine pairs) to approximate
# the annual seasonal shape. Higher K = more flexible curve but more
# parameters. For weekly data with period 52, K can range from 1 (pure
# sine wave) to 26 (maximally flexible, equivalent to seasonal dummies).
# K=5 is a good balance: captures the peaked winter shape and asymmetric
# ascent/descent without overfitting to noise.
#
# A trend() term is included to capture any slow drift in baseline ILI
# levels across the 5-year study period.
#
# WHY THIS MODEL CLASS MATTERS FOR ENSEMBLING:
# TSLM errors have a fundamentally different structure than state-space
# models because it can't adapt to recent data. If ILI is unusually high
# for a given calendar week, TSLM doesn't know -- it predicts the same
# seasonal value regardless. This "rigidity" means it makes DIFFERENT
# mistakes from adaptive models, which is exactly what drives ensemble
# error decorrelation.
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/tslm_fourier.R")
#   result <- run_forecast(config)
#
###############################################################################

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "tslm_fourier",   # 12 chars, within 15-char limit

  # ---- Paths ----
  hub_path = "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026",

  # ---- Simulation settings ----
  n_sims = 1000,

  # Fourier regression needs enough data to estimate the coefficients
  # reliably. With 10 Fourier terms + trend + intercept = 12 parameters,
  # we need well over 12 observations. One full year (52 weeks) is the
  # practical minimum for a meaningful fit.
  min_train_weeks = 52,

  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,

  # ---- Reproducibility ----
  base_seed = 12345,

  # ---- Transform ----
  # No external Box-Cox preprocessing. The log transform is embedded
  # in the model formula (like ets_log) to stabilize variance. Using
  # the +0.001 offset to handle zero-valued summer observations.
  #
  # We use log rather than Box-Cox here because TSLM with Fourier terms
  # models the seasonal pattern in the TRANSFORMED space. A log transform
  # means the Fourier curve describes multiplicative seasonality in the
  # original space (seasonality scales with level), which matches ILI
  # behavior where winter peaks have larger swings than summer troughs.
  transform = "none",

  # ---- Model specification ----
  # TSLM() with trend() and fourier(K=5):
  #   - trend() adds a linear time trend (captures slow drift)
  #   - fourier(K=5) adds 5 pairs of sin/cos terms at annual frequencies
  #     (harmonics 1-5 of the 52-week period)
  #
  # The log(observation + 0.001) transform:
  #   - Stabilizes variance (larger level = larger residuals)
  #   - Makes the Fourier seasonal curve multiplicative
  #   - The +0.001 prevents log(0) for zero-valued summer weeks
  #   - fable back-transforms automatically when generating forecasts
  #
  # lambda argument is ignored (NA) since transform = "none".
  model_fn = function(lambda) {
    TSLM(log(observation + 0.001) ~ trend() + fourier(K = 5))
  },

  # ---- Simulation method ----
  # Bootstrap residuals for prediction intervals. TSLM residuals will
  # contain all the patterns that the fixed Fourier curve can't capture
  # (year-to-year differences in peak timing/height, unusual events).
  # Resampling these residuals produces intervals that reflect this
  # irreducible uncertainty.
  simulate_method = "bootstrap",

  # ---- Custom simulation function ----
  custom_simulate_fn = NULL,

  # ---- Metadata ----
  metadata = list(
    team_name    = "Reger",
    model_name   = "Regression with Fourier Seasonality and Bootstrap",
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
      "Linear regression with Fourier seasonal terms (K=5), log ",
      "transform, and bootstrapped residual intervals via fable."
    ),
    methods_long        = paste0(
      "Uses fable's TSLM() with a linear trend and 5 Fourier harmonic ",
      "pairs (10 sine/cosine terms) to approximate the annual ILI ",
      "seasonal cycle as a smooth curve. Log(x+0.001) transform ",
      "stabilizes variance and makes seasonality multiplicative. ",
      "Prediction intervals from 1000 bootstrapped residual paths ",
      "via fable generate(bootstrap=TRUE). Empirical quantiles at ",
      "23 levels. Values floored at zero for non-negativity."
    ),
    ensemble_of_models      = FALSE,
    ensemble_of_hub_models  = FALSE
  )
)
