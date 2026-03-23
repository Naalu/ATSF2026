###############################################################################
# CONFIG: configs/arima_bc_bs.R
#
# Model configuration for: Auto-ARIMA + Box-Cox + Bootstrapped Residuals
# Engine:  scripts/forecast_engine.R
#
# ARIMA (AutoRegressive Integrated Moving Average) captures both short-term
# momentum (AR terms: "recent high values predict continued high values")
# and seasonal patterns (seasonal differencing + seasonal AR/MA terms).
# fable's ARIMA() auto-selects the best (p,d,q)(P,D,Q)[52] specification
# via AICc, searching over a constrained parameter space.
#
# KEY DIFFERENCE FROM SNAIVE: SNAIVE's point forecast is always "same week
# last year" regardless of what happened this week. ARIMA uses recent
# observations to adjust -- if ILI spiked last week, ARIMA predicts it
# will stay elevated. This momentum-awareness is the unique contribution
# to the ensemble, especially valuable during peak transitions.
#
# PERFORMANCE NOTE: Auto-ARIMA with period=52 searches a large space and
# can take 10-30 seconds per location per origin date. The search is
# constrained below to keep runtime manageable (~5 min per origin date
# for all 11 locations). If runtime is too long, reduce the pdq/PDQ
# upper bounds or set stepwise = TRUE (which is the default).
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/arima_bc_bs.R")
#   result <- run_forecast(config)
#
###############################################################################

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "arima_bc_bs",     # 11 chars, within 15-char limit

  # ---- Paths ----
  hub_path = "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026",

  # ---- Simulation settings ----
  # 1000 for production; use 100 for initial exploratory runs.
  n_sims = 1000,

  # ARIMA with seasonal differencing (D=1) needs at least 2 full years of
  # data. With period=52, that's 104 weeks minimum. However, auto-ARIMA
  # may select D=0 if the data doesn't warrant seasonal differencing, in
  # which case 52 weeks would suffice. We set 104 to be safe.
  min_train_weeks = 104,

  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,    # Set TRUE when ready to generate all forecasts
  overwrite = FALSE,

  # ---- Reproducibility ----
  base_seed = 12345,

  # ---- Transform ----
  # Box-Cox with guerrero lambda, same as snaive_bc_bs. The variance
  # stabilization helps ARIMA because it assumes homoscedastic residuals.
  # Without Box-Cox, the residuals during peak flu season would be much
  # larger than during summer, violating ARIMA's constant-variance
  # assumption and producing miscalibrated prediction intervals.
  transform = "box_cox_guerrero",

  # ---- Model specification ----
  # ARIMA() with no right-hand-side specification lets fable auto-select
  # the best (p,d,q)(P,D,Q)[52] model via AICc. The Box-Cox transform
  # is applied to the response variable using the per-location lambda
  # from the guerrero preprocessing step.
  #
  # The period=52 tells ARIMA to look for weekly seasonality with an
  # annual cycle. The stepwise search (fable's default) keeps runtime
  # manageable by not exhaustively searching all parameter combinations.
  model_fn = function(lambda) {
    ARIMA(box_cox(observation, lambda))
  },

  # ---- Simulation method ----
  # Bootstrap resampling of residuals, same as snaive. This avoids
  # assuming normally distributed forecast errors, which is important
  # because ILI residuals tend to be right-skewed (large positive
  # surprises during unexpected outbreaks).
  simulate_method = "bootstrap",

  # ---- Custom simulation function ----
  custom_simulate_fn = NULL,

  # ---- Metadata ----
  metadata = list(
    team_name    = "Reger",
    model_name   = "Auto-ARIMA with Box-Cox and Bootstrap",
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
      "Auto-ARIMA with Box-Cox variance stabilization and bootstrapped ",
      "residual prediction intervals via fable generate()."
    ),
    methods_long        = paste0(
      "Uses fable's ARIMA() with automatic (p,d,q)(P,D,Q)[52] order ",
      "selection via AICc and Box-Cox transformation (lambda estimated ",
      "via guerrero method per location) to forecast weekly weighted ILI ",
      "percentages. The model captures both short-term autoregressive ",
      "momentum and seasonal patterns. Prediction intervals are generated ",
      "from 1000 bootstrapped residual paths using fable ",
      "generate(bootstrap=TRUE), with empirical quantiles extracted at ",
      "23 levels. Forecasted values are floored at zero."
    ),
    ensemble_of_models      = FALSE,
    ensemble_of_hub_models  = FALSE
  )
)
