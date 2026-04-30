###############################################################################
# CONFIG: configs/theta_bc_bs.R
#
# Model configuration for: Theta Method + Box-Cox + Bootstrapped Residuals
# Engine:  scripts/forecast_engine.R
#
# The Theta method won the M3 forecasting competition (Assimakopoulos &
# Nikolopoulos, 2000). It works by decomposing the series into two
# "theta lines" -- modified versions of the original series with different
# curvatures -- and combining their forecasts. Hyndman & Billah (2003)
# showed that the standard Theta method is equivalent to simple exponential
# smoothing (SES) with drift applied to the seasonally adjusted series.
#
# For our purposes, Theta is interesting as an ensemble candidate because:
#   - It's a hybrid that blends local level (SES) with global trend (drift)
#   - It tends to produce conservative, well-calibrated forecasts
#   - Its error structure may differ from both pure ARIMA and pure ETS
#
# EXPECTED CORRELATION: Theta is theoretically related to ETS(A,Ad,N) --
# additive damped trend with no seasonality. If our ETS model selects a
# similar specification, the two will be highly correlated and one should
# be dropped from the ensemble. The empirical correlation analysis in
# Phase 2 will determine this.
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/theta_bc_bs.R")
#   result <- run_forecast(config)
#
###############################################################################

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "theta_bc_bs",    # 11 chars, within 15-char limit

  # ---- Paths ----
  hub_path = normalizePath(".", mustWork = FALSE),

  # ---- Simulation settings ----
  n_sims = 1000,

  # Theta needs enough data for seasonal decomposition + SES fitting.
  # One full seasonal cycle is the practical minimum.
  min_train_weeks = 52,

  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,

  # ---- Reproducibility ----
  base_seed = 12345,

  # ---- Transform ----
  # Box-Cox for variance stabilization, same as other models. Theta
  # assumes additive structure internally, so stabilizing variance
  # before applying it improves the additivity assumption.
  transform = "box_cox_guerrero",

  # ---- Model specification ----
  # THETA() with no arguments uses the standard theta method.
  # fable's implementation handles seasonal adjustment internally
  # via STL decomposition before applying the theta lines.
  #
  # Box-Cox is applied to the response to stabilize variance.
  model_fn = function(lambda) {
    THETA(box_cox(observation, lambda))
  },

  # ---- Simulation method ----
  # Bootstrap residuals for distribution-free prediction intervals.
  simulate_method = "bootstrap",

  # ---- Custom simulation function ----
  custom_simulate_fn = NULL,

  # ---- Metadata ----
  metadata = list(
    team_name    = "Reger",
    model_name   = "Theta Method with Box-Cox and Bootstrap",
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
      "Theta method (M3 competition winner) with Box-Cox stabilization ",
      "and bootstrapped prediction intervals via fable generate()."
    ),
    methods_long        = paste0(
      "Uses fable's THETA() method, which decomposes the series into ",
      "theta lines and combines SES with drift on the seasonally ",
      "adjusted series. Box-Cox transformation (lambda via guerrero ",
      "per location) stabilizes variance. Prediction intervals from ",
      "1000 bootstrapped residual paths via fable ",
      "generate(bootstrap=TRUE). Empirical quantiles at 23 levels. ",
      "Values floored at zero for non-negativity."
    ),
    ensemble_of_models      = FALSE,
    ensemble_of_hub_models  = FALSE
  )
)
