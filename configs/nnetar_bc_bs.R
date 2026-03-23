###############################################################################
# CONFIG: configs/nnetar_bc_bs.R
#
# Model configuration for: NNETAR + Box-Cox + Bootstrapped Residuals
# Engine:  scripts/forecast_engine.R
#
# NNETAR (Neural Network AutoRegression) is a feed-forward neural network
# with lagged values of the time series as inputs and a single hidden layer.
# fable fits 20 networks with random starting weights and averages them to
# reduce sensitivity to initialization.
#
# For seasonal data, the model is NNAR(p, P, k)[m] -- analogous to
# ARIMA(p,0,0)(P,0,0)[m] but with nonlinear activation functions instead
# of linear combinations. This means NNETAR can capture:
#   - Asymmetric peak shapes (steep ascent, gradual descent)
#   - Nonlinear interactions between recent values and seasonal lags
#   - Regime-dependent dynamics (behaves differently at low vs high ILI)
#
# These are patterns that ALL linear models (ARIMA, ETS, SNAIVE, TSLM)
# cannot capture by construction. This nonlinearity is NNETAR's unique
# contribution to the ensemble.
#
# CAVEATS:
#   - NNETAR is slower than linear models (~5-15 sec per location)
#   - It can overfit with limited training data
#   - Multi-step forecasts are generated recursively (feeding predictions
#     back as inputs), so errors compound more than with direct methods
#   - Prediction intervals come from simulation, not closed-form formulas
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/nnetar_bc_bs.R")
#   result <- run_forecast(config)
#
###############################################################################

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "nnetar_bc_bs",   # 12 chars, within 15-char limit

  # ---- Paths ----
  hub_path = "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026",

  # ---- Simulation settings ----
  # 1000 for production. NNETAR simulation is slower than linear models
  # because each path requires recursive multi-step prediction through
  # the network. For initial exploratory runs, use 500.
  n_sims = 1000,

  # NNETAR with seasonal lags needs at least 2 full seasonal cycles to
  # have enough lagged inputs for training. With period 52, that's 104.
  min_train_weeks = 104,

  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,

  # ---- Reproducibility ----
  # IMPORTANT: NNETAR uses random initialization for its 20 networks.
  # The per-date seed ensures reproducibility, but the 20-network
  # averaging already provides some variance reduction.
  base_seed = 12345,

  # ---- Transform ----
  # Box-Cox stabilizes variance before feeding into the network. This
  # helps because neural networks are sensitive to the scale of inputs.
  # The scale_inputs = TRUE default in NNETAR() further normalizes the
  # lagged inputs to [0, 1] range.
  transform = "box_cox_guerrero",

  # ---- Model specification ----
  # NNETAR() with no explicit p/P/k lets fable auto-select:
  #   p = optimal number of non-seasonal lags (like AR order)
  #   P = number of seasonal lags (typically 1)
  #   k = hidden layer size (default: (p + P + 1) / 2, rounded)
  #
  # The Box-Cox transform is applied to the response. fable handles
  # the back-transformation of generated paths automatically.
  #
  # n_networks is left at the default of 20 for stability.
  model_fn = function(lambda) {
    NNETAR(box_cox(observation, lambda))
  },

  # ---- Simulation method ----
  # NNETAR always uses simulation for prediction intervals (there's no
  # closed-form distribution). Setting "bootstrap" tells generate() to
  # resample from fitted residuals when generating future paths. This
  # is the standard approach for NNETAR uncertainty quantification.
  simulate_method = "bootstrap",

  # ---- Custom simulation function ----
  custom_simulate_fn = NULL,

  # ---- Metadata ----
  metadata = list(
    team_name    = "Reger",
    model_name   = "Neural Network Autoregression with Box-Cox and Bootstrap",
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
      "NNETAR neural network autoregression with Box-Cox stabilization ",
      "and bootstrapped residual simulation via fable generate()."
    ),
    methods_long        = paste0(
      "Uses fable's NNETAR() feed-forward neural network with ",
      "auto-selected lagged inputs (seasonal and non-seasonal) and ",
      "Box-Cox transformation (lambda via guerrero per location). ",
      "Twenty networks with random starting weights are averaged. ",
      "Prediction intervals from 1000 bootstrapped residual simulation ",
      "paths using fable generate(bootstrap=TRUE). Empirical quantiles ",
      "at 23 levels. Values floored at zero for non-negativity."
    ),
    ensemble_of_models      = FALSE,
    ensemble_of_hub_models  = FALSE
  )
)
