###############################################################################
# FILE: configs/nnetar_log.R
#
# PURPOSE:
#   Component model config: NNETAR (Neural Network AutoRegression) with
#   log(observation + 0.001) transform and bootstrapped residual prediction
#   intervals.
#
#   This is one of the Phase 1 transform variants
#   (FLUCAST_IMPLEMENTATION_PLAN.md §4.2 Tier 2 / Checkpoint 1C). It's a
#   single-line variant of the existing nnetar_bc_bs.R config — the model
#   family and engine path are identical; only the variance-stabilizing
#   transform changes from per-location Box-Cox to uniform log+offset.
#
# ROLE IN CANDIDATE POOL:
#   The candidate pool already contains nnetar_bc_bs (NNETAR + Box-Cox),
#   which was the best individual performer in the proposal-era analysis
#   (MAE 0.360, WIS 0.241). This variant tests whether a simpler, uniform
#   log transform yields meaningfully different forecasts when fed into
#   the same nonlinear model architecture.
#
#   WHY MIGHT THIS DIFFER FROM nnetar_bc_bs?
#   Neural networks are sensitive to the scale of their inputs. NNETAR
#   internally normalizes lagged inputs to [0, 1] (its scale_inputs default),
#   which partly mitigates the choice of pre-transform. But the
#   normalization is computed from the training window, so the choice of
#   pre-transform still affects which extreme values are clipped and how
#   the network's nonlinearity engages with the data. Box-Cox produces
#   location-specific dynamic ranges; log produces a uniform multiplicative
#   structure across locations. The two could behave noticeably differently
#   on regions with very different baseline wILI levels (e.g. Region 6).
#
#   Twin-pair caveat: see arima_log.R header. If correlation > 0.93 with
#   nnetar_bc_bs, the twin rule keeps one.
#
# WHY THE +0.001 OFFSET?
#   See configs/arima_log.R or configs/ets_log.R for full rationale.
#   ~10 zero-valued summer observations require a small positive offset
#   to avoid log(0) = -Inf.
#
# INPUTS / OUTPUTS / RUNTIME:
#   See configs/nnetar_bc_bs.R header.
#   Estimated runtime: ~5-10 min per origin date for all 11 locations
#   (NNETAR's recursive multi-step simulation dominates; transform choice
#   is irrelevant to runtime).
#
# HOW TO RUN (validation only, Phase 1):
#   source("scripts/forecast_engine.R")
#   source("configs/nnetar_log.R")
#   source("scripts/run_nnetar_log_validation.R")
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md §4.2 Tier 2, Checkpoint 1C
#   - configs/nnetar_bc_bs.R (companion config — Box-Cox variant)
#   - configs/ets_log.R (precedent for log+offset pattern)
###############################################################################

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "nnetar_log",      # 10 chars, within 15-char limit

  # ---- Paths ----
  hub_path = normalizePath(".", mustWork = FALSE),

  # ---- Simulation settings ----
  n_sims = 1000L,

  # NNETAR with seasonal lags needs at least 2 full cycles to populate
  # the lagged-input vector — same constraint as nnetar_bc_bs.
  min_train_weeks = 104L,

  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,

  # ---- Reproducibility ----
  # NNETAR has its own RNG (initializes 20 networks with random weights
  # which are then averaged). The engine's per-date set.seed() drives all
  # downstream randomness deterministically.
  base_seed = 12345L,

  # ---- Transform ----
  # "none" — fable parses the log+offset symbolically inside model_fn.
  # The engine's preprocess step does nothing.
  transform = "none",

  # ---- Model specification ----
  # NNETAR(log(observation + 0.001)):
  #   - log(...)        variance stabilization
  #   - + 0.001         offset for zero-valued observations
  #   - NNETAR(...)     auto-selects p (non-seasonal lags), P (seasonal lags),
  #                     and k (hidden layer size); fits 20 networks with
  #                     random initialization and averages
  model_fn = function(lambda) {
    NNETAR(log(observation + 0.001))
  },

  # ---- Simulation method ----
  # NNETAR has no closed-form prediction interval; uncertainty must come
  # from simulation regardless. "bootstrap" tells fable's generate() to
  # resample fitted residuals when producing simulated future paths,
  # rather than assuming a parametric residual distribution.
  simulate_method = "bootstrap",

  # ---- Custom simulation function ----
  custom_simulate_fn = NULL,

  # ---- Metadata ----
  metadata = list(
    team_name  = "Reger",
    model_name = "Neural Network Autoregression with Log Transform and Bootstrap",

    model_contributors = list(
      list(
        name        = "Karl Reger",
        affiliation = "Northern Arizona University",
        email       = "kcr28@nau.edu"
      )
    ),

    license          = "CC-BY-4.0",
    designated_model = FALSE,
    data_inputs      = "ILI%",

    methods = paste0(
      "NNETAR neural-network autoregression with log(x + 0.001) variance ",
      "stabilization and bootstrapped residual simulation via fable ",
      "generate()."
    ),

    methods_long = paste0(
      "Uses fable's NNETAR() feed-forward neural network with auto-selected ",
      "lagged inputs (seasonal and non-seasonal) and a log(x + 0.001) ",
      "transformation. The 0.001 offset prevents log(0) for the rare ",
      "zero-valued summer observations. Twenty networks with random ",
      "starting weights are fit and averaged to reduce initialization ",
      "sensitivity. Prediction intervals are generated from 1000 ",
      "bootstrapped residual paths via fable generate(bootstrap=TRUE) — ",
      "NNETAR has no closed-form predictive distribution, so simulation ",
      "is the only path to quantiles. Empirical quantiles at 23 levels, ",
      "floored at zero. Companion candidate to nnetar_bc_bs: same ",
      "architecture, simpler transform. Phase 2 correlation analysis ",
      "determines whether both contribute or are twins."
    ),

    ensemble_of_models     = FALSE,
    ensemble_of_hub_models = FALSE
  )
)
