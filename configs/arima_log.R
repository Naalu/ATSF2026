###############################################################################
# FILE: configs/arima_log.R
#
# PURPOSE:
#   Component model config: Auto-ARIMA with log(observation + 0.001) transform
#   and bootstrapped residual prediction intervals.
#
#   This is one of the Phase 1 transform variants
#   (FLUCAST_IMPLEMENTATION_PLAN.md §4.2 Tier 2 / Checkpoint 1C). It's a
#   single-line variant of the existing arima_bc_bs.R config — the model
#   family and engine path are identical; only the variance-stabilizing
#   transform changes from Box-Cox (per-location lambda via guerrero) to
#   log(x + 0.001) (a single hand-picked transform applied uniformly).
#
# ROLE IN CANDIDATE POOL:
#   The candidate pool already contains arima_bc_bs (ARIMA + Box-Cox). This
#   variant tests whether the simpler log transform — a special case of
#   Box-Cox with lambda=0 — produces materially different forecasts. If the
#   correlation analysis in Phase 2 shows arima_log and arima_bc_bs are
#   "twins" (correlation > 0.93), the twin-pair rule will keep only one.
#   Including both candidates in Phase 1 lets the data decide.
#
#   WHY LOG INSTEAD OF BOX-COX?
#   - Computationally cheaper (no per-location guerrero estimation)
#   - Forces the same transformation across all locations, which can
#     improve interpretability when comparing locations
#   - When the underlying process is multiplicative (which ILI roughly is),
#     log is the theoretically-motivated choice; Box-Cox just empirically
#     finds something close to log when lambda ≈ 0
#
# WHY THE +0.001 OFFSET?
#   Approximately 10 wILI observations in the training data are exactly 0
#   (summer weeks in low-incidence locations). log(0) = -Inf, which crashes
#   the fable model fit. Adding 0.001 (one part in a thousand of a typical
#   1-3% wILI value) is well below measurement noise and preserves
#   monotonicity. fable's formula parser handles log(x + c) symbolically
#   and inverts it correctly during back-transformation (exp(y) - c).
#
# INPUTS / OUTPUTS / RUNTIME:
#   See configs/arima_bc_bs.R header — identical to that config.
#   Estimated runtime: ~3-5 min per origin date for all 11 locations
#   (auto-ARIMA stepwise search dominates; transform choice is irrelevant
#   to runtime).
#
# HOW TO RUN (validation only, Phase 1):
#   source("scripts/forecast_engine.R")
#   source("configs/arima_log.R")
#   source("scripts/run_arima_log_validation.R")
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md §4.2 Tier 2, Checkpoint 1C
#   - configs/arima_bc_bs.R (base config this is patterned on)
#   - configs/ets_log.R (precedent for the log+offset pattern)
###############################################################################

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "arima_log",       # 9 chars, well within hub's 15-char limit

  # ---- Paths ----
  hub_path = normalizePath(".", mustWork = FALSE),

  # ---- Simulation settings ----
  # 1000 sims matches arima_bc_bs and the rest of the bootstrap-based
  # fable models. Keeping all candidates at the same n_sims avoids any
  # apples-to-oranges comparison in the Phase 2 correlation analysis.
  n_sims = 1000L,

  # ---- Training window requirement ----
  # Same as arima_bc_bs: ARIMA with seasonal differencing (D=1) needs at
  # least 2 full annual cycles. With weekly data and period=52, that's
  # 104 weeks. Auto-ARIMA may select D=0 in which case 52 would suffice;
  # 104 is the defensive floor.
  min_train_weeks = 104L,

  # ---- Run control ----
  # All four Phase 1 candidates use the same control settings — driven by
  # per-model run_*_validation.R drivers, not by the engine's auto-loop.
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,

  # ---- Reproducibility ----
  base_seed = 12345L,

  # ---- Transform ----
  # "none" — the engine's preprocess step does nothing. The transform is
  # baked into model_fn via fable's formula syntax, which means fable
  # tracks it symbolically and inverts it during forecast back-transformation.
  transform = "none",

  # ---- Model specification ----
  # ARIMA(log(observation + 0.001)):
  #   - log(...)        the variance-stabilizing transform
  #   - + 0.001         offset to handle the 10 zero-valued obs in training data
  #   - ARIMA(...)      auto-selects (p,d,q)(P,D,Q)[52] via AICc stepwise search
  # The lambda argument is unused (transform = "none" passes NA), but the
  # engine's fable path always passes lambda to model_fn so the signature
  # must accept it.
  model_fn = function(lambda) {
    ARIMA(log(observation + 0.001))
  },

  # ---- Simulation method ----
  # Bootstrap residuals — same as arima_bc_bs. Does not assume Gaussian
  # forecast errors, which matters because log-transformed ILI residuals
  # are still mildly right-skewed (large positive surprises during peaks).
  simulate_method = "bootstrap",

  # ---- Custom simulation function ----
  custom_simulate_fn = NULL,

  # ---- Metadata ----
  metadata = list(
    team_name  = "Reger",
    model_name = "Auto-ARIMA with Log Transform and Bootstrap",

    model_contributors = list(
      list(
        name        = "Karl Reger",
        affiliation = "Northern Arizona University",
        email       = "kcr28@nau.edu"
      )
    ),

    license          = "CC-BY-4.0",
    designated_model = FALSE,        # Candidate, not the team submission
    data_inputs      = "ILI%",

    # <200 chars per hub schema.
    methods = paste0(
      "Auto-ARIMA with log(x + 0.001) variance stabilization and ",
      "bootstrapped residual prediction intervals via fable generate()."
    ),

    methods_long = paste0(
      "Uses fable's ARIMA() with automatic (p,d,q)(P,D,Q)[52] order ",
      "selection via AICc and a log(x + 0.001) transformation to forecast ",
      "weekly weighted ILI percentages. The 0.001 offset prevents log(0) ",
      "for the rare zero-valued summer observations in the training data. ",
      "Fable parses the formula symbolically and back-transforms forecast ",
      "paths automatically. Prediction intervals are generated from 1000 ",
      "bootstrapped residual paths via fable generate(bootstrap=TRUE), ",
      "with empirical quantiles extracted at 23 levels and floored at zero. ",
      "Companion candidate to arima_bc_bs: same model family, simpler ",
      "transform (log = Box-Cox with lambda=0). Phase 2 correlation analysis ",
      "will determine whether both add value to the ensemble or are twins."
    ),

    ensemble_of_models     = FALSE,
    ensemble_of_hub_models = FALSE
  )
)
