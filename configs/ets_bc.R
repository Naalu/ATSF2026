###############################################################################
# FILE: configs/ets_bc.R
#
# PURPOSE:
#   Component model config: ETS state-space model with Box-Cox variance
#   stabilization (per-location lambda via guerrero) and bootstrapped
#   residual prediction intervals.
#
#   This is one of the Phase 1 transform variants
#   (FLUCAST_IMPLEMENTATION_PLAN.md §4.2 Tier 2 / Checkpoint 1C). It's a
#   single-line variant of the existing ets_log.R config — the model
#   family and engine path are identical; only the variance-stabilizing
#   transform changes from log(x + 0.001) (uniform across locations) to
#   Box-Cox with per-location guerrero lambda.
#
# ROLE IN CANDIDATE POOL:
#   The candidate pool already contains ets_log (ETS + log+offset). This
#   variant tests whether per-location Box-Cox lambdas — which adapt the
#   transform to each location's variance-mean relationship — produce
#   forecasts that are sufficiently different to add ensemble value.
#
#   WHY MIGHT BOX-COX HELP?
#   - Some locations have a near-multiplicative structure (lambda ≈ 0,
#     nearly identical to log)
#   - Other locations may have intermediate variance structure (lambda
#     between 0 and 1) where Box-Cox finds a better stabilization than log
#   - For Region 6 (typically high baseline wILI) Box-Cox with lambda > 0
#     might preserve more of the absolute-scale dynamics
#
#   If correlation analysis flags ets_bc and ets_log as twins (>0.93),
#   the twin-pair rule will keep one. Including both in Phase 1 lets data
#   decide rather than us guessing.
#
# INPUTS / OUTPUTS / RUNTIME:
#   See configs/ets_log.R header — identical to that config.
#   Estimated runtime: ~1-2 min per origin date for all 11 locations
#   (ETS fits much faster than ARIMA's stepwise search).
#
# HOW TO RUN (validation only, Phase 1):
#   source("scripts/forecast_engine.R")
#   source("configs/ets_bc.R")
#   source("scripts/run_ets_bc_validation.R")
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md §4.2 Tier 2, Checkpoint 1C
#   - configs/ets_log.R (companion config — log+offset variant)
#   - configs/snaive_bc_bs.R (precedent for the box_cox_guerrero pattern)
###############################################################################

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "ets_bc",          # 6 chars, well within 15-char limit

  # ---- Paths ----
  hub_path = "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026",

  # ---- Simulation settings ----
  n_sims = 1000L,

  # ETS doesn't have ARIMA's seasonal-differencing requirement, but for
  # the seasonal components (additive or multiplicative) to be estimated
  # reliably we want at least 2 full annual cycles. 104 matches the rest
  # of the candidate pool — easier to compare apples-to-apples in Phase 2.
  min_train_weeks = 104L,

  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,

  # ---- Reproducibility ----
  base_seed = 12345L,

  # ---- Transform ----
  # box_cox_guerrero: the engine's preprocess_fn estimates a separate lambda
  # per location using the guerrero method, then passes it to model_fn.
  # This is the same pipeline used by snaive_bc_bs, arima_bc_bs, nnetar_bc_bs.
  transform = "box_cox_guerrero",

  # ---- Model specification ----
  # ETS(box_cox(observation, lambda)):
  #   - box_cox(...)    fable's symbolic Box-Cox transform; fable tracks
  #                     the lambda and inverts on back-transformation
  #   - ETS(...)        auto-selects error/trend/seasonal components
  #                     (e.g., ETS(M,N,A) or ETS(A,A,N)) via AICc
  # The lambda argument here IS used — it comes from the engine's
  # per-location guerrero estimation step.
  model_fn = function(lambda) {
    ETS(box_cox(observation, lambda))
  },

  # ---- Simulation method ----
  # Bootstrap residuals. ETS has a closed-form parametric prediction
  # interval too, but to keep all candidate models on the same uncertainty
  # mechanism we use bootstrap consistently. This makes the diversity
  # claims in Phase 2 cleaner — any difference between models is
  # attributable to model structure, not to interval methodology.
  simulate_method = "bootstrap",

  # ---- Custom simulation function ----
  custom_simulate_fn = NULL,

  # ---- Metadata ----
  metadata = list(
    team_name  = "Reger",
    model_name = "ETS with Box-Cox and Bootstrap",

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
      "ETS state-space model with Box-Cox variance stabilization (per-",
      "location guerrero lambda) and bootstrapped residual prediction ",
      "intervals via fable generate()."
    ),

    methods_long = paste0(
      "Uses fable's ETS() with automatic error/trend/seasonal component ",
      "selection (minimizing AICc) and Box-Cox transformation (lambda ",
      "estimated via guerrero method per location) to forecast weekly ",
      "weighted ILI percentages. The model is fit independently per ",
      "location for each origin date. Prediction intervals are generated ",
      "from 1000 bootstrapped residual paths via fable ",
      "generate(bootstrap=TRUE), with empirical quantiles extracted at ",
      "23 levels and floored at zero. Companion candidate to ets_log: same ",
      "model family, location-adaptive transform (Box-Cox per location vs. ",
      "uniform log+offset). Phase 2 correlation analysis determines whether ",
      "both add value or are twins."
    ),

    ensemble_of_models     = FALSE,
    ensemble_of_hub_models = FALSE
  )
)
