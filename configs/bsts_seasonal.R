###############################################################################
# FILE: configs/bsts_seasonal.R
#
# PURPOSE:
#   Component model config for the FluSight ILI Sandbox Hub: a Bayesian
#   Structural Time Series (BSTS) model with annual seasonality.
#
#   This is one of the new candidate models added in Phase 1 of the final
#   submission plan (FLUCAST_IMPLEMENTATION_PLAN.md, §4.2 / Checkpoint 1A).
#   Its role in the candidate pool is to provide *genuinely Bayesian*
#   uncertainty quantification (MCMC-based) — distinct from all the existing
#   fable models, which produce uncertainty either parametrically (ETS,
#   ARIMA) or via residual bootstrapping (SNAIVE+BC+BS, etc.). Adding a
#   different uncertainty mechanism is what makes bsts a candidate for
#   ensemble diversity rather than just another point forecaster.
#
# MODEL STRUCTURE:
#   y_t = mu_t + tau_t + epsilon_t
#     mu_t   = local linear trend  (level + slope, both Gaussian random walks)
#     tau_t  = seasonal component  (52-week annual cycle, sum-to-zero constraint)
#     eps_t  ~ N(0, sigma^2)
#
#   This is the classic structural decomposition (Harvey 1989, Scott & Varian
#   2014). It's *additive* and works on the raw wILI scale — bsts does not
#   apply a Box-Cox or log transform internally. We accept the (mild)
#   misspecification of Gaussian errors on right-skewed positive data
#   because the resulting MCMC uncertainty is a different shape from the
#   Box-Cox + bootstrap envelopes our other models produce. That difference
#   is what we want for the ensemble.
#
# INPUTS:
#   - wili_tsibble (passed by engine): tsibble[location, target_end_date,
#     observation], filtered per origin date to target_end_date <= origin_date
#   - hub-config/tasks.json (read by engine to enumerate origin dates)
#
# OUTPUTS:
#   - One CSV per origin date in:
#       model-output/KReger-bsts_seasonal/<origin_date>-KReger-bsts_seasonal.csv
#     Each CSV: 1012 rows (4 horizons * 11 locations * 23 quantiles), columns
#     per the FluSight hub spec.
#
# HOW TO RUN (validation dates only — Phase 1):
#   # From the project root in an R session:
#   source("scripts/forecast_engine.R")
#   source("configs/bsts_seasonal.R")
#   # The driver script handles validation-only filtering:
#   source("scripts/run_bsts_validation.R")
#
# HOW TO RUN (smoke test — 1 date):
#   source("scripts/forecast_engine.R")
#   source("configs/bsts_seasonal.R")
#   result <- run_forecast(config)              # initializes only (run_loop=FALSE)
#   result$forecast(as.Date("2015-10-24"))      # produces one CSV
#
# DEPENDENCIES:
#   - bsts (install.packages("bsts"))  -- pulls in Boom, BoomSpikeSlab
#   - The engine + helpers: scripts/forecast_engine.R, scripts/forecast_helpers.R
#
# RUNTIME EXPECTATIONS:
#   ~30-90 seconds per origin date with these settings (1100 iterations,
#   11 locations sequentially). Validation set = 56 dates -> ~30-90 min total.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md §1.3, §4.2 (Tier 1), Checkpoint 1A
#   - bsts CRAN page: https://cran.r-project.org/package=bsts
#   - Scott, S. L., & Varian, H. R. (2014). "Predicting the present with
#     Bayesian structural time series." IJMM&O 5(1).
###############################################################################


# ============================================================================
# DEPENDENCY CHECK
# ============================================================================
# We don't `library(bsts)` because that would spam the user's namespace with
# a couple hundred functions. Instead we verify it's installed and call its
# functions with `bsts::` qualifiers later.
if (!requireNamespace("bsts", quietly = TRUE)) {
  stop(
    "Package 'bsts' is required but not installed. Install with:\n",
    "  install.packages(\"bsts\")\n",
    "Note: this also installs Boom and BoomSpikeSlab as dependencies."
  )
}


# ============================================================================
# BSTS HYPERPARAMETERS
#
# These live as top-level constants (not on the config list) so that the
# closure inside custom_simulate_fn captures them via R's lexical scoping.
# Putting them here, separate from the rest of the config, makes them easy
# to find and tune in isolation when we revisit MCMC budget vs. runtime.
#
# Naming: ALL_CAPS to flag "constant — don't reassign mid-script".
# ============================================================================

# MCMC iteration count.
# Must satisfy: BSTS_NITER >= BSTS_BURN + n_sims (where n_sims = config$n_sims).
# We pick the minimum that satisfies the constraint to keep runtime down on
# the first pass; increase if posterior draws look too coarse.
BSTS_NITER <- 1100L

# Burn-in iterations to discard. bsts's own default is ~10% of niter, so 100
# out of 1100 matches that convention exactly. The first ~100 draws are the
# sampler finding the typical set of the posterior; keeping them would bias
# the predictive distribution toward the (arbitrary) initial state.
BSTS_BURN <- 100L

# Length of the seasonal period. Weekly wILI -> 52-week annual cycle.
# (Yes, calendar years are 52.18 weeks. The hub uses MMWR weeks, which are
# strictly 52- or 53-week years; bsts handles 53-week years as a 1-week
# rounding error. Acceptable for our application.)
BSTS_NSEASONS <- 52L

# MCMC progress reporting cadence. 0 = silent (we have 56 dates x 11
# locations to fit; verbose progress would be unbearable).
BSTS_PING <- 0L


# ============================================================================
# CONFIG LIST
#
# Schema enforced by the engine (see scripts/forecast_engine.R):
#   team_abbr, model_abbr, hub_path, n_sims, min_train_weeks, dev_mode,
#   run_loop, overwrite, base_seed, transform, model_fn, simulate_method,
#   custom_simulate_fn
# ============================================================================

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "bsts_seasonal",   # 13 chars, within hub's 15-char model_abbr limit

  # ---- Paths ----
  # Absolute path to the local clone of the ATSF2026 repo. Edit if you move
  # the repo or run on another machine. Once the model is final, this gets
  # parameterized for the README_submission reproducibility step (Phase 5A).
  hub_path = normalizePath(".", mustWork = FALSE),

  # ---- Simulation settings ----
  # n_sims is the number of posterior predictive draws we keep for each
  # (location, horizon) combination. Must equal BSTS_NITER - BSTS_BURN.
  n_sims = 1000L,

  # ---- Training window requirement ----
  # bsts with nseasons=52 needs at least 2 full annual cycles for the
  # seasonal component's prior to be informative. We have ~12 years of
  # pre-hub data, so 104 is a defensive floor that won't ever bite in
  # practice but protects against accidentally calling bsts on stub data.
  min_train_weeks = 104L,

  # ---- Run control ----
  # dev_mode = FALSE: don't subset dates inside the engine — we control
  # which dates run via the driver script (run_bsts_validation.R).
  dev_mode  = FALSE,

  # run_loop = FALSE: don't auto-loop. The driver script invokes
  # result$forecast(d) for each validation date manually. This lets us
  # constrain to validation-only without touching the engine.
  run_loop  = FALSE,

  # overwrite = FALSE: skip dates that already have CSVs. Safe to rerun
  # the driver after a partial completion — only the missing dates redo.
  overwrite = FALSE,

  # ---- Reproducibility ----
  # base_seed: the engine sets set.seed(base_seed + as.integer(origin_date))
  # before each date's forecast, giving deterministic per-date draws.
  # bsts has its own internal RNG too — we draw bsts's seed from R's RNG
  # inside custom_simulate_fn so the whole pipeline is reproducible.
  base_seed = 12345L,

  # ---- Transform ----
  # "none" because bsts works on the raw scale. The engine's preprocess
  # step won't apply Box-Cox or log; bsts's Gaussian innovations are an
  # intentional design choice to differentiate this model from the
  # transform-based candidates in the pool.
  transform = "none",

  # ---- Model spec ----
  # Custom path: model_fn = NULL signals to the engine that we're NOT a
  # fable model. The engine then routes simulation through custom_simulate_fn.
  model_fn        = NULL,
  simulate_method = "custom",

  #' Generate posterior predictive samples for one origin date.
  #'
  #' Called once per origin date by forecast_engine.R's
  #' generate_single_forecast(). Inside, we loop over the 11 hub locations,
  #' fit an independent BSTS model to each location's series, draw posterior
  #' predictive samples for h weeks ahead, and return a long-format tibble
  #' that the engine converts to quantiles and writes to CSV.
  #'
  #' @param train_data tsibble. Columns: location (chr), target_end_date
  #'   (Date), observation (dbl). Already filtered by the engine to
  #'   target_end_date <= origin_date. Includes all 11 hub locations.
  #' @param h integer. Forecast horizon in weeks (the hub spec asks for 4).
  #' @param n_sims integer. Number of posterior draws to return per
  #'   (location, horizon) cell. Must match BSTS_NITER - BSTS_BURN.
  #' @param preprocess_result list or NULL. Output of config$preprocess_fn
  #'   if defined; here we don't define one, so this is NULL. Argument
  #'   exists for engine-interface consistency with other custom configs
  #'   (e.g., bma_ensemble.R) that DO use a preprocess_fn.
  #' @return tibble with columns: location (chr), target_end_date (Date),
  #'   .sim (dbl). One row per (location, horizon, draw). Total rows =
  #'   length(unique(train_data$location)) * h * n_sims.
  #' @examples
  #'   # Not callable directly — invoked by the engine. To exercise:
  #'   # result <- run_forecast(config)
  #'   # result$forecast(as.Date("2015-10-24"))
  custom_simulate_fn = function(train_data, h, n_sims, preprocess_result) {

    # ------------------------------------------------------------------------
    # Sanity check: the engine should pass exactly n_sims = BSTS_NITER - BSTS_BURN.
    # If somebody bumps n_sims in the config without bumping BSTS_NITER, this
    # catches the mismatch before we waste ~30 sec/loc fitting bsts.
    # ------------------------------------------------------------------------
    expected_draws <- BSTS_NITER - BSTS_BURN
    if (n_sims != expected_draws) {
      stop(sprintf(
        "n_sims (%d) != BSTS_NITER - BSTS_BURN (%d - %d = %d). Update either ",
        n_sims, BSTS_NITER, BSTS_BURN, expected_draws),
        "BSTS_NITER at top of bsts_seasonal.R or n_sims in config to match."
      )
    }

    # ------------------------------------------------------------------------
    # Identify the locations we need to forecast. The hub expects all 11 in
    # every CSV, so we expect length(locations) == 11. If a location is
    # missing from train_data the engine's downstream validation will fail
    # loudly — no need to check here.
    # ------------------------------------------------------------------------
    locations <- unique(train_data$location)

    # ------------------------------------------------------------------------
    # Loop over locations. bsts is a univariate model; there is no native
    # panel/multi-series support, so each location gets its own fit and
    # MCMC chain. purrr::map keeps things tidy and lets us drop in
    # furrr::future_map later for parallelism if 30-90 sec/date is too slow.
    # ------------------------------------------------------------------------
    location_results <- purrr::map(locations, function(loc) {

      # Pull this location's series, sorted ascending by date. arrange() is
      # defensive — train_data should already be sorted, but bsts requires
      # a strictly chronological vector and silently misbehaves otherwise.
      loc_data <- train_data |>
        dplyr::filter(location == loc) |>
        dplyr::arrange(target_end_date)

      # The raw wILI vector. bsts works on plain numeric vectors — no need
      # to preserve the tsibble structure inside the model.
      y <- loc_data$observation

      # Anchor for forecast date arithmetic. The engine guarantees the most
      # recent training observation is on the origin_date Saturday (the
      # hub's "as-of" semantics), so last_date == origin_date and our
      # forecast horizons land on origin_date + 7, +14, +21, +28.
      last_date <- max(loc_data$target_end_date)

      # ----------------------------------------------------------------------
      # State specification.
      #
      # bsts builds the state spec by chaining Add* functions onto a list.
      # Each call appends one component:
      #   AddLocalLinearTrend: level + slope, both Gaussian random walks.
      #     Captures the slowly-drifting baseline of ILI.
      #   AddSeasonal:        sum-to-zero seasonal dummies, one per nseasons.
      #     Captures the strong annual cycle (winter peak, summer trough).
      #
      # Both functions inspect `y` to set sane priors on the variance
      # hyperparameters (Inverse-Gamma scaled by sd(y)). That's why we pass
      # `y` to both of them rather than relying on defaults.
      # ----------------------------------------------------------------------
      ss <- list()
      ss <- bsts::AddLocalLinearTrend(ss, y)
      ss <- bsts::AddSeasonal(ss, y, nseasons = BSTS_NSEASONS)

      # ----------------------------------------------------------------------
      # bsts has its own RNG path independent of base R's. We sample one
      # integer from R's RNG (which the engine has already seeded
      # deterministically per date) to act as bsts's seed. This way:
      #   - rerunning the same date gives identical posterior draws
      #   - different dates get different (but deterministic) bsts seeds
      # ----------------------------------------------------------------------
      bsts_seed <- as.integer(stats::runif(1, 1, .Machine$integer.max))

      # ----------------------------------------------------------------------
      # Fit the model. Wrap in tryCatch so a single-location failure gives
      # us a clear, attributable error rather than a cryptic stack trace.
      # ----------------------------------------------------------------------
      fit <- tryCatch(
        bsts::bsts(
          formula = y,
          state.specification = ss,
          niter = BSTS_NITER,
          ping  = BSTS_PING,
          seed  = bsts_seed
        ),
        error = function(e) {
          stop(sprintf(
            "bsts fit failed for location='%s', last_date=%s, n_obs=%d: %s",
            loc, format(last_date), length(y), conditionMessage(e)
          ), call. = FALSE)
        }
      )

      # ----------------------------------------------------------------------
      # Posterior predictive draws for the next h weeks.
      #
      # predict.bsts(burn = ...) discards the first `burn` MCMC iterations
      # internally — we don't need to do it ourselves. The returned object
      # has $distribution: a matrix with rows = retained draws, cols = horizon.
      # ----------------------------------------------------------------------
      pred <- predict(
        fit,
        horizon   = h,
        burn      = BSTS_BURN,
        # quantiles= argument here only affects the printed summary,
        # not what's in $distribution. Pass anything reasonable.
        quantiles = c(0.025, 0.975)
      )

      # Defensive shape check. If bsts ever changes its convention, fail
      # loudly here instead of producing silently-wrong forecasts.
      stopifnot(
        "bsts returned wrong number of draws" =
          nrow(pred$distribution) == n_sims,
        "bsts returned wrong horizon" =
          ncol(pred$distribution) == h
      )

      # ----------------------------------------------------------------------
      # Build the forecast date sequence: last_date + 7, +14, ..., +h*7.
      # All ATSF2026 hub dates are Saturdays, so weekly stepping is exact.
      # ----------------------------------------------------------------------
      forecast_dates <- last_date + (1:h) * 7L

      # ----------------------------------------------------------------------
      # Reshape the [n_sims x h] matrix to long format: one row per
      # (sim, horizon).
      #
      # as.vector(M) traverses an R matrix in column-major order, so the
      # first n_sims entries correspond to horizon 1, the next n_sims to
      # horizon 2, etc. rep(forecast_dates, each = n_sims) constructs the
      # date column to match: [d1 x n_sims, d2 x n_sims, d3 x n_sims, d4 x n_sims].
      # ----------------------------------------------------------------------
      tibble::tibble(
        location        = loc,
        target_end_date = rep(forecast_dates, each = n_sims),
        .sim            = as.vector(pred$distribution)
      )
    })

    # Combine all 11 location tibbles into one long-format result. The engine
    # takes it from here: groups by (location, target_end_date), computes
    # the 23 hub quantile levels via compute_quantiles_from_sims(), formats,
    # validates, writes CSV.
    dplyr::bind_rows(location_results)
  },

  # ============================================================================
  # METADATA
  #
  # Drives auto-generation of model-metadata/KReger-bsts_seasonal.yml. The
  # engine's validate_config() step requires this list, and write_metadata()
  # serializes it to the hub-spec YAML format (per model-metadata/README.md
  # and hub-config/model-metadata-schema.json upstream).
  #
  # Field references:
  #   team_abbr / model_abbr           — top-level on `config`, not here
  #   team_name                        — full project team name
  #   model_name                       — human-readable model name
  #   model_contributors               — list of contributor info
  #   license                          — CC-BY-4.0 (matches other KReger models)
  #   designated_model                 — FALSE: this is a CANDIDATE component,
  #                                      not the team's primary submission. The
  #                                      ensemble (KReger-bma_ensemble) will be
  #                                      designated_model = TRUE when built.
  #   data_inputs                      — what the model consumes
  #   methods (<200 chars)             — short summary
  #   methods_long                     — detailed description
  #   ensemble_of_models               — FALSE: single model
  #   ensemble_of_hub_models           — FALSE: doesn't combine other hub models
  # ============================================================================
  metadata = list(
    team_name  = "Reger",
    model_name = "Bayesian Structural Time Series with Annual Seasonality",

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

    # Short methods string. Hub schema enforces <200 chars.
    methods = paste0(
      "Bayesian Structural Time Series with local linear trend and 52-week ",
      "seasonal component. MCMC posterior predictive samples yield ",
      "empirical quantiles."
    ),

    # Detailed methods description. No length limit, but keep it focused
    # on what's distinctive about this model in the candidate pool.
    methods_long = paste0(
      "Univariate state-space model fit independently per location using the ",
      "bsts R package (Scott & Varian 2014). State specification: ",
      "AddLocalLinearTrend (level + slope, both Gaussian random walks) plus ",
      "AddSeasonal with nseasons=52 (annual cycle, sum-to-zero constraint). ",
      "MCMC sampler runs 1100 iterations with 100 burn-in, yielding 1000 ",
      "posterior predictive draws per (location, horizon). Forecasts produced ",
      "via predict.bsts() on the raw wILI scale (no Box-Cox transform). The ",
      "23 hub quantile levels are computed empirically from the posterior ",
      "draws and floored at zero. Distinct from the other KReger candidates: ",
      "uncertainty arises from MCMC over the latent state rather than from ",
      "parametric ARIMA/ETS innovations or bootstrap residual resampling."
    ),

    ensemble_of_models     = FALSE,
    ensemble_of_hub_models = FALSE
  )
)
