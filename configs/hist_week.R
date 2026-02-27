###############################################################################
# CONFIG: configs/hist_week.R
#
# Model configuration for: Historical Week Climatology (non-fable)
# Engine:  scripts/forecast_engine.R
#
# A simple climatological baseline that uses NO statistical model. For each
# forecast horizon, it identifies the target week-of-year and resamples
# (with replacement) from all historical observations of that same week
# across all available years. This produces simulation paths that reflect
# the historical distribution of ILI activity for a given calendar week.
#
# This is a pure tidyverse implementation — no fable dependency — and
# serves as a reference example of the custom_simulate_fn escape hatch.
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/hist_week.R")
#   result <- run_forecast(config)
#
###############################################################################

config <- list(
  
  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "hist_week",
  
  # ---- Paths ----
  hub_path = "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026",
  
  # ---- Simulation settings ----
  # 1000 resampled values per location-horizon gives smooth quantile estimates.
  n_sims = 1000,
  
  # Need at least 2 years of data so that each week-of-year has multiple
  # historical observations to resample from.
  min_train_weeks = 104,
  
  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,
  
  # ---- Reproducibility ----
  base_seed = 12345,
  
  # ---- Transform ----
  # No preprocessing needed — resampling operates on raw observations.
  transform = "none",
  
  # ---- Model specification ----
  # NULL signals the engine to use custom_simulate_fn instead of the
  # fable model() |> generate() pipeline.
  model_fn = NULL,
  
  # ---- Simulation method ----
  simulate_method = "custom",
  
  # ---- Custom simulation function ----
  # For each location and each forecast horizon (1-4 weeks ahead):
  #   1. Compute the target_end_date (origin + horizon * 7)
  #   2. Find its ISO week number
  #   3. Pull all historical observations for that week-of-year (same location)
  #   4. Resample n_sims values with replacement
  #
  # Returns a tibble matching the engine's expected schema:
  #   location (chr), target_end_date (Date), .sim (dbl)
  custom_simulate_fn = function(train_data, h, n_sims, preprocess_result) {
    
    # train_data is a tsibble filtered to <= origin_date by the engine.
    # Convert to tibble for simpler manipulation.
    train_df <- tibble::as_tibble(train_data)
    
    # The origin_date is the max date in the training data
    origin_date <- max(train_df$target_end_date)
    
    # Target dates for each horizon
    target_dates <- origin_date + (1:h) * 7
    
    # Pre-compute ISO week for every historical observation
    train_df <- train_df |>
      dplyr::mutate(iso_week = lubridate::isoweek(target_end_date))
    
    locations <- unique(train_df$location)
    
    # Build simulation paths for all locations x horizons
    purrr::map_dfr(locations, function(loc) {
      loc_data <- train_df |> dplyr::filter(location == loc)
      
      purrr::map_dfr(seq_along(target_dates), function(horizon_idx) {
        target_date <- target_dates[horizon_idx]
        target_week <- lubridate::isoweek(target_date)
        
        # All historical observations for this week-of-year
        week_pool <- loc_data |>
          dplyr::filter(iso_week == target_week) |>
          dplyr::pull(observation)
        
        # If fewer than 2 observations for this week, fall back to
        # neighboring weeks (+/- 1) to avoid degenerate resampling
        if (length(week_pool) < 2) {
          neighbor_weeks <- c(target_week - 1, target_week, target_week + 1) %% 53
          neighbor_weeks[neighbor_weeks == 0] <- 53
          week_pool <- loc_data |>
            dplyr::filter(iso_week %in% neighbor_weeks) |>
            dplyr::pull(observation)
        }
        
        # Resample with replacement
        sim_values <- sample(week_pool, size = n_sims, replace = TRUE)
        
        tibble::tibble(
          location        = loc,
          target_end_date = target_date,
          .sim            = sim_values
        )
      })
    })
  },
  
  # ---- Metadata ----
  metadata = list(
    team_name    = "Reger",
    model_name   = "Historical Week Climatology",
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
      "Climatological baseline: resamples historical observations from the ",
      "same week-of-year to produce forecast distributions."
    ),
    methods_long        = paste0(
      "A non-parametric climatological baseline that uses no statistical ",
      "model. For each forecast horizon, the target week-of-year is ",
      "identified and 1000 values are resampled with replacement from all ",
      "historical observations of that same ISO week across available years. ",
      "This produces simulation paths reflecting the historical distribution ",
      "of ILI activity for a given calendar week. Empirical quantiles are ",
      "extracted at 23 levels (0.01 to 0.99). Forecast values are floored ",
      "at zero. Implemented in pure tidyverse with no fable dependency, ",
      "serving as a reference example of a custom (non-fable) model."
    ),
    ensemble_of_models      = FALSE,
    ensemble_of_hub_models  = FALSE
  )
)