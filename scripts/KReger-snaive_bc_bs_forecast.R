###############################################################################
# FILE: scripts/KReger-snaive_bc_bs_forecast.R
#
# PURPOSE:
#   Model-specific forecasting pipeline for the KReger-snaive_bc_bs submission
#   to the FluSight ILI Sandbox Hub (ATSF2026).
#
# MODEL:
#   SNAIVE + Box-Cox + Bootstrapped Residuals
#   fable specification: SNAIVE(box_cox(observation, lambda) ~ lag('year'))
#
#   - Box-Cox transform stabilizes variance (high winter peaks vs low summer)
#   - SNAIVE with yearly lag uses same-week-last-year as the point forecast
#   - Bootstrapped residuals via generate(bootstrap = TRUE) produce simulated
#     future paths from which we compute empirical quantiles
#
#   This is the same approach as the instructor's atsf-meanbs model, but
#   substituting SNAIVE for MEAN.
#
# ARCHITECTURE:
#   This script sources forecast_helpers.R for all model-agnostic utilities
#   (data loading, formatting, validation, I/O). The ONLY model-specific
#   logic is in this file:
#     1. Lambda estimation via guerrero() — this section (Milestone 2)
#     2. generate_single_forecast() — the fit/generate/quantile function
#        (will be added in Milestone 3)
#     3. The loop over all 139 origin dates (will be added in Milestone 5)
#
#   To implement a DIFFERENT model in the future, copy this file, change
#   the model() call inside generate_single_forecast(), adjust preprocessing
#   (e.g., lambda might not be needed), and update naming constants.
#   Everything else (formatting, validation, writing) is inherited from
#   forecast_helpers.R.
#
# HOW TO RUN:
#   1. Open RStudio
#   2. Set working directory to the ATSF2026 repo root:
#      setwd("/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026")
#   3. source("scripts/KReger-snaive_bc_bs_forecast.R")
#
# INPUTS:
#   - target-data/time-series.csv (wILI time series)
#   - hub-config/tasks.json (origin dates)
#
# OUTPUTS:
#   - 139 CSV files in model-output/KReger-snaive_bc_bs/
#   - Console progress messages
#
# DESIGN DOC REFERENCES:
#   - Section 5: Baseline Modeling Approach (model spec, lambda strategy)
#   - Section 6: End-to-End Workflow (steps 0-5)
#   - Section 7.4: Extensibility Architecture
#   - Section 10: Milestones 2-5
#
###############################################################################


# =============================================================================
# LOAD PACKAGES
#
# We load fpp3 (which brings in fable, fabletools, tsibble, feasts) and
# tidyverse (which brings in dplyr, readr, ggplot2, purrr, etc.).
# These are the same packages used in w5d1-solutions.qmd.
# =============================================================================

library(fpp3)
library(tidyverse)


# =============================================================================
# CONFIGURATION — Model-specific constants
#
# All configurable values are gathered here at the top so you can see
# everything that might need changing in one place. Constants that are
# shared across ALL models (quantile levels, locations, etc.) live in
# forecast_helpers.R instead.
# =============================================================================

# --- Team and model identifiers ---
# These must match the directory name and filename convention exactly.
# Reference: flucast-design-document.md Section 4.1
TEAM_MODEL <- "KReger-snaive_bc_bs"

# --- Path to the local ATSF2026 repo clone ---
# This is the root of the hub directory containing hub-config/, model-output/,
# target-data/, etc. All other paths are constructed relative to this.
HUB_PATH <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"

# --- Derived paths (constructed from HUB_PATH and TEAM_MODEL) ---
# These follow the directory layout in design doc Section 7.3.
DATA_PATH    <- file.path(HUB_PATH, "target-data", "time-series.csv")
TASKS_PATH   <- file.path(HUB_PATH, "hub-config", "tasks.json")
OUTPUT_DIR   <- file.path(HUB_PATH, "model-output", TEAM_MODEL)

# --- Bootstrap replicate count ---
# 1000 replicates for the production run. This gives at least 10 observations
# in each tail for the 0.01 and 0.99 quantiles (the most extreme levels),
# producing stable estimates.
#
# For DEVELOPMENT/TESTING: change to 100 for ~10x faster runs.
# For PRODUCTION: use 1000 (the value below).
#
# Reference: flucast-design-document.md Section 5.3
N_BOOTSTRAP <- 1000   # <-- Change to 100 for faster test runs

# --- Dev mode flag ---
# When TRUE, only processes the first 5 origin dates (for quick testing).
# When FALSE, processes all 139 origin dates (production run).
# Tip: use DEV_MODE = TRUE with N_BOOTSTRAP = 100 for a fast sanity check
# (~30 seconds total), then switch both for production.
DEV_MODE <- FALSE     # <-- Set to TRUE for testing, FALSE for production

# --- Run forecast loop flag ---
# When TRUE, the main forecast loop (Milestone 5) executes automatically
# when the script is source()'d. When FALSE, only the setup (Steps 1-3)
# runs, and you can call generate_single_forecast() manually for testing.
RUN_FORECAST_LOOP <- TRUE  # <-- Set to TRUE to run the full pipeline


# =============================================================================
# SOURCE SHARED HELPERS
#
# forecast_helpers.R provides all model-agnostic functions:
#   load_and_prepare_tsibble(), get_origin_dates(),
#   compute_quantiles_from_sims(), format_for_hub(),
#   validate_forecast_df(), write_hub_csv(), plot_quantile_fan()
#
# It also provides the shared constants:
#   DATE_COL, QUANTILE_LEVELS, HUB_LOCATIONS, HUB_HORIZONS, EXPECTED_ROWS
#
# Reference: flucast-design-document.md Section 7.4
# =============================================================================

source(file.path(HUB_PATH, "scripts", "forecast_helpers.R"))


# =============================================================================
# STEP 1: LOAD AND PREPARE DATA
#
# Read the wILI target time series and convert to a tsibble.
# This is the ground truth data that models are trained on.
#
# Reference: flucast-design-document.md Section 6, Steps 1-2
#            w5d1-solutions.qmd Exercise 1
# =============================================================================

cat("=== STEP 1: Loading target data ===\n")
wili_tsibble <- load_and_prepare_tsibble(DATA_PATH)

# Quick summary to confirm data loaded correctly
cat(paste0("  Rows: ", nrow(wili_tsibble), "\n"))
cat(paste0("  Date range: ",
           min(wili_tsibble$target_end_date), " to ",
           max(wili_tsibble$target_end_date), "\n"))
cat(paste0("  Locations: ",
           length(unique(wili_tsibble$location)), "\n\n"))


# =============================================================================
# STEP 2: LOAD ORIGIN DATES
#
# Parse the full list of origin dates from hub-config/tasks.json.
# These are the 139 Saturdays we need to generate forecasts for.
#
# Reference: flucast-design-document.md Section 6, Step 4
# =============================================================================

cat("=== STEP 2: Loading origin dates ===\n")
origin_dates <- get_origin_dates(TASKS_PATH)
cat("\n")


# =============================================================================
# STEP 3: ESTIMATE BOX-COX LAMBDA (per location)
#
# WHY BOX-COX?
#   ILI percentages have strongly heteroskedastic variance: during winter
#   flu peaks, the values spike to 5-10% with high variability, while during
#   summer they hover near 1% with low variability. If we fit a model on the
#   raw values, the residuals during peak season will dominate, producing
#   bootstrap intervals that are too wide in summer and too narrow in winter.
#
#   The Box-Cox transform y' = (y^lambda - 1) / lambda (for lambda != 0)
#   or y' = log(y) (for lambda = 0) stabilizes this variance so that
#   the residuals are more homoskedastic. fable's SNAIVE() applies the
#   transform before fitting and automatically back-transforms when
#   generating forecasts.
#
# WHY GUERRERO?
#   The guerrero() method (Guerrero, 1993) automatically selects the
#   optimal lambda that minimizes the coefficient of variation of
#   subseries (seasonal or non-overlapping time windows). It's the
#   standard approach in fpp3 — see w3d2-solutions.qmd for examples.
#
# WHY ESTIMATE ONCE (not per origin date)?
#   - Efficiency: computing guerrero 139 times × 11 locations = 1,529
#     extra feature computations, with minimal benefit for a baseline model.
#   - Stability: using data before the first origin date means the lambda
#     is fixed and reproducible regardless of which origin date is being
#     processed. This is the strategy recommended in design doc Section 5.2.
#
# HOW IT WORKS:
#   1. Filter the tsibble to keep only data BEFORE the first origin date.
#      This is the "pre-forecast" historical data.
#   2. For each location, compute features(observation, guerrero) which
#      returns a tibble with one row per location and a lambda_guerrero column.
#   3. Store as a named vector: names = location names, values = lambda.
#
# Reference: flucast-design-document.md Section 5.2
#            w3d2-solutions.qmd — features(GDP, features = guerrero)
#            04_transformations.pdf — "Transformations can stabilize variance"
# =============================================================================

cat("=== STEP 3: Estimating Box-Cox lambdas via guerrero ===\n")

# --- Filter to pre-forecast data ---
# We use data strictly BEFORE the first origin date. The first origin date
# is the earliest Saturday we forecast from, so all data up to (but not
# including) that date is "historical" for lambda estimation purposes.
#
# We actually use <= (first_origin - 1 day) to exclude the first origin date
# itself, since on that date we'd be making a forecast, not training.
first_origin <- min(origin_dates)
cat(paste0("  Using data before first origin date: ", first_origin, "\n"))

pre_forecast_data <- wili_tsibble |>
  filter(target_end_date < first_origin)

cat(paste0("  Pre-forecast rows: ", nrow(pre_forecast_data), "\n"))
cat(paste0("  Pre-forecast date range: ",
           min(pre_forecast_data$target_end_date), " to ",
           max(pre_forecast_data$target_end_date), "\n"))

# --- Compute guerrero lambda for each location ---
# features() is a fabletools function that computes summary statistics
# ("features") for each time series in a tsibble. Here we compute just
# one feature: the guerrero lambda.
#
# The result is a tibble with columns: location, target, lambda_guerrero.
# We extract location and lambda_guerrero, then create a named vector.
#
# Pattern from: w3d2-solutions.qmd
#   us_economy |> features(GDP, features = guerrero)
lambda_table <- pre_forecast_data |>
  features(observation, features = guerrero)

# Display the full lambda table for visual inspection
cat("\n  Guerrero lambda estimates by location:\n")
print(lambda_table, n = Inf)

# --- Create a named vector for easy lookup ---
# This lets us do: lambdas["US National"] to get the lambda for that location.
# The model-fitting step (Milestone 3) will use this to pass the correct
# lambda to SNAIVE(box_cox(observation, lambda)) for each location.
lambdas <- lambda_table$lambda_guerrero
names(lambdas) <- lambda_table$location

# --- Sanity checks ---

# Check 1: We should have exactly 11 lambdas (one per location)
if (length(lambdas) != length(HUB_LOCATIONS)) {
  warning(
    "Expected ", length(HUB_LOCATIONS), " lambdas but got ", length(lambdas),
    ". Missing locations: ",
    paste(setdiff(HUB_LOCATIONS, names(lambdas)), collapse = ", ")
  )
}

# Check 2: Lambdas should be in a "reasonable" range
# The Box-Cox lambda typically falls between -2 and 2. Values outside
# this range suggest something unusual about the data (extreme skew,
# outliers, or numerical issues in the guerrero optimization).
#
# If any lambda is extreme, we warn but don't fail — the user can
# investigate and decide whether to use log1p() as a fallback for
# that location (as suggested in design doc Section 5.2).
extreme_lambdas <- lambdas[abs(lambdas) > 2]
if (length(extreme_lambdas) > 0) {
  warning(
    "Some lambdas are outside [-2, 2] — consider using log1p() as fallback:\n",
    paste0("  ", names(extreme_lambdas), ": ", round(extreme_lambdas, 4),
           collapse = "\n")
  )
} else {
  cat("\n  All lambdas are in the reasonable range [-2, 2].\n")
}

# Check 3: If any lambda is very close to 0, note it
# Lambda near 0 means the Box-Cox transform is approximately log().
# This is fine and common for right-skewed data like ILI percentages.
near_zero <- lambdas[abs(lambdas) < 0.05]
if (length(near_zero) > 0) {
  cat(paste0("  Note: ", length(near_zero),
             " location(s) have lambda near zero (≈ log transform): ",
             paste(names(near_zero), collapse = ", "), "\n"))
}

cat("\n  Lambda estimation complete.\n")
cat(paste0("  Range: [", round(min(lambdas), 4), ", ",
           round(max(lambdas), 4), "]\n\n"))


# =============================================================================
# STEP 4: SINGLE-DATE FORECAST FUNCTION (Milestone 3)
#
# This is the core function that produces a complete hub-formatted forecast
# for ONE origin date across all 11 locations. It implements Steps 5a–5f
# from the design doc (Section 6).
#
# The function is designed to be called in a loop over all 139 origin dates
# (Milestone 5), but can also be called standalone for testing.
#
# Reference: flucast-design-document.md Section 6, Steps 5a–5f
#            flucast-implementation-prompt.md — Milestone 3 spec
# =============================================================================

#' Generate a complete hub-formatted forecast for a single origin date
#'
#' Runs the full pipeline: filter training data → fit SNAIVE with Box-Cox →
#' generate bootstrap paths → compute quantiles → format for hub → validate →
#' optionally write CSV. Returns the formatted dataframe.
#'
#' @param origin_date A single Date — the Saturday to forecast from.
#' @param wili_tsibble The full wILI tsibble (all dates, all locations).
#' @param lambdas A named numeric vector: names = location strings,
#'   values = guerrero lambda for each location.
#' @param quantile_levels Numeric vector of quantile probabilities.
#'   Defaults to the hub's 23 required levels.
#' @param n_bootstrap Integer. Number of bootstrap replicates for generate().
#'   Use 100 for dev/testing, 1000 for production.
#' @param team_model Character string (e.g., "KReger-snaive_bc_bs").
#' @param output_dir Character string. Path to write the CSV file.
#' @param write_file Logical. If TRUE, writes the CSV. If FALSE, returns
#'   the dataframe without writing (useful for inspection).
#'
#' @return A tibble with 1,012 rows × 8 columns in hub format, or NULL
#'   if the forecast could not be generated (insufficient data, model error).
generate_single_forecast <- function(origin_date,
                                     wili_tsibble,
                                     lambdas,
                                     quantile_levels = QUANTILE_LEVELS,
                                     n_bootstrap = N_BOOTSTRAP,
                                     team_model = TEAM_MODEL,
                                     output_dir = OUTPUT_DIR,
                                     write_file = TRUE) {
  
  # =========================================================================
  # STEP 5a: FILTER TRAINING DATA
  #
  # Keep only data where target_end_date <= origin_date. This simulates
  # standing on the origin Saturday and looking backward — we can see
  # all data up through and including that date, but nothing after.
  #
  # The <= (not <) is deliberate: the origin_date Saturday is the last
  # day of the reporting epiweek, so data for that week is available
  # when we make the forecast.
  #
  # Pattern from: w5d1-solutions.qmd Exercise 2
  #   wili_tsibble |> filter(target_end_date < '2016-12-05')
  # =========================================================================
  
  train_data <- wili_tsibble |>
    filter(target_end_date <= origin_date)
  
  # --- Check for minimum data requirement ---
  # SNAIVE with lag('year') needs at least 52 weeks (1 full year) of history
  # to produce the seasonal naive forecast. If any location has fewer than
  # 52 weeks, the model will error.
  min_weeks <- train_data |>
    as_tibble() |>
    group_by(location) |>
    summarise(n = n(), .groups = "drop") |>
    pull(n) |>
    min()
  
  if (min_weeks < 52) {
    warning(
      "Skipping ", origin_date, ": only ", min_weeks,
      " weeks of training data (need >= 52 for SNAIVE yearly lag)"
    )
    return(NULL)
  }
  
  # =========================================================================
  # STEP 5b + 5c: FIT MODEL AND GENERATE BOOTSTRAP PATHS (per location)
  #
  # WHY A PER-LOCATION LOOP?
  #   We originally tried passing lambda as a column in the tsibble and
  #   fitting all 11 locations in a single model() call. This failed:
  #   fable's box_cox() stores the lambda value during model(), but
  #   generate() creates new future rows that don't carry the lambda
  #   column. Result: the back-transformation silently fails and we get
  #   forecasts stuck in Box-Cox space (~0.5) instead of real ILI space
  #   (~2-4%). See the Milestone 3 bug fix discussion.
  #
  #   The fix: loop over locations with purrr::map(), passing each
  #   location's lambda as a SCALAR to box_cox(). Fable correctly stores
  #   scalar lambdas in the model object and back-transforms properly.
  #
  # WHAT THE MODEL DOES:
  #   - box_cox(observation, loc_lambda): applies the Box-Cox variance-
  #     stabilizing transform before fitting. When loc_lambda is a scalar,
  #     fable stores it in the model and automatically back-transforms
  #     when generating forecasts.
  #   - SNAIVE(... ~ lag('year')): the seasonal naive model uses the value
  #     from the same week one year ago (52 weeks back) as the point
  #     forecast. This captures the seasonal pattern of flu — "next week's
  #     ILI will look like the same week last year."
  #   - generate(bootstrap = TRUE): resamples from fitted residuals to
  #     create simulated future paths. Each replicate is one possible
  #     future trajectory of wILI over the next 4 weeks.
  #
  # Reference: flucast-design-document.md Section 5.1
  #            w5d1-solutions.qmd Exercise 5
  #            06_forecast_submissions.pdf — "Getting quantile forecasts"
  # =========================================================================
  
  # Get the unique locations in the training data
  locations <- unique(train_data$location)
  
  # Loop over each location: fit model with scalar lambda, generate sims.
  # purrr::map() returns a list of tibbles; bind_rows() combines them.
  sims <- tryCatch(
    {
      purrr::map(locations, function(loc) {
        
        # --- Get this location's scalar lambda ---
        # lambdas is a named vector: lambdas["US National"] returns a scalar
        loc_lambda <- lambdas[loc]
        
        # --- Filter training data to this location only ---
        loc_train <- train_data |>
          filter(location == loc)
        
        # --- Fit the model ---
        # Because loc_lambda is a scalar, fable stores it in the model
        # object and correctly back-transforms during generate().
        loc_fit <- loc_train |>
          model(
            snaive_bc_bs = SNAIVE(
              box_cox(observation, loc_lambda) ~ lag("year")
            )
          )
        
        # --- Generate bootstrap paths ---
        # h = 4: forecast 4 weeks ahead (horizons 1-4)
        # times = n_bootstrap: number of simulated paths (100 dev, 1000 prod)
        # bootstrap = TRUE: resample from residuals (not parametric)
        #
        # Output: a tsibble with columns:
        #   location, target, .model, target_end_date, .rep, .sim
        # Rows: 4 horizons × n_bootstrap replicates
        loc_sims <- loc_fit |>
          generate(h = 4, times = n_bootstrap, bootstrap = TRUE)
        
        # Return as a plain tibble for easier binding across locations
        as_tibble(loc_sims)
        
      }) |>
        # Combine all 11 location results into one dataframe
        bind_rows()
    },
    error = function(e) {
      warning("Fit/generate failed for ", origin_date, ": ", e$message)
      NULL
    }
  )
  
  # If fitting or generation failed, bail out gracefully
  if (is.null(sims)) return(NULL)
  
  # --- Check for NA simulations ---
  # If there are gaps in the time series or numerical issues, some .sim
  # values might be NA. This would produce NA quantiles downstream.
  # Also check for negative .sim values — the Box-Cox back-transform can
  # occasionally produce negatives when bootstrap paths go extreme.
  # Reference: flucast-design-document.md Section 9.4 — edge case checks
  na_count <- sum(is.na(sims$.sim))
  if (na_count > 0) {
    warning(
      origin_date, ": ", na_count, " NA values in bootstrap simulations. ",
      "These will be excluded from quantile computation."
    )
  }
  
  # =========================================================================
  # STEP 5d: COMPUTE QUANTILES
  #
  # Group the simulated paths by location and target_end_date, then compute
  # the empirical quantile at each of the 23 required levels. This converts
  # n_bootstrap simulated values into a probabilistic forecast described
  # by its quantile function.
  #
  # Uses compute_quantiles_from_sims() from forecast_helpers.R, which also
  # applies pmax(value, 0) to enforce non-negativity.
  #
  # Reference: flucast-design-document.md Section 5.4, steps 3-4
  # =========================================================================
  
  quantile_df <- compute_quantiles_from_sims(sims, quantile_levels)
  
  # =========================================================================
  # STEP 5e: FORMAT FOR HUB
  #
  # Add the remaining hub columns (origin_date, target, horizon, output_type),
  # enforce constraints, and select exactly the 8 required columns.
  #
  # Uses format_for_hub() from forecast_helpers.R.
  #
  # Reference: flucast-design-document.md Section 4.2
  # =========================================================================
  
  hub_df <- format_for_hub(quantile_df, origin_date, team_model)
  
  # =========================================================================
  # STEP 5e (continued): VALIDATE
  #
  # Run all 14+ programmatic quality checks BEFORE writing the file.
  # If any check fails, validate_forecast_df() will stop() with a
  # descriptive error message.
  #
  # Uses validate_forecast_df() from forecast_helpers.R.
  #
  # Reference: flucast-design-document.md Section 9.1
  # =========================================================================
  
  validate_forecast_df(hub_df, origin_date, quantile_levels)
  
  # =========================================================================
  # STEP 5f: WRITE CSV
  #
  # Write the formatted dataframe to a CSV with the filename pattern:
  #   YYYY-MM-DD-KReger-snaive_bc_bs.csv
  #
  # Only writes if write_file = TRUE (default). Setting write_file = FALSE
  # is useful for testing — you get the dataframe back without creating files.
  #
  # Uses write_hub_csv() from forecast_helpers.R.
  #
  # Reference: model-output/README.md — filename format
  # =========================================================================
  
  if (write_file) {
    filepath <- write_hub_csv(hub_df, origin_date, output_dir, team_model)
    cat(paste0("  Wrote: ", basename(filepath), "\n"))
  }
  
  # Return the formatted dataframe (useful for inspection and plotting)
  return(hub_df)
}


# =============================================================================
# STEP 5: MAIN FORECAST LOOP (Milestone 5)
#
# Generates forecasts for ALL origin dates (or just the first 5 in DEV_MODE).
# This is the production pipeline that creates the 139 CSV files for
# submission to the FluSight ILI Sandbox Hub.
#
# Features:
#   - Progress reporting with ETA
#   - Skip logic: resumes interrupted runs (skips existing CSVs)
#   - Error handling: logs failures, continues to next date
#   - Summary report at the end
#
# Reference: flucast-design-document.md Section 6, Step 5
#            flucast-design-document.md Milestone 5
# =============================================================================

if (RUN_FORECAST_LOOP) {
  
  cat("\n=== STEP 5: Running forecast loop ===\n")
  
  # --- Determine which dates to process ---
  # In DEV_MODE, only process the first 5 dates for a quick sanity check.
  # In production mode, process all 139 origin dates.
  if (DEV_MODE) {
    dates_to_process <- origin_dates[1:min(5, length(origin_dates))]
    cat("  DEV_MODE = TRUE: processing first", length(dates_to_process),
        "dates only\n")
    cat("  N_BOOTSTRAP =", N_BOOTSTRAP, "\n\n")
  } else {
    dates_to_process <- origin_dates
    cat("  PRODUCTION MODE: processing all", length(dates_to_process),
        "dates\n")
    cat("  N_BOOTSTRAP =", N_BOOTSTRAP, "\n\n")
  }
  
  # --- Create output directory if it doesn't exist ---
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
    cat("  Created output directory:", OUTPUT_DIR, "\n")
  }
  
  # --- Initialize tracking variables ---
  n_total     <- length(dates_to_process)   # Total dates to process
  n_success   <- 0                          # Counter: successful forecasts
  n_skipped   <- 0                          # Counter: skipped (CSV exists)
  n_failed    <- 0                          # Counter: failed forecasts
  failed_log  <- character(0)               # Stores "date: error message"
  loop_start  <- Sys.time()                 # Wall clock start time
  
  # --- Main loop ---
  for (i in seq_along(dates_to_process)) {
    
    this_date <- dates_to_process[i]
    
    # --- Skip logic ---
    # If the CSV already exists for this date, skip it. This allows
    # resuming an interrupted run without regenerating existing files.
    # To force regeneration, delete the file first.
    expected_filename <- paste0(this_date, "-", TEAM_MODEL, ".csv")
    expected_filepath <- file.path(OUTPUT_DIR, expected_filename)
    
    if (file.exists(expected_filepath)) {
      cat(sprintf("  [%3d/%d] %s — SKIPPED (file exists)\n",
                  i, n_total, this_date))
      n_skipped <- n_skipped + 1
      next
    }
    
    # --- Progress reporting ---
    # Show current date, elapsed time, and estimated time remaining.
    elapsed <- as.numeric(difftime(Sys.time(), loop_start, units = "secs"))
    # Only estimate ETA after the first completed forecast
    dates_done <- n_success + n_failed
    if (dates_done > 0) {
      avg_per_date <- elapsed / dates_done
      remaining    <- (n_total - i) * avg_per_date
      eta_str      <- sprintf(" | ETA: %.0f min", remaining / 60)
    } else {
      eta_str <- ""
    }
    
    cat(sprintf("  [%3d/%d] %s%s\n", i, n_total, this_date, eta_str))
    
    # --- Generate forecast with error handling ---
    # suppressWarnings silences the routine Box-Cox back-transform warnings
    # from fable (e.g., "1 simulation replaced with NA"). These are harmless
    # and would otherwise flood the console with 22 warnings × 139 dates.
    result <- tryCatch(
      {
        suppressWarnings(
          generate_single_forecast(
            origin_date  = this_date,
            wili_tsibble = wili_tsibble,
            lambdas      = lambdas,
            n_bootstrap  = N_BOOTSTRAP,
            team_model   = TEAM_MODEL,
            output_dir   = OUTPUT_DIR,
            write_file   = TRUE
          )
        )
      },
      error = function(e) {
        # Log the error but don't stop the loop
        msg <- paste0(this_date, ": ", e$message)
        failed_log <<- c(failed_log, msg)
        cat("    ✗ ERROR:", e$message, "\n")
        NULL
      }
    )
    
    # --- Update counters ---
    if (!is.null(result)) {
      n_success <- n_success + 1
    } else {
      # Only count as failed if we didn't already log it in tryCatch
      # (generate_single_forecast returns NULL for data issues too)
      if (length(failed_log) == 0 ||
          !grepl(as.character(this_date), tail(failed_log, 1))) {
        n_failed <- n_failed + 1
        failed_log <- c(failed_log, paste0(this_date, ": returned NULL"))
      } else {
        n_failed <- n_failed + 1
      }
    }
  }
  
  # --- Summary report ---
  loop_end     <- Sys.time()
  total_elapsed <- as.numeric(difftime(loop_end, loop_start, units = "secs"))
  
  cat("\n", strrep("=", 60), "\n")
  cat("  FORECAST LOOP COMPLETE\n")
  cat(strrep("=", 60), "\n\n")
  cat(sprintf("  Total dates:     %d\n", n_total))
  cat(sprintf("  Successful:      %d\n", n_success))
  cat(sprintf("  Skipped:         %d (already existed)\n", n_skipped))
  cat(sprintf("  Failed:          %d\n", n_failed))
  cat(sprintf("  Total time:      %.1f minutes (%.0f seconds)\n",
              total_elapsed / 60, total_elapsed))
  
  if (n_success > 0) {
    cat(sprintf("  Avg time/date:   %.1f seconds\n",
                total_elapsed / (n_success + n_failed)))
  }
  
  # --- Report failures ---
  if (n_failed > 0) {
    cat("\n  FAILED DATES:\n")
    for (msg in failed_log) {
      cat("    ✗", msg, "\n")
    }
  }
  
  # --- File count check ---
  # Verify the number of CSVs in the output directory matches expectations.
  csv_files <- list.files(OUTPUT_DIR, pattern = "\\.csv$")
  cat(sprintf("\n  Files in output directory: %d\n", length(csv_files)))
  cat(sprintf("  Expected (all dates):     %d\n", length(origin_dates)))
  
  if (length(csv_files) == length(origin_dates)) {
    cat("  ✓ File count matches!\n")
  } else {
    cat(sprintf("  ⚠ Mismatch: %d files vs %d expected dates\n",
                length(csv_files), length(origin_dates)))
  }
  
  # =========================================================================
  # POST-LOOP SPOT CHECKS
  #
  # Validate a sample of 5 CSVs spread across the date range:
  # first, last, and 3 evenly-spaced middle dates.
  # This catches systematic issues that might only appear in certain seasons.
  # =========================================================================
  
  cat("\n  Running spot-check validation on 5 sample files...\n\n")
  
  # Pick 5 representative dates: first, 25th percentile, median, 75th, last
  sample_indices <- unique(round(
    quantile(seq_along(origin_dates), probs = c(0, 0.25, 0.5, 0.75, 1))
  ))
  sample_dates <- origin_dates[sample_indices]
  
  for (j in seq_along(sample_dates)) {
    spot_date <- sample_dates[j]
    spot_file <- paste0(spot_date, "-", TEAM_MODEL, ".csv")
    spot_path <- file.path(OUTPUT_DIR, spot_file)
    
    if (file.exists(spot_path)) {
      # Read it back and run our internal validation
      spot_df <- readr::read_csv(spot_path, show_col_types = FALSE)
      cat(sprintf("  Spot-checking %s...\n", spot_file))
      
      tryCatch(
        {
          validate_forecast_df(spot_df, spot_date, QUANTILE_LEVELS)
        },
        error = function(e) {
          cat("    ✗ VALIDATION FAILED:", e$message, "\n")
        }
      )
    } else {
      cat(sprintf("  Spot-check: %s not found (skipped)\n", spot_file))
    }
  }
  
  cat("\n  Spot checks complete.\n")
  cat("  To run full hubValidations check on a specific file:\n")
  cat('  validate_submission(hub_path = HUB_PATH,\n')
  cat('    file_path = "KReger-snaive_bc_bs/YYYY-MM-DD-KReger-snaive_bc_bs.csv")\n')
  
} else {
  cat("\n=== STEP 5: Forecast loop SKIPPED (RUN_FORECAST_LOOP = FALSE) ===\n")
  cat("  Set RUN_FORECAST_LOOP <- TRUE at the top of the script to run.\n")
  cat("  Or call generate_single_forecast() manually for individual dates.\n")
}


###############################################################################
# MILESTONE 3 — TEST SECTION
#
# Run these commands INTERACTIVELY after source()-ing the script to verify
# the single-date pipeline works end-to-end.
#
# Test date: 2016-12-03 — this matches the origin date used in the
# w5d1-solutions.qmd exercises, so it's a known-good time period with
# plenty of training data and we're in the middle of flu season (good
# for seeing non-trivial forecasts).
###############################################################################

# --- TEST: Generate a single forecast ---
# Uncomment the block below and run it in your console.
# Use n_bootstrap = 100 for a fast test (~10 seconds), not the full 1000.
#
# cat("\n=== MILESTONE 3 TEST: Single-date forecast ===\n")
# test_result <- generate_single_forecast(
#   origin_date   = as.Date("2016-12-03"),
#   wili_tsibble  = wili_tsibble,
#   lambdas       = lambdas,
#   n_bootstrap   = 100,        # Fast test; use 1000 for production
#   write_file    = TRUE         # Writes the CSV so we can validate it
# )
#
# # --- Check dimensions ---
# cat("\nResult dimensions:", nrow(test_result), "rows x",
#     ncol(test_result), "cols\n")
# cat("Expected: 1012 rows x 8 cols\n\n")
#
# # --- Inspect the data ---
# cat("First rows:\n")
# print(head(test_result, 6))
# cat("\nLast rows:\n")
# print(tail(test_result, 6))
#
# # --- Hub validation using hubValidations ---
# # This is the OFFICIAL validator that GitHub Actions will run on your PR.
# # It checks everything: filename format, column names, value ranges, etc.
# library(hubValidations)
# hub_validation_result <- validate_submission(
#   hub_path  = HUB_PATH,
#   file_path = "KReger-snaive_bc_bs/2016-12-03-KReger-snaive_bc_bs.csv"
# )
# print(hub_validation_result)
#
# # --- Spot-check plot ---
# # This fan chart shows the forecast overlaid on the historical series.
# # For SNAIVE, the median (blue line) should be very close to the
# # same-week-last-year values. The fan should widen slightly with horizon.
# p <- plot_quantile_fan(test_result, wili_tsibble,
#                        as.Date("2016-12-03"), "US National")
# print(p)
#
# # --- Clean up the test file ---
# # Remove it so it doesn't get included in the final submission by accident.
# # (Or leave it if you want to inspect it.)
# # file.remove(file.path(OUTPUT_DIR, "2016-12-03-KReger-snaive_bc_bs.csv"))