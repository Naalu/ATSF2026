###############################################################################
# FILE: scripts/KReger-snaive_bc_bs_forecast.R
#
# Model-specific forecasting pipeline for the KReger-snaive_bc_bs submission
# to the FluSight ILI Sandbox Hub (ATSF2026).
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
# ARCHITECTURE:
#   This script sources forecast_helpers.R for all model-agnostic utilities.
#   The only model-specific logic is here:
#     1. Lambda estimation via guerrero()
#     2. generate_single_forecast() — the fit/generate/quantile pipeline
#     3. The loop over all 139 origin dates
#
#   To implement a different model, copy this file, change the model() call
#   inside generate_single_forecast(), and update naming constants. Everything
#   else (formatting, validation, writing) is inherited from forecast_helpers.R.
#
# HOW TO RUN:
#   1. Set working directory to the ATSF2026 repo root
#   2. source("scripts/KReger-snaive_bc_bs_forecast.R")
#
# INPUTS:
#   - target-data/time-series.csv (wILI time series)
#   - hub-config/tasks.json (origin dates)
#
# OUTPUTS:
#   - 139 CSV files in model-output/KReger-snaive_bc_bs/
#
###############################################################################


# =============================================================================
# LOAD PACKAGES
# =============================================================================

# fpp3 brings in fable, fabletools, tsibble, feasts
# tidyverse brings in dplyr, readr, ggplot2, purrr, etc.
library(fpp3)
library(tidyverse)


# =============================================================================
# CONFIGURATION
#
# All configurable values in one place. Shared constants (quantile levels,
# locations, etc.) live in forecast_helpers.R.
# =============================================================================

# Team and model identifiers (must match directory/filename convention)
TEAM_MODEL <- "KReger-snaive_bc_bs"

# Path to the local ATSF2026 repo clone. All other paths derive from this.
HUB_PATH <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"

# Derived paths
DATA_PATH    <- file.path(HUB_PATH, "target-data", "time-series.csv")
TASKS_PATH   <- file.path(HUB_PATH, "hub-config", "tasks.json")
OUTPUT_DIR   <- file.path(HUB_PATH, "model-output", TEAM_MODEL)

# Bootstrap replicate count.
# 1000 gives ≥10 observations in each tail for the 0.01/0.99 quantiles.
# Use 100 for fast dev/testing (~10x speedup).
N_BOOTSTRAP <- 1000

# When TRUE, only processes the first 5 origin dates (quick testing).
DEV_MODE <- FALSE

# When TRUE, the forecast loop runs on source(). When FALSE, only setup
# (Steps 1-3) runs, and you can call generate_single_forecast() manually.
RUN_FORECAST_LOOP <- FALSE


# =============================================================================
# SOURCE SHARED HELPERS
# =============================================================================

source(file.path(HUB_PATH, "scripts", "forecast_helpers.R"))


# =============================================================================
# STEP 1: LOAD AND PREPARE DATA
# =============================================================================

cat("=== STEP 1: Loading target data ===\n")
wili_tsibble <- load_and_prepare_tsibble(DATA_PATH)

cat(paste0("  Rows: ", nrow(wili_tsibble), "\n"))
cat(paste0("  Date range: ",
           min(wili_tsibble$target_end_date), " to ",
           max(wili_tsibble$target_end_date), "\n"))
cat(paste0("  Locations: ",
           length(unique(wili_tsibble$location)), "\n\n"))


# =============================================================================
# STEP 2: LOAD ORIGIN DATES
# =============================================================================

cat("=== STEP 2: Loading origin dates ===\n")
origin_dates <- get_origin_dates(TASKS_PATH)
cat("\n")


# =============================================================================
# STEP 3: ESTIMATE BOX-COX LAMBDA (per location)
#
# WHY BOX-COX?
#   ILI percentages have strongly heteroskedastic variance: winter peaks
#   spike to 5-10% with high variability, while summer hovers near 1% with
#   low variability. Box-Cox stabilizes this so that bootstrap residuals
#   produce properly calibrated intervals across all seasons.
#
# WHY GUERRERO?
#   The guerrero() method (Guerrero, 1993) automatically selects the lambda
#   that minimizes the coefficient of variation of subseries. It's the
#   standard approach in fpp3.
#
# WHY ESTIMATE ONCE (not per origin date)?
#   Efficiency (avoids 139 × 11 = 1,529 redundant computations) and
#   stability (fixed lambda is reproducible regardless of origin date).
#   We use data strictly before the first origin date.
# =============================================================================

cat("=== STEP 3: Estimating Box-Cox lambdas via guerrero ===\n")

# Use data before the first origin date for lambda estimation.
# Strict < excludes the first origin date itself (on that date we forecast).
first_origin <- min(origin_dates)
cat(paste0("  Using data before first origin date: ", first_origin, "\n"))

pre_forecast_data <- wili_tsibble |>
  filter(target_end_date < first_origin)

cat(paste0("  Pre-forecast rows: ", nrow(pre_forecast_data), "\n"))
cat(paste0("  Pre-forecast date range: ",
           min(pre_forecast_data$target_end_date), " to ",
           max(pre_forecast_data$target_end_date), "\n"))

# features() computes guerrero lambda for each location's time series.
# Result: tibble with columns location, target, lambda_guerrero.
lambda_table <- pre_forecast_data |>
  features(observation, features = guerrero)

cat("\n  Guerrero lambda estimates by location:\n")
print(lambda_table, n = Inf)

# Named vector for easy lookup: lambdas["US National"] returns a scalar
lambdas <- lambda_table$lambda_guerrero
names(lambdas) <- lambda_table$location

# Sanity check: correct count
if (length(lambdas) != length(HUB_LOCATIONS)) {
  warning(
    "Expected ", length(HUB_LOCATIONS), " lambdas but got ", length(lambdas),
    ". Missing: ", paste(setdiff(HUB_LOCATIONS, names(lambdas)), collapse = ", ")
  )
}

# Flag extreme lambdas (outside [-2, 2] suggests data issues)
extreme_lambdas <- lambdas[abs(lambdas) > 2]
if (length(extreme_lambdas) > 0) {
  warning(
    "Some lambdas outside [-2, 2] — consider log1p() fallback:\n",
    paste0("  ", names(extreme_lambdas), ": ", round(extreme_lambdas, 4),
           collapse = "\n")
  )
} else {
  cat("\n  All lambdas in reasonable range [-2, 2].\n")
}

# Note lambdas near zero (≈ log transform, common for right-skewed ILI data)
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
# STEP 4: SINGLE-DATE FORECAST FUNCTION
# =============================================================================

#' Generate a complete hub-formatted forecast for a single origin date
#'
#' Pipeline: filter training data → fit SNAIVE with Box-Cox → generate
#' bootstrap paths → compute quantiles → format → validate → write CSV.
#'
#' @param origin_date A single Date — the Saturday to forecast from.
#' @param wili_tsibble The full wILI tsibble (all dates, all locations).
#' @param lambdas Named numeric vector: location names → guerrero lambdas.
#' @param quantile_levels Numeric vector of quantile probabilities (default: 23 hub levels).
#' @param n_bootstrap Number of bootstrap replicates (100 for testing, 1000 for production).
#' @param team_model Character string (e.g., "KReger-snaive_bc_bs").
#' @param output_dir Path to write the CSV file.
#' @param write_file If TRUE, writes CSV. If FALSE, returns dataframe only.
#' @return A 1,012-row × 8-column tibble in hub format, or NULL on failure.
generate_single_forecast <- function(origin_date,
                                     wili_tsibble,
                                     lambdas,
                                     quantile_levels = QUANTILE_LEVELS,
                                     n_bootstrap = N_BOOTSTRAP,
                                     team_model = TEAM_MODEL,
                                     output_dir = OUTPUT_DIR,
                                     write_file = TRUE) {
  
  # --- Filter training data ---
  # Use <= because the origin_date Saturday is the last day of the reporting
  # epiweek, so data for that week is available when we make the forecast.
  train_data <- wili_tsibble |>
    filter(target_end_date <= origin_date)
  
  # SNAIVE with lag('year') needs ≥52 weeks of history
  min_weeks <- train_data |>
    as_tibble() |>
    group_by(location) |>
    summarise(n = n(), .groups = "drop") |>
    pull(n) |>
    min()
  
  if (min_weeks < 52) {
    warning(
      "Skipping ", origin_date, ": only ", min_weeks,
      " weeks of data (need >= 52 for SNAIVE yearly lag)"
    )
    return(NULL)
  }
  
  # --- Fit model and generate bootstrap paths (per location) ---
  # We loop over locations because fable's box_cox() requires a scalar lambda.
  # Passing lambda as a column fails: generate() creates future rows that
  # don't carry the lambda column, so back-transformation silently breaks.
  # With a scalar, fable stores lambda in the model object and back-transforms
  # correctly during generate().
  locations <- unique(train_data$location)
  
  sims <- tryCatch(
    {
      purrr::map(locations, function(loc) {
        loc_lambda <- lambdas[loc]
        loc_train  <- train_data |> filter(location == loc)
        
        # box_cox() stabilizes variance; SNAIVE uses same-week-last-year
        loc_fit <- loc_train |>
          model(
            snaive_bc_bs = SNAIVE(
              box_cox(observation, loc_lambda) ~ lag("year")
            )
          )
        
        # Bootstrap: resample residuals to create simulated future paths
        loc_sims <- loc_fit |>
          generate(h = 4, times = n_bootstrap, bootstrap = TRUE)
        
        as_tibble(loc_sims)
        
      }) |>
        bind_rows()
    },
    error = function(e) {
      warning("Fit/generate failed for ", origin_date, ": ", e$message)
      NULL
    }
  )
  
  if (is.null(sims)) return(NULL)
  
  # Warn about NAs (can occur from time series gaps or numerical issues)
  na_count <- sum(is.na(sims$.sim))
  if (na_count > 0) {
    warning(origin_date, ": ", na_count, " NA values in simulations — ",
            "excluded from quantile computation.")
  }
  
  # --- Compute quantiles, format, validate, write ---
  quantile_df <- compute_quantiles_from_sims(sims, quantile_levels)
  hub_df      <- format_for_hub(quantile_df, origin_date, team_model)
  validate_forecast_df(hub_df, origin_date, quantile_levels)
  
  if (write_file) {
    filepath <- write_hub_csv(hub_df, origin_date, output_dir, team_model)
    cat(paste0("  Wrote: ", basename(filepath), "\n"))
  }
  
  return(hub_df)
}


# =============================================================================
# STEP 5: MAIN FORECAST LOOP
#
# Generates forecasts for all origin dates (or first 5 in DEV_MODE).
# Features: progress/ETA, skip-if-exists, per-date error handling, summary.
# =============================================================================

if (RUN_FORECAST_LOOP) {
  
  cat("\n=== STEP 5: Running forecast loop ===\n")
  
  # Determine which dates to process
  if (DEV_MODE) {
    dates_to_process <- origin_dates[1:min(5, length(origin_dates))]
    cat("  DEV_MODE: processing first", length(dates_to_process), "dates\n")
  } else {
    dates_to_process <- origin_dates
    cat("  PRODUCTION: processing all", length(dates_to_process), "dates\n")
  }
  cat("  N_BOOTSTRAP =", N_BOOTSTRAP, "\n\n")
  
  if (!dir.exists(OUTPUT_DIR)) {
    dir.create(OUTPUT_DIR, recursive = TRUE)
    cat("  Created output directory:", OUTPUT_DIR, "\n")
  }
  
  # Tracking variables
  n_total      <- length(dates_to_process)
  n_success    <- 0
  n_skipped    <- 0
  n_failed     <- 0
  failed_log   <- character(0)
  loop_start   <- Sys.time()
  process_time <- 0  # cumulative time for non-skipped dates (for accurate ETA)
  
  for (i in seq_along(dates_to_process)) {
    
    this_date <- dates_to_process[i]
    
    # Skip if CSV already exists (allows resuming interrupted runs)
    expected_filepath <- file.path(
      OUTPUT_DIR, paste0(this_date, "-", TEAM_MODEL, ".csv")
    )
    
    if (file.exists(expected_filepath)) {
      cat(sprintf("  [%3d/%d] %s — SKIPPED (exists)\n", i, n_total, this_date))
      n_skipped <- n_skipped + 1
      next
    }
    
    # Progress with ETA (based on actual processing time, not wall clock)
    n_processed <- n_success + n_failed
    if (n_processed > 0) {
      avg_secs    <- process_time / n_processed
      remaining   <- (n_total - n_skipped - n_processed - 1) * avg_secs
      eta_str     <- sprintf(" | ETA: %.0f min", remaining / 60)
    } else {
      eta_str <- ""
    }
    cat(sprintf("  [%3d/%d] %s%s\n", i, n_total, this_date, eta_str))
    
    # Generate forecast with error handling.
    # suppressWarnings silences routine fable Box-Cox back-transform messages
    # that would otherwise flood the console (~22 warnings × 139 dates).
    date_start <- Sys.time()
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
        cat("    ✗ ERROR:", e$message, "\n")
        NULL
      }
    )
    process_time <- process_time + as.numeric(difftime(Sys.time(), date_start, units = "secs"))
    
    # Update counters. NULL means failure (either tryCatch error or
    # generate_single_forecast returning NULL for insufficient data).
    if (!is.null(result)) {
      n_success <- n_success + 1
    } else {
      n_failed <- n_failed + 1
      failed_log <- c(failed_log, as.character(this_date))
    }
  }
  
  # --- Summary report ---
  total_elapsed <- as.numeric(difftime(Sys.time(), loop_start, units = "secs"))
  
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
    cat(sprintf("  Avg time/date:   %.1f seconds\n", process_time / (n_success + n_failed)))
  }
  
  if (n_failed > 0) {
    cat("\n  FAILED DATES:\n")
    for (msg in failed_log) cat("    ✗", msg, "\n")
  }
  
  # File count check
  csv_files <- list.files(OUTPUT_DIR, pattern = "\\.csv$")
  cat(sprintf("\n  Files in output directory: %d\n", length(csv_files)))
  cat(sprintf("  Expected (all dates):     %d\n", length(origin_dates)))
  
  if (length(csv_files) == length(origin_dates)) {
    cat("  ✓ File count matches!\n")
  } else {
    cat(sprintf("  ⚠ Mismatch: %d files vs %d expected dates\n",
                length(csv_files), length(origin_dates)))
  }
  
  # --- Post-loop spot checks ---
  # Validate 5 CSVs spread across the date range to catch seasonal issues.
  cat("\n  Running spot-check validation on 5 sample files...\n\n")
  
  sample_indices <- unique(round(
    quantile(seq_along(origin_dates), probs = c(0, 0.25, 0.5, 0.75, 1))
  ))
  sample_dates <- origin_dates[sample_indices]
  
  for (j in seq_along(sample_dates)) {
    spot_date <- sample_dates[j]
    spot_file <- paste0(spot_date, "-", TEAM_MODEL, ".csv")
    spot_path <- file.path(OUTPUT_DIR, spot_file)
    
    if (file.exists(spot_path)) {
      spot_df <- readr::read_csv(spot_path, show_col_types = FALSE)
      cat(sprintf("  Spot-checking %s...\n", spot_file))
      tryCatch(
        validate_forecast_df(spot_df, spot_date, QUANTILE_LEVELS),
        error = function(e) cat("    ✗ VALIDATION FAILED:", e$message, "\n")
      )
    } else {
      cat(sprintf("  Spot-check: %s not found (skipped)\n", spot_file))
    }
  }
  
  cat("\n  Spot checks complete.\n")
  cat("  For full hub validation on a specific file:\n")
  cat('  hubValidations::validate_submission(hub_path = HUB_PATH,\n')
  cat('    file_path = "KReger-snaive_bc_bs/YYYY-MM-DD-KReger-snaive_bc_bs.csv")\n')
  
} else {
  cat("\n=== STEP 5: Forecast loop SKIPPED (RUN_FORECAST_LOOP = FALSE) ===\n")
  cat("  Set RUN_FORECAST_LOOP <- TRUE to run, or call generate_single_forecast() manually.\n")
}


###############################################################################
# INTERACTIVE TESTING
#
# Uncomment and run after source()-ing the script to verify the single-date
# pipeline end-to-end. Test date 2016-12-03 is mid-flu-season with plenty
# of training data.
###############################################################################

# test_result <- generate_single_forecast(
#   origin_date   = as.Date("2016-12-03"),
#   wili_tsibble  = wili_tsibble,
#   lambdas       = lambdas,
#   n_bootstrap   = 100,
#   write_file    = TRUE
# )
#
# cat("Result:", nrow(test_result), "rows x", ncol(test_result), "cols\n")
# print(head(test_result, 6))
#
# # Official hub validator (what GitHub Actions runs on your PR)
# library(hubValidations)
# validate_submission(
#   hub_path  = HUB_PATH,
#   file_path = "KReger-snaive_bc_bs/2016-12-03-KReger-snaive_bc_bs.csv"
# )
#
# # Fan chart for visual sanity check
# p <- plot_quantile_fan(test_result, wili_tsibble,
#                        as.Date("2016-12-03"), "US National")
# print(p)