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
TEAM_MODEL <- "KReger-snaive_bc_bs"

# --- Path to the local ATSF2026 repo clone ---
# This is the root of the hub directory containing hub-config/, model-output/,
# target-data/, etc. All other paths are constructed relative to this.
HUB_PATH <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"

# --- Derived paths (constructed from HUB_PATH and TEAM_MODEL) ---
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
N_BOOTSTRAP <- 1000   # <-- Change to 100 for faster test runs


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
# =============================================================================

source(file.path(HUB_PATH, "scripts", "forecast_helpers.R"))


# =============================================================================
# STEP 1: LOAD AND PREPARE DATA
#
# Read the wILI target time series and convert to a tsibble.
# This is the ground truth data that models are trained on.
#
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
#     processed.
#
# HOW IT WORKS:
#   1. Filter the tsibble to keep only data BEFORE the first origin date.
#      This is the "pre-forecast" historical data.
#   2. For each location, compute features(observation, guerrero) which
#      returns a tibble with one row per location and a lambda_guerrero column.
#   3. Store as a named vector: names = location names, values = lambda.
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
# that location.
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
# (MILESTONE 3 will go here: generate_single_forecast() function)
# =============================================================================

# =============================================================================
# (MILESTONE 5 will go here: loop over all origin dates)
# =============================================================================


###############################################################################
# VERIFICATION SECTION
#
# Uncomment and run these blocks interactively to verify the setup is correct.
# These are NOT part of the production pipeline — they're diagnostic tools.
###############################################################################

# --- VERIFY 1: Print lambda summary ---
# Uncomment these lines to inspect lambdas after running the script:
#
# cat("\nLambda summary:\n")
# print(tibble(location = names(lambdas), lambda = round(lambdas, 4)))

# --- VERIFY 2: Visual check — raw vs Box-Cox-transformed series ---
# This plot compares the raw observation series with the Box-Cox-transformed
# version for one location. After transformation, the seasonal peaks and
# troughs should have more uniform amplitude (stabilized variance).
#
# Uncomment and run:
#
# library(cowplot)
#
# # Pick a location to inspect (US National shows the clearest seasonal pattern)
# check_loc <- "US National"
# check_lambda <- lambdas[check_loc]
#
# # Raw series: note how winter peaks vary in height across years
# p_raw <- pre_forecast_data |>
#   filter(location == check_loc) |>
#   autoplot(observation) +
#   labs(title = paste0(check_loc, " — Raw"),
#        subtitle = "Note: peak heights vary across seasons",
#        y = "wILI %") +
#   theme_minimal()
#
# # Box-Cox transformed: peaks should be more uniform
# p_bc <- pre_forecast_data |>
#   filter(location == check_loc) |>
#   autoplot(box_cox(observation, check_lambda)) +
#   labs(title = paste0(check_loc, " — Box-Cox (λ=", round(check_lambda, 3), ")"),
#        subtitle = "Peak heights should be more uniform",
#        y = "Transformed wILI") +
#   theme_minimal()
#
# # Display side by side
# plot_grid(p_raw, p_bc, ncol = 1)

# --- VERIFY 3: Confirm fable can use box_cox with our lambda ---
# This is a quick smoke test: fit SNAIVE with Box-Cox on one location
# for a single origin date, just to confirm the model() call works.
#
# Uncomment and run:
#
# test_data <- wili_tsibble |>
#   filter(location == "US National", target_end_date < as.Date("2016-12-03"))
# test_lambda <- lambdas["US National"]
#
# test_fit <- test_data |>
#   model(snaive_bc = SNAIVE(box_cox(observation, test_lambda) ~ lag("year")))
#
# cat("Model fit successful. Summary:\n")
# print(test_fit)
#
# # Generate a small number of bootstrap paths to verify generate() works
# test_sims <- test_fit |>
#   generate(h = 4, times = 10, bootstrap = TRUE)
#
# cat("\nGenerate output dimensions: ", dim(test_sims), "\n")
# cat("Columns: ", paste(names(test_sims), collapse = ", "), "\n")
# print(head(test_sims))