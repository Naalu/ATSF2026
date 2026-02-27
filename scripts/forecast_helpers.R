###############################################################################
# FILE: scripts/forecast_helpers.R
#
# PURPOSE:
#   Shared, MODEL-AGNOSTIC utility functions for the FluSight ILI Sandbox Hub
#   forecast submission project. This file is designed to be reused across
#   multiple models (SNAIVE, ETS, ARIMA, etc.) — the model-specific logic
#   lives in separate driver scripts that source() this file.
#
# DESIGN DOC REFERENCE:
#   flucast-design-document.md, Section 7.4 (Extensibility Architecture)
#
# FUNCTIONS PROVIDED:
#   1. load_and_prepare_tsibble()    — Load CSV, create tsibble
#   2. get_origin_dates()            — Parse origin dates from tasks.json
#   3. compute_quantiles_from_sims() — Empirical quantiles from bootstrap sims
#   4. format_for_hub()              — Shape data into the 8-column hub format
#   5. validate_forecast_df()        — Run ALL Section 9.1 quality checks
#   6. write_hub_csv()               — Write CSV with correct filename pattern
#   7. plot_quantile_fan()           — Spot-check visualization
#
# HOW TO USE:
#   In your model-specific script (e.g., KReger-snaive_bc_bs_forecast.R):
#     source("scripts/forecast_helpers.R")
#   Then call these functions as needed. See each function's roxygen block
#   for parameters, return values, and examples.
#
# DEPENDENCIES:
#   tidyverse, tsibble, fable, fabletools, jsonlite, lubridate, ggplot2
#   (All loaded via library(fpp3); library(tidyverse) in your driver script)
#
###############################################################################


# =============================================================================
# CONSTANTS — Hub-specific values used across multiple functions
#
# These are defined here as single-source-of-truth so that if the hub spec
# ever changes, you only need to update ONE place.
# =============================================================================

# ---------------------------------------------------------------------------
# DATE COLUMN NAME
#
# RESOLVED from Open Question #1: The validated CSVs in this hub use
# "origin_date", NOT "reference_date" (which is what the README template says).
# The hubValidations validator checks for "origin_date", and the instructor's
# get_step_ahead_model_output.R writes "origin_date".
#
# If a future hub uses a different name, change ONLY this constant.
# Reference: flucast-design-document.md Section 11, OQ#1
# ---------------------------------------------------------------------------
DATE_COL <- "origin_date"

# ---------------------------------------------------------------------------
# 23 REQUIRED QUANTILE LEVELS
#
# These are the exact quantile probabilities the hub expects in every CSV.
# Defined in hub-config/tasks.json under output_type.quantile.output_type_id.
# The seq() in the middle generates 0.05, 0.10, ..., 0.95 (19 values),
# and we prepend 0.01, 0.025 and append 0.975, 0.99 for the tails.
# Total: 2 + 19 + 2 = 23 levels.
#
# Reference: flucast-design-document.md Section 4.2
# ---------------------------------------------------------------------------
QUANTILE_LEVELS <- round(c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99), 3)

# ---------------------------------------------------------------------------
# 11 REQUIRED LOCATIONS
#
# The hub covers the US national level plus 10 HHS Regions.
# These are the string values that must appear in the "location" column.
#
# Reference: hub-config/tasks.json — location.optional
# ---------------------------------------------------------------------------
HUB_LOCATIONS <- c(
  "US National",
  paste("HHS Region", 1:10)
)

# ---------------------------------------------------------------------------
# HORIZONS: 1 through 4 weeks ahead
# Reference: hub-config/tasks.json — horizon.optional
# ---------------------------------------------------------------------------
HUB_HORIZONS <- 1:4

# ---------------------------------------------------------------------------
# EXPECTED ROW COUNT per CSV file
# 11 locations × 4 horizons × 23 quantile levels = 1,012 rows
# Reference: flucast-design-document.md Section 4.3
# ---------------------------------------------------------------------------
EXPECTED_ROWS <- length(HUB_LOCATIONS) * length(HUB_HORIZONS) * length(QUANTILE_LEVELS)

# ---------------------------------------------------------------------------
# EXPECTED COLUMN NAMES (in the order the hub examples use)
# Reference: model-output/delphi-epicast/2015-10-31-delphi-epicast.csv
# ---------------------------------------------------------------------------
EXPECTED_COLS <- c(
  DATE_COL, "location", "target", "horizon",
  "target_end_date", "output_type", "output_type_id", "value"
)


###############################################################################
# FUNCTION 1: load_and_prepare_tsibble
###############################################################################

#' Load target data CSV and convert to a tsibble
#'
#' Reads the wILI percentage time series from a CSV file (local path or URL)
#' and converts it into a tsibble suitable for use with fable models. The
#' tsibble is indexed by target_end_date (weekly Saturdays) and keyed by
#' location and target.
#'
#' This follows the exact pattern from w5d1-solutions.qmd Exercises 1–2.
#'
#' @param data_path Character string. Path to time-series.csv (local file
#'   or a URL like the GitHub raw URL).
#'
#' @return A tsibble with columns: location, target_end_date, target,
#'   observation. Index = target_end_date [1W], Key = location, target.
#'
#' @examples
#' wili <- load_and_prepare_tsibble("target-data/time-series.csv")
#' wili <- load_and_prepare_tsibble(
#'   "https://raw.githubusercontent.com/Naalu/ATSF2026/refs/heads/main/target-data/time-series.csv"
#' )
load_and_prepare_tsibble <- function(data_path) {
  
  # --- Input validation ---
  # Make sure we got a non-empty string
  stopifnot(
    "data_path must be a single character string" =
      is.character(data_path) && length(data_path) == 1
  )
  
  # --- Read the CSV ---
  # show_col_types = FALSE suppresses the readr column specification message
  # that would otherwise clutter the console during batch runs
  wili_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  
  # --- Verify expected columns exist ---
  required_cols <- c("location", "target_end_date", "target", "observation")
  missing_cols <- setdiff(required_cols, names(wili_raw))
  if (length(missing_cols) > 0) {
    stop(
      "Target data CSV is missing expected columns: ",
      paste(missing_cols, collapse = ", "),
      "\nFound columns: ", paste(names(wili_raw), collapse = ", ")
    )
  }
  
  # --- Convert to tsibble ---
  # index = target_end_date: this is the time axis (weekly Saturdays)
  # key = c(location, target): each combination of location × target
  #   is a separate time series. For this hub there's only one target
  #   ("ili perc"), but including it in the key makes the tsibble
  #   compatible with hubs that have multiple targets.
  #
  # Pattern from: w5d1-solutions.qmd Exercise 1
  wili_tsibble <- wili_raw |>
    tsibble::as_tsibble(index = target_end_date, key = c(location, target))
  
  # --- Verify location count ---
  n_locs <- length(unique(wili_tsibble$location))
  if (n_locs != length(HUB_LOCATIONS)) {
    warning(
      "Expected ", length(HUB_LOCATIONS), " locations but found ", n_locs,
      ". Locations: ", paste(unique(wili_tsibble$location), collapse = ", ")
    )
  }
  
  return(wili_tsibble)
}


###############################################################################
# FUNCTION 2: get_origin_dates
###############################################################################

#' Parse origin dates from the hub's tasks.json configuration
#'
#' Reads hub-config/tasks.json and extracts the vector of all origin dates
#' that forecasts should be submitted for. Returns them as sorted Date objects.
#'
#' Uses simplifyVector = FALSE to handle the nested JSON structure correctly
#' (without this, jsonlite flattens nested arrays into atomic vectors which
#' breaks $ navigation).
#'
#' @param tasks_json_path Character string. Path to hub-config/tasks.json.
#'
#' @return A sorted Date vector of all origin dates (all Saturdays).
#'   For this hub: 139 dates from 2015-10-17 to 2020-02-29.
#'
#' @examples
#' dates <- get_origin_dates("hub-config/tasks.json")
#' length(dates)  # 139
get_origin_dates <- function(tasks_json_path) {
  
  # --- Input validation ---
  stopifnot(
    "tasks_json_path must be a single character string" =
      is.character(tasks_json_path) && length(tasks_json_path) == 1
  )
  
  # --- Parse JSON ---
  # simplifyVector = FALSE is CRITICAL here. Without it, jsonlite collapses
  # deeply nested lists into atomic vectors, causing "$ operator is invalid
  # for atomic vectors" errors when navigating the structure.
  tasks <- jsonlite::fromJSON(tasks_json_path, simplifyVector = FALSE)
  
  # --- Navigate to the origin dates ---
  # Structure: tasks -> rounds[[1]] -> model_tasks[[1]] -> task_ids ->
  #            origin_date -> optional (list of date strings)
  task_ids <- tasks$rounds[[1]]$model_tasks[[1]]$task_ids
  
  # Try "origin_date" first, then "reference_date" as fallback
  # (different hubs may use different names for this field)
  raw_dates <- unlist(task_ids$origin_date$optional)
  if (is.null(raw_dates)) {
    raw_dates <- unlist(task_ids$reference_date$optional)
  }
  if (is.null(raw_dates)) {
    stop(
      "Could not find origin dates in tasks.json. ",
      "Available task_id fields: ", paste(names(task_ids), collapse = ", ")
    )
  }
  
  # --- Convert to Date and sort ---
  origin_dates <- sort(as.Date(raw_dates))
  
  # --- Sanity checks ---
  # All origin dates should be Saturdays (wday: 1=Sun, 7=Sat)
  if (!all(lubridate::wday(origin_dates) == 7)) {
    warning("Not all origin dates are Saturdays — check tasks.json")
  }
  
  cat(paste0("  Loaded ", length(origin_dates), " origin dates (",
             min(origin_dates), " to ", max(origin_dates), ")\n"))
  
  return(origin_dates)
}


###############################################################################
# FUNCTION 3: compute_quantiles_from_sims
###############################################################################

#' Compute empirical quantiles from bootstrap simulation paths
#'
#' Takes the output of fable's generate(bootstrap = TRUE) — a tsibble with
#' a .sim column containing simulated future values — and computes the
#' requested quantile levels for each location × target_end_date group.
#' Floors all values at zero to satisfy the hub's non-negativity constraint.
#'
#' This implements Steps 5d–5e from flucast-design-document.md Section 6.
#'
#' @param sims_tsibble A tsibble produced by fable::generate(), containing
#'   at minimum columns: location, target_end_date, .sim
#' @param quantile_levels Numeric vector of quantile probabilities to compute.
#'   Defaults to the hub's 23 required levels (QUANTILE_LEVELS constant).
#'
#' @return A tibble with columns: location, target_end_date, output_type_id
#'   (the quantile probability), and value (the quantile estimate, >= 0).
#'   Rows: n_locations × n_horizons × length(quantile_levels).
#'
#' @examples
#' # After: sims <- model_fit |> generate(h = 4, times = 1000, bootstrap = TRUE)
#' # quantile_df <- compute_quantiles_from_sims(sims)
compute_quantiles_from_sims <- function(sims_tsibble,
                                        quantile_levels = QUANTILE_LEVELS) {
  
  # --- Input validation ---
  if (!".sim" %in% names(sims_tsibble)) {
    stop("sims_tsibble must contain a '.sim' column (output of generate())")
  }
  if (!"location" %in% names(sims_tsibble)) {
    stop("sims_tsibble must contain a 'location' column")
  }
  if (!"target_end_date" %in% names(sims_tsibble)) {
    stop("sims_tsibble must contain a 'target_end_date' column")
  }
  
  # --- Compute quantiles ---
  # Group by location and target_end_date (each group has `times` simulated
  # values). For each group, compute quantile() at each requested level.
  #
  # The reframe() function is used instead of summarise() because we're
  # returning multiple rows per group (one row per quantile level).
  quantile_df <- sims_tsibble |>
    tibble::as_tibble() |>
    dplyr::group_by(location, target_end_date) |>
    dplyr::reframe(
      # output_type_id: the quantile probability (e.g., 0.01, 0.025, ...)
      # round() fixes IEEE 754 floating-point drift from seq(), e.g.,
      # seq(0.05, 0.95, 0.05) produces 0.15000000000000002 instead of 0.15.
      # 3 decimal places covers our finest granularity (0.025, 0.975).
      output_type_id = round(quantile_levels, 3),
      # value: the empirical quantile of the simulated values at each level
      value = quantile(.sim, probs = quantile_levels, na.rm = TRUE)
    ) |>
    dplyr::ungroup()
  
  # --- Floor at zero ---
  # The hub requires all values to be non-negative (hub-config/tasks.json:
  # value.minimum = 0). ILI percentages can't be negative in reality, but
  # the Box-Cox back-transformation might produce small negatives in the
  # lower tails. pmax() clamps these to 0.
  # Reference: flucast-design-document.md Section 5.4, step 4
  quantile_df <- quantile_df |>
    dplyr::mutate(value = pmax(value, 0))
  
  return(quantile_df)
}


###############################################################################
# FUNCTION 4: format_for_hub
###############################################################################

#' Format quantile forecast data into the hub's 8-column CSV schema
#'
#' Takes a quantile dataframe (from compute_quantiles_from_sims) and an
#' origin date, and produces a dataframe with exactly the 8 columns required
#' by the FluSight ILI Sandbox Hub, in the correct order.
#'
#' Adds: origin_date, target, horizon, output_type columns.
#' Computes: horizon = (target_end_date - origin_date) / 7
#' Enforces: value >= 0, correct column names and types.
#'
#' This implements Step 5e from flucast-design-document.md Section 6.
#'
#' @param quantile_df A tibble with columns: location, target_end_date,
#'   output_type_id, value. Typically from compute_quantiles_from_sims().
#' @param origin_date A single Date object (the forecast origin Saturday).
#' @param team_model Character string like "KReger-snaive_bc_bs". Used only
#'   for logging/error messages, not written to the CSV.
#'
#' @return A tibble with exactly 8 columns in hub order:
#'   origin_date, location, target, horizon, target_end_date,
#'   output_type, output_type_id, value.
#'   Rows: should be 1,012 (11 locations × 4 horizons × 23 quantiles).
#'
#' @examples
#' # hub_df <- format_for_hub(quantile_df, as.Date("2016-12-03"), "KReger-snaive_bc_bs")
format_for_hub <- function(quantile_df, origin_date, team_model = "KReger-snaive_bc_bs") {
  
  # --- Input validation ---
  stopifnot(
    "origin_date must be a single Date" =
      inherits(origin_date, "Date") && length(origin_date) == 1,
    "quantile_df must have columns: location, target_end_date, output_type_id, value" =
      all(c("location", "target_end_date", "output_type_id", "value") %in%
            names(quantile_df))
  )
  
  hub_df <- quantile_df |>
    dplyr::mutate(
      # --- Add the origin date column ---
      # This is the same date as the filename date. Every row gets the same value.
      # We use !!DATE_COL via rlang to set the column name dynamically based on
      # the DATE_COL constant (resolved as "origin_date" for this hub).
      !!DATE_COL := origin_date,
      
      # --- Hardcode the target ---
      # This hub tracks only one target: "ili perc" (weighted ILI percentage).
      target = "ili perc",
      
      # --- Compute the horizon ---
      # Horizon is the number of WEEKS between the origin date and the
      # target_end_date. Since both are Saturdays spaced in whole weeks,
      # dividing the day difference by 7 gives an integer.
      horizon = as.integer(
        as.numeric(difftime(target_end_date, origin_date, units = "days")) / 7
      ),
      
      # --- Hardcode the output type ---
      # All rows are quantile forecasts (not point forecasts or pmf).
      output_type = "quantile",
      
      # --- Ensure non-negativity (belt AND suspenders) ---
      # This should already be handled by compute_quantiles_from_sims(),
      # but we enforce it again here as a safety net.
      value = pmax(value, 0)
    ) |>
    # --- Select exactly the 8 required columns in the standard order ---
    # "No additional columns are allowed" — model-output/README.md
    dplyr::select(
      dplyr::all_of(EXPECTED_COLS)
    )
  
  return(hub_df)
}


###############################################################################
# FUNCTION 5: validate_forecast_df
###############################################################################

#' Run all programmatic quality checks on a formatted forecast dataframe
#'
#' Implements EVERY check from flucast-design-document.md Section 9.1.
#' Each check prints a pass/fail indicator to the console. On any failure,
#' prints the specific violation details (which rows, what values) and then
#' stops execution with an error.
#'
#' This is meant to run BEFORE writing the CSV, catching errors early.
#' For final validation, use hubValidations::validate_submission() on the
#' written file.
#'
#' @param df A tibble formatted by format_for_hub() — should have exactly
#'   8 columns and 1,012 rows.
#' @param origin_date The Date that should match the origin_date column
#'   and the filename date.
#' @param quantile_levels Numeric vector of expected quantile levels.
#'   Defaults to QUANTILE_LEVELS.
#'
#' @return TRUE invisibly if all checks pass. Stops with an error on the
#'   first failure.
#'
#' @examples
#' # validate_forecast_df(hub_df, as.Date("2016-12-03"))
validate_forecast_df <- function(df, origin_date,
                                 quantile_levels = QUANTILE_LEVELS) {
  
  # Internal helper to print pass/fail and optionally stop
  check <- function(condition, description, detail = NULL) {
    if (condition) {
      cat(paste0("    \u2713 ", description, "\n"))
    } else {
      cat(paste0("    \u2717 ", description, "\n"))
      if (!is.null(detail)) cat(paste0("      Detail: ", detail, "\n"))
      stop("Validation failed: ", description, call. = FALSE)
    }
  }
  
  cat(paste0("  Validating forecast for ", origin_date, "...\n"))
  
  # -------------------------------------------------------------------------
  # CHECK 1: Column names — exactly the 8 expected, no extras
  # Reference: Section 9.1 — "Column names" and "No extra columns"
  # -------------------------------------------------------------------------
  check(
    identical(sort(names(df)), sort(EXPECTED_COLS)),
    "Column names match expected",
    detail = paste0(
      "Expected: ", paste(sort(EXPECTED_COLS), collapse = ", "),
      " | Got: ", paste(sort(names(df)), collapse = ", ")
    )
  )
  
  check(
    ncol(df) == 8,
    paste0("Column count is 8 (got ", ncol(df), ")")
  )
  
  # -------------------------------------------------------------------------
  # CHECK 2: Row count — should be exactly 1,012
  # Reference: Section 9.1 — "Row count"
  # -------------------------------------------------------------------------
  check(
    nrow(df) == EXPECTED_ROWS,
    paste0("Row count is ", EXPECTED_ROWS, " (got ", nrow(df), ")")
  )
  
  # -------------------------------------------------------------------------
  # CHECK 3: Non-negative values — all value >= 0
  # Reference: Section 9.1 — "Non-negative values"
  # -------------------------------------------------------------------------
  neg_count <- sum(df$value < 0, na.rm = TRUE)
  check(
    neg_count == 0,
    "All values are non-negative",
    detail = if (neg_count > 0) paste0(neg_count, " negative values found")
  )
  
  # -------------------------------------------------------------------------
  # CHECK 4: No NA values — critical for hub submission
  # -------------------------------------------------------------------------
  na_count <- sum(is.na(df$value))
  check(
    na_count == 0,
    "No NA values",
    detail = if (na_count > 0) paste0(na_count, " NA values found")
  )
  
  # -------------------------------------------------------------------------
  # CHECK 5: Quantile count per location-horizon group — should be 23
  # Reference: Section 9.1 — "Quantile count per group"
  # -------------------------------------------------------------------------
  q_counts <- df |>
    dplyr::group_by(location, horizon) |>
    dplyr::summarise(n_q = dplyr::n_distinct(output_type_id), .groups = "drop")
  
  check(
    all(q_counts$n_q == length(quantile_levels)),
    paste0("Each location-horizon group has ", length(quantile_levels), " quantiles"),
    detail = if (!all(q_counts$n_q == length(quantile_levels))) {
      bad <- q_counts |> dplyr::filter(n_q != length(quantile_levels))
      paste0("Mismatched groups: ", nrow(bad))
    }
  )
  
  # -------------------------------------------------------------------------
  # CHECK 6: Quantile levels match the expected set exactly
  # Reference: Section 9.1 — "Quantile levels match"
  # -------------------------------------------------------------------------
  actual_levels <- sort(unique(df$output_type_id))
  extra <- setdiff(actual_levels, quantile_levels)
  missing_q <- setdiff(quantile_levels, actual_levels)
  check(
    length(extra) == 0 && length(missing_q) == 0,
    "Quantile levels match expected set",
    detail = if (length(extra) > 0 || length(missing_q) > 0) {
      paste0("Extra: ", paste(extra, collapse = ", "),
             " | Missing: ", paste(missing_q, collapse = ", "))
    }
  )
  
  # -------------------------------------------------------------------------
  # CHECK 7: Quantile monotonicity — within each location-horizon group,
  # values must be non-decreasing as output_type_id increases.
  # Reference: Section 9.1 — "Quantile monotonicity"
  # -------------------------------------------------------------------------
  mono_check <- df |>
    dplyr::group_by(location, horizon) |>
    dplyr::arrange(output_type_id, .by_group = TRUE) |>
    dplyr::summarise(
      # Check that each successive value is >= the previous value.
      # diff(value) gives the change between consecutive quantiles;
      # all should be >= 0 for monotonicity.
      is_monotonic = all(diff(value) >= -1e-10),
      .groups = "drop"
    )
  
  check(
    all(mono_check$is_monotonic),
    "Quantile values are monotonically non-decreasing",
    detail = if (!all(mono_check$is_monotonic)) {
      bad <- mono_check |> dplyr::filter(!is_monotonic)
      paste0(nrow(bad), " groups violate monotonicity")
    }
  )
  
  # -------------------------------------------------------------------------
  # CHECK 8: Horizon range — should be exactly 1, 2, 3, 4
  # Reference: Section 9.1 — "Horizon range"
  # -------------------------------------------------------------------------
  actual_horizons <- sort(unique(df$horizon))
  check(
    identical(actual_horizons, as.integer(HUB_HORIZONS)),
    paste0("Horizons are {1,2,3,4} (got {", paste(actual_horizons, collapse = ","), "})"),
    detail = paste0("Expected: ", paste(HUB_HORIZONS, collapse = ","))
  )
  
  # -------------------------------------------------------------------------
  # CHECK 9: Date arithmetic — target_end_date == origin_date + horizon * 7
  # Reference: Section 9.1 — "Date arithmetic"
  # -------------------------------------------------------------------------
  date_check <- df |>
    dplyr::mutate(
      expected_ted = .data[[DATE_COL]] + horizon * 7,
      date_ok = (target_end_date == expected_ted)
    )
  
  check(
    all(date_check$date_ok),
    "target_end_date == origin_date + horizon * 7",
    detail = if (!all(date_check$date_ok)) {
      bad <- date_check |> dplyr::filter(!date_ok)
      paste0(nrow(bad), " rows have incorrect target_end_date")
    }
  )
  
  # -------------------------------------------------------------------------
  # CHECK 10: origin_date matches the expected date parameter
  # Reference: Section 9.1 — "reference_date matches filename"
  # -------------------------------------------------------------------------
  unique_dates <- unique(df[[DATE_COL]])
  check(
    length(unique_dates) == 1 && unique_dates == origin_date,
    paste0("origin_date column matches expected (", origin_date, ")"),
    detail = paste0("Found: ", paste(unique_dates, collapse = ", "))
  )
  
  # -------------------------------------------------------------------------
  # CHECK 11: origin_date is a Saturday
  # Reference: Section 9.1 — "reference_date is Saturday"
  # -------------------------------------------------------------------------
  check(
    lubridate::wday(origin_date) == 7,
    paste0("origin_date (", origin_date, ") is a Saturday"),
    detail = paste0("Day of week: ", lubridate::wday(origin_date, label = TRUE))
  )
  
  # -------------------------------------------------------------------------
  # CHECK 12: Location set — must be exactly the 11 expected locations
  # Reference: Section 9.1 — "Location set"
  # -------------------------------------------------------------------------
  actual_locs <- sort(unique(df$location))
  check(
    setequal(actual_locs, HUB_LOCATIONS),
    "Location set matches all 11 expected locations",
    detail = if (!setequal(actual_locs, HUB_LOCATIONS)) {
      paste0("Missing: ", paste(setdiff(HUB_LOCATIONS, actual_locs), collapse = ", "),
             " | Extra: ", paste(setdiff(actual_locs, HUB_LOCATIONS), collapse = ", "))
    }
  )
  
  # -------------------------------------------------------------------------
  # CHECK 13: target value — must be "ili perc"
  # Reference: Section 9.1 — "target value"
  # -------------------------------------------------------------------------
  check(
    identical(unique(df$target), "ili perc"),
    'target is "ili perc"',
    detail = paste0("Found: ", paste(unique(df$target), collapse = ", "))
  )
  
  # -------------------------------------------------------------------------
  # CHECK 14: output_type value — must be "quantile"
  # Reference: Section 9.1 — "output_type value"
  # -------------------------------------------------------------------------
  check(
    identical(unique(df$output_type), "quantile"),
    'output_type is "quantile"',
    detail = paste0("Found: ", paste(unique(df$output_type), collapse = ", "))
  )
  
  cat("  All checks passed!\n")
  invisible(TRUE)
}


###############################################################################
# FUNCTION 6: write_hub_csv
###############################################################################

#' Write a forecast dataframe to CSV with the correct hub filename convention
#'
#' Constructs the filename as YYYY-MM-DD-team-model.csv and writes the CSV
#' to the specified output directory. Creates the directory if it doesn't exist.
#'
#' @param df A tibble formatted by format_for_hub() — should have exactly
#'   8 columns and 1,012 rows.
#' @param origin_date A single Date object for this forecast.
#' @param output_dir Character string. Path to the output directory
#'   (e.g., "model-output/KReger-snaive_bc_bs").
#' @param team_model Character string (e.g., "KReger-snaive_bc_bs"). Used
#'   to construct the filename.
#'
#' @return The full file path of the written CSV (invisibly).
#'
#' @examples
#' # write_hub_csv(hub_df, as.Date("2016-12-03"),
#' #              "model-output/KReger-snaive_bc_bs", "KReger-snaive_bc_bs")
write_hub_csv <- function(df, origin_date, output_dir, team_model) {
  
  # --- Create the output directory if it doesn't exist yet ---
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # --- Construct the filename ---
  # Pattern: YYYY-MM-DD-team-model.csv
  # Example: 2016-12-03-KReger-snaive_bc_bs.csv
  # Reference: model-output/README.md — "Forecasts" section
  filename <- paste0(origin_date, "-", team_model, ".csv")
  filepath <- file.path(output_dir, filename)
  
  # --- Write the CSV ---
  # readr::write_csv() writes without row names and uses consistent
  # date formatting (YYYY-MM-DD). This matches the hub's expected format.
  readr::write_csv(df, filepath)
  
  return(invisible(filepath))
}


###############################################################################
# FUNCTION 7: plot_quantile_fan
###############################################################################

#' Generate a spot-check fan chart for a single location and origin date
#'
#' Plots the historical time series, the forecast quantile fan, and the
#' actual observed values (if available) for visual QA. Useful for quickly
#' checking that forecasts look reasonable (seasonal shape, widening
#' intervals, no negatives, etc.).
#'
#' @param forecast_df A tibble formatted by format_for_hub() for ONE
#'   origin date (1,012 rows covering all locations).
#' @param actuals_tsibble The full wILI tsibble (for plotting history
#'   and future actuals to compare against).
#' @param origin_date A Date object — the forecast origin.
#' @param location Character string — which location to plot
#'   (e.g., "US National").
#'
#' @return A ggplot2 object.
#'
#' @examples
#' # p <- plot_quantile_fan(hub_df, wili_tsibble, as.Date("2016-12-03"), "US National")
#' # print(p)
plot_quantile_fan <- function(forecast_df, actuals_tsibble,
                              origin_date, location = "US National") {
  
  # --- Filter forecast to the requested location ---
  fc_loc <- forecast_df |>
    dplyr::filter(location == !!location)
  
  if (nrow(fc_loc) == 0) {
    stop("No forecast data found for location '", location, "'")
  }
  
  # --- Pivot quantiles wider for plotting ---
  # We need columns for specific quantile levels to draw the fan.
  # Extract a few key quantile bands: 1-99%, 10-90%, 25-75%, and median
  fc_wide <- fc_loc |>
    dplyr::select(target_end_date, output_type_id, value) |>
    tidyr::pivot_wider(names_from = output_type_id, values_from = value,
                       names_prefix = "q")
  
  # --- Filter actuals to a window around the origin date ---
  # Show ~6 months before and ~2 months after for context
  window_start <- origin_date - lubridate::weeks(26)
  window_end   <- origin_date + lubridate::weeks(8)
  
  actuals_window <- actuals_tsibble |>
    tibble::as_tibble() |>
    dplyr::filter(
      location == !!location,
      target_end_date >= window_start,
      target_end_date <= window_end
    )
  
  # --- Build the plot ---
  p <- ggplot2::ggplot() +
    # 1-99% prediction interval (lightest band)
    ggplot2::geom_ribbon(
      data = fc_wide,
      ggplot2::aes(x = target_end_date, ymin = q0.01, ymax = q0.99),
      fill = "steelblue", alpha = 0.15
    ) +
    # 10-90% prediction interval
    ggplot2::geom_ribbon(
      data = fc_wide,
      ggplot2::aes(x = target_end_date, ymin = q0.1, ymax = q0.9),
      fill = "steelblue", alpha = 0.25
    ) +
    # 25-75% prediction interval (darkest band)
    ggplot2::geom_ribbon(
      data = fc_wide,
      ggplot2::aes(x = target_end_date, ymin = q0.25, ymax = q0.75),
      fill = "steelblue", alpha = 0.4
    ) +
    # Median forecast line
    ggplot2::geom_line(
      data = fc_wide,
      ggplot2::aes(x = target_end_date, y = q0.5),
      color = "steelblue", linewidth = 1
    ) +
    # Historical + future actuals (black line and dots)
    ggplot2::geom_line(
      data = actuals_window,
      ggplot2::aes(x = target_end_date, y = observation),
      color = "black", linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = actuals_window |>
        dplyr::filter(target_end_date > origin_date),
      ggplot2::aes(x = target_end_date, y = observation),
      color = "red", size = 2
    ) +
    # Vertical dashed line at the origin date
    ggplot2::geom_vline(
      xintercept = as.numeric(origin_date),
      linetype = "dashed", color = "gray40"
    ) +
    # Labels and theme
    ggplot2::labs(
      title = paste0("Forecast fan chart: ", location),
      subtitle = paste0("Origin date: ", origin_date),
      x = "Date", y = "wILI %"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}