###############################################################################
# FILE: scripts/forecast_helpers.R
#
# Shared, model-agnostic utility functions for the FluSight ILI Sandbox Hub.
# These are sourced automatically by forecast_engine.R and used by all model
# configurations. You should never need to edit this file to add a new model.
#
# Functions:
#   load_and_prepare_tsibble()    — Load CSV, create tsibble
#   get_origin_dates()            — Parse origin dates from tasks.json
#   compute_quantiles_from_sims() — Empirical quantiles from bootstrap sims
#   format_for_hub()              — Shape data into the 8-column hub format
#   validate_forecast_df()        — Run all quality checks
#   write_hub_csv()               — Write CSV with correct filename pattern
#   plot_quantile_fan()           — Spot-check visualization
#
# Constants:
#   DATE_COL, QUANTILE_LEVELS, HUB_LOCATIONS, HUB_HORIZONS,
#   EXPECTED_ROWS, EXPECTED_COLS
#
# Dependencies: tidyverse, tsibble, fable, fabletools, jsonlite, lubridate
#
###############################################################################


# =============================================================================
# CONSTANTS — Hub-specific values (single source of truth)
# =============================================================================

# The validated CSVs in this hub use "origin_date" (not "reference_date").
# Confirmed via hubValidations and existing model-output CSVs.
DATE_COL <- "origin_date"

# 23 required quantile levels from hub-config/tasks.json.
# round(3) fixes IEEE 754 drift from seq() (e.g., 0.15000000000000002).
QUANTILE_LEVELS <- round(c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99), 3)

# 11 locations: US national + 10 HHS Regions (hub-config/tasks.json)
HUB_LOCATIONS <- c("US National", paste("HHS Region", 1:10))

# Horizons: 1 through 4 weeks ahead
HUB_HORIZONS <- 1:4

# Expected row count: 11 locations × 4 horizons × 23 quantiles = 1,012
EXPECTED_ROWS <- length(HUB_LOCATIONS) * length(HUB_HORIZONS) * length(QUANTILE_LEVELS)

# Expected column names in hub submission order
EXPECTED_COLS <- c(
  DATE_COL, "location", "target", "horizon",
  "target_end_date", "output_type", "output_type_id", "value"
)


###############################################################################
#' Load target data CSV and convert to a tsibble
#'
#' Reads the wILI percentage time series from a CSV file (local path or URL)
#' and converts it into a tsibble indexed by target_end_date (weekly Saturdays)
#' and keyed by location and target.
#'
#' @param data_path Character string. Path to time-series.csv.
#' @return A tsibble with columns: location, target_end_date, target, observation.
#'
#' @examples
#' wili <- load_and_prepare_tsibble("target-data/time-series.csv")
###############################################################################
load_and_prepare_tsibble <- function(data_path) {
  
  stopifnot(
    "data_path must be a single character string" =
      is.character(data_path) && length(data_path) == 1
  )
  
  wili_raw <- readr::read_csv(data_path, show_col_types = FALSE)
  
  # Verify expected columns exist
  required_cols <- c("location", "target_end_date", "target", "observation")
  missing_cols <- setdiff(required_cols, names(wili_raw))
  if (length(missing_cols) > 0) {
    stop(
      "Target data CSV is missing columns: ",
      paste(missing_cols, collapse = ", "),
      "\nFound: ", paste(names(wili_raw), collapse = ", ")
    )
  }
  
  # Convert to tsibble. Including "target" in the key makes this compatible
  # with hubs that track multiple targets (this hub only has "ili perc").
  wili_tsibble <- wili_raw |>
    tsibble::as_tsibble(index = target_end_date, key = c(location, target))
  
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
#' Parse origin dates from the hub's tasks.json configuration
#'
#' Reads hub-config/tasks.json and extracts the sorted Date vector of all
#' origin dates that forecasts should be submitted for.
#'
#' @param tasks_json_path Character string. Path to hub-config/tasks.json.
#' @return A sorted Date vector of origin dates (all Saturdays). 139 for this hub.
#'
#' @examples
#' dates <- get_origin_dates("hub-config/tasks.json")
###############################################################################
get_origin_dates <- function(tasks_json_path) {
  
  stopifnot(
    "tasks_json_path must be a single character string" =
      is.character(tasks_json_path) && length(tasks_json_path) == 1
  )
  
  # simplifyVector = FALSE prevents jsonlite from flattening nested lists
  # into atomic vectors, which would break $ navigation of the structure.
  tasks <- jsonlite::fromJSON(tasks_json_path, simplifyVector = FALSE)
  
  # Navigate: rounds[[1]] -> model_tasks[[1]] -> task_ids -> origin_date
  task_ids <- tasks$rounds[[1]]$model_tasks[[1]]$task_ids
  
  # Try "origin_date" first, fall back to "reference_date" for other hubs
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
  
  origin_dates <- sort(as.Date(raw_dates))
  
  if (!all(lubridate::wday(origin_dates) == 7)) {
    warning("Not all origin dates are Saturdays — check tasks.json")
  }
  
  cat(paste0("  Loaded ", length(origin_dates), " origin dates (",
             min(origin_dates), " to ", max(origin_dates), ")\n"))
  
  return(origin_dates)
}


###############################################################################
#' Compute empirical quantiles from bootstrap simulation paths
#'
#' Takes fable's generate(bootstrap = TRUE) output and computes requested
#' quantile levels for each location × target_end_date group. Floors all
#' values at zero to satisfy the hub's non-negativity constraint.
#'
#' @param sims_tsibble A tsibble from fable::generate() with a .sim column.
#' @param quantile_levels Numeric vector of quantile probabilities.
#'   Defaults to the hub's 23 required levels.
#' @return A tibble with columns: location, target_end_date, output_type_id, value.
#'
#' @examples
#' # quantile_df <- compute_quantiles_from_sims(sims)
###############################################################################
compute_quantiles_from_sims <- function(sims_tsibble,
                                        quantile_levels = QUANTILE_LEVELS) {
  
  # Validate required columns
  required <- c(".sim", "location", "target_end_date")
  missing <- setdiff(required, names(sims_tsibble))
  if (length(missing) > 0) {
    stop("sims_tsibble is missing required columns: ", paste(missing, collapse = ", "))
  }
  
  # Group by location × date, compute quantile at each level.
  # reframe() returns multiple rows per group (one per quantile level).
  quantile_df <- sims_tsibble |>
    tibble::as_tibble() |>
    dplyr::group_by(location, target_end_date) |>
    dplyr::reframe(
      output_type_id = round(quantile_levels, 3),
      value = quantile(.sim, probs = quantile_levels, na.rm = TRUE)
    ) |>
    dplyr::ungroup()
  
  # Floor at zero: ILI% can't be negative, but Box-Cox back-transformation
  # can produce small negatives in the lower tails.
  quantile_df <- quantile_df |>
    dplyr::mutate(value = pmax(value, 0))
  
  return(quantile_df)
}


###############################################################################
#' Format quantile forecast data into the hub's 8-column CSV schema
#'
#' Adds origin_date, target, horizon, and output_type columns. Computes
#' horizon = (target_end_date - origin_date) / 7. Selects exactly the 8
#' required columns in submission order.
#'
#' @param quantile_df A tibble from compute_quantiles_from_sims().
#' @param origin_date A single Date (the forecast origin Saturday).
#' @param team_model Character string (e.g., "KReger-snaive_bc_bs"). Required.
#' @return A tibble with 8 columns in hub order, typically 1,012 rows.
#'
#' @examples
#' # hub_df <- format_for_hub(quantile_df, as.Date("2016-12-03"), "KReger-snaive_bc_bs")
###############################################################################
format_for_hub <- function(quantile_df, origin_date, team_model) {
  
  stopifnot(
    "origin_date must be a single Date" =
      inherits(origin_date, "Date") && length(origin_date) == 1,
    "team_model must be a non-empty character string" =
      is.character(team_model) && length(team_model) == 1 && nchar(team_model) > 0,
    "quantile_df must have columns: location, target_end_date, output_type_id, value" =
      all(c("location", "target_end_date", "output_type_id", "value") %in%
            names(quantile_df))
  )
  
  hub_df <- quantile_df |>
    dplyr::mutate(
      # Set the date column name dynamically via rlang (resolves to "origin_date")
      !!DATE_COL := origin_date,
      target = "ili perc",
      # Both dates are Saturdays, so day difference / 7 = integer horizon
      horizon = as.integer(
        as.numeric(difftime(target_end_date, origin_date, units = "days")) / 7
      ),
      output_type = "quantile",
      # Safety net: enforce non-negativity even if already floored upstream
      value = pmax(value, 0)
    ) |>
    dplyr::select(dplyr::all_of(EXPECTED_COLS))
  
  return(hub_df)
}


###############################################################################
#' Run all programmatic quality checks on a formatted forecast dataframe
#'
#' Each check prints a pass/fail indicator. On failure, prints violation
#' details and stops with an error. Run this before writing the CSV.
#' For final validation, use hubValidations::validate_submission() on the file.
#'
#' @param df A tibble from format_for_hub() (8 columns, 1,012 rows).
#' @param origin_date The Date that should match the origin_date column.
#' @param quantile_levels Expected quantile levels. Defaults to QUANTILE_LEVELS.
#' @return TRUE invisibly if all checks pass.
#'
#' @examples
#' # validate_forecast_df(hub_df, as.Date("2016-12-03"))
###############################################################################
validate_forecast_df <- function(df, origin_date,
                                 quantile_levels = QUANTILE_LEVELS) {
  
  # Helper: print pass/fail, stop on failure
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
  
  # CHECK 1: Column names
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
  
  # CHECK 2: Row count
  check(
    nrow(df) == EXPECTED_ROWS,
    paste0("Row count is ", EXPECTED_ROWS, " (got ", nrow(df), ")")
  )
  
  # CHECK 3: Non-negative values
  neg_count <- sum(df$value < 0, na.rm = TRUE)
  check(
    neg_count == 0,
    "All values are non-negative",
    detail = if (neg_count > 0) paste0(neg_count, " negative values found")
  )
  
  # CHECK 4: No NAs
  na_count <- sum(is.na(df$value))
  check(
    na_count == 0,
    "No NA values",
    detail = if (na_count > 0) paste0(na_count, " NA values found")
  )
  
  # CHECK 5: 23 quantiles per location-horizon group
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
  
  # CHECK 6: Quantile levels match expected set
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
  
  # CHECK 7: Quantile monotonicity within each location-horizon group
  mono_check <- df |>
    dplyr::group_by(location, horizon) |>
    dplyr::arrange(output_type_id, .by_group = TRUE) |>
    dplyr::summarise(
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
  
  # CHECK 8: Horizons are {1,2,3,4}
  actual_horizons <- sort(unique(as.integer(df$horizon)))
  check(
    identical(actual_horizons, as.integer(HUB_HORIZONS)),
    paste0("Horizons are {1,2,3,4} (got {", paste(actual_horizons, collapse = ","), "})"),
    detail = paste0("Expected: ", paste(HUB_HORIZONS, collapse = ","))
  )
  
  # CHECK 9: Date arithmetic — target_end_date == origin_date + horizon * 7
  date_check <- df |>
    dplyr::mutate(
      expected_ted = .data[[DATE_COL]] + horizon * 7,
      date_ok = (target_end_date == expected_ted)
    )
  
  check(
    all(date_check$date_ok),
    "target_end_date == origin_date + horizon * 7",
    detail = if (!all(date_check$date_ok)) {
      paste0(sum(!date_check$date_ok), " rows have incorrect target_end_date")
    }
  )
  
  # CHECK 10: origin_date column matches expected date
  unique_dates <- unique(df[[DATE_COL]])
  check(
    length(unique_dates) == 1 && unique_dates == origin_date,
    paste0("origin_date column matches expected (", origin_date, ")"),
    detail = paste0("Found: ", paste(unique_dates, collapse = ", "))
  )
  
  # CHECK 11: origin_date is a Saturday (wday: 7 = Sat)
  check(
    lubridate::wday(origin_date) == 7,
    paste0("origin_date (", origin_date, ") is a Saturday"),
    detail = paste0("Day of week: ", lubridate::wday(origin_date, label = TRUE))
  )
  
  # CHECK 12: Location set matches all 11 expected
  actual_locs <- sort(unique(df$location))
  check(
    setequal(actual_locs, HUB_LOCATIONS),
    "Location set matches all 11 expected locations",
    detail = if (!setequal(actual_locs, HUB_LOCATIONS)) {
      paste0("Missing: ", paste(setdiff(HUB_LOCATIONS, actual_locs), collapse = ", "),
             " | Extra: ", paste(setdiff(actual_locs, HUB_LOCATIONS), collapse = ", "))
    }
  )
  
  # CHECK 13: target is "ili perc"
  check(
    identical(unique(df$target), "ili perc"),
    'target is "ili perc"',
    detail = paste0("Found: ", paste(unique(df$target), collapse = ", "))
  )
  
  # CHECK 14: output_type is "quantile"
  check(
    identical(unique(df$output_type), "quantile"),
    'output_type is "quantile"',
    detail = paste0("Found: ", paste(unique(df$output_type), collapse = ", "))
  )
  
  cat("  All checks passed!\n")
  invisible(TRUE)
}


###############################################################################
#' Write a forecast dataframe to CSV with the hub filename convention
#'
#' Filename pattern: YYYY-MM-DD-team-model.csv. Creates the output
#' directory if it doesn't exist.
#'
#' @param df A tibble from format_for_hub().
#' @param origin_date A single Date.
#' @param output_dir Path to output directory.
#' @param team_model Character string (e.g., "KReger-snaive_bc_bs").
#' @return The full file path (invisibly).
###############################################################################
write_hub_csv <- function(df, origin_date, output_dir, team_model) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  filename <- paste0(origin_date, "-", team_model, ".csv")
  filepath <- file.path(output_dir, filename)
  readr::write_csv(df, filepath)
  
  return(invisible(filepath))
}


###############################################################################
#' Generate a spot-check fan chart for a single location and origin date
#'
#' Plots historical actuals, the forecast quantile fan (1-99%, 10-90%,
#' 25-75%, median), and future actuals (red dots) for visual QA.
#'
#' @param forecast_df A tibble from format_for_hub() for one origin date.
#' @param actuals_tsibble The full wILI tsibble.
#' @param origin_date A Date — the forecast origin.
#' @param location Which location to plot (default: "US National").
#' @return A ggplot2 object.
#'
#' @examples
#' # p <- plot_quantile_fan(hub_df, wili_tsibble, as.Date("2016-12-03"), "US National")
###############################################################################
plot_quantile_fan <- function(forecast_df, actuals_tsibble,
                              origin_date, location = "US National") {
  
  fc_loc <- forecast_df |>
    dplyr::filter(location == !!location)
  
  if (nrow(fc_loc) == 0) {
    stop("No forecast data found for location '", location, "'")
  }
  
  # Pivot quantiles to wide format for ribbon plotting
  fc_wide <- fc_loc |>
    dplyr::select(target_end_date, output_type_id, value) |>
    tidyr::pivot_wider(names_from = output_type_id, values_from = value,
                       names_prefix = "q")
  
  # Show ~6 months before and ~2 months after the origin date
  window_start <- origin_date - lubridate::weeks(26)
  window_end   <- origin_date + lubridate::weeks(8)
  
  actuals_window <- actuals_tsibble |>
    tibble::as_tibble() |>
    dplyr::filter(
      location == !!location,
      target_end_date >= window_start,
      target_end_date <= window_end
    )
  
  p <- ggplot2::ggplot() +
    # Prediction interval bands (lightest to darkest)
    ggplot2::geom_ribbon(
      data = fc_wide,
      ggplot2::aes(x = target_end_date, ymin = q0.01, ymax = q0.99),
      fill = "steelblue", alpha = 0.15
    ) +
    ggplot2::geom_ribbon(
      data = fc_wide,
      ggplot2::aes(x = target_end_date, ymin = q0.1, ymax = q0.9),
      fill = "steelblue", alpha = 0.25
    ) +
    ggplot2::geom_ribbon(
      data = fc_wide,
      ggplot2::aes(x = target_end_date, ymin = q0.25, ymax = q0.75),
      fill = "steelblue", alpha = 0.4
    ) +
    # Median forecast
    ggplot2::geom_line(
      data = fc_wide,
      ggplot2::aes(x = target_end_date, y = q0.5),
      color = "steelblue", linewidth = 1
    ) +
    # Historical actuals (black line)
    ggplot2::geom_line(
      data = actuals_window,
      ggplot2::aes(x = target_end_date, y = observation),
      color = "black", linewidth = 0.5
    ) +
    # Future actuals (red dots — what we're trying to predict)
    ggplot2::geom_point(
      data = actuals_window |>
        dplyr::filter(target_end_date > origin_date),
      ggplot2::aes(x = target_end_date, y = observation),
      color = "red", size = 2
    ) +
    ggplot2::geom_vline(
      xintercept = as.numeric(origin_date),
      linetype = "dashed", color = "gray40"
    ) +
    ggplot2::labs(
      title = paste0("Forecast fan chart: ", location),
      subtitle = paste0("Origin date: ", origin_date),
      x = "Date", y = "wILI %"
    ) +
    ggplot2::theme_minimal()
  
  return(p)
}