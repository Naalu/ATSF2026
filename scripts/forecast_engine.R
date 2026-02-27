###############################################################################
# FILE: scripts/forecast_engine.R
#
# Generic forecasting engine for the FluSight ILI Sandbox Hub (ATSF2026).
# Takes a model configuration list and runs the complete pipeline:
#   load data → preprocess → fit/simulate → quantiles → format → validate → write
#
# This file is model-agnostic. All model-specific logic lives in config files
# (e.g., configs/snaive_bc_bs.R). The engine handles two paths:
#   1. Fable path — config provides model_fn; engine runs per-location
#      model() |> generate() pipeline
#   2. Custom path — config provides custom_simulate_fn; engine delegates
#      simulation entirely to the user's function
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/snaive_bc_bs.R")    # defines `config`
#   result <- run_forecast(config)
#   run_diagnostics(result)              # optional: visual QA plots
#
# Provides:
#   run_forecast(config)     — run the full pipeline, return result object
#   run_diagnostics(result)  — produce diagnostic plots from a result object
#
# Dependencies:
#   fpp3, tidyverse (loaded automatically)
#   patchwork (required only for run_diagnostics)
#
###############################################################################


###############################################################################
#' Run the full forecast pipeline for a given model configuration
#'
#' @param config A named list defining the model. See configs/snaive_bc_bs.R
#'   for the full specification and examples.
#' @return A list (invisibly) containing everything needed for diagnostics
#'   and ad-hoc forecasts: config, team_model, output_dir, wili_tsibble,
#'   origin_dates, preprocess, and a callable forecast() function.
###############################################################################
run_forecast <- function(config) {
  
  # =========================================================================
  # VALIDATE CONFIG
  # =========================================================================
  
  required <- c("team_abbr", "model_abbr", "hub_path", "n_sims",
                "min_train_weeks", "dev_mode", "run_loop",
                "transform", "simulate_method")
  missing <- setdiff(required, names(config))
  if (length(missing) > 0) {
    stop("Config missing required fields: ", paste(missing, collapse = ", "))
  }
  
  if (is.null(config$model_fn) && is.null(config$custom_simulate_fn)) {
    stop("Config must provide either model_fn (fable) or custom_simulate_fn (custom)")
  }
  
  # =========================================================================
  # SETUP
  # =========================================================================
  
  # Derived identifiers and paths
  team_model <- paste0(config$team_abbr, "-", config$model_abbr)
  data_path  <- file.path(config$hub_path, "target-data", "time-series.csv")
  tasks_path <- file.path(config$hub_path, "hub-config", "tasks.json")
  output_dir <- file.path(config$hub_path, "model-output", team_model)
  
  cat(sprintf("=== Forecast Engine: %s ===\n\n", team_model))
  
  # Load packages
  if (!requireNamespace("fpp3", quietly = TRUE)) stop("Package 'fpp3' required")
  if (!requireNamespace("tidyverse", quietly = TRUE)) stop("Package 'tidyverse' required")
  suppressPackageStartupMessages({
    library(fpp3)
    library(tidyverse)
  })
  
  # Source shared helper functions (constants + utilities)
  helpers_path <- file.path(config$hub_path, "scripts", "forecast_helpers.R")
  if (!file.exists(helpers_path)) {
    stop("forecast_helpers.R not found at: ", helpers_path)
  }
  source(helpers_path)
  
  # =========================================================================
  # STEP 1: LOAD DATA
  # =========================================================================
  
  cat("--- Loading target data ---\n")
  wili_tsibble <- load_and_prepare_tsibble(data_path)
  cat(sprintf("  %d rows, %s to %s, %d locations\n\n",
              nrow(wili_tsibble),
              min(wili_tsibble$target_end_date),
              max(wili_tsibble$target_end_date),
              length(unique(wili_tsibble$location))))
  
  # =========================================================================
  # STEP 2: LOAD ORIGIN DATES
  # =========================================================================
  
  cat("--- Loading origin dates ---\n")
  origin_dates <- get_origin_dates(tasks_path)
  cat("\n")
  
  # =========================================================================
  # STEP 3: PREPROCESSING
  #
  # For Box-Cox models, estimate per-location lambda via guerrero() using
  # data strictly before the first origin date (so lambda is fixed and
  # reproducible regardless of which origin date is being processed).
  # =========================================================================
  
  preprocess_result <- list()
  
  if (config$transform == "box_cox_guerrero") {
    cat("--- Estimating Box-Cox lambdas via guerrero ---\n")
    first_origin <- min(origin_dates)
    
    pre_data <- wili_tsibble |>
      dplyr::filter(target_end_date < first_origin)
    
    lambda_table <- pre_data |>
      fabletools::features(observation, features = guerrero)
    
    lambdas <- stats::setNames(lambda_table$lambda_guerrero, lambda_table$location)
    
    # Consolidated sanity checks
    if (length(lambdas) != length(HUB_LOCATIONS)) {
      warning("Expected ", length(HUB_LOCATIONS), " lambdas, got ", length(lambdas))
    }
    extreme <- lambdas[abs(lambdas) > 2]
    if (length(extreme) > 0) {
      warning("Lambdas outside [-2, 2]: ",
              paste(names(extreme), round(extreme, 3), sep = "=", collapse = ", "))
    }
    
    cat(sprintf("  %d lambdas, range [%.4f, %.4f]\n\n",
                length(lambdas), min(lambdas), max(lambdas)))
    preprocess_result$lambdas <- lambdas
  }
  
  # =========================================================================
  # BUILD THE SIMULATION FUNCTION
  #
  # The engine constructs an internal function that takes training data for
  # one origin date and returns simulated future paths for all locations.
  # This function closes over the config and preprocess results.
  #
  # Two paths:
  #   (a) Fable: engine loops over locations, calls config$model_fn(lambda)
  #       to get a model spec, then model() |> generate().
  #   (b) Custom: engine delegates entirely to config$custom_simulate_fn.
  # =========================================================================
  
  if (!is.null(config$custom_simulate_fn)) {
    
    # ---- Custom path ----
    simulate_one_date <- function(train_data, h, n_sims) {
      config$custom_simulate_fn(train_data, h, n_sims, preprocess_result)
    }
    
  } else {
    
    # ---- Fable path ----
    simulate_one_date <- function(train_data, h, n_sims) {
      locations <- unique(train_data$location)
      is_bootstrap <- identical(config$simulate_method, "bootstrap")
      
      purrr::map(locations, function(loc) {
        loc_train <- train_data |> dplyr::filter(location == loc)
        
        # Get per-location lambda (NA if no Box-Cox transform)
        loc_lambda <- if (!is.null(preprocess_result$lambdas)) {
          preprocess_result$lambdas[loc]
        } else {
          NA_real_
        }
        
        # Build model spec from config's model_fn, then fit
        spec <- config$model_fn(loc_lambda)
        loc_fit <- loc_train |> fabletools::model(m = spec)
        
        # Generate simulated future paths
        loc_sims <- loc_fit |>
          fabletools::generate(h = h, times = n_sims, bootstrap = is_bootstrap)
        
        tibble::as_tibble(loc_sims)
      }) |>
        dplyr::bind_rows()
    }
  }
  
  # =========================================================================
  # SINGLE-DATE FORECAST FUNCTION
  #
  # Produces a complete hub-formatted forecast for one origin date:
  #   filter training data → simulate → quantiles → format → validate → write
  #
  # This closure captures wili_tsibble, config, preprocess_result, etc.
  # from the enclosing run_forecast() scope. It's also attached to the
  # return object so the user can call result$forecast(date) for ad-hoc use.
  # =========================================================================
  
  generate_single_forecast <- function(origin_date, write_file = TRUE) {
    
    # Filter training data: use <= because the origin_date Saturday is the
    # last day of the epiweek, so that week's data is available.
    train_data <- wili_tsibble |>
      dplyr::filter(target_end_date <= origin_date)
    
    # Check minimum training data
    min_weeks <- train_data |>
      tibble::as_tibble() |>
      dplyr::group_by(location) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::pull(n) |>
      min()
    
    if (min_weeks < config$min_train_weeks) {
      warning("Skipping ", origin_date, ": only ", min_weeks,
              " weeks (need ", config$min_train_weeks, ")")
      return(NULL)
    }
    
    # Simulate future paths
    sims <- tryCatch(
      simulate_one_date(train_data, h = length(HUB_HORIZONS), n_sims = config$n_sims),
      error = function(e) {
        warning("Simulation failed for ", origin_date, ": ", e$message)
        NULL
      }
    )
    if (is.null(sims)) return(NULL)
    
    # Warn about NAs (can arise from time series gaps or numerical issues)
    na_count <- sum(is.na(sims$.sim))
    if (na_count > 0) {
      warning(origin_date, ": ", na_count, " NA values in simulations")
    }
    
    # Quantiles → format → validate → write
    quantile_df <- compute_quantiles_from_sims(sims)
    hub_df      <- format_for_hub(quantile_df, origin_date, team_model)
    validate_forecast_df(hub_df, origin_date)
    
    if (write_file) {
      filepath <- write_hub_csv(hub_df, origin_date, output_dir, team_model)
      cat(paste0("  Wrote: ", basename(filepath), "\n"))
    }
    
    return(hub_df)
  }
  
  # =========================================================================
  # STEP 5: MAIN FORECAST LOOP
  #
  # Iterates over all origin dates (or first 5 in dev_mode). Features:
  #   - Progress with ETA (based on actual processing time, not wall clock)
  #   - Skip-if-exists (for resuming interrupted runs)
  #   - Per-date error handling (failures logged, loop continues)
  #   - Post-loop summary and spot-check validation
  # =========================================================================
  
  if (config$run_loop) {
    
    cat("--- Running forecast loop ---\n")
    
    if (config$dev_mode) {
      dates_to_process <- origin_dates[1:min(5, length(origin_dates))]
      cat("  DEV_MODE: first", length(dates_to_process), "dates\n")
    } else {
      dates_to_process <- origin_dates
      cat("  PRODUCTION: all", length(dates_to_process), "dates\n")
    }
    cat("  n_sims =", config$n_sims, "\n\n")
    
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Tracking
    n_total      <- length(dates_to_process)
    n_success    <- 0
    n_skipped    <- 0
    n_failed     <- 0
    failed_log   <- character(0)
    loop_start   <- Sys.time()
    process_time <- 0  # cumulative time for non-skipped dates
    
    for (i in seq_along(dates_to_process)) {
      this_date <- dates_to_process[i]
      expected_path <- file.path(output_dir,
                                 paste0(this_date, "-", team_model, ".csv"))
      
      # Skip existing files (allows resuming interrupted runs)
      if (file.exists(expected_path)) {
        cat(sprintf("  [%3d/%d] %s — SKIPPED\n", i, n_total, this_date))
        n_skipped <- n_skipped + 1
        next
      }
      
      # Progress with ETA
      n_processed <- n_success + n_failed
      eta_str <- if (n_processed > 0) {
        remaining <- (n_total - n_skipped - n_processed - 1) *
          (process_time / n_processed)
        sprintf(" | ETA: %.0f min", remaining / 60)
      } else ""
      cat(sprintf("  [%3d/%d] %s%s\n", i, n_total, this_date, eta_str))
      
      # Generate forecast. suppressWarnings silences routine fable Box-Cox
      # back-transform messages that would flood the console.
      date_start <- Sys.time()
      result_df <- tryCatch(
        suppressWarnings(generate_single_forecast(this_date, write_file = TRUE)),
        error = function(e) {
          cat("    ✗ ERROR:", e$message, "\n")
          NULL
        }
      )
      process_time <- process_time +
        as.numeric(difftime(Sys.time(), date_start, units = "secs"))
      
      if (!is.null(result_df)) {
        n_success <- n_success + 1
      } else {
        n_failed <- n_failed + 1
        failed_log <- c(failed_log, as.character(this_date))
      }
    }
    
    # ---- Summary ----
    total_elapsed <- as.numeric(difftime(Sys.time(), loop_start, units = "secs"))
    cat(sprintf("\n%s\n  FORECAST LOOP COMPLETE\n%s\n\n",
                strrep("=", 60), strrep("=", 60)))
    cat(sprintf("  Total: %d | Success: %d | Skipped: %d | Failed: %d\n",
                n_total, n_success, n_skipped, n_failed))
    cat(sprintf("  Time: %.1f min (%.0f sec)\n",
                total_elapsed / 60, total_elapsed))
    if (n_success > 0) {
      cat(sprintf("  Avg: %.1f sec/date\n",
                  process_time / (n_success + n_failed)))
    }
    if (n_failed > 0) {
      cat("\n  FAILED:\n")
      for (d in failed_log) cat("    ✗", d, "\n")
    }
    
    # File count
    csv_files <- list.files(output_dir, pattern = "\\.csv$")
    cat(sprintf("\n  Files: %d / %d expected\n",
                length(csv_files), length(origin_dates)))
    
    # ---- Spot-check validation on 5 sample files ----
    cat("\n  Spot-checking 5 sample files...\n")
    sample_idx <- unique(round(
      quantile(seq_along(origin_dates), probs = c(0, 0.25, 0.5, 0.75, 1))
    ))
    for (j in sample_idx) {
      d <- origin_dates[j]
      f <- file.path(output_dir, paste0(d, "-", team_model, ".csv"))
      if (file.exists(f)) {
        df <- readr::read_csv(f, show_col_types = FALSE)
        cat(sprintf("    %s...\n", basename(f)))
        tryCatch(
          validate_forecast_df(df, d),
          error = function(e) cat("      ✗", e$message, "\n")
        )
      }
    }
    cat("\n")
    
  } else {
    cat("--- Forecast loop SKIPPED (run_loop = FALSE) ---\n")
    cat("  Set run_loop = TRUE in config, or use: result$forecast(date)\n\n")
  }
  
  # =========================================================================
  # RETURN RESULT OBJECT
  #
  # Contains everything needed for diagnostics, ad-hoc forecasts, and
  # downstream analysis. The forecast() function is a closure that still
  # has access to wili_tsibble, lambdas, etc.
  # =========================================================================
  
  result <- list(
    config       = config,
    team_model   = team_model,
    output_dir   = output_dir,
    wili_tsibble = wili_tsibble,
    origin_dates = origin_dates,
    preprocess   = preprocess_result,
    forecast     = generate_single_forecast
  )
  
  cat(sprintf("=== Engine complete: %s ===\n", team_model))
  cat("  Use result$forecast(as.Date(\"YYYY-MM-DD\")) for ad-hoc forecasts\n")
  cat("  Use run_diagnostics(result) for visual QA\n")
  
  invisible(result)
}


###############################################################################
#' Produce diagnostic plots from a forecast engine result
#'
#' Generates three types of plots:
#'   (a) Fan charts for 5 representative origin dates × 3 locations
#'   (b) Sanity dashboard: horizon-1 median vs actuals + 90% PI width
#'   (c) Lambda summary (only if transform = "box_cox_guerrero")
#'
#' @param result The list returned by run_forecast().
#' @param qa_dates Optional Date vector of origin dates to plot fan charts for.
#'   Defaults to 5 dates spanning different flu season phases.
#' @param qa_locations Optional character vector of locations for fan charts.
#'   Defaults to c("US National", "HHS Region 4", "HHS Region 10").
#' @param dashboard_loc Location for the sanity dashboard (default: "US National").
###############################################################################
run_diagnostics <- function(result,
                            qa_dates = NULL,
                            qa_locations = c("US National",
                                             "HHS Region 4",
                                             "HHS Region 10"),
                            dashboard_loc = "US National") {
  
  # patchwork is only needed for diagnostics, not forecasting
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' required for diagnostics: install.packages(\"patchwork\")")
  }
  library(patchwork)
  
  # Unpack result for readability
  cfg          <- result$config
  team_model   <- result$team_model
  output_dir   <- result$output_dir
  wili_tsibble <- result$wili_tsibble
  origin_dates <- result$origin_dates
  
  cat(sprintf("\n=== Diagnostics for %s ===\n", team_model))
  
  # Default QA dates span different flu season phases
  if (is.null(qa_dates)) {
    qa_dates <- as.Date(c(
      "2015-11-07",  # Early season
      "2017-01-14",  # Near winter peak
      "2018-03-10",  # Descending from peak
      "2019-12-28",  # Holiday period
      "2020-02-15"   # Late in last season
    ))
  }
  
  # =========================================================================
  # PLOT (a): FAN CHARTS
  #
  # For each QA date, read the CSV and plot fan charts for each location.
  # What to look for:
  #   - Red dots (actuals) mostly within the fan
  #   - Fan widens with horizon (h1 narrower than h4)
  #   - No negative y values
  #   - Larger upward spread during peak season (Box-Cox effect)
  # =========================================================================
  
  cat("\n--- Plot (a): Fan charts ---\n")
  
  for (i in seq_along(qa_dates)) {
    d <- qa_dates[i]
    csv_file <- file.path(output_dir,
                          paste0(d, "-", team_model, ".csv"))
    
    if (!file.exists(csv_file)) {
      cat("  Skipping", as.character(d), "— file not found\n")
      next
    }
    
    fc_df <- readr::read_csv(csv_file, show_col_types = FALSE)
    
    plots <- list()
    for (loc in qa_locations) {
      p <- plot_quantile_fan(fc_df, wili_tsibble, d, location = loc) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10),
                       plot.subtitle = ggplot2::element_text(size = 8))
      plots <- c(plots, list(p))
    }
    
    combined <- patchwork::wrap_plots(plots, ncol = 1) +
      patchwork::plot_annotation(
        title = paste0(team_model, " — Fan charts for ", d),
        subtitle = "Blue = forecast fan (1-99%, 10-90%, 25-75%, median) | Red dots = actuals",
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold"))
      )
    print(combined)
    cat("  Plotted:", as.character(d), "\n")
  }
  
  # =========================================================================
  # PLOT (b): SANITY DASHBOARD — median vs actuals over time
  #
  # Reads all CSVs and extracts horizon-1 data for one location.
  # Panel 1: median forecast + 90% PI vs observed actuals
  # Panel 2: 90% PI width over time (should peak during winter)
  # =========================================================================
  
  cat("\n--- Plot (b): Sanity dashboard ---\n")
  
  dashboard_data <- purrr::map_dfr(origin_dates, function(d) {
    csv_file <- file.path(output_dir, paste0(d, "-", team_model, ".csv"))
    if (!file.exists(csv_file)) return(NULL)
    readr::read_csv(csv_file, show_col_types = FALSE) |>
      dplyr::filter(location == dashboard_loc, horizon == 1) |>
      dplyr::select(origin_date, target_end_date, output_type_id, value)
  })
  
  if (nrow(dashboard_data) == 0) {
    cat("  No CSV files found — skipping dashboard\n")
  } else {
    dash <- dashboard_data |>
      dplyr::filter(output_type_id %in% c(0.05, 0.5, 0.95)) |>
      tidyr::pivot_wider(names_from = output_type_id, values_from = value,
                         names_prefix = "q") |>
      dplyr::mutate(pi_width = q0.95 - q0.05)
    
    actuals <- wili_tsibble |>
      tibble::as_tibble() |>
      dplyr::filter(location == dashboard_loc) |>
      dplyr::select(target_end_date, observation)
    
    dash <- dash |> dplyr::left_join(actuals, by = "target_end_date")
    
    p_median <- ggplot2::ggplot(dash, ggplot2::aes(x = target_end_date)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q0.05, ymax = q0.95),
                           fill = "steelblue", alpha = 0.2) +
      ggplot2::geom_line(ggplot2::aes(y = observation),
                         color = "black", linewidth = 0.4) +
      ggplot2::geom_line(ggplot2::aes(y = q0.5),
                         color = "steelblue", linewidth = 0.6) +
      ggplot2::labs(
        title = paste0("Horizon-1 median vs actuals: ", dashboard_loc),
        subtitle = "Blue = median | Black = actual | Shaded = 90% PI",
        x = "Date", y = "wILI %"
      ) +
      ggplot2::theme_minimal()
    
    p_width <- ggplot2::ggplot(dash,
                               ggplot2::aes(x = target_end_date, y = pi_width)) +
      ggplot2::geom_line(color = "darkorange", linewidth = 0.6) +
      ggplot2::labs(
        title = "90% prediction interval width (horizon 1)",
        subtitle = "Should be wider during peak flu season (Jan-Feb)",
        x = "Date", y = "Interval width"
      ) +
      ggplot2::theme_minimal()
    
    dashboard <- p_median / p_width +
      patchwork::plot_annotation(
        title = paste0(team_model, " — Sanity Dashboard"),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold"))
      )
    print(dashboard)
    
    # Numeric summary
    cat(sprintf("  Median forecast range: %.2f to %.2f wILI%%\n",
                min(dash$q0.5, na.rm = TRUE),
                max(dash$q0.5, na.rm = TRUE)))
    cat(sprintf("  Actual range: %.2f to %.2f wILI%%\n",
                min(dash$observation, na.rm = TRUE),
                max(dash$observation, na.rm = TRUE)))
    cat(sprintf("  90%% PI width range: %.2f to %.2f\n",
                min(dash$pi_width, na.rm = TRUE),
                max(dash$pi_width, na.rm = TRUE)))
  }
  
  # =========================================================================
  # PLOT (c): LAMBDA SUMMARY (only for Box-Cox models)
  # =========================================================================
  
  if (!is.null(result$preprocess$lambdas)) {
    cat("\n--- Plot (c): Lambda summary ---\n")
    lambdas <- result$preprocess$lambdas
    
    lambda_df <- tibble::tibble(
      location = names(lambdas),
      lambda   = unname(lambdas)
    ) |>
      dplyr::arrange(lambda) |>
      dplyr::mutate(location = factor(location, levels = location))
    
    p_lambda <- ggplot2::ggplot(lambda_df,
                                ggplot2::aes(x = location, y = lambda)) +
      ggplot2::geom_col(fill = "steelblue", alpha = 0.7) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", lambda)),
                         hjust = -0.1, size = 3) +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste0(team_model, " — Guerrero Box-Cox lambdas"),
        subtitle = paste0("\u03bb \u2248 0 \u2192 log | ",
                          "\u03bb \u2248 1 \u2192 no transform | ",
                          "\u03bb < 0 \u2192 stronger stabilization"),
        x = NULL, y = "Lambda (\u03bb)"
      ) +
      ggplot2::scale_y_continuous(
        limits = c(min(lambda_df$lambda) - 0.05,
                   max(lambda_df$lambda) + 0.1)
      ) +
      ggplot2::theme_minimal()
    
    print(p_lambda)
    cat(sprintf("  Range: [%.4f, %.4f], median: %.4f\n",
                min(lambdas), max(lambdas), median(lambdas)))
  }
  
  cat(sprintf("\n=== Diagnostics complete for %s ===\n", team_model))
  invisible(NULL)
}