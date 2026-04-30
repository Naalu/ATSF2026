###############################################################################
# CONFIG: configs/bma_ensemble.R
#
# Bayesian Model Averaging Ensemble
# Engine:  scripts/forecast_engine.R
#
# Combines 5 component models using adaptive weights derived from
# recent rolling WIS scores. For each origin date:
#   1. Read precomputed forecast CSVs from the 5 component models
#   2. Compute rolling WIS for each over the last k weeks
#   3. Convert to weights via softmax (lower WIS = higher weight)
#   4. Produce weighted quantile combination
#
# REQUIRES: All 5 component models must be run first.
#
# Usage:
#   source("scripts/forecast_engine.R")
#   source("configs/bma_ensemble.R")
#   config$run_loop <- TRUE
#   result <- run_forecast(config)
#
###############################################################################

config <- list(

  # ---- Identity ----
  team_abbr  = "KReger",
  model_abbr = "bma_ensemble",

  # ---- Paths ----
  hub_path = normalizePath(".", mustWork = FALSE),

  # ---- Simulation settings ----
  n_sims = 1000,
  min_train_weeks = 52,

  # ---- Run control ----
  dev_mode  = FALSE,
  run_loop  = FALSE,
  overwrite = FALSE,
  base_seed = 12345,

  # ---- Transform ----
  transform = "none",

  # ---- Model specification ----
  model_fn = NULL,
  simulate_method = "custom",

  # ---- BMA-specific settings ----
  bma = list(
    components  = c("snaive_bc_bs", "ets_log", "arima_bc_bs",
                    "nnetar_bc_bs", "hist_week"),
    k           = 8,       # rolling window (number of prior origin dates)
    temperature = 1.0      # softmax temp (1=standard, Inf=equal weights)
  ),

  # ---- Custom simulation function ----
  custom_simulate_fn = function(train_data, h, n_sims, preprocess_result) {

    bma_cfg    <- config$bma
    hub_path   <- config$hub_path
    team_abbr  <- config$team_abbr
    components <- bma_cfg$components
    k          <- bma_cfg$k
    temp       <- bma_cfg$temperature

    # Determine origin date from training data
    origin_date <- max(train_data$target_end_date)

    # Load actuals for WIS calculation
    actuals_tbl <- tibble::as_tibble(train_data) |>
      dplyr::select(location, target_end_date, actual = observation)

    # Load component forecasts for THIS origin date
    component_forecasts <- list()
    for (comp in components) {
      comp_id  <- paste0(team_abbr, "-", comp)
      csv_path <- file.path(hub_path, "model-output", comp_id,
                            paste0(origin_date, "-", comp_id, ".csv"))
      if (!file.exists(csv_path)) next
      df <- readr::read_csv(csv_path, show_col_types = FALSE) |>
        dplyr::mutate(component = comp)
      component_forecasts[[comp]] <- df
    }

    if (length(component_forecasts) == 0) {
      stop("BMA: no component forecasts found for ", origin_date)
    }

    all_components <- dplyr::bind_rows(component_forecasts)

    # Get all valid origin dates to find the prior window
    tasks_path <- file.path(hub_path, "hub-config", "tasks.json")
    tasks <- jsonlite::fromJSON(tasks_path, simplifyVector = FALSE)
    all_dates <- sort(as.Date(unlist(
      tasks$rounds[[1]]$model_tasks[[1]]$task_ids$origin_date$optional
    )))
    prior_dates <- all_dates[all_dates < origin_date]

    # Compute weights
    if (length(prior_dates) < k) {
      # Burn-in: equal weights
      n_comp <- length(component_forecasts)
      weights <- stats::setNames(rep(1 / n_comp, n_comp),
                                 names(component_forecasts))
    } else {
      eval_dates <- utils::tail(prior_dates, k)

      comp_wis <- sapply(names(component_forecasts), function(comp) {
        comp_id <- paste0(team_abbr, "-", comp)
        wis_values <- c()

        for (d in eval_dates) {
          csv_path <- file.path(hub_path, "model-output", comp_id,
                                paste0(d, "-", comp_id, ".csv"))
          if (!file.exists(csv_path)) next

          fc <- readr::read_csv(csv_path, show_col_types = FALSE) |>
            dplyr::filter(horizon == 1) |>
            dplyr::inner_join(actuals_tbl,
                              by = c("location", "target_end_date"))

          if (nrow(fc) == 0) next

          wis_val <- fc |>
            dplyr::mutate(
              tau = output_type_id,
              pl = dplyr::if_else(
                actual >= value,
                (actual - value) * tau,
                (value - actual) * (1 - tau)
              )
            ) |>
            dplyr::summarise(wis = mean(2 * pl, na.rm = TRUE)) |>
            dplyr::pull(wis)

          wis_values <- c(wis_values, wis_val)
        }

        if (length(wis_values) == 0) return(Inf)
        mean(wis_values, na.rm = TRUE)
      })

      # Softmax: lower WIS = higher weight
      neg_wis <- -comp_wis / temp
      neg_wis <- neg_wis - max(neg_wis)
      exp_wis <- exp(neg_wis)
      weights <- exp_wis / sum(exp_wis)
      names(weights) <- names(component_forecasts)
    }

    # Weighted quantile combination
    weighted_df <- all_components |>
      dplyr::mutate(weight = weights[component]) |>
      dplyr::group_by(location, target_end_date, output_type_id) |>
      dplyr::summarise(
        combined_value = sum(value * weight, na.rm = TRUE),
        .groups = "drop"
      )

    # Return as synthetic simulation paths
    # Replicate each combined value so quantile() returns it exactly
    sims_out <- weighted_df |>
      dplyr::select(location, target_end_date, .sim = combined_value) |>
      dplyr::slice(rep(seq_len(dplyr::n()), each = ceiling(n_sims / 23))) |>
      dplyr::group_by(location, target_end_date) |>
      dplyr::slice_head(n = n_sims) |>
      dplyr::ungroup()

    return(sims_out)
  },

  # ---- Metadata ----
  metadata = list(
    team_name    = "Reger",
    model_name   = "Bayesian Model Averaging Ensemble",
    model_contributors = list(
      list(
        name        = "Karl Reger",
        affiliation = "Northern Arizona University",
        email       = "kcr28@nau.edu"
      )
    ),
    license             = "CC-BY-4.0",
    designated_model    = TRUE,
    data_inputs         = "ILI%",
    methods             = paste0(
      "BMA ensemble of 5 models with adaptive softmax weights based ",
      "on rolling WIS over the preceding 8 origin dates."
    ),
    methods_long        = paste0(
      "Bayesian Model Averaging ensemble combining SNAIVE+Box-Cox, ",
      "ETS+log, ARIMA+Box-Cox, NNETAR+Box-Cox, and historical week ",
      "climatology. Weights are softmax of negative rolling WIS over ",
      "8 prior origin dates (posterior model probabilities). During ",
      "burn-in (first 8 dates), equal weights. Produces weighted ",
      "linear combination of component quantile forecasts. Values ",
      "floored at zero for non-negativity."
    ),
    ensemble_of_models      = TRUE,
    ensemble_of_hub_models  = FALSE
  )
)
