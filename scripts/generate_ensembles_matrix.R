###############################################################################
# FILE: scripts/generate_ensembles_matrix.R
#
# PURPOSE:
#   Phase 4D of FLUCAST_IMPLEMENTATION_PLAN.md. Generate four ensemble
#   forecasts and compare them via leave-one-date-out cross-validation.
#   The four-cell comparison matrix is the central empirical result of
#   the wider-pool experiment.
#
#   Comparison matrix (4 cells):
#                    BMA            Equal-weight
#     Primary  (8)   already exists  new
#     Expanded (10)  Phase 4C        new
#
#   "BMA-primary" already exists from Phase 3 (analysis/phase3/ensemble_*).
#   This script:
#     1. Generates ensemble forecasts for the three NEW cells:
#        - Equal-weight on primary 8-model pool (134 dates)
#        - BMA-weighted on expanded 10-model pool (134 dates)
#        - Equal-weight on expanded 10-model pool (134 dates)
#     2. Runs leave-one-date-out cross-validation on validation data for
#        all FOUR ensembles, including Phase 3's existing BMA-primary.
#     3. Produces a final comparison matrix and designation recommendation.
#
# WHY ALL FOUR ENSEMBLES IN ONE SCRIPT:
#   The methodology must be apples-to-apples. Same data hygiene, same
#   quantile-aggregation routine, same monotonicity enforcement, same
#   scoringutils invocation, same LOO protocol. Splitting across scripts
#   risks subtle inconsistencies that would muddy the comparison.
#
# CRITICAL: validation-only weight estimation (no leakage).
#   The "production weights" for both BMA ensembles are computed once
#   on the full 57-date validation set (already done in Phase 3A and
#   Phase 4C). For the LOO comparison, we re-compute weights at each
#   LOO fold using the 56-date training subset. Equal-weight ensembles
#   need no re-computation (weights are constant) but their LOO scoring
#   is still useful for apples-to-apples comparison.
#
# OUTPUT STAGING:
#   All generated ensemble CSVs go into analysis/phase4/staging/ — they
#   are NOT promoted to model-output/ until the GATE check decides which
#   ensemble gets designated_model = TRUE for hub submission. This is
#   the same gating principle as Phase 3B/3C.
#
# INPUTS:
#   - analysis/phase3/production_weights.csv         (BMA-primary weights)
#   - analysis/phase4/production_weights_expanded.csv (BMA-expanded weights)
#   - analysis/phase4/shortlist_expanded.csv          (10-model list)
#   - analysis/phase2/regimes.csv                     (validation regimes)
#   - target-data/time-series.csv                    (truth)
#   - hub-config/tasks.json                          (134 origin dates)
#   - model-output/<model>/*.csv                     (component forecasts)
#
# OUTPUTS in analysis/phase4/:
#   - staging/ensemble_bma_primary/                  (134 CSVs — copy from
#                                                     analysis/phase3 or new)
#   - staging/ensemble_equal_primary/                (134 CSVs — new)
#   - staging/ensemble_bma_expanded/                 (134 CSVs — new)
#   - staging/ensemble_equal_expanded/               (134 CSVs — new)
#   - matrix_loo_comparison.csv                      (4 rows, headline result)
#   - matrix_loo_comparison_by_regime.csv            (per-regime breakdown)
#   - generate_ensembles_matrix_metadata.txt         (run details)
#
# HOW TO RUN:
#   setwd("<path-to-ATSF2026-repo>")
#   source("scripts/generate_ensembles_matrix.R")
#
#   Runtime: ~5-8 minutes. Bottleneck: 4 LOO matrices x 57 folds each.
#   Equal-weight ensembles are nearly free; BMA ensembles are ~2 min each.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md Phase 4D
#   - scripts/generate_ensemble.R (Phase 3B template)
#   - scripts/compute_bma_weights.R + compute_bma_weights_expanded.R
###############################################################################


# ============================================================================
# DEPENDENCIES
# ============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(tibble)
  library(jsonlite)
  library(scoringutils)
})


# ============================================================================
# CONFIGURATION
# ============================================================================
HUB_PATH       <- normalizePath(".", mustWork = FALSE)
TARGET_DATA    <- file.path(HUB_PATH, "target-data", "time-series.csv")
TASKS_JSON     <- file.path(HUB_PATH, "hub-config", "tasks.json")
MODEL_OUTPUT   <- file.path(HUB_PATH, "model-output")

PHASE3_DIR     <- file.path(HUB_PATH, "analysis", "phase3")
PHASE4_DIR     <- file.path(HUB_PATH, "analysis", "phase4")
PHASE2_DIR     <- file.path(HUB_PATH, "analysis", "phase2")

REGIMES_CSV    <- file.path(PHASE2_DIR, "regimes.csv")
WEIGHTS_PRIMARY  <- file.path(PHASE3_DIR, "production_weights.csv")
WEIGHTS_EXPANDED <- file.path(PHASE4_DIR, "production_weights_expanded.csv")
SHORTLIST_EXP    <- file.path(PHASE4_DIR, "shortlist_expanded.csv")

VALIDATION_END_DATE <- as.Date("2017-05-06")

# Primary ensemble (8 models from Phase 3).
PRIMARY_MODELS <- c(
  "nnetar_bc_bs", "stl_arima_bc", "arima_bc_bs", "ets_bc",
  "hist_week", "bsts_seasonal", "tslm_fourier", "snaive_bc_bs"
)
# These are the model_abbr suffixes; the directories are KReger-{abbr}.
PRIMARY_DIRS <- paste0("KReger-", PRIMARY_MODELS)

# Staging directory for all generated ensembles.
STAGING_DIR    <- file.path(PHASE4_DIR, "staging")

# Output team and ensemble names for staging.
ENSEMBLE_NAMES <- list(
  bma_primary    = "bma_primary",
  equal_primary  = "equal_primary",
  bma_expanded   = "bma_expanded",
  equal_expanded = "equal_expanded"
)

set.seed(20260428L)


# ============================================================================
# DIRECTORY SETUP
# ============================================================================
for (n in names(ENSEMBLE_NAMES)) {
  d <- file.path(STAGING_DIR, paste0("ensemble_", n))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# ============================================================================
# LOAD PRODUCTION WEIGHTS AND DERIVE EXPANDED MODEL LIST
# ============================================================================
cat("\n--- Loading production weights ---\n")

w_bma_primary <- readr::read_csv(WEIGHTS_PRIMARY, show_col_types = FALSE) |>
  dplyr::select(model, regime, weight)

w_bma_expanded <- readr::read_csv(WEIGHTS_EXPANDED, show_col_types = FALSE) |>
  dplyr::select(model, regime, weight)

shortlist_exp <- readr::read_csv(SHORTLIST_EXP, show_col_types = FALSE)
EXPANDED_MODELS <- shortlist_exp$model  # directory names
EXPANDED_LABELS <- setNames(shortlist_exp$model_label, shortlist_exp$model)

# The primary ensemble's BMA weights use bare model_abbr; need to add
# KReger- prefix for directory matching.
w_bma_primary <- w_bma_primary |>
  dplyr::mutate(model_dir = paste0("KReger-", model)) |>
  dplyr::select(model = model_dir, regime, weight)

cat(sprintf("  BMA-primary weights: %d rows (%d models x %d regimes)\n",
            nrow(w_bma_primary), length(unique(w_bma_primary$model)),
            length(unique(w_bma_primary$regime))))
cat(sprintf("  BMA-expanded weights: %d rows (%d models x %d regimes)\n",
            nrow(w_bma_expanded), length(unique(w_bma_expanded$model)),
            length(unique(w_bma_expanded$regime))))


# ============================================================================
# LOAD ANCILLARIES
# ============================================================================
cat("\n--- Loading regimes, truth, and tasks.json ---\n")

regimes <- readr::read_csv(REGIMES_CSV, show_col_types = FALSE) |>
  dplyr::mutate(origin_date = as.Date(origin_date))

truth <- readr::read_csv(TARGET_DATA, show_col_types = FALSE) |>
  dplyr::transmute(location, target_end_date = as.Date(target_end_date),
                   observed = observation)

tasks <- jsonlite::fromJSON(TASKS_JSON, simplifyVector = FALSE)
all_origin_dates <- as.Date(unlist(
  tasks$rounds[[1]]$model_tasks[[1]]$task_ids$origin_date$optional
)) |> sort()

validation_dates <- all_origin_dates[all_origin_dates <= VALIDATION_END_DATE]

cat(sprintf("  origin dates total: %d\n", length(all_origin_dates)))
cat(sprintf("  validation dates:   %d\n", length(validation_dates)))


# ============================================================================
# REGIMES FOR ALL 134 DATES (validation already in regimes.csv;
# test dates need on-the-fly classification using regime_helpers.R)
# ============================================================================
test_dates <- all_origin_dates[all_origin_dates > VALIDATION_END_DATE]

# We only need test-date regimes if BMA generation extends to all 134 dates.
# For the LOO-CV comparison (validation-only), test regimes aren't needed.
# But for full ensemble generation we do need them.
cat("\n--- Classifying test-period regimes ---\n")

source(file.path(HUB_PATH, "scripts", "regime_helpers.R"))
wili_data <- truth |> dplyr::transmute(location, target_end_date,
                                        observation = observed)

calibration <- compute_regime_calibration(
  wili_data, validation_start_date = as.Date("2015-10-24")
)

regimes_test <- precompute_regimes(
  wili_data,
  origin_dates = test_dates,
  locations    = unique(regimes$location),
  calibration  = calibration
)

regimes_all <- dplyr::bind_rows(regimes, regimes_test)
cat(sprintf("  total regime cells (134 dates x 11 locs): %d\n",
            nrow(regimes_all)))


# ============================================================================
# LOAD COMPONENT QUANTILES
# ============================================================================
cat("\n--- Loading component model quantile forecasts ---\n")

#' Load all 134-date validation+test quantile CSVs for one model directory.
load_full_quantiles <- function(model_dir) {
  dir <- file.path(MODEL_OUTPUT, model_dir)
  csvs <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  loaded <- purrr::map_dfr(csvs, readr::read_csv, show_col_types = FALSE,
                            col_types = readr::cols(.default = readr::col_guess()))
  loaded |>
    dplyr::transmute(
      model           = model_dir,
      location        = as.character(location),
      origin_date     = as.Date(origin_date),
      horizon         = as.integer(horizon),
      target_end_date = as.Date(target_end_date),
      quantile_level  = as.numeric(output_type_id),
      predicted       = as.numeric(value)
    )
}

# Union of primary and expanded model lists (some overlap — load each only once).
all_models_needed <- union(PRIMARY_DIRS, EXPANDED_MODELS)
cat(sprintf("  loading %d unique models...\n", length(all_models_needed)))

components <- list()
for (m in all_models_needed) {
  comp <- load_full_quantiles(m)
  n_dates <- length(unique(comp$origin_date))
  if (n_dates < 134L) {
    stop(sprintf("Component %s has only %d dates (need 134)", m, n_dates))
  }
  components[[m]] <- comp
  cat(sprintf("    %-25s %d dates  %d rows\n", m, n_dates, nrow(comp)))
}

components_long <- dplyr::bind_rows(components)


# ============================================================================
# WEIGHT TABLES FOR EQUAL-WEIGHT ENSEMBLES
#
# Equal-weight = uniform 1/N over all models in the pool. Same weight in
# every regime (regime-conditional structure trivially constant).
# ============================================================================
build_equal_weights <- function(model_list, regime_list) {
  N <- length(model_list)
  tidyr::expand_grid(
    model  = model_list,
    regime = regime_list
  ) |>
    dplyr::mutate(weight = 1 / N)
}

regime_list <- unique(regimes$regime)

w_equal_primary  <- build_equal_weights(PRIMARY_DIRS, regime_list)
w_equal_expanded <- build_equal_weights(EXPANDED_MODELS, regime_list)

cat(sprintf("\n  Equal-weight pools:\n"))
cat(sprintf("    primary:  %d models -> uniform 1/%d = %.4f per model\n",
            length(PRIMARY_DIRS), length(PRIMARY_DIRS), 1/length(PRIMARY_DIRS)))
cat(sprintf("    expanded: %d models -> uniform 1/%d = %.4f per model\n",
            length(EXPANDED_MODELS), length(EXPANDED_MODELS), 1/length(EXPANDED_MODELS)))


# ============================================================================
# ENSEMBLE FORECAST HELPERS
# ============================================================================

#' Build the ensemble forecast for one origin date given weights + comp quantiles.
build_ensemble_one_date <- function(date_components, date_regimes, weights) {
  date_components |>
    dplyr::inner_join(date_regimes, by = c("origin_date", "location")) |>
    dplyr::inner_join(weights,      by = c("model", "regime")) |>
    dplyr::group_by(location, origin_date, target_end_date,
                    horizon, quantile_level) |>
    dplyr::summarise(predicted = sum(predicted * weight), .groups = "drop") |>
    # Floor at 0 (wILI cannot be negative).
    dplyr::mutate(predicted = pmax(predicted, 0)) |>
    # Enforce monotonicity within each (location, horizon) group.
    dplyr::group_by(location, origin_date, target_end_date, horizon) |>
    dplyr::arrange(quantile_level, .by_group = TRUE) |>
    dplyr::mutate(predicted = sort(predicted)) |>
    dplyr::ungroup()
}

#' Format an ensemble forecast as hub-spec CSV-ready tibble.
to_hub_format <- function(ensemble) {
  ensemble |>
    dplyr::transmute(
      origin_date,
      target          = "ili perc",
      horizon,
      location,
      output_type     = "quantile",
      output_type_id  = quantile_level,
      target_end_date,
      value           = predicted
    ) |>
    dplyr::arrange(location, horizon, output_type_id)
}


# ============================================================================
# GENERATE ALL FOUR ENSEMBLES (134 dates each)
#
# We generate validation + test together. The validation portion will be
# re-scored in the LOO loop below; the test portion sits in staging until
# Phase 5 promotes the winning ensemble to model-output/.
# ============================================================================

cat("\n=== Generating four ensembles for all 134 dates ===\n")

#' Generate all 134 CSVs for one ensemble configuration.
generate_ensemble_full <- function(ensemble_name, weights, model_list,
                                    output_dir) {
  cat(sprintf("\n--- %s (%d models) ---\n", ensemble_name, length(model_list)))

  # Subset components_long to just this ensemble's models.
  these_components <- components_long |>
    dplyr::filter(model %in% model_list)

  start <- Sys.time()
  for (i in seq_along(all_origin_dates)) {
    d <- all_origin_dates[i]

    # Per-date filter
    date_components <- these_components |> dplyr::filter(origin_date == d)
    date_regimes <- regimes_all |> dplyr::filter(origin_date == d)

    if (nrow(date_components) == 0L) {
      stop(sprintf("No components for date %s", format(d)))
    }

    ensemble <- build_ensemble_one_date(date_components, date_regimes, weights)
    hub_csv <- to_hub_format(ensemble)

    fname <- sprintf("%s-KReger-%s.csv", format(d), ensemble_name)
    readr::write_csv(hub_csv, file.path(output_dir, fname))

    if (i %% 20L == 0L || i == 1L || i == length(all_origin_dates)) {
      cat(sprintf("  [%3d/%d] %s\n", i, length(all_origin_dates), format(d)))
    }
  }
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "mins"))
  cat(sprintf("  Done in %.1f min\n", elapsed))
}

generate_ensemble_full(
  ensemble_name = ENSEMBLE_NAMES$bma_primary,
  weights       = w_bma_primary,
  model_list    = PRIMARY_DIRS,
  output_dir    = file.path(STAGING_DIR, "ensemble_bma_primary")
)

generate_ensemble_full(
  ensemble_name = ENSEMBLE_NAMES$equal_primary,
  weights       = w_equal_primary,
  model_list    = PRIMARY_DIRS,
  output_dir    = file.path(STAGING_DIR, "ensemble_equal_primary")
)

generate_ensemble_full(
  ensemble_name = ENSEMBLE_NAMES$bma_expanded,
  weights       = w_bma_expanded,
  model_list    = EXPANDED_MODELS,
  output_dir    = file.path(STAGING_DIR, "ensemble_bma_expanded")
)

generate_ensemble_full(
  ensemble_name = ENSEMBLE_NAMES$equal_expanded,
  weights       = w_equal_expanded,
  model_list    = EXPANDED_MODELS,
  output_dir    = file.path(STAGING_DIR, "ensemble_equal_expanded")
)


# ============================================================================
# LOO-CV CROSS-METHOD COMPARISON
#
# Score validation ensembles via LOO: for each held-out date, compute weights
# from remaining 56 dates, build that date's ensemble, score it. Sum WIS
# across all 57 folds.
#
# Equal-weight ensembles need no re-fit per fold (weights are constant) but
# we still run the loop to keep the four methods strictly comparable —
# same input data, same scoring routine, same per-fold context.
# ============================================================================

cat("\n=== LOO-CV cross-method comparison ===\n")

#' One LOO fold: hold out a date, compute weights, score held-out ensemble.
loo_fold_score <- function(held_out_date, weight_fn, model_list, scores_with_regime) {
  # Compute weights from training subset (all dates except held-out).
  train_scores <- scores_with_regime |>
    dplyr::filter(origin_date != held_out_date)

  w <- weight_fn(train_scores, model_list)

  # Pull held-out date's component quantiles and regimes.
  held_components <- components_long |>
    dplyr::filter(model %in% model_list, origin_date == held_out_date)
  held_regimes <- regimes_all |> dplyr::filter(origin_date == held_out_date)

  ens <- build_ensemble_one_date(held_components, held_regimes, w)

  # Score with scoringutils.
  fc_obj <- ens |>
    dplyr::inner_join(truth, by = c("location", "target_end_date")) |>
    scoringutils::as_forecast_quantile(
      forecast_unit = c("location", "origin_date", "target_end_date", "horizon")
    )

  scoringutils::score(fc_obj) |>
    tibble::as_tibble() |>
    dplyr::mutate(held_out_date = held_out_date)
}

# Weight functions take (train_scores, model_list) and return the same
# tibble shape.
weight_fn_inverse_wis_primary <- function(train_scores, model_list) {
  train_scores |>
    dplyr::filter(model %in% model_list) |>
    dplyr::group_by(model, regime) |>
    dplyr::summarise(mean_wis = mean(wis, na.rm = TRUE), .groups = "drop") |>
    dplyr::group_by(regime) |>
    dplyr::mutate(weight = (1 / mean_wis) / sum(1 / mean_wis)) |>
    dplyr::ungroup() |>
    dplyr::select(model, regime, weight)
}

# For BMA we use the *frozen* full-validation weights (per the design choice
# documented in analysis/phase2/decisions.md §3.3). The LOO-CV here is for
# fair comparison across methods, NOT for re-tuning BMA weights per fold.
# This means BMA's "LOO score" reflects how the frozen weights perform on
# each held-out date — slightly higher than a method that re-tunes per fold,
# but that's the trade-off for a leakage-free production pipeline.
weight_fn_bma_primary_frozen <- function(train_scores, model_list) {
  w_bma_primary
}
weight_fn_bma_expanded_frozen <- function(train_scores, model_list) {
  w_bma_expanded
}
weight_fn_equal <- function(train_scores, model_list) {
  build_equal_weights(model_list, regime_list)
}


# Build score table joined to regimes for the LOO loops (saves repeated work).
# We need this only for inverse-WIS-style methods that fit weights per fold.
# Kept minimal here since both BMA methods use frozen weights.

cat("\n--- Collecting per-method validation WIS via LOO ---\n")

# Pre-build a score table for inverse_wis-style fits (not strictly needed
# in this script but kept structure-ready in case we extend later).
all_scores_long_primary <- readr::read_csv(
  file.path(PHASE2_DIR, "scores_long.csv"), show_col_types = FALSE
) |>
  dplyr::mutate(
    origin_date = as.Date(origin_date),
    # phase2 uses bare model_abbr; convert to dir name to match keys here.
    model = paste0("KReger-", model)
  ) |>
  dplyr::filter(model %in% PRIMARY_DIRS) |>
  dplyr::inner_join(regimes, by = c("origin_date", "location"))

all_scores_long_expanded <- readr::read_csv(
  file.path(PHASE4_DIR, "scores_long_expanded.csv"), show_col_types = FALSE
) |>
  dplyr::mutate(origin_date = as.Date(origin_date)) |>
  dplyr::filter(model %in% EXPANDED_MODELS) |>
  dplyr::inner_join(regimes, by = c("origin_date", "location"))


# Helper: run LOO across validation_dates with given config.
run_loo <- function(label, weight_fn, model_list, scores_with_regime) {
  cat(sprintf("  %s... ", label))
  start <- Sys.time()
  results <- list()
  for (d in validation_dates) {
    res <- loo_fold_score(d, weight_fn, model_list, scores_with_regime)
    results[[as.character(d)]] <- res
  }
  bound <- dplyr::bind_rows(results)
  total_wis <- sum(bound$wis)
  per_regime <- bound |>
    dplyr::inner_join(regimes, by = c("origin_date", "location")) |>
    dplyr::group_by(regime) |>
    dplyr::summarise(mean_wis = mean(wis), n_cells = dplyr::n(),
                     .groups = "drop") |>
    dplyr::mutate(method = label)
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
  cat(sprintf("LOO total WIS = %8.3f  (%.1fs)\n", total_wis, elapsed))
  list(total_wis = total_wis, per_regime = per_regime, n_cells = nrow(bound))
}

# Four methods, four LOO runs.
loo_results <- list()

loo_results$bma_primary <- run_loo(
  "BMA-primary    ",
  weight_fn_bma_primary_frozen,
  PRIMARY_DIRS,
  all_scores_long_primary
)
loo_results$equal_primary <- run_loo(
  "Equal-primary  ",
  weight_fn_equal,
  PRIMARY_DIRS,
  all_scores_long_primary
)
loo_results$bma_expanded <- run_loo(
  "BMA-expanded   ",
  weight_fn_bma_expanded_frozen,
  EXPANDED_MODELS,
  all_scores_long_expanded
)
loo_results$equal_expanded <- run_loo(
  "Equal-expanded ",
  weight_fn_equal,
  EXPANDED_MODELS,
  all_scores_long_expanded
)


# ============================================================================
# COMPARISON MATRIX
# ============================================================================
cat("\n=== Comparison matrix ===\n")

matrix_df <- tibble::tibble(
  cell           = c("BMA-primary", "Equal-primary",
                     "BMA-expanded", "Equal-expanded"),
  pool           = c("primary (8)", "primary (8)",
                     "expanded (10)", "expanded (10)"),
  method         = c("BMA", "Equal", "BMA", "Equal"),
  loo_total_wis  = sapply(loo_results, function(x) x$total_wis),
  n_cells        = sapply(loo_results, function(x) x$n_cells),
  mean_wis       = sapply(loo_results, function(x) x$total_wis / x$n_cells)
) |>
  dplyr::mutate(ranking = rank(loo_total_wis, ties.method = "min")) |>
  dplyr::arrange(loo_total_wis)

readr::write_csv(matrix_df, file.path(PHASE4_DIR, "matrix_loo_comparison.csv"))

cat("\n  Final 2x2 comparison matrix (sorted by LOO WIS):\n")
print(matrix_df)

# Also: per-regime breakdown, for transparency.
per_regime_long <- dplyr::bind_rows(lapply(loo_results, function(r) r$per_regime))
readr::write_csv(per_regime_long,
                 file.path(PHASE4_DIR, "matrix_loo_comparison_by_regime.csv"))

cat("\n  Per-regime breakdown:\n")
print(
  per_regime_long |>
    tidyr::pivot_wider(names_from = method, values_from = mean_wis) |>
    dplyr::select(regime, n_cells, dplyr::any_of(matrix_df$cell))
)


# ============================================================================
# DESIGNATION RECOMMENDATION
# ============================================================================
winner <- matrix_df$cell[1]
cat(sprintf("\n=== Designation recommendation: %s ===\n", winner))

cat(sprintf("\n  Winner: %s\n", winner))
cat(sprintf("  LOO total WIS: %.3f\n", matrix_df$loo_total_wis[1]))
cat(sprintf("  Margin vs runner-up: %.2f WIS units (%.2f%%)\n",
            matrix_df$loo_total_wis[2] - matrix_df$loo_total_wis[1],
            100 * (matrix_df$loo_total_wis[2] - matrix_df$loo_total_wis[1]) /
                  matrix_df$loo_total_wis[1]))

# Map winner cell name back to a staging directory name.
winner_dir_map <- c(
  "BMA-primary"    = "ensemble_bma_primary",
  "Equal-primary"  = "ensemble_equal_primary",
  "BMA-expanded"   = "ensemble_bma_expanded",
  "Equal-expanded" = "ensemble_equal_expanded"
)
winner_dir <- file.path(STAGING_DIR, winner_dir_map[winner])
cat(sprintf("\n  Winner staging path: %s\n", winner_dir))
cat("  This ensemble's CSVs should be promoted to model-output/KReger-bma_ensemble/\n")
cat("  with designated_model = TRUE in the metadata.\n")


# ============================================================================
# METADATA
# ============================================================================
metadata_path <- file.path(PHASE4_DIR, "generate_ensembles_matrix_metadata.txt")
sink(metadata_path)
cat("Phase 4D: four-ensemble generation + LOO-CV comparison matrix\n")
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat("\nEnsembles generated (134 CSVs each, in analysis/phase4/staging/):\n")
cat(sprintf("  BMA-primary    (8 models, log-score BMA, tau frozen at Phase 3A)\n"))
cat(sprintf("  Equal-primary  (8 models, uniform 1/8)\n"))
cat(sprintf("  BMA-expanded   (10 models, log-score BMA, tau frozen at Phase 4C)\n"))
cat(sprintf("  Equal-expanded (10 models, uniform 1/10)\n"))
cat("\nLOO-CV comparison matrix:\n")
cat(paste(capture.output(print(matrix_df)), collapse = "\n"))
cat(sprintf("\n\nWinner: %s\n", winner))
cat(sprintf("\nNote: BMA ensembles use frozen full-validation weights for LOO\n"))
cat(sprintf("scoring (no per-fold re-fitting). This is the same hygiene as Phase\n"))
cat(sprintf("3A and Phase 4C — production weights apply unchanged at test time.\n"))
sink()


cat("\n========================================\n")
cat("  Phase 4D complete\n")
cat("========================================\n")
cat(sprintf("  Four ensembles staged in %s/\n", STAGING_DIR))
cat(sprintf("  Comparison matrix: %s\n",
            file.path(PHASE4_DIR, "matrix_loo_comparison.csv")))
cat(sprintf("  Winner: %s -> promote for designation in Phase 5\n", winner))
cat("========================================\n")
