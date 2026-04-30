###############################################################################
# FILE: scripts/generate_ensemble.R
#
# PURPOSE:
#   Phase 3B+3C of FLUCAST_IMPLEMENTATION_PLAN.md. Generate the BMA ensemble
#   forecasts for all 134 origin dates by linear-pooling the 8 candidate
#   models' quantile values using the frozen per-(model, regime) weights
#   from Phase 3A. Then validate the ensemble's WIS on the validation set
#   only — test-set scoring deferred to Phase 4.
#
#   Output is written to a *staging* directory:
#     analysis/phase3/ensemble_validation/  (57 CSVs, validation period)
#     analysis/phase3/ensemble_test/        (77 CSVs, test period)
#
#   The staging structure separates "candidate ensemble output" from
#   "production hub submission." Only after Phase 3C diagnostics pass do
#   we promote these CSVs into model-output/KReger-bma_ensemble/ via a
#   separate script (scripts/promote_ensemble.R).
#
# CONSTRUCTION ALGORITHM:
#   For each (origin_date, location, target_end_date, horizon, quantile_level):
#     1. Look up the regime for (origin_date, location) from regimes.csv
#     2. Look up the 8 weights for that regime from production_weights.csv
#     3. Linear-pool: ensemble_q = sum_k w_k(regime) * component_k_q
#     4. Floor at 0 (wILI cannot be negative)
#
#   After computing all quantile values for one (origin_date, location,
#   horizon), enforce monotonicity by sorting the 23 quantile values
#   ascending. Linear pooling of monotonic distributions is mathematically
#   guaranteed to produce monotonic results, but floating-point arithmetic
#   can introduce micro-violations. The post-hoc sort handles those without
#   any statistical consequence (we're sorting values that should already
#   be sorted; sorts only fix numerical noise).
#
# VALIDATION-ONLY DIAGNOSTIC:
#   After generating all 134 CSVs, score the 57 validation-period
#   ensemble forecasts using scoringutils. Compare to:
#     - The best individual component (nnetar_bc_bs at WIS 0.284)
#     - Per-regime breakdowns
#     - The Phase 3A LOO-CV WIS estimate (638 / 2508 = 0.254 per cell)
#
#   This diagnostic is the GATE for Phase 3B output. If the ensemble
#   does not beat its best component on validation, we have a bug or
#   design flaw and should NOT promote to production.
#
# INPUTS:
#   - analysis/phase3/production_weights.csv  (8 models x 3 regimes)
#   - analysis/phase2/regimes.csv             (627 cells, validation only)
#   - hub-config/tasks.json                   (134 origin dates)
#   - target-data/time-series.csv             (truth for diagnostic)
#   - model-output/KReger-{model}/*.csv       (8 candidates' forecasts)
#
# OUTPUTS in analysis/phase3/:
#   - ensemble_validation/                    (57 CSVs, hub format)
#   - ensemble_test/                          (77 CSVs, hub format)
#   - ensemble_validation_diagnostics.csv     (per-regime WIS comparison)
#   - generate_ensemble_metadata.txt          (run details)
#
# HOW TO RUN:
#   setwd("/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026")
#   source("scripts/generate_ensemble.R")
#
#   Runtime: ~2-3 minutes. Bottleneck is the 134-date generation loop;
#   each date involves a join+summarise across 8 models x 11 locations
#   x 4 horizons x 23 quantiles = 8,096 rows.
#
# CRITICAL: regime classification for test-period dates.
#   The regimes.csv file contains regime labels for the 57 validation
#   dates only. For the 77 test dates, we need to classify on-the-fly
#   using the same calibration that was fit on pre-2015 data. This
#   script sources scripts/regime_helpers.R, recomputes the calibration
#   from pre-2015 wILI, and classifies each test date as it's processed.
#   No leakage: the calibration uses no data from 2015 onward.
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md Phase 3B + 3C
#   - analysis/phase2/decisions.md §3.1-3.3
#   - scripts/compute_bma_weights.R (produces production_weights.csv)
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
  library(jsonlite)     # for tasks.json parsing
  library(scoringutils) # for the validation diagnostic
})


# ============================================================================
# CONFIGURATION
# ============================================================================
HUB_PATH      <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"
MODEL_OUTPUT  <- file.path(HUB_PATH, "model-output")
TARGET_DATA   <- file.path(HUB_PATH, "target-data", "time-series.csv")
TASKS_JSON    <- file.path(HUB_PATH, "hub-config", "tasks.json")
WEIGHTS_CSV   <- file.path(HUB_PATH, "analysis", "phase3", "production_weights.csv")
REGIMES_CSV   <- file.path(HUB_PATH, "analysis", "phase2", "regimes.csv")

PHASE3_DIR    <- file.path(HUB_PATH, "analysis", "phase3")
VAL_DIR       <- file.path(PHASE3_DIR, "ensemble_validation")
TEST_DIR      <- file.path(PHASE3_DIR, "ensemble_test")

# Validation/test cutoff (matches every other script in this project).
VALIDATION_END_DATE <- as.Date("2017-05-06")

# 8-model shortlist (must match scripts/compute_bma_weights.R).
ENSEMBLE_MODELS <- c(
  "nnetar_bc_bs",
  "stl_arima_bc",
  "arima_bc_bs",
  "ets_bc",
  "hist_week",
  "bsts_seasonal",
  "tslm_fourier",
  "snaive_bc_bs"
)

# Output naming for the staged ensemble. Final hub model_abbr.
ENSEMBLE_TEAM  <- "KReger"
ENSEMBLE_NAME  <- "bma_ensemble"


# ============================================================================
# DIRECTORY SETUP
# ============================================================================
for (d in c(VAL_DIR, TEST_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}


# ============================================================================
# LOAD WEIGHTS, REGIMES (VALIDATION), AND DEPENDENCIES
# ============================================================================
cat("\n--- Loading weights and validation regimes ---\n")

weights <- readr::read_csv(WEIGHTS_CSV, show_col_types = FALSE) |>
  dplyr::select(model, regime, weight)

# Sanity: weights should cover all 8 models in 3 regimes (24 rows total),
# and sum to 1 within each regime.
stopifnot(
  "Weights table not 24 rows" = nrow(weights) == 24L,
  "Weights cover wrong model set" =
    setequal(unique(weights$model), ENSEMBLE_MODELS)
)
weight_sums <- weights |>
  dplyr::group_by(regime) |>
  dplyr::summarise(s = sum(weight), .groups = "drop")
stopifnot(
  "Per-regime weight sums not 1" = all(abs(weight_sums$s - 1) < 1e-9)
)

cat(sprintf("  Loaded weights: %d rows, %d regimes\n",
            nrow(weights), length(unique(weights$regime))))

# Validation-period regimes already on disk.
regimes_validation <- readr::read_csv(REGIMES_CSV, show_col_types = FALSE) |>
  dplyr::mutate(origin_date = as.Date(origin_date))

cat(sprintf("  Validation regimes: %d cells\n", nrow(regimes_validation)))


# ============================================================================
# REGIME CLASSIFIER FOR TEST DATES
#
# Test-period (post-2017-05-06) regimes are not in regimes.csv. We need
# to classify each test (origin_date, location) cell using the same
# calibration that was fit on pre-2015 data. Source the helper, fit
# calibration once, classify test cells once.
#
# Critical: the calibration uses validation_start_date = 2015-10-24, so
# nothing past that date enters calibration. This matches Phase 2's
# methodology. No leakage.
# ============================================================================
cat("\n--- Classifying test-period regimes ---\n")

source(file.path(HUB_PATH, "scripts", "regime_helpers.R"))

# Load full wili history for classification.
wili_data <- readr::read_csv(TARGET_DATA, show_col_types = FALSE) |>
  dplyr::transmute(location, target_end_date = as.Date(target_end_date),
                   observation)

# Recompute calibration (fast, ~30 sec).
calibration <- compute_regime_calibration(
  wili_data, validation_start_date = as.Date("2015-10-24")
)

cat(sprintf("  Calibration recomputed (slope_z_threshold=%.3f)\n",
            calibration$slope_z_threshold))


# ============================================================================
# LOAD ALL ORIGIN DATES FROM tasks.json
# ============================================================================
cat("\n--- Loading 134 origin dates from tasks.json ---\n")

tasks <- jsonlite::fromJSON(TASKS_JSON, simplifyVector = FALSE)
all_origin_dates <- as.Date(unlist(
  tasks$rounds[[1]]$model_tasks[[1]]$task_ids$origin_date$optional
))
all_origin_dates <- sort(all_origin_dates)

stopifnot(
  "Expected 134 origin dates" = length(all_origin_dates) == 134L
)
cat(sprintf("  %d origin dates (%s to %s)\n",
            length(all_origin_dates),
            format(min(all_origin_dates)),
            format(max(all_origin_dates))))

validation_dates <- all_origin_dates[all_origin_dates <= VALIDATION_END_DATE]
test_dates       <- all_origin_dates[all_origin_dates >  VALIDATION_END_DATE]

cat(sprintf("  Validation: %d dates  Test: %d dates\n",
            length(validation_dates), length(test_dates)))


# ============================================================================
# CLASSIFY TEST DATE REGIMES
# ============================================================================
all_locations <- unique(regimes_validation$location)

regimes_test <- precompute_regimes(
  wili_data,
  origin_dates = test_dates,
  locations    = all_locations,
  calibration  = calibration
)

# Combine validation + test regimes into one lookup. Both are needed for
# the generation loop.
regimes_all <- dplyr::bind_rows(regimes_validation, regimes_test)
stopifnot(
  "Expected 134*11 = 1474 regime rows" = nrow(regimes_all) == 134L * 11L
)
cat(sprintf("  Test regimes: %d cells\n", nrow(regimes_test)))


# ============================================================================
# LOAD COMPONENT QUANTILE FORECASTS
#
# All 8 candidates' forecasts for all 134 dates. The Phase 1+ work
# produced CSVs only for validation (57 dates) for nnetar_log/etc., but
# the Phase-3B-eligible 8 models all have full 134-date coverage (they're
# the original 7 plus bsts and ets_bc, all of which we generated for the
# full hub before twin-pair pruning).
#
# BUT WAIT: ets_bc and bsts_seasonal have validation-only output. We need
# to run them on the test set too before this script can complete. If
# they're missing test CSVs we'll fail loudly here.
# ============================================================================
cat("\n--- Loading 8 components' forecasts (all 134 dates) ---\n")

#' Load all CSVs for one model and concatenate.
load_component <- function(model_abbr) {
  dir <- file.path(MODEL_OUTPUT, paste0("KReger-", model_abbr))
  csvs <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)

  loaded <- purrr::map_dfr(csvs, readr::read_csv, show_col_types = FALSE)

  loaded |>
    dplyr::transmute(
      model           = model_abbr,
      location        = location,
      origin_date     = as.Date(origin_date),
      horizon         = as.integer(horizon),
      target_end_date = as.Date(target_end_date),
      quantile_level  = as.numeric(output_type_id),
      predicted       = value
    )
}

components <- list()
for (m in ENSEMBLE_MODELS) {
  comp <- load_component(m)
  n_dates <- length(unique(comp$origin_date))
  components[[m]] <- comp
  cat(sprintf("  %-15s %5d dates  %d rows\n", m, n_dates, nrow(comp)))
  if (n_dates < 134L) {
    stop(sprintf(
      "Component %s has only %d dates; need all 134. Run validation/test ",
      m, n_dates),
      "generation for this model first."
    )
  }
}

components_long <- dplyr::bind_rows(components)
cat(sprintf("\n  Total component rows: %d\n", nrow(components_long)))


# ============================================================================
# ENSEMBLE GENERATION
# ============================================================================

#' Build the ensemble forecast for one origin date and write CSV.
#'
#' @param origin_date Date.
#' @param output_dir Character. Directory for the CSV.
#' @return Tibble of the ensemble forecast (long quantile format).
generate_one_date <- function(origin_date, output_dir) {

  # Pull this date's component forecasts and regimes.
  date_components <- components_long |>
    dplyr::filter(origin_date == !!origin_date)

  date_regimes <- regimes_all |>
    dplyr::filter(origin_date == !!origin_date)

  if (nrow(date_components) == 0L) {
    stop(sprintf("No component forecasts found for origin_date = %s",
                 format(origin_date)))
  }

  # Build the ensemble: join components -> regimes -> weights, then
  # weighted-sum the quantile values within each
  # (location, target_end_date, horizon, quantile_level) cell.
  ensemble <- date_components |>
    dplyr::inner_join(date_regimes, by = c("origin_date", "location")) |>
    dplyr::inner_join(weights,      by = c("model", "regime")) |>
    dplyr::group_by(location, origin_date, target_end_date,
                    horizon, quantile_level) |>
    dplyr::summarise(
      predicted = sum(predicted * weight),
      .groups   = "drop"
    ) |>
    # Floor at 0 (wILI cannot be negative).
    dplyr::mutate(predicted = pmax(predicted, 0))

  # Enforce monotonicity within each (location, horizon) group. Linear
  # pooling of monotonic distributions is mathematically monotonic, but
  # floating-point can introduce micro-violations.
  ensemble <- ensemble |>
    dplyr::group_by(location, origin_date, target_end_date, horizon) |>
    dplyr::arrange(quantile_level, .by_group = TRUE) |>
    dplyr::mutate(predicted = sort(predicted)) |>
    dplyr::ungroup()

  # Format to hub spec. Hub columns: origin_date, target, horizon, location,
  # output_type, output_type_id, target_end_date, value.
  hub_format <- ensemble |>
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

  # Write CSV. Filename: YYYY-MM-DD-Team-model.csv per hub convention.
  file_path <- file.path(
    output_dir,
    sprintf("%s-%s-%s.csv", format(origin_date), ENSEMBLE_TEAM, ENSEMBLE_NAME)
  )
  readr::write_csv(hub_format, file_path)

  hub_format
}


# ============================================================================
# MAIN GENERATION LOOP
# ============================================================================
cat("\n--- Generating validation ensemble forecasts (57 dates) ---\n")
val_start <- Sys.time()
for (i in seq_along(validation_dates)) {
  d <- validation_dates[i]
  if (i %% 10L == 0L || i == 1L || i == length(validation_dates)) {
    cat(sprintf("  [%2d/%d] %s\n", i, length(validation_dates), format(d)))
  }
  generate_one_date(d, VAL_DIR)
}
val_elapsed <- as.numeric(difftime(Sys.time(), val_start, units = "mins"))
cat(sprintf("  Done in %.1f min\n", val_elapsed))

cat("\n--- Generating test ensemble forecasts (77 dates) ---\n")
test_start <- Sys.time()
for (i in seq_along(test_dates)) {
  d <- test_dates[i]
  if (i %% 10L == 0L || i == 1L || i == length(test_dates)) {
    cat(sprintf("  [%2d/%d] %s\n", i, length(test_dates), format(d)))
  }
  generate_one_date(d, TEST_DIR)
}
test_elapsed <- as.numeric(difftime(Sys.time(), test_start, units = "mins"))
cat(sprintf("  Done in %.1f min\n", test_elapsed))


# ============================================================================
# VALIDATION DIAGNOSTIC
#
# Score the 57 validation-period ensemble CSVs and compare to:
#   - The best individual component (per Phase 2's per_model_summary)
#   - Per-regime breakdowns
#   - The LOO-CV estimate from Phase 3A (638 total / 2508 cells = 0.254)
# ============================================================================
cat("\n--- Validation diagnostic ---\n")

# Read all 57 validation ensemble CSVs back in.
val_csvs <- list.files(VAL_DIR, pattern = "\\.csv$", full.names = TRUE)
ensemble_val <- purrr::map_dfr(val_csvs, readr::read_csv, show_col_types = FALSE) |>
  dplyr::transmute(
    location, origin_date = as.Date(origin_date),
    horizon = as.integer(horizon),
    target_end_date = as.Date(target_end_date),
    quantile_level  = as.numeric(output_type_id),
    predicted = value
  )

truth <- wili_data |>
  dplyr::transmute(location, target_end_date, observed = observation)

# Score with scoringutils.
fc_obj <- ensemble_val |>
  dplyr::inner_join(truth, by = c("location", "target_end_date")) |>
  scoringutils::as_forecast_quantile(
    forecast_unit = c("location", "origin_date", "target_end_date", "horizon")
  )

ens_scores <- scoringutils::score(fc_obj) |>
  tibble::as_tibble()

# Add regime to scored rows.
ens_scores <- ens_scores |>
  dplyr::inner_join(regimes_validation, by = c("origin_date", "location"))

overall_wis <- mean(ens_scores$wis)
overall_total_wis <- sum(ens_scores$wis)

cat(sprintf("  Ensemble overall mean WIS: %.4f  (total: %.2f over %d cells)\n",
            overall_wis, overall_total_wis, nrow(ens_scores)))
cat(sprintf("  Phase 3A LOO-CV estimate:  %.4f  (total: 638.4)\n", 638.4 / 2508))
cat(sprintf("  Best component (nnetar):   0.2840\n"))

# Per-regime breakdown.
per_regime <- ens_scores |>
  dplyr::group_by(regime) |>
  dplyr::summarise(
    mean_wis = mean(wis),
    n_cells  = dplyr::n(),
    .groups  = "drop"
  ) |>
  dplyr::arrange(regime)

cat("\n  Per-regime breakdown:\n")
print(per_regime)


# ============================================================================
# COMPARE PER-MODEL: ENSEMBLE VS BEST COMPONENT BY REGIME
#
# For each regime, also pull the best component's mean WIS in that regime.
# This gives us the per-regime "ensemble lift" diagnostic.
# ============================================================================
cat("\n  Per-regime comparison vs best component:\n")
scores_long_path <- file.path(HUB_PATH, "analysis", "phase2", "scores_long.csv")
scores_long <- readr::read_csv(scores_long_path, show_col_types = FALSE) |>
  dplyr::mutate(origin_date = as.Date(origin_date)) |>
  dplyr::filter(model %in% ENSEMBLE_MODELS) |>
  dplyr::inner_join(regimes_validation, by = c("origin_date", "location"))

per_model_per_regime <- scores_long |>
  dplyr::group_by(model, regime) |>
  dplyr::summarise(mean_wis = mean(wis), .groups = "drop")

best_per_regime <- per_model_per_regime |>
  dplyr::group_by(regime) |>
  dplyr::slice_min(mean_wis, n = 1) |>
  dplyr::ungroup() |>
  dplyr::transmute(regime, best_component = model, best_wis = mean_wis)

comparison <- per_regime |>
  dplyr::inner_join(best_per_regime, by = "regime") |>
  dplyr::mutate(
    ensemble_lift_pct = round(100 * (best_wis - mean_wis) / best_wis, 2)
  )

print(comparison)

# Verdict: is the ensemble beating its best component?
overall_best_wis <- min(per_model_per_regime$mean_wis)
if (overall_wis < overall_best_wis) {
  cat(sprintf("\n  GATE PASSED: Ensemble (%.4f) beats best single-regime component (%.4f).\n",
              overall_wis, overall_best_wis))
  cat("  Cleared to promote staged CSVs to model-output/KReger-bma_ensemble/.\n")
} else {
  cat(sprintf("\n  GATE FAILED: Ensemble WIS %.4f >= best component %.4f.\n",
              overall_wis, overall_best_wis))
  cat("  DO NOT promote. Investigate before continuing.\n")
}

# Save the diagnostic.
diag_path <- file.path(PHASE3_DIR, "ensemble_validation_diagnostics.csv")
readr::write_csv(comparison, diag_path)
cat(sprintf("\n  Wrote: %s\n", diag_path))


# ============================================================================
# METADATA
# ============================================================================
metadata_path <- file.path(PHASE3_DIR, "generate_ensemble_metadata.txt")
sink(metadata_path)
cat("Phase 3B+3C: Ensemble forecast generation + validation diagnostic\n")
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("\nEnsemble models (%d):\n", length(ENSEMBLE_MODELS)))
cat(paste("  -", ENSEMBLE_MODELS, collapse = "\n"))
cat(sprintf("\n\nDates generated:\n"))
cat(sprintf("  Validation: %d (%s to %s)\n",
            length(validation_dates),
            format(min(validation_dates)),
            format(max(validation_dates))))
cat(sprintf("  Test:       %d (%s to %s)\n",
            length(test_dates),
            format(min(test_dates)),
            format(max(test_dates))))
cat(sprintf("\nOutput directories:\n"))
cat(sprintf("  Validation: %s\n", VAL_DIR))
cat(sprintf("  Test:       %s\n", TEST_DIR))
cat(sprintf("\nValidation diagnostic:\n"))
cat(sprintf("  Mean WIS:        %.4f\n", overall_wis))
cat(sprintf("  Total WIS:       %.2f\n", overall_total_wis))
cat(sprintf("  Best component:  %.4f\n", overall_best_wis))
cat(paste(capture.output(print(comparison)), collapse = "\n"))
sink()


cat("\n========================================\n")
cat("  Phase 3B+3C complete\n")
cat("========================================\n")
cat(sprintf("  Ensemble CSVs in staging:\n"))
cat(sprintf("    %s/  (validation, 57)\n", VAL_DIR))
cat(sprintf("    %s/  (test, 77)\n",       TEST_DIR))
cat(sprintf("\n  Diagnostic outputs:\n"))
cat(sprintf("    %s\n", diag_path))
cat(sprintf("    %s\n", metadata_path))
cat("========================================\n")
