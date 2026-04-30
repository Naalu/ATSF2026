###############################################################################
# FILE: scripts/score_candidates.R
#
# PURPOSE:
#   Phase 2 of FLUCAST_IMPLEMENTATION_PLAN.md. Score every candidate model on
#   the validation-period forecasts (57 origin dates, seasons 1-2 only),
#   summarize per model, and compute pairwise correlation matrices to feed
#   the Phase 3 ensemble-selection step.
#
#   This script combines the work of plan checkpoints 2A (per-forecast
#   scoring) and 2B (per-model summary + correlations) because they share
#   inputs, dependencies, and output directory. Splitting them into two
#   scripts would be a forced separation.
#
# CRITICAL: Validation-data only.
#   The 57 validation origin dates (2015-10-24 through 2017-05-06) are the
#   only forecasts used here. Test-set seasons (3-5, dates from 2017-10-28
#   onward) MUST remain unseen during model selection — using them would
#   contaminate the test-set evaluation later. The script verifies this
#   constraint and aborts if any test-period forecast slips through.
#
# INPUTS:
#   - All 11 candidate models' CSVs in:
#       model-output/KReger-{model_abbr}/*-KReger-{model_abbr}.csv
#     11 models × 57 dates = 627 CSVs total. Each is hub-format
#     (8 columns, 1012 rows: 4 horizons × 11 locations × 23 quantiles).
#   - Ground truth: target-data/time-series.csv
#       Columns: location, target_end_date, observation
#   - Hub config: hub-config/tasks.json (used to enumerate origin dates)
#
# OUTPUTS:
#   analysis/phase2/scores_long.csv
#     One row per (model, location, horizon, origin_date). 11 × 11 × 4 × 57
#     = 27,588 rows. Columns:
#       model, location, origin_date, horizon, target_end_date,
#       observed, point_forecast (median), wis, wis_overprediction,
#       wis_underprediction, wis_dispersion, ae_median,
#       cov_50 (logical), cov_95 (logical)
#
#   analysis/phase2/per_model_summary.csv
#     One row per model. Columns:
#       model, mae, wis, wis_over, wis_under, wis_disp,
#       cov_50_rate, cov_95_rate, n_forecasts
#
#   analysis/phase2/correlation_point_error.csv
#     11x11 symmetric matrix of average per-(location, horizon) Pearson
#     correlations of point-forecast errors (observed - median).
#
#   analysis/phase2/correlation_wis.csv
#     Same shape, correlations of WIS values across (location, horizon, date).
#
#   analysis/phase2/scores_metadata.txt
#     Run timestamp, R + package versions, sanity counts.
#
# HOW TO RUN:
#   # From project root:
#   setwd("/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026")
#   source("scripts/score_candidates.R")
#
#   Runtime: ~1-2 minutes. Bottleneck is reading 627 CSVs; scoring itself
#   is fast (scoringutils is data.table-backed).
#
# REFERENCES:
#   - FLUCAST_IMPLEMENTATION_PLAN.md §3.1 (validation-only constraint),
#     Checkpoints 2A + 2B
#   - scoringutils 2.x API (current as of 2024+):
#     https://epiforecasts.io/scoringutils/
#   - Bracher et al. (2021) "Evaluating epidemic forecasts in an interval
#     format." PLoS Comp Biol 17(2). [WIS definition]
###############################################################################


# ============================================================================
# DEPENDENCIES
# ============================================================================
# Load (don't just check) — we use these throughout.
suppressPackageStartupMessages({
  library(dplyr)        # data wrangling
  library(tidyr)        # pivot_longer, pivot_wider
  library(readr)        # CSV I/O
  library(purrr)        # map functions for the per-model loop
  library(scoringutils) # WIS, coverage, CRPS
})


# ============================================================================
# CONFIGURATION
# ============================================================================

# Path constants. Edit if running on a different machine or layout.
HUB_PATH        <- "/Users/chrisreger/Documents/NAU/Grad/Informatics/INF 599 TS/Project/ATSF2026"
TARGET_DATA     <- file.path(HUB_PATH, "target-data", "time-series.csv")
MODEL_OUTPUT    <- file.path(HUB_PATH, "model-output")
ANALYSIS_DIR    <- file.path(HUB_PATH, "analysis", "phase2")

# Validation cutoff — same constant as scripts/run_validation.R.
# Anything strictly after this date is test-set and must be excluded.
VALIDATION_END_DATE <- as.Date("2017-05-06")

# Candidate models — every directory under model-output/ that's a KReger model.
# We list them explicitly rather than scanning the directory so the script's
# behavior is reproducible even if someone drops a stale folder in there.
CANDIDATES <- c(
  # The original 7 (proposal-era)
  "snaive_bc_bs",
  "ets_log",
  "arima_bc_bs",
  "nnetar_bc_bs",
  "tslm_fourier",
  "stl_arima_bc",
  "hist_week",
  # Phase 1 additions
  "bsts_seasonal",
  "ets_bc",
  "arima_log",
  "nnetar_log"
)

# Hub spec: 23 quantile levels and the central-interval coverage targets.
# We compute coverage for the 50% (q0.25/q0.75) and 95% (q0.025/q0.975)
# intervals, which are the two reported in most FluSight evaluation papers.
COV_50_LOWER <- 0.25
COV_50_UPPER <- 0.75
COV_95_LOWER <- 0.025
COV_95_UPPER <- 0.975

# Reproducibility seed. Most operations here are deterministic given inputs,
# but set anyway in case scoringutils' internals do anything stochastic.
set.seed(20260426L)


# ============================================================================
# DIRECTORY SETUP
# ============================================================================
if (!dir.exists(ANALYSIS_DIR)) {
  dir.create(ANALYSIS_DIR, recursive = TRUE)
}
cat(sprintf("Analysis output dir: %s\n", ANALYSIS_DIR))


# ============================================================================
# LOAD GROUND TRUTH
# ============================================================================
cat("\n--- Loading ground truth ---\n")

# The target-data CSV has columns: location, target_end_date, observation.
# We filter to the validation window's needed dates (Saturdays only in the
# target_end_date range matching our forecasts' h=1..4 weeks ahead).
truth <- readr::read_csv(TARGET_DATA, show_col_types = FALSE) |>
  dplyr::select(location, target_end_date, observed = observation)

cat(sprintf("  %d ground-truth rows (%d locations, %s to %s)\n",
            nrow(truth),
            length(unique(truth$location)),
            format(min(truth$target_end_date)),
            format(max(truth$target_end_date))))


# ============================================================================
# LOAD ALL CANDIDATE FORECASTS
#
# Each CSV has columns: origin_date, target, horizon, location,
# output_type, output_type_id, target_end_date, value (8 columns).
# We read them one model at a time, tag with the model name, and bind.
#
# Critically: filter to validation dates only. Even though Phase 1 produced
# only 57 CSVs per model, we double-check by origin_date in case of any
# stray test-set CSV that snuck in (e.g., from an earlier ad-hoc run).
# ============================================================================
cat("\n--- Loading candidate forecasts ---\n")

#' Load all forecasts for one candidate model from disk.
#'
#' @param model_abbr Character. The model_abbr matching the directory name
#'   under model-output/ (without the team prefix).
#' @return A tibble with columns: model, location, origin_date, horizon,
#'   target_end_date, quantile_level, predicted. One row per quantile level
#'   per (location, horizon, origin_date).
load_model_forecasts <- function(model_abbr) {
  dir <- file.path(MODEL_OUTPUT, paste0("KReger-", model_abbr))
  if (!dir.exists(dir)) {
    stop(sprintf("Model directory not found: %s", dir))
  }

  csvs <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csvs) == 0L) {
    stop(sprintf("No CSVs in %s", dir))
  }

  # Read all and bind. show_col_types = FALSE silences the per-file column
  # spec messages; readr is plenty fast for 57 small CSVs per model.
  all_csvs <- purrr::map_dfr(csvs, readr::read_csv, show_col_types = FALSE)

  all_csvs |>
    # Defensive validation-window filter. Any CSV whose origin_date is
    # past the validation cutoff is silently dropped; we report the count
    # at the end.
    dplyr::filter(origin_date <= VALIDATION_END_DATE) |>
    # Rename to scoringutils 2.x conventions and add the model column.
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

# Load all 11 candidates. Tracking per-model row count so we can assert at
# the end that nothing was silently dropped.
forecast_list <- list()
for (m in CANDIDATES) {
  fc <- load_model_forecasts(m)
  forecast_list[[m]] <- fc
  cat(sprintf("  %-15s  %5d rows  (%d unique origin_dates)\n",
              m, nrow(fc), length(unique(fc$origin_date))))
}

# Combine into one big long-format tibble. ~310k rows total
# (11 models × 57 dates × 11 locations × 4 horizons × 23 quantiles).
forecasts <- dplyr::bind_rows(forecast_list)
cat(sprintf("\n  Total forecast rows: %d\n", nrow(forecasts)))

# Sanity check: every model should have exactly 57 origin dates.
date_counts <- forecasts |>
  dplyr::group_by(model) |>
  dplyr::summarise(n_dates = dplyr::n_distinct(origin_date), .groups = "drop")
if (any(date_counts$n_dates != 57L)) {
  stop("Validation-period origin date count mismatch:\n",
       paste(capture.output(print(date_counts)), collapse = "\n"))
}


# ============================================================================
# JOIN GROUND TRUTH
# ============================================================================
cat("\n--- Joining ground truth ---\n")

# Inner join on (location, target_end_date). Any forecast whose target
# date is missing from the truth file gets dropped — but the validation
# window (max h=4 ahead from 2017-05-06 = 2017-06-03) is well inside
# the truth window (extends to 2020-08-29), so this should be a no-op.
forecasts_scored <- forecasts |>
  dplyr::inner_join(truth, by = c("location", "target_end_date"))

dropped <- nrow(forecasts) - nrow(forecasts_scored)
if (dropped > 0L) {
  warning(sprintf(
    "Dropped %d forecast rows during truth join (missing observed values). ",
    dropped),
    "Investigate before trusting downstream analysis."
  )
}
cat(sprintf("  %d rows joined to ground truth (%d dropped)\n",
            nrow(forecasts_scored), dropped))


# ============================================================================
# SCORE WITH scoringutils
#
# scoringutils 2.x workflow:
#   1. as_forecast_quantile(forecast_unit = ...): wrap the data.frame in a
#      forecast object after declaring which columns identify a single
#      "forecast" (here: model + location + origin_date + horizon).
#   2. score(): compute WIS components, ae_median, coverage flags, etc.
#      Output is one row per forecast unit (NOT per quantile-level row of input).
# ============================================================================
cat("\n--- Computing scores via scoringutils ---\n")

# Tell scoringutils which columns identify a single forecast distribution.
# Everything else (predicted, observed, quantile_level) is the data of the
# forecast itself.
forecast_obj <- forecasts_scored |>
  scoringutils::as_forecast_quantile(
    forecast_unit = c("model", "location", "origin_date",
                      "target_end_date", "horizon")
  )

# score() returns a data.table with one row per forecast unit and columns:
#   wis, overprediction, underprediction, dispersion, bias,
#   interval_coverage_50, interval_coverage_90, ae_median
# (Plus any forecast_unit columns we declared above.)
#
# Note scoringutils' default coverage levels are 50 and 90, NOT 50 and 95.
# To get 95% coverage, we have to compute it manually below from the raw
# quantile data — scoringutils' get_coverage() exposes per-quantile coverage
# but it's cleaner to just compute the 95% interval flag ourselves.
scores <- scoringutils::score(forecast_obj) |>
  tibble::as_tibble()

cat(sprintf("  Scored %d forecast units\n", nrow(scores)))


# ============================================================================
# COMPUTE 95% COVERAGE FLAG MANUALLY
#
# scoringutils gives us interval_coverage_50 and interval_coverage_90 by
# default. We want 95% (the FluSight standard) and 50%. Compute 95% from
# the q0.025 and q0.975 quantiles directly, then merge back.
# ============================================================================
cov_95 <- forecasts_scored |>
  dplyr::filter(quantile_level %in% c(COV_95_LOWER, COV_95_UPPER)) |>
  tidyr::pivot_wider(
    id_cols     = c(model, location, origin_date, target_end_date, horizon, observed),
    names_from  = quantile_level,
    values_from = predicted,
    names_prefix = "q"
  ) |>
  dplyr::mutate(
    cov_95 = (observed >= .data[[paste0("q", COV_95_LOWER)]]) &
             (observed <= .data[[paste0("q", COV_95_UPPER)]])
  ) |>
  dplyr::select(model, location, origin_date, target_end_date, horizon, cov_95)


# ============================================================================
# BUILD THE LONG SCORES TABLE (artifact 1)
# ============================================================================
cat("\n--- Building scores_long.csv ---\n")

# Pull point forecasts (medians) for the table.
medians <- forecasts_scored |>
  dplyr::filter(abs(quantile_level - 0.5) < 1e-9) |>
  dplyr::select(model, location, origin_date, target_end_date, horizon,
                point_forecast = predicted)

scores_long <- scores |>
  dplyr::transmute(
    model, location, origin_date, target_end_date, horizon,
    wis,
    wis_overprediction  = overprediction,
    wis_underprediction = underprediction,
    wis_dispersion      = dispersion,
    ae_median,
    # scoringutils calls it interval_coverage_50; rename for clarity.
    cov_50 = interval_coverage_50
  ) |>
  dplyr::left_join(medians,
                   by = c("model", "location", "origin_date",
                          "target_end_date", "horizon")) |>
  dplyr::left_join(cov_95,
                   by = c("model", "location", "origin_date",
                          "target_end_date", "horizon")) |>
  # Also pull observed back in for human inspection.
  dplyr::left_join(
    forecasts_scored |>
      dplyr::distinct(model, location, origin_date, target_end_date,
                      horizon, observed),
    by = c("model", "location", "origin_date", "target_end_date", "horizon")
  ) |>
  # Tidy column order
  dplyr::select(model, location, origin_date, horizon, target_end_date,
                observed, point_forecast,
                wis, wis_overprediction, wis_underprediction, wis_dispersion,
                ae_median,
                cov_50, cov_95)

# Expected row count: 11 models × 11 locations × 4 horizons × 57 dates = 27,588
expected_rows <- length(CANDIDATES) * 11L * 4L * 57L
if (nrow(scores_long) != expected_rows) {
  warning(sprintf("scores_long has %d rows; expected %d. Investigate.",
                  nrow(scores_long), expected_rows))
}
cat(sprintf("  %d rows (expected %d)\n", nrow(scores_long), expected_rows))

readr::write_csv(scores_long, file.path(ANALYSIS_DIR, "scores_long.csv"))
cat(sprintf("  Wrote: %s\n", file.path(ANALYSIS_DIR, "scores_long.csv")))


# ============================================================================
# PER-MODEL SUMMARY (artifact 2)
# ============================================================================
cat("\n--- Building per_model_summary.csv ---\n")

# Aggregate across all (location, horizon, date) to one row per model.
per_model <- scores_long |>
  dplyr::group_by(model) |>
  dplyr::summarise(
    mae          = mean(ae_median, na.rm = TRUE),
    wis          = mean(wis, na.rm = TRUE),
    wis_over     = mean(wis_overprediction, na.rm = TRUE),
    wis_under    = mean(wis_underprediction, na.rm = TRUE),
    wis_disp     = mean(wis_dispersion, na.rm = TRUE),
    cov_50_rate  = mean(cov_50, na.rm = TRUE),
    cov_95_rate  = mean(cov_95, na.rm = TRUE),
    n_forecasts  = dplyr::n(),
    .groups      = "drop"
  ) |>
  # Sort by WIS ascending so the best model is at the top.
  dplyr::arrange(wis)

readr::write_csv(per_model, file.path(ANALYSIS_DIR, "per_model_summary.csv"))
cat(sprintf("  Wrote: %s\n", file.path(ANALYSIS_DIR, "per_model_summary.csv")))

# Echo to console — this is the main slide-worthy table.
cat("\n  Per-model summary (sorted by WIS):\n")
print(per_model, n = Inf)


# ============================================================================
# CORRELATION MATRICES (artifact 3 + 4)
#
# Two correlation matrices are computed:
#
#   (1) Point-error correlations.
#       For each (location, horizon) pair, compute Pearson correlations
#       across the 57 origin dates of (observed - point_forecast). Then
#       average the resulting correlation matrices across all 44
#       (location, horizon) cells to get a single 11x11 correlation matrix.
#       This matches the proposal-era methodology.
#
#   (2) WIS-based correlations.
#       Same procedure but on the WIS values themselves: do models tend to
#       have a bad day on the same dates?
#
# Both matrices are symmetric with 1.0 on the diagonal.
# ============================================================================
cat("\n--- Computing correlation matrices ---\n")

#' Compute average per-(location, horizon) correlation matrix for a metric.
#'
#' @param data Tibble with columns: model, location, horizon, origin_date,
#'   and the metric column.
#' @param metric_col Character name of the metric column to correlate
#'   (e.g., "point_error" or "wis").
#' @return A model x model correlation matrix (numeric matrix).
compute_avg_correlation <- function(data, metric_col) {
  # Pivot to wide format: one column per model, one row per
  # (location, horizon, origin_date) triple.
  wide <- data |>
    dplyr::select(model, location, horizon, origin_date,
                  value = dplyr::all_of(metric_col)) |>
    tidyr::pivot_wider(names_from = model, values_from = value)

  # For each (location, horizon) group of 57 rows, compute the model x model
  # correlation matrix. Then average across groups.
  groups <- wide |>
    dplyr::group_by(location, horizon) |>
    dplyr::group_split()

  cor_matrices <- purrr::map(groups, function(g) {
    # g is one (location, horizon) cell. Extract just the model columns.
    mat <- g |> dplyr::select(dplyr::all_of(CANDIDATES)) |> as.matrix()
    # use = "complete.obs" drops any row with any NA; with our data this
    # should never trigger but is defensive against unexpected gaps.
    stats::cor(mat, use = "complete.obs", method = "pearson")
  })

  # Element-wise mean across the 44 matrices.
  Reduce(`+`, cor_matrices) / length(cor_matrices)
}

# Augment scores_long with point_error column.
with_point_error <- scores_long |>
  dplyr::mutate(point_error = observed - point_forecast)

cor_point <- compute_avg_correlation(with_point_error, "point_error")
cor_wis   <- compute_avg_correlation(scores_long,      "wis")

# Round to 3 decimals for human readability — full precision is lost
# noise anyway.
cor_point <- round(cor_point, 3)
cor_wis   <- round(cor_wis,   3)

# Write as CSV with row + column names = model abbreviations.
write_cor_csv <- function(cor_mat, path) {
  out <- as.data.frame(cor_mat) |>
    tibble::rownames_to_column("model")
  readr::write_csv(out, path)
}

write_cor_csv(cor_point, file.path(ANALYSIS_DIR, "correlation_point_error.csv"))
write_cor_csv(cor_wis,   file.path(ANALYSIS_DIR, "correlation_wis.csv"))

cat("\n  Point-error correlation matrix (rounded to 3 dp):\n")
print(cor_point)
cat("\n  WIS correlation matrix (rounded to 3 dp):\n")
print(cor_wis)

# Print the off-diagonal extremes — these are the candidates for the
# twin-pair rule and the most-decorrelated pairs that anchor diversity.
off_diag <- function(m) {
  m[upper.tri(m)] |>
    setNames(apply(which(upper.tri(m), arr.ind = TRUE), 1, function(idx) {
      paste(rownames(m)[idx[1]], colnames(m)[idx[2]], sep = " <-> ")
    }))
}

cat("\n  Top 5 highest point-error correlations (twin candidates):\n")
print(sort(off_diag(cor_point), decreasing = TRUE)[1:5])

cat("\n  Top 5 lowest point-error correlations (most diverse pairs):\n")
print(sort(off_diag(cor_point), decreasing = FALSE)[1:5])


# ============================================================================
# METADATA
# ============================================================================
metadata_path <- file.path(ANALYSIS_DIR, "scores_metadata.txt")
sink(metadata_path)
cat(sprintf("Phase 2 candidate scoring run\n"))
cat(sprintf("Generated: %s\n", format(Sys.time())))
cat(sprintf("R version: %s\n", R.version.string))
cat(sprintf("scoringutils: %s\n", packageVersion("scoringutils")))
cat(sprintf("dplyr: %s\n", packageVersion("dplyr")))
cat(sprintf("\nCandidate models (%d):\n", length(CANDIDATES)))
cat(paste("  -", CANDIDATES, collapse = "\n"))
cat(sprintf("\n\nValidation cutoff: %s\n", VALIDATION_END_DATE))
cat(sprintf("Validation origin dates: %d\n", 57L))
cat(sprintf("Total forecast rows scored: %d\n", nrow(scores_long)))
sink()
cat(sprintf("\n  Wrote: %s\n", metadata_path))


cat("\n========================================\n")
cat("  Phase 2 scoring complete\n")
cat("========================================\n")
cat(sprintf("  4 artifacts in %s\n", ANALYSIS_DIR))
cat("    - scores_long.csv\n")
cat("    - per_model_summary.csv\n")
cat("    - correlation_point_error.csv\n")
cat("    - correlation_wis.csv\n")
cat("    - scores_metadata.txt\n")
cat("========================================\n")
