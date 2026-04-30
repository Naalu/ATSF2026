# =============================================================================
# scripts/score_test_set.R
#
# Phase 5D — Test-set evaluation
#
# Scores all four ensemble configurations (BMA-primary, BMA-expanded,
# Equal-primary, Equal-expanded) on the test window (Seasons 3-5, 77 origin
# dates from 2017-10-28 to 2020-02-29) using frozen validation-period weights.
#
# Compares against:
#   - Validation 4-cell matrix from Phase 4D (in the report's Table 1)
#   - Per-regime breakdown comparable to the report's Table 2
#
# Outputs (saved to analysis/phase4/):
#   - test_evaluation_summary.csv:    4-cell matrix on test
#   - test_evaluation_by_regime.csv:  per-regime mean WIS
#   - test_evaluation_report.md:      readable markdown summary
#
# Run from project root:
#   Rscript scripts/score_test_set.R
# Or in RStudio: setwd() to project root, then source("scripts/score_test_set.R")
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(scoringutils)
})

# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

# Test window: Seasons 3-5
TEST_START <- as.Date("2017-10-28")
TEST_END   <- as.Date("2020-02-29")

# Validation reference numbers from Phase 4D / report Table 1
VALIDATION_REF <- tibble::tribble(
  ~ensemble,        ~val_total_wis, ~val_mean_wis,
  "bma_expanded",            614.0,         0.245,
  "bma_primary",             630.3,         0.251,
  "equal_primary",           732.3,         0.292,
  "equal_expanded",          750.6,         0.299
)

# Where each ensemble's forecasts live
ENSEMBLE_DIRS <- list(
  bma_primary    = "analysis/phase4/staging/ensemble_bma_primary",
  bma_expanded   = "analysis/phase4/staging/ensemble_bma_expanded",
  equal_primary  = "analysis/phase4/staging/ensemble_equal_primary",
  equal_expanded = "analysis/phase4/staging/ensemble_equal_expanded"
)

# Output paths
OUT_DIR  <- "analysis/phase4"
OUT_SUMMARY <- file.path(OUT_DIR, "test_evaluation_summary.csv")
OUT_REGIME  <- file.path(OUT_DIR, "test_evaluation_by_regime.csv")
OUT_REPORT  <- file.path(OUT_DIR, "test_evaluation_report.md")

# -----------------------------------------------------------------------------
# 1. Load truth from oracle-output.csv
# -----------------------------------------------------------------------------

cat("\n=== Loading truth data ===\n")

truth <- readr::read_csv("target-data/oracle-output.csv",
                         show_col_types = FALSE) |>
  dplyr::filter(target == "ili perc") |>
  dplyr::transmute(
    location,
    target_end_date = as.Date(target_end_date),
    observed        = oracle_value
  ) |>
  # In case oracle has duplicate rows for the same (location, target_end_date)
  dplyr::distinct(location, target_end_date, .keep_all = TRUE)

cat(sprintf("  Loaded %d (location, target_end_date) truth rows\n", nrow(truth)))
cat(sprintf("  Date range: %s to %s\n",
            min(truth$target_end_date), max(truth$target_end_date)))

# -----------------------------------------------------------------------------
# 2. Function: load test-period forecasts from one ensemble's staging dir
# -----------------------------------------------------------------------------

#' Load all CSVs from a staging directory and filter to test period
#'
#' @param staging_dir Path to staging directory containing per-origin-date CSVs
#' @return Long-format tibble with one row per (origin_date, location, horizon,
#'   quantile_level) — restricted to test-period origin dates.
load_test_forecasts <- function(staging_dir) {
  csv_files <- list.files(staging_dir, pattern = "\\.csv$", full.names = TRUE)
  
  # Filter filenames to test period before reading (avoid loading validation CSVs)
  origin_dates <- as.Date(stringr::str_extract(basename(csv_files),
                                                "^\\d{4}-\\d{2}-\\d{2}"))
  test_files <- csv_files[origin_dates >= TEST_START & origin_dates <= TEST_END]
  
  if (length(test_files) == 0L) {
    stop("No test-period CSVs found in ", staging_dir)
  }
  
  # Read and combine
  purrr::map_dfr(test_files, function(f) {
    origin_date <- as.Date(stringr::str_extract(basename(f),
                                                 "^\\d{4}-\\d{2}-\\d{2}"))
    readr::read_csv(f, show_col_types = FALSE) |>
      dplyr::mutate(origin_date = origin_date)
  }) |>
    dplyr::filter(output_type == "quantile") |>
    dplyr::transmute(
      origin_date,
      location,
      horizon,
      target_end_date = as.Date(target_end_date),
      quantile_level  = as.numeric(output_type_id),
      predicted       = value
    )
}

# -----------------------------------------------------------------------------
# 3. Function: score one ensemble against truth
# -----------------------------------------------------------------------------

#' Score forecasts against truth using scoringutils
#'
#' @param forecasts Long-format forecast tibble from load_test_forecasts()
#' @param ensemble_label Character label for tagging output
#' @return Tibble with one row per (origin_date, location, horizon) plus scores
score_ensemble <- function(forecasts, ensemble_label) {
  # Join with truth by (location, target_end_date)
  combined <- forecasts |>
    dplyr::left_join(truth, by = c("location", "target_end_date"))
  
  # Report and drop rows missing truth (forecasts past available data)
  n_total   <- nrow(combined)
  n_missing <- sum(is.na(combined$observed))
  if (n_missing > 0L) {
    cat(sprintf("  %s: dropped %d/%d rows (%.1f%%) with missing truth\n",
                ensemble_label, n_missing, n_total,
                100 * n_missing / n_total))
    combined <- dplyr::filter(combined, !is.na(observed))
  }
  
  # Convert to scoringutils forecast object and score
  forecast_obj <- combined |>
    scoringutils::as_forecast_quantile(
      forecast_unit  = c("origin_date", "location", "horizon", "target_end_date"),
      observed       = "observed",
      predicted      = "predicted",
      quantile_level = "quantile_level"
    )
  
  scores <- scoringutils::score(forecast_obj)
  
  scores |> dplyr::mutate(ensemble = ensemble_label)
}

# -----------------------------------------------------------------------------
# 4. Apply to all 4 ensembles
# -----------------------------------------------------------------------------

cat("\n=== Loading and scoring 4 ensembles on test window ===\n")

all_scores <- purrr::imap_dfr(ENSEMBLE_DIRS, function(dir, label) {
  cat(sprintf("\n%s:\n", label))
  fc <- load_test_forecasts(dir)
  cat(sprintf("  Loaded %d (origin_date, location, horizon, quantile) rows\n",
              nrow(fc)))
  cat(sprintf("  Origin dates: %d unique\n",
              length(unique(fc$origin_date))))
  score_ensemble(fc, label)
})

# -----------------------------------------------------------------------------
# 5. Summary: 4-cell comparison matrix
# -----------------------------------------------------------------------------

summary_table <- all_scores |>
  dplyr::group_by(ensemble) |>
  dplyr::summarize(
    test_total_wis        = sum(wis),
    test_mean_wis         = mean(wis),
    mean_overprediction   = mean(overprediction),
    mean_underprediction  = mean(underprediction),
    mean_dispersion       = mean(dispersion),
    n_cells               = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::left_join(VALIDATION_REF, by = "ensemble") |>
  dplyr::mutate(
    delta_mean_wis        = test_mean_wis - val_mean_wis,
    pct_change            = 100 * (test_mean_wis - val_mean_wis) / val_mean_wis
  ) |>
  dplyr::arrange(test_mean_wis)

cat("\n=== Test-set 4-cell comparison ===\n\n")
summary_table |>
  dplyr::select(ensemble, n_cells, test_total_wis, test_mean_wis,
                val_mean_wis, delta_mean_wis, pct_change) |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3))) |>
  print(n = Inf)

readr::write_csv(summary_table, OUT_SUMMARY)
cat(sprintf("\nSaved: %s\n", OUT_SUMMARY))

# -----------------------------------------------------------------------------
# 6. Per-regime breakdown
# -----------------------------------------------------------------------------

cat("\n=== Per-regime breakdown ===\n")

# Load regime classifications
regimes_path <- "analysis/phase2/regimes.csv"
regimes <- readr::read_csv(regimes_path, show_col_types = FALSE)

# Detect schema. Standard expected columns: origin_date, location, regime.
expected_cols <- c("origin_date", "location", "regime")
missing_cols  <- setdiff(expected_cols, names(regimes))
if (length(missing_cols) > 0L) {
  cat(sprintf("WARNING: regimes.csv missing columns: %s\n",
              paste(missing_cols, collapse = ", ")))
  cat(sprintf("Available columns: %s\n",
              paste(names(regimes), collapse = ", ")))
  cat("Skipping per-regime breakdown — fix the regimes file or rename columns.\n")
} else {
  regimes_test <- regimes |>
    dplyr::mutate(origin_date = as.Date(origin_date)) |>
    dplyr::filter(origin_date >= TEST_START, origin_date <= TEST_END) |>
    dplyr::select(dplyr::all_of(expected_cols))
  
  cat(sprintf("Loaded %d regime classifications for test period\n",
              nrow(regimes_test)))
  cat("Test regime distribution:\n")
  print(table(regimes_test$regime))
  
  by_regime <- all_scores |>
    dplyr::left_join(regimes_test,
                     by = c("origin_date", "location")) |>
    dplyr::filter(!is.na(regime)) |>
    dplyr::group_by(ensemble, regime) |>
    dplyr::summarize(
      mean_wis = mean(wis),
      n_cells  = dplyr::n(),
      .groups  = "drop"
    )
  
  # Pivot wider for the report's Table 2 format
  regime_table <- by_regime |>
    dplyr::select(ensemble, regime, mean_wis) |>
    tidyr::pivot_wider(names_from = ensemble, values_from = mean_wis) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 3)))
  
  cat("\nPer-regime mean WIS (test set):\n\n")
  print(regime_table)
  
  readr::write_csv(by_regime, OUT_REGIME)
  cat(sprintf("\nSaved: %s\n", OUT_REGIME))
}

# -----------------------------------------------------------------------------
# 7. Markdown summary
# -----------------------------------------------------------------------------

cat("\n=== Writing markdown report ===\n")

md_lines <- c(
  "# Test-Set Evaluation: 4-Cell Comparison",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  sprintf("Test window: %s to %s (Seasons 3-5)", TEST_START, TEST_END),
  sprintf("Total origin dates: %d",
          length(unique(all_scores$origin_date))),
  "",
  "## Validation vs Test: Mean WIS",
  "",
  "| Ensemble | Validation mean WIS | Test mean WIS | Δ | % change |",
  "|---|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(summary_table))) {
  row <- summary_table[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %s | %.3f | %.3f | %+.3f | %+.1f%% |",
    row$ensemble,
    row$val_mean_wis,
    row$test_mean_wis,
    row$delta_mean_wis,
    row$pct_change
  ))
}

md_lines <- c(md_lines,
              "",
              "## Test-set 4-cell matrix",
              "",
              "| Cell | Pool | Method | Test total WIS | Test mean WIS |",
              "|---|---|---|---:|---:|")

cell_meta <- tibble::tribble(
  ~ensemble,       ~pool,           ~method,
  "bma_expanded",   "expanded (10)", "BMA",
  "bma_primary",    "primary (8)",   "BMA",
  "equal_primary",  "primary (8)",   "Equal",
  "equal_expanded", "expanded (10)", "Equal"
)

ranked <- summary_table |>
  dplyr::left_join(cell_meta, by = "ensemble") |>
  dplyr::arrange(test_mean_wis)

for (i in seq_len(nrow(ranked))) {
  row <- ranked[i, ]
  md_lines <- c(md_lines, sprintf(
    "| %s | %s | %s | %.1f | %.3f |",
    row$ensemble, row$pool, row$method,
    row$test_total_wis, row$test_mean_wis
  ))
}

writeLines(md_lines, OUT_REPORT)
cat(sprintf("Saved: %s\n", OUT_REPORT))

cat("\n=== Done ===\n")
