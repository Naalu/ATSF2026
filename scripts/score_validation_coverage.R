# =============================================================================
# scripts/score_validation_coverage.R
#
# Phase 5E — Validation-period prediction interval coverage
#
# Computes 50% and 95% prediction interval coverage for all four ensembles
# (BMA-primary, BMA-expanded, Equal-primary, Equal-expanded) on the
# 57-date validation window.
#
# Implementation note: uses an inner_join on (forecast cell, lower quantile)
# vs (forecast cell, upper quantile) rather than pivot_wider. This avoids
# float-precision issues with quantile_level values like 0.025 that may
# pivot to slightly-mismatched column names depending on input precision.
#
# Outputs:
#   - analysis/phase4/validation_coverage.csv
#
# Run from project root:
#   Rscript scripts/score_validation_coverage.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

VAL_START <- as.Date("2015-10-24")
VAL_END   <- as.Date("2017-05-06")

ENSEMBLE_DIRS <- list(
  bma_primary    = "analysis/phase4/staging/ensemble_bma_primary",
  bma_expanded   = "analysis/phase4/staging/ensemble_bma_expanded",
  equal_primary  = "analysis/phase4/staging/ensemble_equal_primary",
  equal_expanded = "analysis/phase4/staging/ensemble_equal_expanded"
)

OUT_CSV <- "analysis/phase4/validation_coverage.csv"

# -----------------------------------------------------------------------------
# Load truth
# -----------------------------------------------------------------------------

cat("=== Loading truth ===\n")
truth <- readr::read_csv("target-data/oracle-output.csv",
                         show_col_types = FALSE) |>
  dplyr::filter(target == "ili perc") |>
  dplyr::transmute(
    location,
    target_end_date = as.Date(target_end_date),
    observed        = oracle_value
  ) |>
  dplyr::distinct(location, target_end_date, .keep_all = TRUE)

# -----------------------------------------------------------------------------
# Load validation forecasts from one staging directory
# -----------------------------------------------------------------------------

#' Load validation-period quantile forecasts.
load_validation_forecasts <- function(staging_dir) {
  csv_files <- list.files(staging_dir, pattern = "\\.csv$", full.names = TRUE)
  origin_dates <- as.Date(stringr::str_extract(basename(csv_files),
                                                "^\\d{4}-\\d{2}-\\d{2}"))
  val_files <- csv_files[origin_dates >= VAL_START & origin_dates <= VAL_END]
  
  purrr::map_dfr(val_files, function(f) {
    od <- as.Date(stringr::str_extract(basename(f), "^\\d{4}-\\d{2}-\\d{2}"))
    readr::read_csv(f, show_col_types = FALSE) |>
      dplyr::mutate(origin_date = od)
  }) |>
    dplyr::filter(output_type == "quantile") |>
    dplyr::transmute(
      origin_date,
      location        = as.character(location),
      horizon         = as.integer(horizon),
      target_end_date = as.Date(target_end_date),
      quantile_level  = as.numeric(output_type_id),
      predicted       = as.numeric(value)
    )
}

# -----------------------------------------------------------------------------
# Compute coverage by join (avoids pivot precision issues)
# -----------------------------------------------------------------------------

#' Coverage for a given prediction interval level.
#'
#' Filters forecasts to the lower and upper quantiles, joins them on the
#' forecast cell key, then tests containment of the observed value.
#'
#' Uses near-equality (within 1e-6) on quantile_level to handle any
#' float-precision quirks in the input CSVs.
#'
#' @param forecasts Forecasts joined with truth (must have `observed`)
#' @param level     0.50 or 0.95
#' @return numeric scalar in [0,1]
compute_coverage <- function(forecasts, level) {
  alpha <- (1 - level) / 2
  q_lo  <- alpha
  q_hi  <- 1 - alpha
  
  cell_keys <- c("origin_date", "location", "horizon", "target_end_date")
  
  lower <- forecasts |>
    dplyr::filter(abs(quantile_level - q_lo) < 1e-6) |>
    dplyr::select(dplyr::all_of(cell_keys), observed,
                  q_lower = predicted)
  
  upper <- forecasts |>
    dplyr::filter(abs(quantile_level - q_hi) < 1e-6) |>
    dplyr::select(dplyr::all_of(cell_keys),
                  q_upper = predicted)
  
  joined <- lower |>
    dplyr::inner_join(upper, by = cell_keys) |>
    dplyr::filter(!is.na(observed), !is.na(q_lower), !is.na(q_upper))
  
  if (nrow(joined) == 0L) {
    return(NA_real_)
  }
  
  mean(joined$observed >= joined$q_lower & joined$observed <= joined$q_upper)
}

# -----------------------------------------------------------------------------
# Score all 4 ensembles
# -----------------------------------------------------------------------------

cat("\n=== Computing PI coverage per ensemble ===\n")

results <- purrr::imap_dfr(ENSEMBLE_DIRS, function(dir, label) {
  cat(sprintf("\n%s:\n", label))
  
  fc <- load_validation_forecasts(dir) |>
    dplyr::left_join(truth, by = c("location", "target_end_date"))
  
  # Diagnostic: what unique quantile levels are present?
  qls <- sort(unique(fc$quantile_level))
  cat(sprintf("  unique quantile levels: %d (range %.4f to %.4f)\n",
              length(qls), min(qls), max(qls)))
  
  cov50 <- compute_coverage(fc, level = 0.50)
  cov95 <- compute_coverage(fc, level = 0.95)
  
  cat(sprintf("  50%% PI coverage: %.1f%%\n", cov50 * 100))
  cat(sprintf("  95%% PI coverage: %.1f%%\n", cov95 * 100))
  
  tibble::tibble(
    ensemble       = label,
    pi_coverage_50 = cov50 * 100,
    pi_coverage_95 = cov95 * 100
  )
})

cat("\n=== Final Table ===\n\n")
results |>
  dplyr::mutate(dplyr::across(where(is.numeric), ~round(., 1))) |>
  print()

readr::write_csv(results, OUT_CSV)
cat(sprintf("\nSaved: %s\n", OUT_CSV))
